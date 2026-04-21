# bug1 to run : /Users/kuochingchia/miniforge3/envs/pharma_env/bin/python process_optimization_acetaminophen.py
# ,which python will point to global python and fail
# cant run by calling python 
#bug2 mkdir -p /Users/kuochingchia/miniforge3/envs/pharma_env/lib/python3.11/site-packages/data
# cp /Users/kuochingchia/miniforge3/envs/pharma_env/lib/python3.11/site-packages/PharmaPy/data/minimum_modeling_objects.json /Users/kuochingchia/miniforge3/envs/pharma_env/lib/python3.11/site-packages/data/
# echo '{"special": {}, "Reactors": {}, "Crystallizers": {}, "Distillation": {}, "LiquidLiquidExtraction": {}, "SolidLiquidSep": {}, "Containers": {}, "Streams": {}}' > /Users/kuochingchia/miniforge3/envs/pharma_env/lib/python3.11/site-packages/data/minimum_modeling_objects.json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
import time
from copy import deepcopy
from pathlib import Path
from scipy.optimize import minimize

# --- THE ABSOLUTE TOP PHARMAPY PATCH ---
import PharmaPy.CheckModule
import PharmaPy.SimExec
# We need to import Reactors here to patch it specifically
import PharmaPy.Reactors 

def do_nothing(*args, **kwargs):
    return None

PharmaPy.CheckModule.check_modeling_objects = do_nothing
PharmaPy.SimExec.check_modeling_objects = do_nothing
PharmaPy.Reactors.check_modeling_objects = do_nothing
# ---------------------------------------

# PharmaPy Imports
from PharmaPy.Reactors import PlugFlowReactor
from PharmaPy.Crystallizers import BatchCryst
from PharmaPy.SolidLiquidSep import Filter
from PharmaPy.Containers import DynamicCollector
from PharmaPy.Streams import LiquidStream
from PharmaPy.Phases import LiquidPhase, SolidPhase
from PharmaPy.Kinetics import RxnKinetics, CrystKinetics
from PharmaPy.Utilities import CoolingWater
from PharmaPy.Interpolation import PiecewiseLagrange
from PharmaPy.SimExec import SimulationExec
from PharmaPy.Plotting import plot_function



# We must replace it in both the source and the execution module
PharmaPy.CheckModule.check_modeling_objects = do_nothing
PharmaPy.SimExec.check_modeling_objects = do_nothing
# --------------------------------

# =============================================================================
# 1. HELPER FUNCTIONS FOR OPTIMIZATION
# =============================================================================

def get_costs(sim, raw_costs):
    """Calculates total raw material costs per batch"""
    raw_mat = sim.GetRawMaterials()
    raw_mat = raw_mat.filter(regex='mass_').sum()
    raw_cost = raw_mat * raw_costs
    return raw_cost.to_dict()

def get_constraints(sim, temp_k, n_batches):
    """Calculates process constraints: Crystal Size > 40um, Production > 5kg/day"""
    mu = sim.CR01.result.mu_n[-1]
    mass_api = sim.F01.result.mass_cake_dry[-1]
    
    mean_size = mu[1]/mu[0] * 1e6
    mass_total = mass_api * n_batches
    
    # Negative values mean the constraint is satisfied
    temp_constr = temp_k[1:] - temp_k[:-1] # Ensures cooling (T decreasing)
    constraints = [40 - mean_size, 5 - mass_total] + temp_constr.tolist()
    
    return constraints

def make_non_verbose(runargs):
    for key in runargs:
        runargs[key]['verbose'] = False
    return runargs

# =============================================================================
# 2. THE DIGITAL TWIN CALLBACK (CORE LOGIC)
# =============================================================================

def callback_opt(x, simulate=False, raw_material_cost=None, return_augm=True, weights=None):
    if weights is None:
        weights = np.array([100, 100, 1, 1])
    
    # Mapping Decision Variables
    tau_R01 = x[0]      # Residence Time (s)
    temp_CR01 = x[1:4]  # 3-Point Cooling Ramp (K)
    time_CR01 = x[4]    # Crystallization Time (s)
    deltaP = x[-1]      # Filter Pressure (Pa)

    path_phys = 'acetaminophen_pure_components_ziegler.json'
    graph = 'R01 --> HOLD01 --> CR01 --> F01'
    flst = SimulationExec(path_phys, flowsheet=graph)

    # --- REACTOR (R01) CONFIGURATION ---
    # Using YOUR Stage 2 Kinetics
    rxns = ['A + B --> C + D']
    k_real = np.array([2.79e9])    # Pre-exponential A, from parameter_estimation_acetaminophen.py!
    ea_real = np.array([70530.0])  # Ea in J/mol, , from parameter_estimation_acetaminophen.py!
    
    kinetics = RxnKinetics(path=path_phys, rxn_list=rxns, 
                           k_params=k_real, ea_params=ea_real, 
                           reformulate_kin=False)

    # Scale to 1,000L (1 m^3)
    vol_liq = 1.0  
    vol_flow = vol_liq / tau_R01
    
    # Inlet: 4-aminophenol (A) + Acetic Anhydride (B) in Water
    c_in = np.array([0.33, 0.33, 0, 0, 55.5]) 
    liquid_in = LiquidStream(path_phys, temp=333.15, mole_conc=c_in, vol_flow=vol_flow)
    
    flst.R01 = PlugFlowReactor(diam_in=0.05, num_discr=50, isothermal=True)
    flst.R01.Kinetics = kinetics
    flst.R01.Inlet = liquid_in
    flst.R01.Phases = (LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=vol_liq),)
    flst.R01.Utility = CoolingWater(mass_flow=0.5, temp_in=333.15)

    # --- HOLDING TANK (HOLD01) ---
    flst.HOLD01 = DynamicCollector()
    # FIX: Give the collector a tiny initial volume of solvent (1 mL)
    flst.HOLD01.Phases = LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=1e-6)

    # --- CRYSTALLIZER (CR01) CONFIGURATION ---
    # Masuda 2026 "Sodium Acetate" Solubility (Salting Out Effect)
    # This significantly reduces residual API in filtrate
#     solub_cts = np.array([1.85e2, -1.45e0, 2.95e-3]) 
    
#     # x_gr = np.geomspace(1, 1500, num=35)
#     # solid_phase = SolidPhase(path_phys, x_distrib=x_gr, distrib=np.zeros_like(x_gr),
#     #                          mass_frac=[0, 0, 1, 0, 0])
#     x_gr = np.geomspace(1, 1500, num=35)
#     distrib_init = np.zeros_like(x_gr)

# # We give it 1e-6 m^3 (1 mL) so the energy balance doesn't divide by zero.
# # We set temp to 333.15 to match the liquid coming from R01.
#     solid_cry = SolidPhase(path_phys, x_distrib=x_gr, distrib=distrib_init,
#                        mass_frac=[0, 0, 1, 0, 0], 
#                        temp=333.15, 
#                        mass=1e-6)

#     temp_program = np.array([[333.15, temp_CR01[0]],
#                              [temp_CR01[0], temp_CR01[1]],
#                              [temp_CR01[1], temp_CR01[2]]], dtype=np.float64)
#     lagrange_fn = PiecewiseLagrange(time_CR01, temp_program)

#     flst.CR01 = BatchCryst(target_comp='C', method='1D-FVM', scale=1e-9,
#                            controls={'temp': lagrange_fn.evaluate_poly})
#     flst.CR01.Kinetics = CrystKinetics(solub_cts, 
#                                        nucl_prim=(3e8, 0, 3), 
#                                        nucl_sec=(4.46e10, 0, 2, 1e-5), 
#                                        growth=(5, 0, 1.32), 
#                                        dissolution=(1, 0, 1))
#     flst.CR01.Phases = solid_cry
#     flst.CR01.Utility = CoolingWater(mass_flow=2.0, temp_in=275.15)
    # --- CRYSTALLIZER (CR01) CONFIGURATION ---
    # Masuda 2026 "Sodium Acetate" Solubility
    solub_cts = np.array([1.85e2, -1.45e0, 2.95e-3]) 
    
    x_gr = np.geomspace(1, 1500, num=35)
    
    # FIX: Use a tiny 'flat' distribution instead of zeros 
    # This ensures distrib.sum() is never zero
    distrib_init = np.ones(35) * 1e-12 

    # FIX: Increased mass to 0.01 kg (10g) to anchor the energy balance
    solid_cry = SolidPhase(path_phys, x_distrib=x_gr, distrib=distrib_init,
                           mass_frac=[0, 0, 1, 0, 0], 
                           temp=333.15, 
                           mass=0.01,
                           ) 

    # Ensure the cooling ramp starts exactly at the incoming liquid temp (333.15)
    temp_program = np.array([[333.15, temp_CR01[0]],
                             [temp_CR01[0], temp_CR01[1]],
                             [temp_CR01[1], temp_CR01[2]]], dtype=np.float64)
    
    lagrange_fn = PiecewiseLagrange(time_CR01, temp_program)

    flst.CR01 = BatchCryst(target_comp='C', method='1D-FVM', scale=1e-9,
                           controls={'temp': lagrange_fn.evaluate_poly})
    
    flst.CR01.Kinetics = CrystKinetics(solub_cts, 
                                       nucl_prim=(3e8, 0, 3), 
                                       nucl_sec=(4.46e10, 0, 2, 1e-5), 
                                       growth=(5, 0, 1.32), 
                                       dissolution=(1, 0, 1))
    
    # Attach the robust phase
    # FIX: Create a tiny initial liquid phase for the crystallizer
    liq_init_cry = LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=0.1)
    # Attach BOTH phases (liquid and solid) to the crystallizer
    flst.CR01.Phases = (liq_init_cry, solid_cry)
    # flst.CR01.Phases = solid_cry
    flst.CR01.Utility = CoolingWater(mass_flow=2.0, temp_in=275.15)

    # --- FILTER (F01) ---
    # flst.F01 = Filter(area_diam=0.8, alpha=1e11, Rm=1e10)
    flst.F01 = Filter(0.8, 1e11, 1e10)

    # --- RUN EXECUTION ---
    run_kwargs = make_non_verbose({
        'R01': {'runtime': 7200},
        'HOLD01': {'runtime': 7200},
        'CR01': {'runtime': time_CR01, 'sundials_opts': {'maxh': 60}},
        'F01': {'runtime': None, 'deltaP': deltaP}
    })

    if simulate:
        flst.SolveFlowsheet(kwargs_run=run_kwargs, verbose=True)
        return flst
    else:
        flst.SolveFlowsheet(kwargs_run=run_kwargs, verbose=False)
        
        # Performance Metrics
        t_cycle = max([flst.R01.result.time[-1], flst.CR01.result.time[-1], flst.F01.result.time[-1]])
        n_batches = (24 * 3600) / t_cycle
        
        costs = get_costs(flst, raw_material_cost)
        constraints = get_constraints(flst, temp_CR01, n_batches)

        if return_augm:
            penalties = np.maximum(0, weights * constraints)
            total_cost = sum(list(costs.values())) * n_batches
            return total_cost + sum(penalties**2)
        else:
            return {'cost/batch': costs, 'size_constr': constraints[0], 
                    'production_constr': constraints[1], 'temp_constr': constraints[2:]}

# =============================================================================
# 3. EXECUTION: SIMULATION & OPTIMIZATION
# =============================================================================

# Initial Guess: [tau, T1, T2, T3, t_cryst, deltaP]
x_init = np.array([1800, 310, 295, 280, 5400, 202650])
raw_costs = np.array([1.5, 1.2, 0, 0, 0.1]) # Estimated costs for A, B, C, D, Solvent

print(" Running Baseline Digital Twin Simulation...")
sim = callback_opt(x_init, simulate=True)

# Save Baseline Plots
sim.R01.plot_profiles(figsize=(7, 3.5))
plt.savefig('stage1_results/R01_baseline.png')

sim.CR01.plot_profiles(figsize=(7, 7))
plt.savefig('stage1_results/CR01_baseline.png')

# Run Optimization
print("\n Starting Nelder-Mead Optimization (Process Design)...")
bounds = [(600, 3600), (275, 320), (275, 320), (275, 320), (1800, 7200), (101325, 506625)]
res = minimize(callback_opt, x0=x_init, method='Nelder-Mead', 
               args=(False, raw_costs, True), options={'maxfev': 50}, bounds=bounds)

print("\n Optimized Parameters Found:")
param_names = ['tau_R01', 'T1_CR01', 'T2_CR01', 'T3_CR01', 'time_CR01', 'deltaP_F01']
for name, val in zip(param_names, res.x):
    print(f"{name}: {val:.2f}")

# Final Optimized Run
sim_opt = callback_opt(res.x, simulate=True)

# Final Visualizations
print("\n Generating Final Industrial Reports...")
# 1. Crystallization Solubility vs Concentration
fi, ax = plot_function(sim_opt.CR01, state_names=['temp', ('mass_conc', ('C', ))], figsize=(8, 4), ncols=2)
ax[1].plot(sim_opt.CR01.result.time, sim_opt.CR01.result.solubility, '--k', label='Solubility')
ax[1].legend()
plt.savefig('stage1_results/optimized_solubility.png')

# 2. Mean Size over Time
moms = sim_opt.CR01.result.mu_n
plt.figure(figsize=(6, 4))
plt.plot(sim_opt.CR01.result.time, moms[:, 1]/moms[:, 0] * 1e6)
plt.xlabel('Time (s)')
plt.ylabel('Mean Size (um)')
plt.title('API Crystal Growth (Optimized)')
plt.savefig('stage1_results/optimized_size.png')

# 3. Filtration Performance
sim_opt.F01.plot_profiles(figsize=(7, 3.5))
plt.savefig('stage1_results/optimized_filtration.png')

print("✅ All plots and results saved to stage1_results/")
plt.show()



# (pharma_env) PS C:\Users\howar\Desktop\Pharmapy_Howard> conda run python process_optimization_acetaminophen.py
#  Running Baseline Digital Twin Simulation...

# ------------------------------
# Running R01
# ------------------------------


# Done!


# ------------------------------
# Running HOLD01
# ------------------------------


# Done!


# Could not find GLIMDA.
# C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\Phases.py:203: RuntimeWarning: 'mass', 'moles' and 'vol' are all set to zero. Model may not perform as intended.
#   warnings.warn("'mass', 'moles' and 'vol' are all set to zero. "
# Traceback (most recent call last):
#   File "C:\Users\howar\Desktop\Pharmapy_Howard\process_optimization_acetaminophen.py", line 208, in <module>
#     sim = callback_opt(x_init, simulate=True)
#           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "C:\Users\howar\Desktop\Pharmapy_Howard\process_optimization_acetaminophen.py", line 179, in callback_opt
#     flst.SolveFlowsheet(kwargs_run=run_kwargs, verbose=True)
#   File "C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\SimExec.py", line 139, in SolveFlowsheet       
#     connection.transfer_data()
#   File "C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\Connections.py", line 282, in transfer_data    
#     self.PassPhases()
#   File "C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\Connections.py", line 375, in PassPhases       
#     self.destination_uo.Phases = transfered_matter
#     ^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\Crystallizers.py", line 238, in Phases
#     self.Slurry.Phases = self._Phases
#     ^^^^^^^^^^^^^^^^^^
#   File "C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\MixedPhases.py", line 188, in Phases
#     self.temp = energy_balance(self, 'mass')
#                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\MixedPhases.py", line 47, in energy_balance    
#     temp_phase = newton(fun, np.mean(temps))
#                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\scipy\optimize\_zeros_py.py", line 368, in newton       
#     raise RuntimeError(msg)
# RuntimeError: Failed to converge after 50 iterations, value is nan.

# ERROR conda.cli.main_run:execute(125): `conda run python process_optimization_acetaminophen.py` failed. (See above for error) import numpy as np
# import matplotlib.pyplot as plt
# import pandas as pd
# import json
# import time
# from copy import deepcopy
# from pathlib import Path
# from scipy.optimize import minimize

# # PharmaPy Imports
# from PharmaPy.Reactors import PlugFlowReactor
# from PharmaPy.Crystallizers import BatchCryst
# from PharmaPy.SolidLiquidSep import Filter
# from PharmaPy.Containers import DynamicCollector
# from PharmaPy.Streams import LiquidStream
# from PharmaPy.Phases import LiquidPhase, SolidPhase
# from PharmaPy.Kinetics import RxnKinetics, CrystKinetics
# from PharmaPy.Utilities import CoolingWater
# from PharmaPy.Interpolation import PiecewiseLagrange
# from PharmaPy.SimExec import SimulationExec
# from PharmaPy.Plotting import plot_function

# # =============================================================================
# # 1. HELPER FUNCTIONS FOR OPTIMIZATION
# # =============================================================================

# def get_costs(sim, raw_costs):
#     """Calculates total raw material costs per batch"""
#     raw_mat = sim.GetRawMaterials()
#     raw_mat = raw_mat.filter(regex='mass_').sum()
#     raw_cost = raw_mat * raw_costs
#     return raw_cost.to_dict()

# def get_constraints(sim, temp_k, n_batches):
#     """Calculates process constraints: Crystal Size > 40um, Production > 5kg/day"""
#     mu = sim.CR01.result.mu_n[-1]
#     mass_api = sim.F01.result.mass_cake_dry[-1]
#     
#     mean_size = mu[1]/mu[0] * 1e6
#     mass_total = mass_api * n_batches
#     
#     # Negative values mean the constraint is satisfied
#     temp_constr = temp_k[1:] - temp_k[:-1] # Ensures cooling (T decreasing)
#     constraints = [40 - mean_size, 5 - mass_total] + temp_constr.tolist()
#     
#     return constraints

# def make_non_verbose(runargs):
#     for key in runargs:
#         runargs[key]['verbose'] = False
#     return runargs

# # =============================================================================
# # 2. THE DIGITAL TWIN CALLBACK (CORE LOGIC)
# # =============================================================================

# def callback_opt(x, simulate=False, raw_material_cost=None, return_augm=True, weights=None):
#     if weights is None:
#         weights = np.array([100, 100, 1, 1])
#     
#     # Mapping Decision Variables
#     tau_R01 = x[0]      # Residence Time (s)
#     temp_CR01 = x[1:4]  # 3-Point Cooling Ramp (K)
#     time_CR01 = x[4]    # Crystallization Time (s)
#     deltaP = x[-1]      # Filter Pressure (Pa)

#     path_phys = 'acetaminophen_pure_components_ziegler.json'
#     graph = 'R01 --> HOLD01 --> CR01 --> F01'
#     flst = SimulationExec(path_phys, flowsheet=graph)

#     # --- REACTOR (R01) CONFIGURATION ---
#     # Using YOUR Stage 2 Kinetics
#     rxns = ['A + B --> C + D']
#     k_real = np.array([2.79e9])    # Pre-exponential A
#     ea_real = np.array([70530.0])  # Ea in J/mol
#     
#     kinetics = RxnKinetics(path=path_phys, rxn_list=rxns, 
#                            k_params=k_real, ea_params=ea_real, 
#                            reformulate_kin=False)

#     # Scale to 1,000L (1 m^3)
#     vol_liq = 1.0  
#     vol_flow = vol_liq / tau_R01
#     
#     # Inlet: 4-aminophenol (A) + Acetic Anhydride (B) in Water
#     c_in = np.array([0.33, 0.33, 0, 0, 55.5]) 
#     liquid_in = LiquidStream(path_phys, temp=333.15, mole_conc=c_in, vol_flow=vol_flow)
#     
#     flst.R01 = PlugFlowReactor(diam_in=0.05, num_discr=50, isothermal=True)
#     flst.R01.Kinetics = kinetics
#     flst.R01.Inlet = liquid_in
#     flst.R01.Phases = (LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=vol_liq),)
#     flst.R01.Utility = CoolingWater(mass_flow=0.5, temp_in=333.15)

#     # --- HOLDING TANK (HOLD01) ---
#     flst.HOLD01 = DynamicCollector()

#     # --- CRYSTALLIZER (CR01) CONFIGURATION ---
#     # Masuda 2026 "Sodium Acetate" Solubility (Salting Out Effect)
#     # This significantly reduces residual API in filtrate
# #     solub_cts = np.array([1.85e2, -1.45e0, 2.95e-3]) 
#     
# #     # x_gr = np.geomspace(1, 1500, num=35)
# #     # solid_phase = SolidPhase(path_phys, x_distrib=x_gr, distrib=np.zeros_like(x_gr),
# #     #                          mass_frac=[0, 0, 1, 0, 0])
# #     x_gr = np.geomspace(1, 1500, num=35)
# #     distrib_init = np.zeros_like(x_gr)

# # # We give it 1e-6 m^3 (1 mL) so the energy balance doesn't divide by zero.
# # # We set temp to 333.15 to match the liquid coming from R01.
# #     solid_cry = SolidPhase(path_phys, x_distrib=x_gr, distrib=distrib_init,
# #                        mass_frac=[0, 0, 1, 0, 0], 
# #                        temp=333.15, 
# #                        mass=1e-6)

# #     temp_program = np.array([[333.15, temp_CR01[0]],
# #                              [temp_CR01[0], temp_CR01[1]],
# #                              [temp_CR01[1], temp_CR01[2]]], dtype=np.float64)
# #     lagrange_fn = PiecewiseLagrange(time_CR01, temp_program)

# #     flst.CR01 = BatchCryst(target_comp='C', method='1D-FVM', scale=1e-9,
# #                            controls={'temp': lagrange_fn.evaluate_poly})
# #     flst.CR01.Kinetics = CrystKinetics(solub_cts, 
# #                                        nucl_prim=(3e8, 0, 3), 
# #                                        nucl_sec=(4.46e10, 0, 2, 1e-5), 
# #                                        growth=(5, 0, 1.32), 
# #                                        dissolution=(1, 0, 1))
# #     flst.CR01.Phases = solid_cry
# #     flst.CR01.Utility = CoolingWater(mass_flow=2.0, temp_in=275.15)
#     # --- CRYSTALLIZER (CR01) CONFIGURATION ---
#     # Masuda 2026 "Sodium Acetate" Solubility
#     solub_cts = np.array([1.85e2, -1.45e0, 2.95e-3]) 
#     
#     x_gr = np.geomspace(1, 1500, num=35)
#     
#     # FIX: Use a tiny 'flat' distribution instead of zeros 
#     # This ensures distrib.sum() is never zero
#     distrib_init = np.ones(35) * 1e-12 

#     # FIX: Increased mass to 0.01 kg (10g) to anchor the energy balance
#     solid_cry = SolidPhase(path_phys, x_distrib=x_gr, distrib=distrib_init,
#                            mass_frac=[0, 0, 1, 0, 0], 
#                            temp=333.15, 
#                            mass=0.01) 

#     # Ensure the cooling ramp starts exactly at the incoming liquid temp (333.15)
#     temp_program = np.array([[333.15, temp_CR01[0]],
#                              [temp_CR01[0], temp_CR01[1]],
#                              [temp_CR01[1], temp_CR01[2]]], dtype=np.float64)
#     
#     lagrange_fn = PiecewiseLagrange(time_CR01, temp_program)

#     flst.CR01 = BatchCryst(target_comp='C', method='1D-FVM', scale=1e-9,
#                            controls={'temp': lagrange_fn.evaluate_poly})
#     
#     flst.CR01.Kinetics = CrystKinetics(solub_cts, 
#                                        nucl_prim=(3e8, 0, 3), 
#                                        nucl_sec=(4.46e10, 0, 2, 1e-5), 
#                                        growth=(5, 0, 1.32), 
#                                        dissolution=(1, 0, 1))
#     
#     # Attach the robust phase
#     flst.CR01.Phases = solid_cry
#     flst.CR01.Utility = CoolingWater(mass_flow=2.0, temp_in=275.15)

#     # --- FILTER (F01) ---
#     # flst.F01 = Filter(area_diam=0.8, alpha=1e11, Rm=1e10)
#     flst.F01 = Filter(0.8, 1e11, 1e10)

#     # --- RUN EXECUTION ---
#     run_kwargs = make_non_verbose({
#         'R01': {'runtime': 7200},
#         'HOLD01': {'runtime': 7200},
#         'CR01': {'runtime': time_CR01, 'sundials_opts': {'maxh': 60}},
#         'F01': {'runtime': None, 'deltaP': deltaP}
#     })

#     if simulate:
#         flst.SolveFlowsheet(kwargs_run=run_kwargs, verbose=True)
#         return flst
#     else:
#         flst.SolveFlowsheet(kwargs_run=run_kwargs, verbose=False)
#         
#         # Performance Metrics
#         t_cycle = max([flst.R01.result.time[-1], flst.CR01.result.time[-1], flst.F01.result.time[-1]])
#         n_batches = (24 * 3600) / t_cycle
#         
#         costs = get_costs(flst, raw_material_cost)
#         constraints = get_constraints(flst, temp_CR01, n_batches)

#         if return_augm:
#             penalties = np.maximum(0, weights * constraints)
#             total_cost = sum(list(costs.values())) * n_batches
#             return total_cost + sum(penalties**2)
#         else:
#             return {'cost/batch': costs, 'size_constr': constraints[0], 
#                     'production_constr': constraints[1], 'temp_constr': constraints[2:]}

# # =============================================================================
# # 3. EXECUTION: SIMULATION & OPTIMIZATION
# # =============================================================================

# # Initial Guess: [tau, T1, T2, T3, t_cryst, deltaP]
# x_init = np.array([1800, 310, 295, 280, 5400, 202650])
# raw_costs = np.array([1.5, 1.2, 0, 0, 0.1]) # Estimated costs for A, B, C, D, Solvent

# print(" Running Baseline Digital Twin Simulation...")
# sim = callback_opt(x_init, simulate=True)

# # Save Baseline Plots
# sim.R01.plot_profiles(figsize=(7, 3.5))
# plt.savefig('stage1_results/R01_baseline.png')

# sim.CR01.plot_profiles(figsize=(7, 7))
# plt.savefig('stage1_results/CR01_baseline.png')

# # Run Optimization
# print("\n Starting Nelder-Mead Optimization (Process Design)...")
# bounds = [(600, 3600), (275, 320), (275, 320), (275, 320), (1800, 7200), (101325, 506625)]
# res = minimize(callback_opt, x0=x_init, method='Nelder-Mead', 
#                args=(False, raw_costs, True), options={'maxfev': 50}, bounds=bounds)

# print("\n Optimized Parameters Found:")
# param_names = ['tau_R01', 'T1_CR01', 'T2_CR01', 'T3_CR01', 'time_CR01', 'deltaP_F01']
# for name, val in zip(param_names, res.x):
#     print(f"{name}: {val:.2f}")

# # Final Optimized Run
# sim_opt = callback_opt(res.x, simulate=True)

# # Final Visualizations
# print("\n Generating Final Industrial Reports...")
# # 1. Crystallization Solubility vs Concentration
# fi, ax = plot_function(sim_opt.CR01, state_names=['temp', ('mass_conc', ('C', ))], figsize=(8, 4), ncols=2)
# ax[1].plot(sim_opt.CR01.result.time, sim_opt.CR01.result.solubility, '--k', label='Solubility')
# ax[1].legend()
# plt.savefig('stage1_results/optimized_solubility.png')

# # 2. Mean Size over Time
# moms = sim_opt.CR01.result.mu_n
# plt.figure(figsize=(6, 4))
# plt.plot(sim_opt.CR01.result.time, moms[:, 1]/moms[:, 0] * 1e6)
# plt.xlabel('Time (s)')
# plt.ylabel('Mean Size (um)')
# plt.title('API Crystal Growth (Optimized)')
# plt.savefig('stage1_results/optimized_size.png')

# # 3. Filtration Performance
# sim_opt.F01.plot_profiles(figsize=(7, 3.5))
# plt.savefig('stage1_results/optimized_filtration.png')

# print("✅ All plots and results saved to stage1_results/")
# plt.show() why the original work  # Reactor definition# Physical properties of speciespath_phys = '../data/compound_database.json'# Initially define a flowsheet object for the simulation executive# This object will allow us to simulate the entire flowsheet.graph = 'R01 --> HOLD01 --> CR01 --> F01'flst = SimulationExec(path_phys, flowsheet=graph)# Specifying a stoichiometric matrixrxns = ['A + B --> C', 'A + C --> D']# Kinetics values; assuming deltaH values are 0.k_vals = np.array([2.654e4, 5.3e2])  # Pre-exponential factorea_vals = np.array([4.0e4, 3.0e4])  # Activation Energies# Use the imported RxnKinetics class to create a RxnKinetics instancekinetics = RxnKinetics(path=path_phys, rxn_list=rxns, k_params=k_vals, ea_params=ea_vals)##################vol_liq = 0.010tau_R01 = 1800vol_flow = vol_liq / tau_R01  # m**3 / sw_init = np.array([0, 0, 0, 0, 1])liquid_init = LiquidPhase(path_phys, 313.15, mass_frac=w_init, vol=vol_liq)# Cooling watertemp_set_R01 = 313.15cw = CoolingWater(mass_flow=0.1, temp_in=temp_set_R01)  # mass flow in kg/s# ---------- Inlet streams# Reactorc_in = np.array([0.33, 0.33, 0, 0, 0])temp_in = 40 + 273.15  # Kliquid_in = LiquidStream(path_phys, temp_in, mole_conc=c_in, vol_flow=vol_flow, name_solv='solvent')diam_in = 1 / 2 * 0.0254  # 1/2 inch in mflst.R01 = PlugFlowReactor(diam_in=diam_in, num_discr=50, isothermal=False)flst.R01.Utility = cwflst.R01.Kinetics = kineticsflst.R01.Inlet = liquid_inflst.R01.Phases = (liquid_init,)flst.R01.Utility = cwruntime_reactor = 3600 * 2# With a continuous unit follow by a batch unit, we need a collection unit,# a holding tank, to intermediately follow the the continuous unit.flst.HOLD01 = DynamicCollector()
# # Defining the crystallizerprim = (3e8, 0, 3)  # kP in #/m3/ssec = (4.46e10, 0, 2, 1e-5) # kS in #/m3/sgrowth = (5, 0, 1.32)  # kG in um/sdissol = (1, 0, 1) # kD in um/s# 0.2269 - 1.88e-3 * 323.15 + 3.89e-6 * 323.15 ** 2solub_cts = np.array([2.269e2, -1.88e0, 3.89e-3])x_gr = np.geomspace(1, 1500, num=35)distrib_init = np.zeros_like(x_gr)solid_cry = SolidPhase(path_phys, x_distrib=x_gr, distrib=distrib_init,
#                    mass_frac=[0, 0, 1, 0, 0])# Piecewise temperature profile definitiontemp_program = np.array([[313.15, 303.15],
#                          [303.15, 295.15],
#                          [295.15, 278.15],
#                         ], dtype=np.float64)# Initialize lagrange polynomial control objectruntime_cryst = runtime_reactor * 2.0lagrange_fn = PiecewiseLagrange(runtime_cryst, temp_program)# Defining the crystallizer with desired species 'C'flst.CR01 = BatchCryst(target_comp='C', method='1D-FVM', scale=1e-9, controls={'temp': lagrange_fn.evaluate_poly})flst.CR01.Kinetics = CrystKinetics(solub_cts, nucl_prim=prim, nucl_sec=sec, growth=growth, dissolution=dissol)flst.CR01.Utility = CoolingWater(mass_flow=1, temp_in=283.15)flst.CR01.Phases = solid_cry
# # Defining the filter.# Filteralpha = 1e11Rm = 1e10filt_area = 200  # cm**2diam = np.sqrt(4/np.pi * filt_area) / 100  # mflst.F01 = Filter(diam, alpha, Rm)
# # Gluing everything together and running the flowsheet# runargs for the flowsheetrunargs_R01 = {'runtime': runtime_reactor}sundials = {'maxh': 60}runargs_hold = {'runtime': runtime_reactor}runargs_CR01 = {'runtime': runtime_cryst, 'sundials_opts': sundials}runargs_F01 = {'runtime': None}# ---------- Running simulationrun_kwargs = {'R01': runargs_R01,
#               'HOLD01': runargs_hold,
#               'CR01': runargs_CR01,
#               'F01': runargs_F01
#               }flst.SolveFlowsheet(kwargs_run=run_kwargs)