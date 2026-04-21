import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
from pathlib import Path
from scipy.optimize import minimize

# --- 1. MAC/PHARMAPY STABILITY PATCHES ---
import PharmaPy.CheckModule
import PharmaPy.SimExec
import PharmaPy.Reactors 

def do_nothing(*args, **kwargs): return None

PharmaPy.CheckModule.check_modeling_objects = do_nothing
PharmaPy.SimExec.check_modeling_objects = do_nothing
PharmaPy.Reactors.check_modeling_objects = do_nothing

# --- 2. PHARMAPY IMPORTS ---
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

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_costs(sim, raw_costs):
    raw_mat = sim.GetRawMaterials().filter(regex='mass_').sum()
    return (raw_mat * raw_costs).to_dict()

# def get_constraints(sim, temp_k, n_batches):
#     mu = sim.CR01.result.mu_n[-1]
#     mass_api = sim.F01.result.mass_cake_dry[-1]
#     # mean_size = (mu[1]/mu[0]) * 1e6
#     mean_size = np.divide(
#     moms[:,1],
#     moms[:,0],
#     out=np.full_like(moms[:,1], np.nan),
#     where=moms[:,0] > 1e-12) * 1e6
#     mean_size = np.clip(mean_size, 0, 500)
#     mass_total = mass_api * n_batches
#     temp_constr = temp_k[1:] - temp_k[:-1] 

#     return [40 - mean_size, 5 - mass_total] + temp_constr.tolist()
def get_constraints(sim, temp_k, n_batches):
    # 1. Access the final state of the moments
    # mu_n shape is usually (time_points, moment_index)
    moms_final = sim.CR01.result.mu_n[-1]
    
    # 2. Extract mu0 and mu1 for the final time point
    mu0_final = moms_final[0]
    mu1_final = moms_final[1]
    
    # 3. Handle the "No Crystals" or "Solver Failure" case
    # If mu0 is effectively zero, the batch failed. 
    # We return a large positive value for the constraints to trigger a heavy penalty.
    if mu0_final < 1e-10:
        # Returning large values ensures (Constraint > 0), which increases the penalty
        return [100.0, 100.0, 10.0, 10.0] 

    # 4. Calculate Mean Size safely (L1,0 in microns)
    # Use np.divide only on the final scalar values for the constraint
    mean_size = (mu1_final / mu0_final) * 1e6
    
    # Optional: Clip the mean size to prevent extreme values from warping the gradient
    mean_size = np.clip(mean_size, 0, 1000)

    # 5. Mass and Temperature logic
    mass_api = sim.F01.result.mass_cake_dry[-1]
    mass_total = mass_api * n_batches
    
    # Cooling constraint: ensures T_{i+1} <= T_i (Temperature must decrease)
    # temp_constr will be positive if temperature increases, which we want to penalize
    temp_constr = (temp_k[1:] - temp_k[:-1]).tolist() 

    # Return list of values where > 0 means the constraint is violated
    return [
        40 - mean_size,    # Positive if size < 40um
        5 - mass_total,   # Positive if total production < 5kg
        *temp_constr      # Positive if T increases
    ]

# =============================================================================
# CORE DIGITAL TWIN CALLBACK
# =============================================================================

# def callback_opt(x, simulate=False, raw_material_cost=None, return_augm=True, weights=None):
#     if weights is None: weights = np.array([100, 100, 1, 1])
    
#     # Mapping Decision Variables
#     tau_R01, t1, t2, t3, time_CR01, deltaP = x
#     temp_CR01 = np.array([t1, t2, t3])

#     path_phys = 'acetaminophen_pure_components_ziegler.json'
#     flst = SimulationExec(path_phys, flowsheet='R01 --> HOLD01 --> CR01 --> F01')
# def callback_opt(x, simulate=False, raw_material_cost=None, return_augm=True, weights=None):
#     if weights is None: weights = np.array([100, 100, 1, 1])
    
#     # STABILITY FIX 1: Explicitly force all variables to be pure Python floats.
#     # This prevents the "object too deep" error in the NumPy searchsorted function.
#     tau_R01   = float(x[0])
#     t1        = float(x[1])
#     t2        = float(x[2])
#     t3        = float(x[3])
#     time_CR01 = float(x[4])
#     deltaP    = float(x[5])
    
#     temp_CR01 = np.array([t1, t2, t3])
#     path_phys = 'acetaminophen_pure_components_ziegler.json'
#     flst = SimulationExec(path_phys, flowsheet='R01 --> HOLD01 --> CR01 --> F01')

#     # --- 1. REACTOR (R01) ---
#     kinetics = RxnKinetics(path=path_phys, rxn_list=['A + B --> C + D'], 
#                            k_params=[2.79e9], ea_params=[70530.0], 
#                            reformulate_kin=False)

#     vol_liq = 1.0  # 1,000L
#     vol_flow = vol_liq / tau_R01
#     c_in = np.array([0.33, 0.33, 0, 0, 55.5]) 
    
#     liquid_in = LiquidStream(path_phys, temp=333.15, mole_conc=c_in, vol_flow=vol_flow)
#     flst.R01 = PlugFlowReactor(diam_in=0.05, num_discr=50, isothermal=True)
#     flst.R01.Kinetics = kinetics
#     flst.R01.Inlet = liquid_in
#     flst.R01.Phases = (LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=vol_liq),)
#     flst.R01.Utility = CoolingWater(mass_flow=0.5, temp_in=333.15)

#     # --- 2. HOLDING TANK (HOLD01) ---
#     flst.HOLD01 = DynamicCollector()
#     # Anchor volume to prevent 'zero moles' warning
#     flst.HOLD01.Phases = LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=0.01)

#     # # --- 3. CRYSTALLIZER (CR01) ---
#     # solub_cts = np.array([185.0, -1.45, 2.95e-3]) 
#     # x_gr = np.geomspace(1, 1500, num=35)
    
#     # # Fix 1: Time and Temperature Grids for Lagrange
#     # t_grid = np.array([0.0, time_CR01/3.0, 2.0*time_CR01/3.0, time_CR01])
#     # T_grid = np.array([333.15, t1, t2, t3])
#     # lagrange_fn = PiecewiseLagrange(t_grid, T_grid)

#     # flst.CR01 = BatchCryst(target_comp='C', method='1D-FVM', scale=1e-9,
#     #                        controls={'temp': lagrange_fn.evaluate_poly})
    
#     # flst.CR01.Kinetics = CrystKinetics(solub_cts, nucl_prim=(3e8, 0, 3), 
#     #                                    nucl_sec=(4.46e10, 0, 2, 1e-5), 
#     #                                    growth=(5, 0, 1.32), dissolution=(1, 0, 1))

#     # # Fix 2: Use a LIST [ ] instead of a TUPLE ( ) to allow appending
#     # # Fix 3: Use industrial anchors (0.1m3 and 1e-6 mass) for Mac stability
#     # liq_anchor = LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=0.1)
#     # sol_anchor = SolidPhase(path_phys, x_distrib=x_gr, distrib=np.ones(35)*1e-12,
#     #                         mass_frac=[0, 0, 1, 0, 0], temp=333.15, mass=1e-6)

#     # flst.CR01.Phases = [liq_anchor, sol_anchor] 
#     # flst.CR01.Utility = CoolingWater(mass_flow=2.0, temp_in=275.15)
#     # --- 3. CRYSTALLIZER (CR01) ---
#     solub_cts = np.array([185.0, -1.45, 2.95e-3]) 
#     x_gr = np.geomspace(1, 1500, num=35)
#     distrib_init = np.ones(35) * 1e-12 

#     # STABILITY FIX 2: Use the standard Tutorial Format for PiecewiseLagrange.
#     # We create a 2D matrix of [StartTemp, EndTemp] for each segment.
#     # This is the most stable format for the PharmaPy interpolation engine.
#     temp_prog = np.array([
#         [333.15, t1],  # Segment 1: From 60°C to T1
#         [t1, t2],      # Segment 2: From T1 to T2
#         [t2, t3]       # Segment 3: From T2 to T3
#     ], dtype=np.float64)
    
#     lagrange_fn = PiecewiseLagrange(time_CR01, temp_prog)

#     flst.CR01 = BatchCryst(target_comp='C', method='1D-FVM', scale=1e-9,
#                            controls={'temp': lagrange_fn.evaluate_poly})
    
#     flst.CR01.Kinetics = CrystKinetics(solub_cts, nucl_prim=(3e8, 0, 3), 
#                                        nucl_sec=(4.46e10, 0, 2, 1e-5), 
#                                        growth=(5, 0, 1.32), dissolution=(1, 0, 1))

#     # STABILITY FIX 3: Use a LIST for phases (not a tuple)
#     # Use the MacBook Air "Industrial Anchors"
#     liq_anchor = LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=0.1)
#     sol_anchor = SolidPhase(path_phys, x_distrib=x_gr, distrib=distrib_init,
#                             mass_frac=[0, 0, 1, 0, 0], temp=333.15, mass=1e-6)

#     flst.CR01.Phases = [liq_anchor, sol_anchor] 
#     flst.CR01.Utility = CoolingWater(mass_flow=2.0, temp_in=275.15)

#     # --- 4. FILTER (F01) ---
#     flst.F01 = Filter(0.8, 1e11, 1e10)

#     # --- EXECUTION ---
#     run_kwargs = {'R01': {'runtime': 7200}, 'HOLD01': {'runtime': 7200},
#                   'CR01': {'runtime': time_CR01, 'sundials_opts': {'maxh': 60}},
#                   'F01': {'runtime': None, 'deltaP': deltaP}}

#     flst.SolveFlowsheet(kwargs_run=run_kwargs, verbose=simulate)
    
#     if simulate:
#         return flst
#     else:
#         t_cycle = max([flst.R01.result.time[-1], flst.CR01.result.time[-1], flst.F01.result.time[-1]])
#         n_batches = (24 * 3600) / t_cycle
#         costs = get_costs(flst, raw_material_cost)
#         constraints = get_constraints(flst, temp_CR01, n_batches)

#         if return_augm:
#             penalties = np.maximum(0, weights * constraints)
#             return (sum(costs.values()) * n_batches) + sum(penalties**2)
#         return {'cost': costs, 'constraints': constraints}

def callback_opt(x, simulate=False, raw_material_cost=None, return_augm=True, weights=None):
    if weights is None:
        weights = np.array([100, 100, 1, 1])

    tau_R01   = float(x[0])
    t1        = float(x[1])
    t2        = float(x[2])
    t3        = float(x[3])
    time_CR01 = float(x[4])
    deltaP    = float(x[5])

    temp_CR01 = np.array([t1, t2, t3])
    path_phys = 'acetaminophen_pure_components_ziegler.json'
    flst = SimulationExec(path_phys, flowsheet='R01 --> HOLD01 --> CR01 --> F01')

    # --- 1. REACTOR ---
    kinetics = RxnKinetics(
        path=path_phys,
        rxn_list=['A + B --> C + D'],
        k_params=[2.79e9],
        ea_params=[70530.0],
        reformulate_kin=False
    )

    vol_liq = 1.0
    vol_flow = vol_liq / tau_R01
    c_in = np.array([0.33, 0.33, 0, 0, 55.5])

    liquid_in = LiquidStream(path_phys, temp=333.15, mole_conc=c_in, vol_flow=vol_flow)
    flst.R01 = PlugFlowReactor(diam_in=0.05, num_discr=50, isothermal=True)
    flst.R01.Kinetics = kinetics
    flst.R01.Inlet = liquid_in
    flst.R01.Phases = (LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=vol_liq),)
    flst.R01.Utility = CoolingWater(mass_flow=0.5, temp_in=333.15)

    # --- 2. HOLD TANK ---
    flst.HOLD01 = DynamicCollector()
    flst.HOLD01.Phases = LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=0.01)

    # --- 3. CRYSTALLIZER ---
    solub_cts = np.array([185.0, -1.45, 2.95e-3])
    x_gr = np.geomspace(1, 1500, num=35)
    distrib_init = np.ones(35) * 1e-12

    temp_prog = np.array([
        [333.15, t1],
        [t1, t2],
        [t2, t3]
    ], dtype=np.float64)

    lagrange_fn = PiecewiseLagrange(time_CR01, temp_prog)

    flst.CR01 = BatchCryst(
        target_comp='C',
        method='1D-FVM',
        scale=1e-9,
        controls={'temp': lagrange_fn.evaluate_poly}
    )

    flst.CR01.Kinetics = CrystKinetics(
        solub_cts,
        # nucl_prim=(3e8, 0, 3),
        # nucl_sec=(4.46e10, 0, 2, 1e-5),
        #decrease 
        nucl_prim=(1e6, 0, 2),
        nucl_sec=(1e8, 0, 2, 1e-5),
        # growth=(5, 0, 1.32),
        growth=(1.0, 0, 1.0),
        dissolution=(1, 0, 1)
    )

    liq_anchor = LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=0.1)
    sol_anchor = SolidPhase(
        path_phys,
        x_distrib=x_gr,
        distrib=distrib_init,
        mass_frac=[0, 0, 1, 0, 0],
        temp=333.15,
        # mass=1e-6 
        mass=1e-3 #upgrade seed
    )

    flst.CR01.Phases = [liq_anchor, sol_anchor]
    flst.CR01.Utility = CoolingWater(mass_flow=2.0, temp_in=275.15)

    # --- 4. FILTER ---
    flst.F01 = Filter(0.8, 1e11, 1e10)

    run_kwargs = {
        'R01': {'runtime': 7200},
        'HOLD01': {'runtime': 7200},
        'CR01': {'runtime': time_CR01, 'sundials_opts': {'maxh': 60}},
        'F01': {'runtime': None, 'deltaP': deltaP}
    }

    try:
        flst.SolveFlowsheet(kwargs_run=run_kwargs, verbose=simulate)
    except Exception:
        if simulate:
            raise
        return 1e12

    if simulate:
        return flst

    t_cycle = max([
        flst.R01.result.time[-1],
        flst.CR01.result.time[-1],
        flst.F01.result.time[-1]
    ])
    n_batches = (24 * 3600) / t_cycle
    costs = get_costs(flst, raw_material_cost)
    constraints = get_constraints(flst, temp_CR01, n_batches)


    if return_augm:
        penalties = np.maximum(0, weights * constraints)
        return (sum(costs.values()) * n_batches) + sum(penalties**2)

    return {'cost': costs, 'constraints': constraints}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    output_dir = Path.cwd() / "stage1_results"
    output_dir.mkdir(exist_ok=True)

    x_init = np.array([1800, 310, 295, 280, 5400, 202650])
    raw_costs = np.array([1.5, 1.2, 0, 0, 0.1]) 

    print("🚀 Running Baseline Digital Twin...")
    sim = callback_opt(x_init, simulate=True)
    
    # if sim:
    #     sim.R01.plot_profiles(figsize=(7, 3.5))
    #     plt.savefig(output_dir / 'R01_baseline.png')
    #     sim.CR01.plot_profiles(figsize=(7, 7))
    #     plt.savefig(output_dir / 'CR01_baseline.png')
    if sim:
        try:
            sim.R01.plot_profiles(x_vals=np.linspace(0, 1.0, 50), figsize=(7, 3.5))
            plt.savefig(output_dir / 'R01_baseline.png')
            plt.close()
        except Exception as e:
            print(f"R01 plot skipped: {e}")

        try:
            sim.CR01.plot_profiles(figsize=(7, 7))
            plt.savefig(output_dir / 'CR01_baseline.png')
            plt.close()
        except Exception as e:
            print(f"CR01 plot skipped: {e}")


    print("\n⚖️ Starting Nelder-Mead Optimization...")
    res = minimize(callback_opt, x0=x_init, method='Nelder-Mead', 
                   args=(False, raw_costs, True), options={'maxfev': 50})

    print("\n✅ Optimization Complete.")
    for name, val in zip(['tau_R01', 'T1', 'T2', 'T3', 't_cryst', 'deltaP'], res.x):
        print(f"{name}: {val:.2f}")

    # Final Plotting
    sim_opt = callback_opt(res.x, simulate=True)
    moms = sim_opt.CR01.result.mu_n
    plt.figure(figsize=(6, 4))
    plt.plot(sim_opt.CR01.result.time, (moms[:, 1]/moms[:, 0]) * 1e6)
    plt.ylabel('Mean Size (um)'); plt.title('Optimized Crystal Growth')
    plt.savefig(output_dir / 'optimized_size.png')
    plt.show()

# ✅ Optimization Complete.
# tau_R01: 1813.33
# T1: 312.30
# T2: 293.36
# T3: 292.57
# t_cryst: 5440.00
# deltaP: 202180.90

# ------------------------------
# Running R01
# ------------------------------

# Final Run Statistics: --- 

#  Number of steps                                 : 195
#  Number of function evaluations                  : 241
#  Number of Jacobian*vector evaluations           : 366
#  Number of function eval. due to Jacobian eval.  : 241
#  Number of error test failures                   : 0
#  Number of nonlinear iterations                  : 238
#  Number of nonlinear convergence failures        : 1

# Solver options:

#  Solver                   : CVode
#  Linear multistep method  : BDF
#  Nonlinear solver         : Newton
#  Linear solver type       : SPGMR
#  Maximal order            : 5
#  Tolerances (absolute)    : 1e-06
#  Tolerances (relative)    : 1e-06

# Simulation interval    : 0.0 - 7200.0 seconds.
# Elapsed simulation time: 0.023334625002462417 seconds.

# Done!


# ------------------------------
# Running HOLD01
# ------------------------------

# /Users/kuochingchia/miniforge3/envs/pharma_env/lib/python3.11/site-packages/PharmaPy/Phases.py:203: RuntimeWarning: 'mass', 'moles' and 'vol' are all set to zero. Model may not perform as intended.
#   warnings.warn("'mass', 'moles' and 'vol' are all set to zero. "
# Final Run Statistics: --- 

#  Number of steps                                 : 83
#  Number of function evaluations                  : 117
#  Number of Jacobian evaluations                  : 5
#  Number of function eval. due to Jacobian eval.  : 35
#  Number of error test failures                   : 6
#  Number of nonlinear iterations                  : 114
#  Number of nonlinear convergence failures        : 3

# Solver options:

#  Solver                   : CVode
#  Linear multistep method  : BDF
#  Nonlinear solver         : Newton
#  Linear solver type       : DENSE
#  Maximal order            : 5
#  Tolerances (absolute)    : 1e-06
#  Tolerances (relative)    : 1e-06

# Simulation interval    : 0.0 - 7200.0 seconds.
# Elapsed simulation time: 0.009762457997567253 seconds.

# Done!


# ------------------------------
# Running CR01
# ------------------------------

# Final Run Statistics: --- 

#  Number of steps                                 : 181
#  Number of function evaluations                  : 221
#  Number of Jacobian*vector evaluations           : 744
#  Number of function eval. due to Jacobian eval.  : 221
#  Number of error test failures                   : 1
#  Number of nonlinear iterations                  : 218
#  Number of nonlinear convergence failures        : 7

# Solver options:

#  Solver                   : CVode
#  Linear multistep method  : BDF
#  Nonlinear solver         : Newton
#  Linear solver type       : SPGMR
#  Maximal order            : 5
#  Tolerances (absolute)    : 1e-06
#  Tolerances (relative)    : 1e-06

# Simulation interval    : 0.0 - 5440.0 seconds.
# Elapsed simulation time: 0.08197600000130478 seconds.

# Done!


# ------------------------------
# Running F01
# ------------------------------

# Terminating simulation at t = 9.839885 after signal from handle_event.
# Final Run Statistics: --- 

#  Number of steps                                 : 5
#  Number of function evaluations                  : 6
#  Number of Jacobian evaluations                  : 1
#  Number of function eval. due to Jacobian eval.  : 2
#  Number of error test failures                   : 0
#  Number of nonlinear iterations                  : 5
#  Number of nonlinear convergence failures        : 0
#  Number of state function evaluations            : 21
#  Number of state events                          : 1

# Solver options:

#  Solver                   : CVode
#  Linear multistep method  : BDF
#  Nonlinear solver         : Newton
#  Linear solver type       : DENSE
#  Maximal order            : 5
#  Tolerances (absolute)    : 1e-06
#  Tolerances (relative)    : 1e-06

# Simulation interval    : 0.0 - 9.839884946352415 seconds.
# Elapsed simulation time: 5.012499968870543e-05 seconds.

# Done!