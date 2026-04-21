import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
import os
from pathlib import Path
from scipy.optimize import minimize

# --- 1. MAC/PHARMAPY STABILITY PATCHES ---
# Bypasses the missing modeling objects files on ARM64 installations
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
    """Calculates total raw material costs per batch."""
    raw_mat = sim.GetRawMaterials().filter(regex='mass_').sum()
    return (raw_mat * raw_costs).to_dict()

def get_constraints(sim, temp_k, n_batches):
    """Checks process constraints: Crystal Size > 40um, Production > 5kg/day."""
    mu = sim.CR01.result.mu_n[-1]
    mass_api = sim.F01.result.mass_cake_dry[-1]
    
    mean_size = (mu[1]/mu[0]) * 1e6  # Convert to microns
    mass_total = mass_api * n_batches
    
    # Cooling constraint: T must be decreasing (T_next - T_prev < 0)
    temp_constr = temp_k[1:] - temp_k[:-1] 
    return [40 - mean_size, 5 - mass_total] + temp_constr.tolist()

def make_non_verbose(runargs):
    for key in runargs: runargs[key]['verbose'] = False
    return runargs

# =============================================================================
# CORE DIGITAL TWIN CALLBACK
# =============================================================================

def callback_opt(x, simulate=False, raw_material_cost=None, return_augm=True, weights=None):
    if weights is None:
        weights = np.array([100, 100, 1, 1])
    
    # Mapping Decision Variables
    tau_R01, t1, t2, t3, time_CR01, deltaP = x[0], x[1], x[2], x[3], x[4], x[5]
    temp_CR01 = np.array([t1, t2, t3])

    path_phys = 'acetaminophen_pure_components_ziegler.json'
    flst = SimulationExec(path_phys, flowsheet='R01 --> HOLD01 --> CR01 --> F01')

    # --- 1. REACTOR (R01) ---
    # Verified Kinetics: A + B -> C + D
    kinetics = RxnKinetics(path=path_phys, rxn_list=['A + B --> C + D'], 
                           k_params=[2.79e9], ea_params=[70530.0], 
                           reformulate_kin=False)

    vol_liq = 1.0  # 1,000L Scale
    vol_flow = vol_liq / tau_R01
    c_in = np.array([0.33, 0.33, 0, 0, 55.5]) # 0.33M reactants in Water
    
    liquid_in = LiquidStream(path_phys, temp=333.15, mole_conc=c_in, vol_flow=vol_flow)
    flst.R01 = PlugFlowReactor(diam_in=0.05, num_discr=50, isothermal=True)
    flst.R01.Kinetics = kinetics
    flst.R01.Inlet = liquid_in
    flst.R01.Phases = (LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=vol_liq),)
    flst.R01.Utility = CoolingWater(mass_flow=0.5, temp_in=333.15)

    # --- 2. HOLDING TANK (HOLD01) ---
    flst.HOLD01 = DynamicCollector()
    # Mac Stabilizer: Start with 10L solvent to anchor energy balance
    flst.HOLD01.Phases = LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=0.01)

    # --- 3. CRYSTALLIZER (CR01) ---
    solub_cts = np.array([185.0, -1.45, 2.95e-3]) # Masuda 2026
    x_gr = np.geomspace(1, 1500, num=35)
    distrib_init = np.ones(35) * 1e-12 

    # Mac Stabilizers: Anchor with 10L liquid and 1kg solid seed
    liq_anchor = LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=0.01)
    sol_anchor = SolidPhase(path_phys, x_distrib=x_gr, distrib=distrib_init,
                            mass_frac=[0, 0, 1, 0, 0], temp=333.15, mass=1.0)

    temp_prog = np.array([[333.15, t1], [t1, t2], [t2, t3]], dtype=np.float64)
    lagrange_fn = PiecewiseLagrange(time_CR01, temp_prog)

    flst.CR01 = BatchCryst(target_comp='C', method='1D-FVM', scale=1e-9,
                           controls={'temp': lagrange_fn.evaluate_poly})
    flst.CR01.Kinetics = CrystKinetics(solub_cts, nucl_prim=(3e8, 0, 3), 
                                       nucl_sec=(4.46e10, 0, 2, 1e-5), 
                                       growth=(5, 0, 1.32), dissolution=(1, 0, 1))
    flst.CR01.Phases = (liq_anchor, sol_anchor)
    flst.CR01.Utility = CoolingWater(mass_flow=2.0, temp_in=275.15)
    # --- 3. CRYSTALLIZER (CR01) ---
    # solub_cts = np.array([185.0, -1.45, 2.95e-3]) 
    
    # # Ensure the cooling ramp starts at exactly 333.15
    # temp_prog = np.array([[333.15, t1], [t1, t2], [t2, t3]], dtype=np.float64)
    # lagrange_fn = PiecewiseLagrange(time_CR01, temp_prog)

    # flst.CR01 = BatchCryst(target_comp='C', method='1D-FVM', scale=1e-9,
    #                        controls={'temp': lagrange_fn.evaluate_poly})
    
    # flst.CR01.Kinetics = CrystKinetics(solub_cts, nucl_prim=(3e8, 0, 3), 
    #                                    nucl_sec=(4.46e10, 0, 2, 1e-5), 
    #                                    growth=(5, 0, 1.32), dissolution=(1, 0, 1))

    # # MAC STABILITY FIX: Initialize with ONLY the liquid phase.
    # # This prevents the 'MixedPhases' energy balance crash on ARM64.
    # # We use a 1.0 m3 volume to match the incoming reactor flow (reduces stiffness).
    # flst.CR01.Phases = LiquidPhase(path_phys, 333.15, mole_conc=[0,0,0,0,55.5], vol=1.0)
    
    # # Assign Utility AFTER the phases
    # flst.CR01.Utility = CoolingWater(mass_flow=2.0, temp_in=275.15)

    # --- 4. FILTER (F01) ---
    flst.F01 = Filter(0.8, 1e11, 1e10)

    # --- EXECUTION ---
    run_kwargs = make_non_verbose({
        'R01': {'runtime': 7200}, 'HOLD01': {'runtime': 7200},
        'CR01': {'runtime': time_CR01, 'sundials_opts': {'maxh': 60}},
        'F01': {'runtime': None, 'deltaP': deltaP}
    })

    if simulate:
        flst.SolveFlowsheet(kwargs_run=run_kwargs, verbose=True)
        return flst
    else:
        flst.SolveFlowsheet(kwargs_run=run_kwargs, verbose=False)
        t_cycle = max([flst.R01.result.time[-1], flst.CR01.result.time[-1], flst.F01.result.time[-1]])
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

    # Initial Guess: [tau, T1, T2, T3, t_cryst, deltaP]
    x_init = np.array([1800, 310, 295, 280, 5400, 202650])
    raw_costs = np.array([1.5, 1.2, 0, 0, 0.1]) 

    print("🚀 Running Baseline Digital Twin...")
    sim = callback_opt(x_init, simulate=True)
    
    sim.R01.plot_profiles(figsize=(7, 3.5))
    plt.savefig(output_dir / 'R01_baseline.png')
    sim.CR01.plot_profiles(figsize=(7, 7))
    plt.savefig(output_dir / 'CR01_baseline.png')

    print("\n⚖️ Starting Nelder-Mead Optimization...")
    bounds = [(600, 3600), (275, 320), (275, 320), (275, 320), (1800, 7200), (101325, 506625)]
    res = minimize(callback_opt, x0=x_init, method='Nelder-Mead', 
                   args=(False, raw_costs, True), options={'maxfev': 50}, bounds=bounds)

    print("\n✅ Optimization Found Parameters:")
    for name, val in zip(['tau_R01', 'T1', 'T2', 'T3', 't_cryst', 'deltaP'], res.x):
        print(f"{name}: {val:.2f}")

    # Final Plotting
    sim_opt = callback_opt(res.x, simulate=True)
    moms = sim_opt.CR01.result.mu_n
    plt.figure(figsize=(6, 4))
    plt.plot(sim_opt.CR01.result.time, (moms[:, 1]/moms[:, 0]) * 1e6)
    plt.ylabel('Mean Size (um)'); plt.title('Optimized Crystal Growth')
    plt.savefig(output_dir / 'optimized_size.png')
    
    print(f"📊 Results saved to {output_dir}")
    plt.show()