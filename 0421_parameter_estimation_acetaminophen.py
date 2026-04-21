import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
import os
from pathlib import Path
from scipy.stats import linregress

# --- 1. THE MAC/PHARMAPY STABILITY PATCHES ---
import PharmaPy.CheckModule
import PharmaPy.SimExec
import PharmaPy.Reactors 

def do_nothing(*args, **kwargs):
    return None

PharmaPy.CheckModule.check_modeling_objects = do_nothing
PharmaPy.SimExec.check_modeling_objects = do_nothing
PharmaPy.Reactors.check_modeling_objects = do_nothing

# --- 2. PHARMAPY IMPORTS ---
from PharmaPy.Phases import LiquidPhase
from PharmaPy.Kinetics import RxnKinetics
from PharmaPy.Reactors import BatchReactor
from PharmaPy.SimExec import SimulationExec
from PharmaPy.Utilities import CoolingWater

# --- 3. PATH SETUP ---
path_pure = Path.cwd() / "acetaminophen_pure_components_ziegler.json"
path_init = Path.cwd() / "acetaminophen_ziegler_conc_init.json"
output_dir = Path.cwd() / "stage1_results"
output_dir.mkdir(exist_ok=True)

with open(path_init) as f:
    c_init_dict = json.load(f)

# Temperatures to fit: 323K (50C), 343K (70C), 353K (80C)
temps = [323, 343, 353]
results_list = []

# =============================================================================
# STAGE 1: MULTI-TEMPERATURE ISOTHERMAL FITTING
# =============================================================================
print("🚀 Starting Multi-Temperature Stage 1 Fitting...")
print("-" * 50)

for T in temps:
    print(f"🌡️ Processing T = {T}K...")
    
    # Load synthetic lab data for this temp
    csv_name = f"acetaminophen_data_{T}k.csv"
    data = np.genfromtxt(Path.cwd() / csv_name, delimiter=',')
    t_exp = data[:, 0]
    c_exp = data[:, 1:] # Columns for A and C

    # Configure Kinetics (Fitting local k, so Ea is set to 0 initially)
    kin = RxnKinetics(path_pure, rxn_list=['A + B --> C + D'], 
                      k_params=[0.01], ea_params=[0], 
                      reformulate_kin=False)

    # Initial Liquid Phase (Anchor: 1.0 m3)
    liquid = LiquidPhase(path_pure, temp=T, mole_conc=c_init_dict[f'{T}K'], 
                         name_solv='solvent', vol=1.0)

    reactor = BatchReactor()
    reactor.Kinetics = kin
    reactor.Phases = liquid
    # Utility anchor for solver stability
    reactor.Utility = CoolingWater(mass_flow=1.0, temp_in=T)

    # Parameter Estimation Engine
    sim = SimulationExec(path_pure, flowsheet='R01')
    sim.R01 = reactor
    
    # Tracking Species A (Index 0) and C (Index 2)
    sim.SetParamEstimation(x_data=t_exp, y_data=c_exp, 
                           measured_ind=(0, 2), optimize_flags=[True, False])
    sim.ParamInst.bounds = [(1e-6, 5.0)] 

    # Run Optimizer
    params_fitted, covar, info = sim.EstimateParams(verbose=False)
    k_val = params_fitted[0]
    
    results_list.append({
        'Temp_K': T,
        'Inv_Temp': 1/T,
        'k_value': k_val,
        'ln_k': np.log(k_val)
    })
    
    # Save fit visualization
    fig, ax = sim.ParamInst.plot_data_model(figsize=(6, 4))
    plt.title(f"Fit at {T}K")
    fig.savefig(output_dir / f"apap_fit_{T}K.png")
    plt.close(fig)

# Save Stage 1 results for Stage 2
df_results = pd.DataFrame(results_list)
df_results.to_csv(output_dir / "stage1_summary.csv", index=False)
print("\n✅ Stage 1 Complete. Summary saved.")

# =============================================================================
# STAGE 2: ARRHENIUS REGRESSION
# =============================================================================
print("\n🔬 Starting Stage 2: Arrhenius Regression...")

x = df_results['Inv_Temp'].values
y = df_results['ln_k'].values
slope, intercept, r_value, p_value, std_err = linregress(x, y)

R = 8.314  # Gas Constant
Ea = -slope * R
A = np.exp(intercept)

print("-" * 50)
print(f"✅ Activation Energy (Ea): {Ea/1000:.2f} kJ/mol")
print(f"✅ Pre-exponential (A):    {A:.2e} L/(mol*s)")
print(f"✅ R-squared:              {r_value**2:.4f}")
print("-" * 50)

# Final Plot
plt.figure(figsize=(8, 5))
plt.scatter(x * 1000, y, color='red', label='Isothermal Fits')
plt.plot(x * 1000, slope * x + intercept, color='blue', linestyle='--', label='Arrhenius Fit')
plt.xlabel('1000/T [K$^{-1}$]')
plt.ylabel('ln(k)')
plt.title('Global Kinetic Fit (Acetaminophen)')
plt.legend()
plt.savefig(output_dir / "apap_arrhenius_plot.png")
plt.show()


# e-packages/assimulo/lib/__init__.py)
# Could not find ODEPACK functions.
# Could not find RADAR5
# Could not find GLIMDA.
# 🚀 Starting Multi-Temperature Stage 1 Fitting...
# --------------------------------------------------
# 🌡️ Processing T = 323K...
# Optimization time: 3.76e-02 s.
# 🌡️ Processing T = 343K...
# Optimization time: 6.16e-02 s.
# 🌡️ Processing T = 353K...
# Optimization time: 1.30e-01 s.

# ✅ Stage 1 Complete. Summary saved.

# 🔬 Starting Stage 2: Arrhenius Regression...
# --------------------------------------------------
# ✅ Activation Energy (Ea): 70.53 kJ/mol
# ✅ Pre-exponential (A):    2.79e+09 L/(mol*s)
# ✅ R-squared:              0.9885
# --------------------------------------------------