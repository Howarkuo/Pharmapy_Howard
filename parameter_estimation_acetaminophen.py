# #### Standalone simulation Starts ###

# #p-Aminophenol (A) + Acetic Anhydride (B) $\rightarrow$ Paracetamol (C) + Acetic Acid (D)
# # temp=353 (80°C) , Paracetamol synthesis is typically done between 50°C and 95°C
# # solvent: water 
# # Goal: "What would happen if $k=0.1$?"
# # Action: You provide the $k$ value, and the computer gives you the concentration curves ($C_i$).
# # Result: The plots you generated earlier showing $A$ and $B$ dropping while $C$ and $D$ rise. You are simply predicting the future based on a guess.





# from pathlib import Path
# import matplotlib.pyplot as plt   
# import os 
# # Gets the full path to the file in the current directory
# path_pure = Path.cwd() / "acetaminophen_pure_components_ziegler.json"
# filename1 = "physicalstate-acetaminophen_water.png"
# filename2 = "SensitivityMatrix-acetaminophen_water.png"
# # You can even check if it exists before trying to open it
# if path_pure.exists():
#     print("File found!")


# from PharmaPy.Phases import LiquidPhase
# from PharmaPy.Kinetics import RxnKinetics
# from PharmaPy.Reactors import BatchReactor


# reactor = BatchReactor()

# # path_pure = '../data/pure_components_ziegler.json'

# # run simulation in 66 deg 
# liquid = LiquidPhase(path_pure, temp=343, mole_conc=[1, 1, 0, 0, 0], name_solv='solvent')
# #5 distinct entries (A, B, C, D, and solvent)
# # Index	Component	Initial Value in your code	Description
# # 0 A	1	Starting reactant
# # 1	B	0	Intermediate/Product
# # 2	C	0	Intermediate/Product
# # 3	D	0	Final Product
# # 4	solvent	0	Reaction medium
# rxns = ['A + B --> C + D']
# k_vals = [0.1] 
# ea_vals = [0] 
# kin = RxnKinetics(path_pure, rxn_list=rxns, k_params=k_vals, ea_params=ea_vals)

# reactor.Phases = liquid
# reactor.Kinetics = kin

# time, states, sens = reactor.solve_unit(runtime=120, eval_sens=True)  # run for 2 minutes
# # 4. Save Physical State (Concentration, Temp, Heat)
# # reactor.plot_profiles returns a list of figures or handles internal axes
# reactor.plot_profiles(figsize=(12, 4))
# plt.savefig(filename1, dpi=300, bbox_inches='tight')
# print(f"Saved: {filename1}")

# # 5. Save Sensitivity Matrix
# reactor.plot_sens() 
# plt.savefig(filename2, dpi=300, bbox_inches='tight')
# print(f"Saved: {filename2}")

# # Optional: show them on screen as well
# plt.show()




# # (pharma_env) PS C:\Users\howar\Desktop\Pharmapy_Howard> python .\parameter_estimation_acetaminophen.py
# # File found!
# # Could not find GLIMDA.
# # C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\Phases.py:203: RuntimeWarning: 'mass', 'moles' and 'vol' are all set to zero. Model may not perform as intended.
# #   warnings.warn("'mass', 'moles' and 'vol' are all set to zero. "
# # Final Run Statistics: --- 

# #  Number of steps                                    : 93
# #  Number of function evaluations                     : 117
# #  Number of Jacobian evaluations                     : 2
# #  Number of function eval. due to Jacobian eval.     : 0
# #  Number of error test failures                      : 2
# #  Number of nonlinear iterations                     : 114
# #  Number of nonlinear convergence failures           : 0
# #  Number of sensitivity evaluations                  : 117
# #  Number of function eval. due to sensitivity eval.  : 0
# #  Number of sensitivity error test failures          : 0

# # Sensitivity options:

# #  Method                   : SIMULTANEOUS
# #  Difference quotient type : CENTERED
# #  Suppress Sens            : False

# # Solver options:

# #  Solver                   : CVode
# #  Linear multistep method  : BDF
# #  Nonlinear solver         : Newton
# #  Linear solver type       : DENSE
# #  Maximal order            : 5
# #  Tolerances (absolute)    : 1e-06
# #  Tolerances (relative)    : 1e-06

# # Simulation interval    : 0.0 - 120.0 seconds.
# # Elapsed simulation time: 0.016072300000359974 seconds.
# # Saved: physicalstate-acetaminophen_water.png
# # Saved: SensitivityMatrix-acetaminophen_water.png




# #SensitivityMatrix-acetaminophen_water.png
# #  Sensitivity to $k_1$ ($\frac{\partial C_i}{\partial k_1}$)
# # The Positive Humps ($C_C, C_D$): These represent your products, Acetaminophen and Acetic Acid. The positive values mean that as you increase $k_1$, you get more product. The peak around 15–20 seconds is the "Golden Window"—this is the time when your data is most useful for calculating $k_1$.
# # The Negative Humps ($C_A, C_B$): These are your reactants, p-Aminophenol and Acetic Anhydride. As $k_1$ increases, their concentration drops faster (hence the negative sensitivity).


# # Sensitivity to $E_{a,1}$ ($\frac{\partial C_i}{\partial E_{a,1}}$)
# # The Magnitude: Notice the y-axis scale is $1e-5$. This is very small compared to the $k_1$ plot.The Meaning: This tells us that at this specific temperature (343K), changing $E_a$ slightly has almost no effect on the concentration. 

# #### Standalone simulation finished###

# ### Parameter Estimation Procedure
# ### Stage 1: Isothermal Parameter Estimation
# # In Ziegler Demo: the goal was to find parameters for four reactions ($i = 1, 2, 3, 4$) simultaneously by tracking the rise and fall of intermediates B and C.
# # Goal: Estimate : A_i and E_i (i = 1,2,3, 4)
# # time variable: B , C in different Temperature

# # In Acetaminophen Project, the complexity is lower (usually only one main reaction
# # Goal: for a single temperature run (e.g., 353K), your goal is to estimate the Local Rate Constant ($k$) for the acetylation.
# # using the Centered Arrhenius reformulation ($reformulate\_kin=True$), your objective is to find:$\varphi_{1, \text{APAP}}$: The log-rate constant at your reference temperature ($T_{ref}$).
# # fix $\varphi_{2}$ (which represents $E_a$) as a "constant guess" during Stage 1 and only let it move during Stage 2.
# # time variable: A (Reactant (Decreases over time)), C (Product (Increases over time))in different Temperature
# #  Isothermal Strategy:Run 1 (323K / 50°C): Fit apap_data_323K.csv $\rightarrow$ Get $\varphi_{1, 323}
# # $Run 2 (343K / 70°C): Fit apap_data_343K.csv $\rightarrow$ Get $\varphi_{1, 343}
# # $Run 3 (353K / 80°C): Fit apap_data_353K.csv $\rightarrow$ Get $\varphi_{1, 353}$

# #input: We don't have actually lab data, so Gemini pro 3 sythetic data 
# # acetaminophen_data_673k
#     # 0.0, 1.0, 0.0
#     # 2.0, 0.904, 0.092
#     # 5.0, 0.782, 0.216
#     # 10.0, 0.614, 0.395
#     # 15.0, 0.471, 0.523
#     # 20.0, 0.366, 0.625
#     # 25.0, 0.294, 0.720
#     # 35.0, 0.177, 0.825
#     # 45.0, 0.103, 0.894
#     # 55.0, 0.066, 0.928
#     # 75.0, 0.021, 0.973
#     # 95.0, 0.006, 0.991
#     # 115.0, 0.004, 0.991
#     # 145.0, 0.0, 1.0

# #time , A (Reactant (Decreases over time)), C (Product (Increases over time))
# # Parameter: 1000L water as solvent!
# import numpy as np
# import matplotlib.pyplot as plt
# import copy
# import json
# import os
# from pathlib import Path

# from PharmaPy.Phases import LiquidPhase
# from PharmaPy.Kinetics import RxnKinetics
# from PharmaPy.Reactors import BatchReactor
# from PharmaPy.SimExec import SimulationExec
# from PharmaPy.Utilities import CoolingWater # Using the class found in your local files
# # =============================================================================
# # 1. SETUP PATHS AND DATA
# # =============================================================================
# path_pure = Path.cwd() / "acetaminophen_pure_components_ziegler.json"
# path_init = Path.cwd() / "acetaminophen_ziegler_conc_init.json"
# # This is the lab data you will prepare (e.g., [Time, Reactant_A, Product_C])
# path_csv = Path.cwd() / "acetaminophen_data_353k.csv" 

# # Create output folder for Stage 1 results
# output_dir = Path.cwd() / "stage1_results"
# output_dir.mkdir(exist_ok=True)

# # Load initial concentrations
# with open(path_init) as f:
#     c_init_dict = json.load(f)

# # Load your experimental data
# # Expected format: Time(s), Conc_A(mol/L), Conc_C(mol/L)
# data = np.genfromtxt(path_csv, delimiter=',')
# t_exp = data[:, 0]
# c_exp = data[:, 1:]

# # =============================================================================
# # 2. CONFIGURE PHARMAPY OBJECTS
# # =============================================================================
# temp_ref = 353  
# k_guesses = [0.1]    # Standard k guess (avoids log bug)
# ea_guesses = [80000]  # Standard Ea guess

# rxns = ['A + B --> C + D']

# # Kinetics with Centered Arrhenius Reformulation
# kin = RxnKinetics(path_pure, rxn_list=rxns, k_params=k_guesses, 
#                   ea_params=ea_guesses, temp_ref=temp_ref, 
#                   reformulate_kin=False)

# # Initial Liquid Phase
# liquid = LiquidPhase(path_pure, temp=temp_ref, mole_conc=c_init_dict['353K'], 
#                      name_solv='solvent', vol=1.0)

# # Batch Reactor with Cooling Jacket (using the CoolingWater class)
# reactor = BatchReactor()
# reactor.Kinetics = kin
# reactor.Phases = liquid
# reactor.Utility = CoolingWater(vol_flow=0.01, temp_in=temp_ref)

# # =============================================================================
# # 3. PARAMETER ESTIMATION (THE OPTIMIZER)
# # =============================================================================
# # Create Simulation Engine FIRST
# sim = SimulationExec(path_pure, flowsheet='R01')

# # THEN attach the reactor
# sim.R01 = reactor


# # 4. Set the Target (Equation 7: Minimizing WSSE)
# # measured_ind=(0, 2) means your CSV tracks p-Aminophenol (0) and Acetaminophen (2)
# sim.SetParamEstimation(
#     x_data=t_exp, 
#     y_data=c_exp, 
#     measured_ind=(0, 2), 
#     optimize_flags=[True, False]
# )
# # Safety Bounds for phi_1 (ln_k)
# sim.ParamInst.bounds = [(-15, 5)]
# print(f"🚀 Starting Isothermal Optimization for Acetaminophen at {temp_ref}K...")
# params_fitted, covar, info = sim.EstimateParams()

# # =============================================================================
# # 4. RESULTS & VISUALIZATION
# # =============================================================================
# print("\n✅ Optimization Complete!")
# print("Estimated phi_1 (ln_k at T_ref):", params_fitted[0])

# # Save the fit plot
# fig, ax = sim.ParamInst.plot_data_model(figsize=(8, 6))
# plt.title(f"Acetaminophen Fit at {temp_ref}K")
# fig.savefig(output_dir / "apap_fit_353K.png")

# # Save the Sensitivity Matrix (to verify the "Golden Window")
# reactor.plot_sens()
# plt.savefig(output_dir / "apap_sens_353K.png")

# print(f"📊 Results saved to {output_dir}")
# plt.show()


# 🚀 Starting Isothermal Optimization for Acetaminophen at 353K...
# ------------------------------------------------------------
# eval    fun_val    ||step||   gradient   dampening_factor
# ------------------------------------------------------------
# 0       1.583e-01  ---        6.496e-01  7.520e-01
# 1       1.536e-01  8.552e-03  1.473e-01  5.299e-01 
# 2       1.534e-01  1.569e-03  1.936e-02  3.487e-01 
# 3       1.534e-01  2.152e-04  2.911e-03  1.178e-01 
# 4       1.534e-01  3.226e-05  4.360e-04  3.927e-02 
# 5       1.534e-01  4.840e-06  6.357e-05  1.309e-02 
# 6       1.534e-01  7.058e-07  6.357e-05  2.618e-02 
# 7       1.534e-01  7.057e-07  6.357e-05  1.047e-01 
# 8       1.534e-01  7.050e-07  6.357e-05  8.377e-01 
# 9       1.534e-01  6.994e-07  6.357e-05  1.340e+01 
# 10      1.534e-01  6.144e-07  6.357e-05  4.289e+02 
# 11      1.534e-01  1.225e-07  6.357e-05  2.745e+04 
# 12      1.534e-01  2.308e-09  6.357e-05  3.514e+06 
# ------------------------------------------------------------

# Optimization time: 6.35e-01 s.

# ✅ Optimization Complete!
# Estimated phi_1 (ln_k at T_ref): 0.0928289509018834

# ln(k_353) = ln (0.0928) = -2.377
# First data point for Arrhenius Plot (Exponential form)



# ### Multi Temparuture loop
# import numpy as np
# import matplotlib.pyplot as plt
# import json
# import pandas as pd
# from pathlib import Path

# from PharmaPy.Phases import LiquidPhase
# from PharmaPy.Kinetics import RxnKinetics
# from PharmaPy.Reactors import BatchReactor
# from PharmaPy.SimExec import SimulationExec
# from PharmaPy.Utilities import CoolingWater

# # --- SETUP ---
# path_pure = Path.cwd() / "acetaminophen_pure_components_ziegler.json"
# path_init = Path.cwd() / "acetaminophen_ziegler_conc_init.json"
# output_dir = Path.cwd() / "stage1_results"
# output_dir.mkdir(exist_ok=True)

# with open(path_init) as f:
#     c_init_dict = json.load(f)

# # The temperatures we want to process
# temps = [323, 343, 353]
# results_list = []

# print("🚀 Starting Multi-Temperature Stage 1 Fitting...")
# print("-" * 50)

# for T in temps:
#     print(f"🌡️ Processing T = {T}K...")
    
#     # 1. Load specific CSV for this temperature
#     csv_name = f"acetaminophen_data_{T}k.csv"
#     data = np.genfromtxt(Path.cwd() / csv_name, delimiter=',')
#     t_exp = data[:, 0]
#     c_exp = data[:, 1:]

#     # 2. Configure Objects
#     # Note: Using reformulate_kin=False for simpler isothermal k-value extraction
#     kin = RxnKinetics(path_pure, rxn_list=['A + B --> C + D'], 
#                       k_params=[0.01], ea_params=[0], 
#                       reformulate_kin=False)

#     liquid = LiquidPhase(path_pure, temp=T, mole_conc=c_init_dict[f'{T}K'], 
#                          name_solv='solvent', vol=1.0)

#     reactor = BatchReactor()
#     reactor.Kinetics = kin
#     reactor.Phases = liquid
#     reactor.Utility = CoolingWater(vol_flow=0.01, temp_in=T)

#     # 3. Parameter Estimation
#     sim = SimulationExec(path_pure, flowsheet='R01')
#     sim.R01 = reactor
    
#     sim.SetParamEstimation(x_data=t_exp, y_data=c_exp, 
#                            measured_ind=(0, 2), optimize_flags=[True, False])
#     sim.ParamInst.bounds = [(1e-6, 1.0)] 

#     # Run Optimizer
#     params_fitted, covar, info = sim.EstimateParams(verbose=False)
#     k_val = params_fitted[0]
    
#     # 4. Store Results
#     results_list.append({
#         'Temp_K': T,
#         'Inv_Temp': 1/T,
#         'k_value': k_val,
#         'ln_k': np.log(k_val)
#     })
    
#     # Save individual plots
#     fig, ax = sim.ParamInst.plot_data_model(figsize=(6, 4))
#     plt.title(f"Fit at {T}K")
#     fig.savefig(output_dir / f"apap_fit_{T}K.png")
#     plt.close(fig)

# # --- FINAL SUMMARY ---
# df_results = pd.DataFrame(results_list)
# print("\n✅ All Isothermal Fits Complete!")
# print("-" * 50)
# print(df_results[['Temp_K', 'k_value', 'ln_k']])

# # Save results for Stage 2
# df_results.to_csv(output_dir / "stage1_summary.csv", index=False)


# 🚀 Starting Multi-Temperature Stage 1 Fitting...
# --------------------------------------------------
# 🌡️ Processing T = 323K...
# Optimization time: 3.15e-01 s.
# 🌡️ Processing T = 343K...
# Optimization time: 5.76e-01 s.
# 🌡️ Processing T = 353K...
# Optimization time: 1.06e+00 s.

# ✅ All Isothermal Fits Complete!
# --------------------------------------------------
#    Temp_K   k_value      ln_k
# 0     323  0.010491 -4.557219
# 1     343  0.058249 -2.843034
# 2     353  0.092827 -2.37701


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from pathlib import Path

# 1. Load your results from the loop
data_path = Path.cwd() / "stage1_results" / "stage1_summary.csv"
df = pd.read_csv(data_path)

# 2. Linear Regression (ln(k) vs 1/T)
# y = ln(k)
# x = 1/T
x = df['Inv_Temp'].values
y = df['ln_k'].values

slope, intercept, r_value, p_value, std_err = linregress(x, y)

# 3. Extract Arrhenius Parameters
# Slope (m) = -Ea / R  => Ea = -m * R
# Intercept (b) = ln(A) => A = exp(b)
R = 8.314  # Gas Constant (J/mol*K)

Ea = -slope * R
A = np.exp(intercept)

print("🔬 FINAL ARRHENIUS PARAMETERS FOR ACETAMINOPHEN")
print("-" * 50)
print(f"✅ Activation Energy (Ea): {Ea/1000:.2f} kJ/mol")
print(f"✅ Pre-exponential (A):    {A:.2e} L/(mol*s)")
print(f"✅ R-squared (Fit Quality): {r_value**2:.4f}")
print("-" * 50)

# 4. The Arrhenius Plot
plt.figure(figsize=(8, 5))
plt.scatter(x * 1000, y, color='red', label='Stage 1 Isothermal Fits') # Plotting 1000/T for readability
plt.plot(x * 1000, slope * x + intercept, color='blue', linestyle='--', label='Arrhenius Regression')

plt.xlabel('Inverse Temperature ($1000/T$) [K$^{-1}$]')
plt.ylabel('Natural Log of Rate Constant ($\ln k$)')
plt.title('Stage 2: Global Kinetic Fit (Acetaminophen)')
plt.legend()
plt.grid(True, alpha=0.3)

# Save the final result
plt.savefig(Path.cwd() / "stage1_results" / "apap_arrhenius_plot.png", dpi=300)
plt.show()


# 🔬 FINAL ARRHENIUS PARAMETERS FOR ACETAMINOPHEN
# --------------------------------------------------
# ✅ Activation Energy (Ea): 70.53 kJ/mol
# ✅ Pre-exponential (A):    2.79e+09 L/(mol*s)
# ✅ R-squared (Fit Quality): 0.9885
# --------------------------------------------------