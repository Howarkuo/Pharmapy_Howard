#p-Aminophenol (A) + Acetic Anhydride (B) $\rightarrow$ Paracetamol (C) + Acetic Acid (D)
# temp=353 (80°C) , Paracetamol synthesis is typically done between 50°C and 95°C
# solvent: water 


#acetaminophen_pure_components_ziegler.json
# {
#     "A": {
#         "Name": "p-Aminophenol",
#         "Formula": "C6H7NO",
#         "CAS": "123-30-8",
#         "mw": 109.13,
#         "p_vap": [10.2, 2800, -50],
#         "rho_liq": 1130,
#         "cp_liq": [150.5, 0, 0, 0, 0],
#         "t_crit": 815.0,
#         "p_crit": 5370000.0
#     },
#     "B": {
#         "Name": "Acetic Anhydride",
#         "Formula": "C4H6O3",
#         "CAS": "108-24-7",
#         "mw": 102.09,
#         "p_vap": [9.2466, 1427.77, -75.11],
#         "rho_liq": 1082,
#         "cp_liq": [184.1, 0, 0, 0, 0],
#         "t_crit": 569.0,
#         "p_crit": 4000000.0
#     },
#     "C": {
#         "Name": "Paracetamol",
#         "Formula": "C8H9NO2",
#         "CAS": "103-90-2",
#         "mw": 151.16,
#         "p_vap": [11.5, 4500, -80],
#         "rho_liq": 1260,
#         "cp_liq": [210.0, 0, 0, 0, 0],
#         "t_crit": 850.0,
#         "p_crit": 4800000.0
#     },
#     "D": {
#         "Name": "Acetic Acid",
#         "Formula": "C2H4O2",
#         "CAS": "64-19-7",
#         "mw": 60.05,
#         "p_vap": [9.6821, 1642.54, -39.76],
#         "rho_liq": 1049,
#         "cp_liq": [124.3, 0, 0, 0, 0],
#         "t_crit": 592.0,
#         "p_crit": 5786000.0
#     },
#     "solvent": {
#         "Name": "Water",
#         "Formula": "H2O",
#         "CAS": "7732-18-5",
#         "mw": 18.015,
#         "p_vap": [10.1, 1687, -42],
#         "rho_liq": 1000,
#         "cp_liq": [75.3, 0, 0, 0, 0],
#         "t_crit": 647.1,
#         "p_crit": 22064000.0
#     }
# }

# acetaminophen_ziegler_conc_init.json

# {
#     "323K": [
#         1.0,
#         1.0,
#         0.0,
#         0.0,
#         0.0
#     ],
#     "343K": [
#         1.0,
#         1.0,
#         0.0,
#         0.0,
#         0.0
#     ],
#     "353K": [
#         1.0,
#         1.0,
#         0.0,
#         0.0,
#         0.0
#     ],
#     "363K": [
#         1.0,
#         1.0,
#         0.0,
#         0.0,
#         0.0
#     ],
#     "368K": [
#         1.0,
#         1.0,
#         0.0,
#         0.0,
#         0.0
#     ]
# }




from pathlib import Path
import matplotlib.pyplot as plt   
import os 
# Gets the full path to the file in the current directory
path_pure = Path.cwd() / "acetaminophen_pure_components_ziegler.json"
filename1 = "physicalstate-acetaminophen_water.png"
filename2 = "SensitivityMatrix-acetaminophen_water.png"
# You can even check if it exists before trying to open it
if path_pure.exists():
    print("File found!")


from PharmaPy.Phases import LiquidPhase
from PharmaPy.Kinetics import RxnKinetics
from PharmaPy.Reactors import BatchReactor


reactor = BatchReactor()

# path_pure = '../data/pure_components_ziegler.json'

# run simulation in 66 deg 
liquid = LiquidPhase(path_pure, temp=343, mole_conc=[1, 1, 0, 0, 0], name_solv='solvent')
#5 distinct entries (A, B, C, D, and solvent)
# Index	Component	Initial Value in your code	Description
# 0 A	1	Starting reactant
# 1	B	0	Intermediate/Product
# 2	C	0	Intermediate/Product
# 3	D	0	Final Product
# 4	solvent	0	Reaction medium
rxns = ['A + B --> C + D']
k_vals = [0.1] 
ea_vals = [0] 
kin = RxnKinetics(path_pure, rxn_list=rxns, k_params=k_vals, ea_params=ea_vals)

reactor.Phases = liquid
reactor.Kinetics = kin

time, states, sens = reactor.solve_unit(runtime=120, eval_sens=True)  # run for 2 minutes
# 4. Save Physical State (Concentration, Temp, Heat)
# reactor.plot_profiles returns a list of figures or handles internal axes
reactor.plot_profiles(figsize=(12, 4))
plt.savefig(filename1, dpi=300, bbox_inches='tight')
print(f"Saved: {filename1}")

# 5. Save Sensitivity Matrix
reactor.plot_sens() 
plt.savefig(filename2, dpi=300, bbox_inches='tight')
print(f"Saved: {filename2}")

# Optional: show them on screen as well
plt.show()




# (pharma_env) PS C:\Users\howar\Desktop\Pharmapy_Howard> python .\parameter_estimation_acetaminophen.py
# File found!
# Could not find GLIMDA.
# C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\Phases.py:203: RuntimeWarning: 'mass', 'moles' and 'vol' are all set to zero. Model may not perform as intended.
#   warnings.warn("'mass', 'moles' and 'vol' are all set to zero. "
# Final Run Statistics: --- 

#  Number of steps                                    : 93
#  Number of function evaluations                     : 117
#  Number of Jacobian evaluations                     : 2
#  Number of function eval. due to Jacobian eval.     : 0
#  Number of error test failures                      : 2
#  Number of nonlinear iterations                     : 114
#  Number of nonlinear convergence failures           : 0
#  Number of sensitivity evaluations                  : 117
#  Number of function eval. due to sensitivity eval.  : 0
#  Number of sensitivity error test failures          : 0

# Sensitivity options:

#  Method                   : SIMULTANEOUS
#  Difference quotient type : CENTERED
#  Suppress Sens            : False

# Solver options:

#  Solver                   : CVode
#  Linear multistep method  : BDF
#  Nonlinear solver         : Newton
#  Linear solver type       : DENSE
#  Maximal order            : 5
#  Tolerances (absolute)    : 1e-06
#  Tolerances (relative)    : 1e-06

# Simulation interval    : 0.0 - 120.0 seconds.
# Elapsed simulation time: 0.016072300000359974 seconds.
# Saved: physicalstate-acetaminophen_water.png
# Saved: SensitivityMatrix-acetaminophen_water.png




#SensitivityMatrix-acetaminophen_water.png
#  Sensitivity to $k_1$ ($\frac{\partial C_i}{\partial k_1}$)
# The Positive Humps ($C_C, C_D$): These represent your products, Acetaminophen and Acetic Acid. The positive values mean that as you increase $k_1$, you get more product. The peak around 15–20 seconds is the "Golden Window"—this is the time when your data is most useful for calculating $k_1$.
# The Negative Humps ($C_A, C_B$): These are your reactants, p-Aminophenol and Acetic Anhydride. As $k_1$ increases, their concentration drops faster (hence the negative sensitivity).


# Sensitivity to $E_{a,1}$ ($\frac{\partial C_i}{\partial E_{a,1}}$)
# The Magnitude: Notice the y-axis scale is $1e-5$. This is very small compared to the $k_1$ plot.The Meaning: This tells us that at this specific temperature (343K), changing $E_a$ slightly has almost no effect on the concentration. 