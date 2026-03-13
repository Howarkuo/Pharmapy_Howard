# C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\data

# trial: 
# Reactant: hexane     C6H14 and CAS 110-54-3
# solvent: C4H8O and CAS 109-99-9, your solvent is Tetrahydrofuran (THF).

##### 1. the standalone simulation: BatchReactor, LiquidPhase, and RxnKinetics objects.####
##### Simulation with guesses. ########


# q1> what are elements of batch reaction modeling ?
# 2. how the sensitivity calculations in dynamic parameter estimation?

from pathlib import Path
import matplotlib.pyplot as plt   
import os 
# Gets the full path to the file in the current directory
path_pure = Path.cwd() / "pure_components_ziegler.json"
filename1 = "physicalstate-hexane_THF.png"
filename2 = "SensitivityMatrix-hexane_THF.png"
# You can even check if it exists before trying to open it
if path_pure.exists():
    print("File found!")


from PharmaPy.Phases import LiquidPhase
from PharmaPy.Kinetics import RxnKinetics
from PharmaPy.Reactors import BatchReactor


reactor = BatchReactor()

# path_pure = '../data/pure_components_ziegler.json'
liquid = LiquidPhase(path_pure, temp=673, mole_conc=[1, 0, 0, 0, 0], name_solv='solvent')
#5 distinct entries (A, B, C, D, and solvent)
# Index	Component	Initial Value in your code	Description
# 0 A	1	Starting reactant
# 1	B	0	Intermediate/Product
# 2	C	0	Intermediate/Product
# 3	D	0	Final Product
# 4	solvent	0	Reaction medium
rxns = ['A --> B', 'A --> C', 'B --> C', 'B --> D']
k_vals = [0.1] * 4
ea_vals = [0] * 4
kin = RxnKinetics(path_pure, rxn_list=rxns, k_params=k_vals, ea_params=ea_vals)

reactor.Phases = liquid
reactor.Kinetics = kin

time, states, sens = reactor.solve_unit(runtime=120, eval_sens=True)  # run for 2 minutes
# 4. Save Physical State (Concentration, Temp, Heat)
# reactor.plot_profiles returns a list of figures or handles internal axes
reactor.plot_profiles(figsize=(12, 4))
# plt.savefig(filename1, dpi=300, bbox_inches='tight')
print(f"Saved: {filename1}")

# 5. Save Sensitivity Matrix
reactor.plot_sens() 
# plt.savefig(filename2, dpi=300, bbox_inches='tight')
print(f"Saved: {filename2}")

# Optional: show them on screen as well
plt.show()




# C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\Phases.py:203: RuntimeWarning: 'mass', 'moles' and 'vol' are all set to zero. Model may not perform as intended.
#   warnings.warn("'mass', 'moles' and 'vol' are all set to zero. "
# Final Run Statistics: --- 

#  Number of steps                                    : 99
#  Number of function evaluations                     : 124
#  Number of Jacobian evaluations                     : 2
#  Number of function eval. due to Jacobian eval.     : 0
#  Number of error test failures                      : 1
#  Number of nonlinear iterations                     : 121
#  Number of nonlinear convergence failures           : 0
#  Number of sensitivity evaluations                  : 124
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
# Elapsed simulation time: 0.05655210000259103 seconds.
# ------------------------
# PharmaPy result object
# ------------------------

# Fields shown in the tables below can be accessed as result.<field>, e.g. result.mole_conc

# ---------------------------------------
# states      dim   units     index
# ---------------------------------------
# mole_conc   4     mol/L     A, B, C, D
# ---------------------------------------

# ------------------------------
# f(states)   dim   units   index
# ------------------------------
# q_rxn       1     W
# q_ht        1     W
# temp        1     K
# ------------------------------

# Time vector can be accessed as result.time

#######2. isothermal: parameter estimates are found for each isothermal dataset on an individual basis.#####
######  oved from "simulating" to "optimizing."#######
#### Levenberg-Marquardt optimizer.####
# The Physics: Use the ODEs from the BatchReactor (Eq. 2).
# The Chemistry: Use the stoichiometric map (A --> B, etc.) and the Rate Laws (Eq. 3) from your RxnKinetics.
# The Goal: Adjust the "knobs" ($k$ values) until the output of this model matches the experimental dots in your .csv files.

# Stage 1 : Fit $k$ at 673K, 723K, etc.
# Uses the model to find the best $k$ for each temp.
# Linear RegressionPlot $\ln(k)$ vs $1/T$.
# Converts Stage 1 results into $A$ and $E_a$ guesses.
# Stage 2
# Global Fitting.Uses $A$ and $E_a$ seeds to fit all data at once.

# stage1: one temp 673k at a time
# exp(-E_a / RT) is only a constant now 
# goal : estimate rate constant ($k$) for that specific temperature. You repeat this for all available temperatures (673K, 698K, 723K, etc.) to get a list of $k$ values.
# problem : hard to estimate difference between k's since A and E are large number 
# solution : centered arrhenius equation
# \ln(k_i) = \varphi_{1,i} + \exp(\varphi_{2,i}) \cdot \left( \frac{1}{T_{mean}} - \frac{1}{T} \right) \quad \forall i \in \{1, 2, 3, 4\}
# \varphi_{1,i} = \ln(A_i) - \frac{E_i}{RT_{mean}}
# \varphi_{2,i} = \ln \left( \frac{E_i}{R} \right)
# The log of the rate constant at a reference temperature.
# The activation energy divided by the gas constant.
# Why use this?
# It brings the parameters ($\phi_1, \phi_2$) and the Sensitivities ($S$) to a similar order of magnitude. When the values are closer together (e.g., between -5 and 20), the optimizer (the "GPS") can find the minimum error much more accurately without crashing.

# "Divide and Conquer" strategy used in chemical kinetic modeling. Instead of trying to guess everything at once, we solve the problem in pieces to ensure the computer doesn't get "lost" in the math.
import numpy as np
import matplotlib.pyplot as plt

# Let's import one of the datasets
# data = Path.cwd() / "ziegler_data_673K.csv"

data = np.genfromtxt('ziegler_data_673K.csv', delimiter=',')
t_exp = data[:, 0]
# Time column
c_exp = data[:, 1:]
# Concentration columns for B and C

print(t_exp.shape, c_exp.shape)

plt.plot(t_exp, c_exp, 'o', mfc='None')
plt.xlabel('time (s)')
plt.ylabel('$C_j$ (mol/L)')
plt.legend(('B', 'C'))

# Stage 1 to Stage 2 Bridge
#after you get your $\phi_1$ values (the $k$ at each temp) for all 5 temperatures, you will perform a Linear Regression:
# X-axis: $1/T$Y-axis: $\ln(k)$Slope: Gives you the first "real" estimate of $E_a$.Intercept: Gives you the first "real" estimate of $A$.

from PharmaPy.SimExec import SimulationExec
import numpy as np
import copy
import json

# Initial concentration
with open('ziegler_conc_init.json') as f:
    c_init = json.load(f)

sim = SimulationExec(path_pure, flowsheet='R01')

reactor = BatchReactor()

kin = RxnKinetics(path_pure, rxn_list=rxns, k_params=k_vals, ea_params=ea_vals, temp_ref=723, reformulate_kin=True)
kin_seed = copy.deepcopy(kin)
liquid = LiquidPhase(path_pure, mole_conc=c_init['673K'], name_solv='solvent')

reactor.Kinetics = kin
reactor.Phases = liquid

sim.R01 = reactor  # The aggregated instance has to have the same name as the names passed to `sim` through the 'flowsheet' argument

param_bools = [True] * 4 + [False] * 4  # (only k_i's are estimated)
sim.SetParamEstimation(x_data=t_exp, y_data=c_exp, measured_ind=(1, 2), optimize_flags=param_bools)

k_estimated, covar, info_di = sim.EstimateParams()

print(sim.ParamInst.paramest_df)

fig, axis = sim.ParamInst.plot_data_model(figsize=(6, 4))

reactor.plot_sens()

# Simulate the system with the parameter seeds and plot
reactor.reset()
reactor.Kinetics = kin_seed

reactor.solve_unit(time_grid=t_exp, verbose=False)
axis.plot(reactor.result.time, reactor.result.mole_conc[:, [1, 2]], 'k--', alpha=0.5)


from PharmaPy.SimExec import SimulationExec
import numpy as np
import copy
import json

# Initial concentration
with open('ziegler_conc_init.json') as f:
    c_init = json.load(f)
# 1. Create directory if it doesn't exist
output_dir = 'simulation_results'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)



sim = SimulationExec(path_pure, flowsheet='R01')

reactor = BatchReactor()

kin = RxnKinetics(path_pure, rxn_list=rxns, k_params=k_vals, ea_params=ea_vals, temp_ref=723, reformulate_kin=True)
kin_seed = copy.deepcopy(kin)
liquid = LiquidPhase(path_pure, mole_conc=c_init['673K'], name_solv='solvent')

reactor.Kinetics = kin
reactor.Phases = liquid

sim.R01 = reactor  # The aggregated instance has to have the same name as the names passed to `sim` through the 'flowsheet' argument

param_bools = [True] * 4 + [False] * 4  # (only k_i's are estimated)
sim.SetParamEstimation(x_data=t_exp, y_data=c_exp, measured_ind=(1, 2), optimize_flags=param_bools)

k_estimated, covar, info_di = sim.EstimateParams()

print(sim.ParamInst.paramest_df)
# 4. Generate Plot 1: modelfix.png
# This plot shows the model fit after parameter optimization
fig_fit, axis = sim.ParamInst.plot_data_model(figsize=(6, 4))
fig_fit.savefig(os.path.join(output_dir, 'modelfix.png'))

# 5. Generate Plot 2: SensitivityMatrix.png
reactor.plot_sens() 
fig_sens = plt.gcf()  # Capture the sensitivity figure
fig_sens.savefig(os.path.join(output_dir, 'SensitivityMatrix.png'))

# 6. Generate Plot 3: SanityCheck_673k.png
# We reset the reactor, run the original "seed" parameters, 
# and overlay them on the existing fit plot.
reactor.reset()
reactor.Kinetics = kin_seed
reactor.solve_unit(time_grid=t_exp, verbose=False)

# Overlay the seed results as black dashed lines
axis.plot(reactor.result.time, reactor.result.mole_conc[:, [1, 2]], 'k--', alpha=0.5, label='Seed')
axis.legend() # Optional: add legend to distinguish lines
fig_fit.savefig(os.path.join(output_dir, 'SanityCheck_673k.png'))

print(f"\nSuccess! 3 images generated in '{output_dir}':")
print("- modelfix.png")
print("- SensitivityMatrix.png")
print("- SanityCheck_673k.png")


# C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\Phases.py:203: RuntimeWarning: 'mass', 'moles' and 'vol' are all set to zero. Model may not perform as intended.
#   warnings.warn("'mass', 'moles' and 'vol' are all set to zero. "
# ------------------------------------------------------------
# eval    fun_val    ||step||   gradient   dampening_factor
# ------------------------------------------------------------
# 0       2.405e+00  ---        1.191e+00  3.529e-03
# 1       4.506e-01  2.874e+00  3.891e-01  1.246e-03 
# 2       2.045e-02  2.825e+00  4.999e-02  4.154e-04 
# 3       9.994e-03  2.770e-01  3.622e-03  1.385e-04 
# 4       9.900e-03  4.018e-02  1.204e-04  4.615e-05 
# 5       9.900e-03  2.421e-03  2.199e-06  1.538e-05 
# 6       9.900e-03  6.168e-05  5.942e-08  5.128e-06 
# 7       9.900e-03  2.158e-06  5.942e-08  1.026e-05 
# 8       9.900e-03  2.157e-06  5.942e-08  4.103e-05 
# 9       9.900e-03  2.152e-06  5.942e-08  3.282e-04 
# 10      9.900e-03  2.105e-06  5.942e-08  5.251e-03 
# 11      9.900e-03  1.563e-06  5.942e-08  1.680e-01 
# 12      9.900e-03  2.604e-07  5.942e-08  1.075e+01 
# ------------------------------------------------------------

# Optimization time: 2.58e-01 s.
#      obj_fun  \varphi_{1, 1}  \varphi_{1, 2}  \varphi_{1, 3}  \varphi_{1, 4}
# 0   2.405431       -2.302585       -2.302585       -2.302585       -2.302585
# 1   0.450650       -3.623486       -3.600098       -4.333372       -1.461591
# 2   0.020447       -4.400277       -4.694309       -4.615062       -3.931676
# 3   0.009994       -4.136805       -4.723977       -4.567379       -3.996344
# 4   0.009900       -4.165404       -4.725049       -4.595061       -4.001763
# 5   0.009900       -4.165386       -4.724137       -4.597293       -4.001548
# 6   0.009900       -4.165415       -4.724127       -4.597343       -4.001569
# 7   0.009900       -4.165415       -4.724126       -4.597344       -4.001569
# 8   0.009900       -4.165415       -4.724126       -4.597344       -4.001569
# 9   0.009900       -4.165415       -4.724126       -4.597344       -4.001569
# 10  0.009900       -4.165415       -4.724126       -4.597344       -4.001569
# 11  0.009900       -4.165415       -4.724127       -4.597344       -4.001569
# 12  0.009900       -4.165415       -4.724127       -4.597343       -4.001569
# C:\ProgramData\anaconda3\envs\pharma_env\Lib\site-packages\PharmaPy\Phases.py:203: RuntimeWarning: 'mass', 'moles' and 'vol' are all set to zero. Model may not perform as intended.
#   warnings.warn("'mass', 'moles' and 'vol' are all set to zero. "
# ------------------------------------------------------------
# eval    fun_val    ||step||   gradient   dampening_factor
# ------------------------------------------------------------
# 0       2.405e+00  ---        1.191e+00  3.529e-03
# 1       4.506e-01  2.874e+00  3.891e-01  1.246e-03 
# 2       2.045e-02  2.825e+00  4.999e-02  4.154e-04 
# 3       9.994e-03  2.770e-01  3.622e-03  1.385e-04 
# 4       9.900e-03  4.018e-02  1.204e-04  4.615e-05 
# 5       9.900e-03  2.421e-03  2.199e-06  1.538e-05 
# 6       9.900e-03  6.168e-05  5.942e-08  5.128e-06 
# 7       9.900e-03  2.158e-06  5.942e-08  1.026e-05 
# 8       9.900e-03  2.157e-06  5.942e-08  4.103e-05 
# 9       9.900e-03  2.152e-06  5.942e-08  3.282e-04 
# 10      9.900e-03  2.105e-06  5.942e-08  5.251e-03 
# 11      9.900e-03  1.563e-06  5.942e-08  1.680e-01 
# 12      9.900e-03  2.604e-07  5.942e-08  1.075e+01 
# ------------------------------------------------------------

# Optimization time: 2.80e-01 s.
#      obj_fun  \varphi_{1, 1}  \varphi_{1, 2}  \varphi_{1, 3}  \varphi_{1, 4}
# 0   2.405431       -2.302585       -2.302585       -2.302585       -2.302585
# 1   0.450650       -3.623486       -3.600098       -4.333372       -1.461591
# 2   0.020447       -4.400277       -4.694309       -4.615062       -3.931676
# 3   0.009994       -4.136805       -4.723977       -4.567379       -3.996344
# 4   0.009900       -4.165404       -4.725049       -4.595061       -4.001763
# 5   0.009900       -4.165386       -4.724137       -4.597293       -4.001548
# 6   0.009900       -4.165415       -4.724127       -4.597343       -4.001569
# 7   0.009900       -4.165415       -4.724126       -4.597344       -4.001569
# 8   0.009900       -4.165415       -4.724126       -4.597344       -4.001569
# 9   0.009900       -4.165415       -4.724126       -4.597344       -4.001569
# 10  0.009900       -4.165415       -4.724126       -4.597344       -4.001569
# 11  0.009900       -4.165415       -4.724127       -4.597344       -4.001569
# 12  0.009900       -4.165415       -4.724127       -4.597343       -4.001569

# Success! 3 images generated in 'simulation_results':
# - modelfix.png
# - SensitivityMatrix.png
# - SanityCheck_673k.png


#### 3.  parameter estimates are used as seed values in the second stage, where the all the datasets are used simultaneously to obtain representative parameter estimates and confidence intervals.
