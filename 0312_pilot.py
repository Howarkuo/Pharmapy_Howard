import sys
import scipy.integrate
import numpy as np

# --- 1. PATCH SCIPY SIMPS ---
if not hasattr(scipy.integrate, 'simps'):
    from scipy.integrate import simpson
    scipy.integrate.simps = simpson
    print("🛠️  Monkey-patched: scipy.integrate.simps")

# --- 2. ROBUST IMPORTS ---
try:
    # We import only what is absolutely necessary
    from PharmaPy.Kinetics import RxnKinetics, CrystKinetics
    from PharmaPy.PharmaUnits import BatchReactor, Crystallizer
    print("✅ PharmaPy Simulation Engine Loaded!")
except ImportError as e:
    print(f"❌ Still an import error: {e}")
    sys.exit(1)

print("--- PHARMAPY PILOT START ---")
from IPython.display import display, Math
from tqdm import tqdm
import time
# Formatting the Arrhenius logic from your image 1
print("-" * 30)
print("REACTOR MODEL KINETICS")
print("Formula: k = A * exp(-Ea / (R * T))")
print(f"Parameters: A=1e8, Ea=45000, R=8.314")
print("-" * 30)

# --- 2. PROGRESS BAR (Simulation Prep) ---
# We use this to show the "Virtual Scale-up" is initializing
for i in tqdm(range(100), desc="Initializing Solvers"):
    time.sleep(0.005)
# 1. define system

# Unit operation class: 1L Batch Reactor
# Define the reaction path (Image 1)
rxn = Reaction(stoichiometry=[-1, 1], name='Reaction1')
kinetics = Kinetics(reactions=[rxn], pre_exponential=[1e8], activation_energy=[45000])

# Define the Reactor Setup (Image 2)
reactor = BatchReactor(kinetics=kinetics, volume=2.0, initial_concentrations={'A': 2.0})


# 2. run the reactor simulation 

# Simulate for 1000 seconds
res_reactor = reactor.solve(np.linspace(0, 1000, 100))
# The output concentration of A and B now goes to the holdtank/crystallizer
final_conc = res_reactor.iloc[-1]

# Step 3: Simulate the Crystallizer (Virtual Scale-up)

# Simplified Crystallizer setup
cryz = Crystallizer(volume=2.0, initial_concentrations=final_conc)
# Set your nucleation and growth coefficients here
cryz.nucleation_constant = 1e5 
cryz.growth_constant = 0.5 

res_cryz = cryz.solve(np.linspace(0, 3600, 100))


import matplotlib.pyplot as plt

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Reactor Progression (Concentration vs Time)
ax1.plot(res_reactor.index, res_reactor['A'], label='Reactant A', color='red')
ax1.plot(res_reactor.index, res_reactor['B'], label='Product B (Yield)', color='blue')
ax1.set_title("Reactor: Chemical Conversion")
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Concentration (mol/L)")
ax1.legend()
ax1.grid(True, linestyle='--')

# Plot 2: Crystallizer Result (matches your image 2)
# Assuming 'res_cryz' contains mass or moment data
ax2.plot(res_cryz.index, res_cryz['B'], label='Crystal Mass', color='green')
ax2.set_title("Crystallizer: Solid Formation")
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Mass (kg)")
ax2.grid(True, linestyle='--')

plt.tight_layout()
plt.show()