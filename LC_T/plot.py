import matplotlib.pyplot as plt
import pandas as pd

# Load the data, skipping comment lines
df = pd.read_csv('results_Fe_BCC.txt', comment='#', header=None,
                 names=['In-Temp', 'lx', 'ly', 'lz', 'Reported-Temp'])

# Create the plot
fig, ax1 = plt.subplots()

# Plot lattice parameters
ax1.plot(df['In-Temp'], df['lx'], label='lx', marker='o')
ax1.plot(df['In-Temp'], df['ly'], label='ly', marker='s')
ax1.plot(df['In-Temp'], df['lz'], label='lz', marker='^')

ax1.set_xlabel('Input Temperature (K)')
ax1.set_ylabel('Lattice Parameters (Å)')
ax1.set_title('Lattice Parameters vs Temperature')
ax1.grid(True)
ax1.legend(loc='upper left')

# Optional: plot Reported-Temp on a secondary y-axis
ax2 = ax1.twinx()
ax2.plot(df['In-Temp'], df['Reported-Temp'], label='Reported Temp', color='gray', linestyle='--')
ax2.set_ylabel('Reported Energy')
ax2.legend(loc='lower right')

plt.tight_layout()
plt.show()
