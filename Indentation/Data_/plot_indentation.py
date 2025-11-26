#!/usr/bin/env python3
"""
Script to plot indentation data from GoodResult.xlsx

Plots:
- X-axis: Position Z
- Left axis: Force Z(*1.602) [nN]
- Right axis: Hardness [GPa]
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the Excel file
excel_file = 'GoodResult.xlsx'

try:
    df = pd.read_excel(excel_file)
    print("File loaded successfully!")
    print(f"Columns: {df.columns.tolist()}")
    print(f"\nData shape: {df.shape}")
    print(f"\nFirst few rows:")
    print(df.head())
    
except FileNotFoundError:
    print(f"Error: {excel_file} not found!")
    exit(1)
except Exception as e:
    print(f"Error reading file: {e}")
    exit(1)

# Try to identify the columns (case-insensitive matching)
# Common column name variations
position_col = None
position_z_col = None
force_z_col = None
hardness_col = None

# Search for columns (case-insensitive)
for col in df.columns:
    col_lower = str(col).lower()
    if 'position' in col_lower and 'z' not in col_lower and position_col is None:
        position_col = col
    if 'position' in col_lower and 'z' in col_lower:
        position_z_col = col
    if 'force' in col_lower and 'z' in col_lower:
        force_z_col = col
    if 'hardness' in col_lower:
        hardness_col = col

# If columns not found, show available columns
if position_col is None:
    print("\nWarning: 'Position' column not found. Available columns:")
    print(df.columns.tolist())
    # Try to find a column that might be Position
    for col in df.columns:
        col_lower = str(col).lower()
        if 'position' in col_lower and position_z_col is None:
            position_col = col
            break

if position_z_col is None:
    print("\nWarning: 'Position Z' column not found. Available columns:")
    print(df.columns.tolist())

if force_z_col is None:
    print("\nWarning: 'Force Z' column not found. Available columns:")
    print(df.columns.tolist())

if hardness_col is None:
    print("\nWarning: 'Hardness' column not found. Available columns:")
    print(df.columns.tolist())

if position_col is None or position_z_col is None or force_z_col is None or hardness_col is None:
    print("\nError: Could not identify required columns. Please check the Excel file structure.")
    print(f"Found - Position: {position_col}, Position Z: {position_z_col}, Force Z: {force_z_col}, Hardness: {hardness_col}")
    exit(1)

# Extract data
position = df[position_col].values
position_z = df[position_z_col].values
force_z = df[force_z_col].values
hardness = df[hardness_col].values

# Calculate X-axis: (Position Z - 878.902)
x_axis = (1150.014221- position_z)

# Calculate Force Z: Force Z * 1.602 [nN]
force_z_calculated = force_z * 1.602

# Remove any NaN values
mask = ~(np.isnan(x_axis) | np.isnan(force_z_calculated) | np.isnan(hardness))
x_axis = x_axis[mask]
force_z_calculated = force_z_calculated[mask]
hardness = hardness[mask]

# Create the plot with dual y-axes
fig, ax1 = plt.subplots(figsize=(12, 8))

# Left axis: Force Z
color1 = 'tab:blue'
ax1.set_xlabel('Position [Angstrom]', fontsize=12, fontweight='bold')
ax1.set_ylabel('Force [nN]', color=color1, fontsize=12, fontweight='bold')
line1 = ax1.plot(x_axis, force_z_calculated, color=color1, linewidth=2, label='Force Z')
ax1.tick_params(axis='y', labelcolor=color1)
ax1.grid(True, alpha=0.3)

# Right axis: Hardness
ax2 = ax1.twinx()
color2 = 'tab:red'
ax2.set_ylabel('Hardness [GPa]', color=color2, fontsize=12, fontweight='bold')
line2 = ax2.plot(x_axis, hardness, color=color2, linewidth=2, label='Hardness')
ax2.tick_params(axis='y', labelcolor=color2)

# Add title
plt.title('Indentation Analysis: Force Z and Hardness vs Position', fontsize=14, fontweight='bold', pad=20)

# Add legend
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper left', fontsize=10)

# Adjust layout
fig.tight_layout()

# Save the plot
output_file = 'indentation_plot.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nPlot saved as: {output_file}")

# Show the plot
plt.show()

print("\nPlotting complete!")

