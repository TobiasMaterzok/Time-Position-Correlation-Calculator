"""
This script calculates the correlations between properties (e.g., position, time) and data values 
(e.g., energy, pressure, number of hydrogen bonds) for data from molecular simulations (GROMACS, LAMMPS) or any time/position-correlated data.
The script uses Welford's method to compute the standard deviation and bins the data to reduce noise and reveal patterns.
Input file requirements:
- f: Data values file (e.g., time_pressure.xvg): The first column represents a property (e.g., time, position), 
      and the following columns represent data values associated with that property.
      Column units depend on the data being analyzed (e.g., time in picoseconds, pressure in bar, stress in pascals).
- x: Properties file (e.g., time_position.xvg): The first column represents one property (e.g., time), and the 
      following columns represent other properties (e.g., position). 
      Column units depend on the specific properties being analyzed (e.g., time in picoseconds, position in nanometers).
Ensure that the units of the input files are consistent with each other and the specific application.
Ensure that the first column in the data values file and the properties file is of the same property (e.g., either time or position)
The output file will contain correlations between the properties and data values, with the same units as the input files.
"""

from fileio.fileread import *
import sys
import operator
import numpy as np
from iocontrol.options import get_options as getopt

# Define command-line options, types, and default values
options = ['-expl', '-ncol_d', '-ncol_p', '-bw', '-f', '-x', '-o']
types = ['int', 'int', 'int', 'float', 'str', 'str', 'str']
defaults = ['1', 1, 1, 0.33, "time_pressure.xvg", "time_position.xvg", "position_pressure.xvg"]

# -expl: toggle between exclusive property column counting (1) or non-exclusive (0)
# -ncol_d: number of data value columns (e.g., energies, pressures, etc.)
# -ncol_p: number of property columns (e.g., positions, times, etc.)
# -bw: bin width for the histogram
# -f: input file containing data values
# -x: input file containing properties
# -o: output file for the correlations between properties and data values

# Parse command-line arguments
try:
    expl, data_cols, prop_cols, bin_w, fin, fin2, fout = getopt(sys.argv, options, types, defaults)
except Exception as e:
    print("Error parsing command-line arguments:", e)
    sys.exit(1)

if expl != 0 and expl != 1:
    sys.exit("Error: Invalid value for -expl option. Allowed values are 0 or 1.")
if data_cols < 1:
    sys.exit("Error: Invalid value for -ncol_d option. Must be greater than or equal to 1.")
if prop_cols < 1:
    sys.exit("Error: Invalid value for -ncol_p option. Must be greater than or equal to 1.")

if len(sys.argv) < 2:
    print("Usage: EXEC -expl [arg0] -ncol_d [arg1] -ncol_p [arg2] -bw [arg3] -f [INPUT] -x [INPUT2] -o [OUTPUT]")
    print(" -expl [arg0]: Toggle between exclusive property column counting (1) or non-exclusive (0)")
    print(" -ncol_d [arg1]: Number of data value columns (e.g., energies, pressures, etc.)")
    print(" -ncol_p [arg2]: Number of property columns (e.g., positions, times, etc.)")
    print(" -bw [arg3]: Bin width for the histogram")
    print(" -f [INPUT]: Input file containing data values (e.g., time_pressure.xvg)")
    print(" -x [INPUT2]: Input file containing properties (e.g., time_position.xvg)")
    print(" -o [OUTPUT]: Output file for the correlations between properties and data values (e.g., position_pressure.xvg)")
    print("\n Example usage: EXEC -expl 1 -ncol_d 1 -ncol_p 1 -bw 0.33 -f time_pressure.xvg -x time_position.xvg -o position_pressure.xvg")
    sys.exit()

print("\nExlusiv: %d  Columns in data values file: %d  Columns in property file: %d" % (expl, data_cols, prop_cols))
if expl and (data_cols != prop_cols):
    print("If columns in %s are functions of columns in %s [arg1] and [arg2] have to be equal." % (fin, fin2))
    print("Are you sure that your properties in %s are explicit to each individual element or data point?" % (fin))
    sys.exit()


# Function to calculate standard deviation using Welford's online algorithm
# Welford's method is more numerically stable and less prone to round-off errors than naïve approaches
def stddev(size, arr):
    m_new, m_old, s_new, s_old = 0, 0, 0, 0
    m_new = m_old = arr[0]

    for i in range(1, int(size)):
        m_new = m_old + (arr[i] - m_old) / i
        s_new = s_old + (arr[i] - m_old) * (arr[i] - m_new)
        m_old, s_old = m_new, s_new

    return np.sqrt(s_new / (size - 1)) if size > 2 else 0

# Load data from input files (e.g., GROMACS .xvg or LAMMPS output)
try:
    prop_data = np.loadtxt(fin2)
    data = np.loadtxt(fin)
    prop_len = len(prop_data)
except Exception as e:
    print("Error loading input files:", e)
    sys.exit(1)

# Calculate maximum squared value and number of bins for histogram
max_sq = np.max(prop_data[:, 1]) ** 2
nbin = int(np.sqrt(max_sq) * 2 / bin_w)
inv_bin_w = 2.0 / bin_w
prop_count = np.zeros(nbin + 1)
data_count = np.zeros(nbin + 1)
value_arr = np.zeros((nbin + 10, prop_len + 10))


# Generate histogram for data with 1 property (e.g., position, time)
# This can be used to analyze data from simulations, experiments or any time/position-correlated data
# Binning is used to group data into intervals and reduce noise, simplifying the data for further analysis
for i in range(prop_len):
    if(expl):
        for j in range(1, prop_cols + 1):
            idx = int(prop_data[i, 1] * inv_bin_w)
            if(idx >= 0):
                prop_count[idx] += 1
                data_count[idx] += data[i, j]
                value_arr[idx, int(prop_count[idx] - 1)] = data[i, j]
    if not (expl):
        for j in range(1, data_cols + 1):
            idx = int(prop_data[i, 1] * inv_bin_w)
            if(j == 1):
                prop_count[idx] += 1
            data_count[idx] += data[i, j]
            value_arr[idx, int(prop_count[idx] - 1)] = data[i, j]


bin_count = int((nbin + 1) / 2)
res = np.zeros((bin_count, 3))

# Calculate correlations between properties and data values
for idx in range(bin_count):
    r = bin_w * idx
    if idx == 0 and prop_count[idx] >= 0:
        res[idx, 0] = r
        res[idx, 1] = data_count[idx].sum() / prop_count[idx]
        res[idx, 2] = stddev(prop_count[idx], value_arr[idx])
    count = prop_count[2 * idx - 1: 2 * idx + 1].sum()
    elif count > 0:
        data_sum = data_count[2 * idx - 1: 2 * idx + 1].sum()
        res[idx, 0] = r
        res[idx, 1] = data_sum / count
        res[idx, 2] = np.sqrt(np.power(stddev(prop_count[2 * idx - 1], value_arr[2 * idx - 1]), 2) +
                              np.power(stddev(prop_count[2 * idx], value_arr[2 * idx]), 2))

# Save output to a file
# The output file will contain correlations between properties and data values
np.savetxt(fout, res[res[:, 0] != 0], fmt='%8.8f')
