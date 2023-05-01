# Time-Position-Correlation-Calculator

This script calculates the correlations between properties (e.g., position, time) and data values (e.g., energy, pressure, number of hydrogen bonds) for data from molecular simulations (GROMACS, LAMMPS) or any time/position-correlated data. The script uses Welford's method to compute the standard deviation and bins the data to reduce noise and reveal more statistically robust patterns.

Was used for, e.g., computing [change in hydrogen bond number as a function of distance to a surface during non-equilibrium steered molecular dynamics](https://pubs.acs.org/doi/full/10.1021/acsnano.2c08627) or computing the [average dissipated adhesive energy during peel-off a mesoscale spatula model](https://onlinelibrary.wiley.com/doi/full/10.1002/smll.202201674).

## Input File Requirements

- `-f`: Data values file (e.g., time_pressure.xvg): The first column represents a property (e.g., time, position), and the following columns represent data values associated with that property. Column units depend on the data being analyzed (e.g., time in picoseconds, pressure in bar, stress in pascals).
- `-x`: Properties file (e.g., time_position.xvg): The first column represents one property (e.g., time), and the following columns represent other properties (e.g., position). Column units depend on the specific properties being analyzed (e.g., time in picoseconds, position in nanometers).

Ensure that the units of the input files are consistent with each other and the specific application.
Ensure that the first column in the data values file and the properties file is of the same property (e.g., either time or position)
The output file will contain correlations between the properties and data values, with the same units as the input files.

## Command-Line Options

- `-expl`: Toggle between exclusive property column counting (1) or non-exclusive (0)
- `-ncol_d`: Number of data value columns (e.g., energies, pressures, etc.)
- `-ncol_p`: Number of property columns (e.g., positions, times, etc.)
- `-bw`: Bin width for the histogram
- `-f`: Input file containing data values
- `-x`: Input file containing properties
- `-o`: Output file for the correlations between properties and data values

## Usage
```
EXEC -expl [arg0] -ncol_d [arg1] -ncol_p [arg2] -bw [arg3] -f [INPUT] -x [INPUT2] -o [OUTPUT]
```
## Example
```
python time-pos-corr-calc.py -expl 1 -ncol_d 1 -ncol_p 1 -bw 0.33 -f time_pressure.xvg -x time_position.xvg -o position_pressure.xvg
```
