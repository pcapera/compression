# Chemical Reaction Simulation

This project simulates chemical reaction networks of oligomerization from monomers to dimers, with temperature-dependent kinetics and thermodynamics.

## Files Structure

### Core Files
- `parameters.json` - Configuration file containing all simulation parameters in JSON format
- `compression_main.py` - Main simulation script that reads parameters and saves results of species concentrations over time, and of number of species above threshold of detection
- `plot_concentrations.py` - Plotting script that reads simulation results and creates visualizations of concentrations over time
- `plot_number_species.py` - Plotting script that reads simulation results and creates visualizations of number of species above threshold of detection as a function of parameter values
- `concentrations.txt` - File containing simulation results (generated after running simulation) of time and concentration of each monomer and dimer species
- `number_species.txt` - File containing simulation results (generated after running simulation) of number of species above threshold, temperature, reaction time and concentration of compressor

## Usage

### 1. Configure Parameters
Edit `parameters.json` to modify simulation parameters:
- Physical constants (Boltzmann constant, Planck constant, gas constant)
- Thermodynamic parameters (enthalpies, entropies, randomness factors)
- Kinetic parameters (activation energies)
- Initial conditions (concentrations, temperature range)
- Simulation settings (time steps, detection thresholds)

### 2. Run Simulation
```bash
python3 compression_main.py
```
This will:
- Load parameters from `parameters.json`
- Run the simulation, with a loop over the specified temperature, reaction time or compressor concentration (if specified)
- Save all results to `concentrations.txt` and `number_species.txt`
- Print progress, summary statistics and plot concentrations (if indicated)

### 3. Generate Plots
```bash
python3 plot_concentrations.py
python3 plot_number_species.py
```
This will:
- Load results from `concentrations.txt` or `number_species.txt`
- Create various plots showing:
  - Evolution of concentrations of all species over time
  - Number of species depending on temperature, reaction time or compressor initial concentration
- Save all plots in `.svg` format

########################################################################################

## Parameters Description

### Simulation Parameters
- `num_species`: Number of molecular species (default: 12)
- `feeding_0`, `feeding_1`: Feeding conditions for species 0 and 1
- `precision_dt`: Time step precision factor
- `threshold_detection`: Minimum concentration for species detection
- `loop_range`: Temperature range for parameter sweep

### Physical Constants
- `kappa`: Boltzmann constant
- `planck`: Planck constant  
- `R`: Universal gas constant

### Thermodynamic Parameters
- `HM_base`, `SM_base`: Base enthalpy and entropy values
- `randomness_H`, `randomness_Ea`: Random variation factors
- `conversion_factor`: Unit conversion factor (cal to J)

### Initial Conditions
- `M_initial`: Initial monomer concentrations
- `P2_initial`: Initial dimer concentrations
- `IC`: Initial condition scaling factor

## Output Description

The `output.data` file contains:
- `all_simulations`: List of individual simulation results
- `summary`: Aggregated results across all simulations
- `parameters`: Copy of input parameters used
- `timestamp`: When the simulation was run

### Parameters to change

`connectivity`: '0' for isolates, '1' for partially connected, '2' for fully connected
`dominance_thermo`: enthalpy thermodynamic compressor, for thermodynamic dominance 'dominance_thermo'<'enthalpy'
`dominance_kinetic`: activation energy kinetic compressor, for kinetic dominance 'dominance_kinetic'<'activation_energy'
`IC`: initial concentration of compressor
`feeding_0`: feeding of thermodynamic compressor (0 means no, 1 means yes)
`feeding_1`: feeding of kinetic compressor (0 means no, 1 means yes)
`T`: temperature
`time_end`: total time

`plotting`: 1 for plotting, 0 for not plotting
`target_time`: time at which number of species should be measured
`num_species`: number of monomers
`threshold_detection`: minimal concentration to be considered measurable
`randomnes_H`: randomness associated to enthalpy
`randomnes_Ea`: randomness associated to entropy

`loop_start`: starting loop
`loop_end`: ending loop
`loop_increment`: increment loop

`precision_dt`: size of the time increments, it must be <(1*'max_concentration')

`IC_all_monomers`: initial concentration for all monomers
`enthalpy`: enthalpy for all monomers other than the thermodynamic compressor
`entropy`: entropy for all monomers
`activation_energy`: activation free energy for all monomers other than the kinetic compressor
`feeding_quantity`: feeding quantity, if feeding is applied
`IC_eps`: initial concentration for all dimmers
`water`: concentration of water	


#########################################################################################

## Dependencies

```bash
pip install numpy matplotlib pickle json datetime
```

## Notes

- Results are saved in binary format for efficient storage
- All plots are saved as SVG files
- The code handles both thermodynamic and kinetic dominance effects
- Random seed is set for reproducible results
