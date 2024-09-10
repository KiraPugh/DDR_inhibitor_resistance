# DDR_inhibitor_resistance

The code base is an extension of Hamis et al. "Targeting Cellular DNA Damage Responses in Cancer: An In Vitro-Calibrated Agent-Based Model Simulating Monolayer and Spheroid Treatment Responses to ATR-Inhibiting Drugs." Bull Math Biol. 2021 Aug 30;83(10):103. doi: https://doi.org/10.1007/s11538-021-00935-y (GitHub reposiory: https://github.com/SJHamis/DDRinhibitors).


The codes in this repository are used to produce the in silico results using a cell crowding ABM which is implemented in C/C++. The file "main_update" is the main file where all model parameters can be modified. The remaining files are header files:
1. h_allocate_deallocate_update: this file allocates and deallocates memory for certain variables.
2. h_cell_cycle_progression_update: this file tracks cells through the cell cycle and includes the drug effects.
3. h_cell_proliferation_and_death_update: this file includes functions to place cells on the lattice and to put cells into the non-cycling state.
4. h_generate_output_update: this file creates outputs files.
5. h_lattice_setup_update: this file creates the initial cell configurations.

To compile the main file "main_update" into an executable file, run the command:
g++ -std=c++11 -O3 main_update.cpp -o secondrun.o

To execute the compiled file, run the command:
./secondrun.o

In the in silico experiments, we vary 6 things (List 1 in the main manuscript).
%
To vary these 6 items in the C++ code, change the values of the following parameters in the main file:
1. no_circles is the number of clusters in which the cells are seeded.
2. DR_fraction is the initial fraction of drug-resistant cells.
3. DR_type indicates which drug(s) the drug-resistant cells are resistant to (a value of 1, 2, 3 is drug 1, drug 2, and both drugs 1 and 2 resistance respectively).
4. drug_dose_ATRi and drug_dose_PARPi are the doses of drug 1 and drug 2 respectively, in micromole.
5. mu_in and mu_DR_in are the mean doubling times of the drug-sensitive and drug-resistant cells respectively.

Plotting the figures in the manuscript
- Fig. 1b: ODE_calibration
- Fig. 1d: ABM_calibration
- Figs. 3-5: plot_experiment1
- Fig. 6: plot_experiment2
- Fig. 7: plot_experiment3
Note that the latter 3 matlab codes also include plots for the dynamic drug-resistant fraction and total cell count over 310 hours.
