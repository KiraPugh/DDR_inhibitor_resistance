/*************
*
* This file includes functions to allocate and deallocate memory for data arrays.
* Data variables for lattice point-obesrvables are stored in pointers.
* Example: observable[n] denotes the value of a certain observable in lattice point n.
* This is (spatially) a 2D lattice, there are NN=N*N lattice points, (where N=100).
* Lattice point n is located on row n/N (integer division) and column n%N.
*
* Variable details:
* cell_cycle_phase[n] = 0 (no cell), 1 (G1), 2 (S), 3 (SD), 4 (G2/M), 5 (G0), 20 (NC)
*
*************/



#ifndef H_ALLOCATE_DEALLOCATE
#define H_ALLOCATE_DEALLOCATE

int *cell_cycle_phase, *cell_cycle_clock, *cell_cycle_clockSD, *cell_cycle_length, *cell_lives_here, *scheduled_lattice_removal, *cell_resistant_type, *cell_cycle_resistant_phase;
double *oxygen, *oxygen_scaled, *drug_ATRi, *drug_PARPi, *cell_cycle_length_G1, *cell_cycle_length_S, *cell_cycle_length_SD, *extra_G1_time, *extra_S_time, *extra_SD_time, *extra_G2_time;

void Allocate();
void Deallocate();

/** Allocate memory for data arrays */
void Allocate()
{

    cell_cycle_phase = (int*)malloc(NN*sizeof(int*));
    cell_cycle_clock = (int*)malloc(NN*sizeof(int*));
    cell_cycle_clockSD = (int*)malloc(NN*sizeof(int*));
    cell_cycle_length = (int*)malloc(NN*sizeof(int*));
    cell_cycle_length_G1 = (double*)malloc(NN*sizeof(double*));
    cell_cycle_length_S = (double*)malloc(NN*sizeof(double*));
    cell_cycle_length_SD = (double*)malloc(NN*sizeof(double*));
    extra_G1_time = (double*)malloc(NN*sizeof(double*));
    extra_S_time = (double*)malloc(NN*sizeof(double*));
    extra_G2_time = (double*)malloc(NN*sizeof(double*));
    extra_SD_time = (double*)malloc(NN*sizeof(double*));
    cell_lives_here = (int*)malloc(NN*sizeof(int*));
    oxygen_scaled = (double*)malloc(NN*sizeof(double*));
    drug_ATRi = (double*)malloc(NN*sizeof(double*));
    drug_PARPi = (double*)malloc(NN*sizeof(double*));
    scheduled_lattice_removal = (int*)malloc(NN*sizeof(int*));
	cell_resistant_type = (int*)malloc(NN*sizeof(int*));
	cell_cycle_resistant_phase = (int*)malloc(NN*sizeof(int*));
}

/** Deallocate memory for data arrays */
void Deallocate()
{
    free(cell_cycle_phase);
    free(cell_cycle_clock);
    free(cell_cycle_clockSD);
    free(cell_cycle_length);
    free(cell_cycle_length_G1);
    free(cell_cycle_length_S);
    free(cell_cycle_length_SD);
    free(extra_G1_time);
    free(extra_S_time);
    free(extra_G2_time);
    free(extra_SD_time);
    free(cell_lives_here);
    free(oxygen_scaled);
    free(drug_ATRi);
    free(drug_PARPi);
    free(scheduled_lattice_removal);
    free(cell_resistant_type);
	free(cell_cycle_resistant_phase);
}

#endif
