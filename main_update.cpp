/*********************************************************
 *
 * This file contains the main file.
 * The scheduling for drugs ceralasertib (ATRi) and olaparib (PARPi) are set here.
 *
 *********************************************************/

/** includes */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <string>
#include <cstring>
#include <vector>
#include <random>

/* defines */
#define N 100
#define NN 10000
#define dt 0.001
#define dx 0.2
#define dy 0.2
#define hour 1000

double dtodxy = dt/(dx*dx);
#define invitro 1

/* include header files */
#include "h_allocate_deallocate_update.h"
#include "h_cell_proliferation_and_death_update.h"
#include "h_lattice_setup_update.h"
#include "h_cell_cycle_progression_update.h"
#include "h_generate_output_update.h"



/** main function */
int main(int argc, char *argv[])


{
    // set the calibrated parameters
    double drug_dose_ATRi=0;//Drug dose of the ATRi (in micromole)
    double drug_dose_PARPi=0;//Drug dose of the PARPi (in micromole)

    // set the initial setup
    double DR_fraction = 0.3; //fraction of drug resistant cells
    double DR_type = 3; //0 no drug resistant cells, 1 resisrant to the ATRi, 2 resistant to the PARPi, 3 resistant to both ATRi and PARPi
    int no_circles = 482; //number of seeded clusters: choose to be one of the following: 1, 2, 4, 8, 16, 32, 64, 128, 256, 419, 482
    int circle_rad = CircleRadius(no_circles); //the radius of each cluster

    //initial fractions of each compartment
    double frac_G1_in = 0.36; //G1_0 in the manuscript
    double frac_S_in = 0.14; //S_0 in the manuscript
    double frac_SD_in = 0.04; //SD_0 in the manuscript
    double frac_G2M_in = 0.46; //G2/M_0 in the manuscript

    double mu_cells = 482; // mean number of cells (mu_P^*)
    double sigma_cells = 21; // sigma_P^*

    double mu_ds_in=41000; //mean cell cycle doubling time of drug sensitive cells
    double sigma_ds_in=8200;// standard deviation of cell cycle doubling time of drug sensitive cells
    double mu_dr_in=82000; //mean cell cycle doubling time of drug resistant cells
    double sigma_dr_in=8200;// standard deviation of cell cycle doubling time of drug resistant cells

    double mu_G1_in = 0.36; //mean fraction of time spent in G1 (tau_G1 in the manuscipt)
    double mu_S_in = 0.18; //mean fraction of time spent in S (tau_S in the manuscipt)
    double mu_SD_in = 0.04; //mean fraction of time spent in SD (tau_SD in the manuscipt)

    SetInputOptions1(mu_ds_in, sigma_ds_in, mu_dr_in, sigma_dr_in, mu_G1_in, mu_S_in, mu_SD_in);

    /* set the drug parameters */
    double Emax_ATRi_in = 1; // Emax_1 in the manuscript
    double EC50_ATRi_in=0.2579; // EC50_1 in the manuscript
    double hc_ATRi_in = 1.5187; // h_1 in the manuscript
    double Emax_PARPi_in = 0.5609; // Emax_2 in the manuscript
    double EC50_PARPi_in=0.1275; // EC50_2 in the manuscript
    double hc_PARPi_in = 1.0962; // h_2 in the manuscript
    double prob_healthy_in = 0.3558; //p_0 in the manuscript
    double prob_repair_in = 1; //q_0 in the manuscript

    int drug_clear_time=0; //assume no drug clearance

    SetInputOptions2(Emax_ATRi_in, Emax_PARPi_in, hc_ATRi_in, EC50_ATRi_in, hc_PARPi_in, EC50_PARPi_in, prob_healthy_in, prob_repair_in);

    int start_time, end_time, drug_time;
    start_time=0;
    end_time=1000*hour;
    drug_time=0;

    int n_caption=100;
    int no_cells=0;
    int no_alive_cells=0;
    srand(time(0));

    /** Allocate Grid variables [h_allocate_deallocate.h] */
    Allocate();


    // constant oxygen in space and time
    for(int n=0; n<NN; n++)
    {
        oxygen_scaled[n]=100;
    }

    /** Place the first cancer cell(s) on the grid [h_lattice_setup.h] */
    place_cells(no_circles, circle_rad, DR_type, DR_fraction, frac_G1_in, frac_S_in, frac_SD_in, frac_G2M_in);

    int init_number_of_cells = InitialCellCount(mu_cells,sigma_cells);
    printf("initial number of cells %d \n", init_number_of_cells);

    /** Time loop */
    for(int t=0; t <= end_time; t++)
    {
        /** Progress cells in cell-cycle [h_cell_cyle.h] */

        Compute_CellCycle(t);
        CheckResistantPhase(t);

        if( t>0 && start_time==0 )
        {
            no_cells=CountCells(t);


            if(no_cells>=init_number_of_cells)
            {
                start_time=t;
                drug_time=t;
                drug_clear_time=drug_time+500*hour; //no drug elimination
                end_time=start_time+310*hour;
                printf("Time: %d is start time.\n Time: %d is end time. drug_clear_time %d \n ", start_time, end_time, drug_clear_time );
                //CountCellsDetailed(t);
            }
        }

        if(start_time>0)
        {

            if(t==drug_time)
            {
                // Drugs are constant in space and time
                for(int n=0; n<NN; n++)
                {
                    drug_ATRi[n]=drug_dose_ATRi;
                    drug_PARPi[n]=drug_dose_PARPi;
                }
            }


            if(t>=drug_time )
            {
                if( (t-start_time)%(1*hour)==0 )
                {

                    no_cells=CountCells(t);
                    no_alive_cells = CountCyclingRepairableCells(t);

                    CreateOutputScalarData(t, (int)mu_ds_in, (int)sigma_ds_in, (int)mu_dr_in, (int)sigma_dr_in,    (int)(drug_dose_ATRi*100), (int)(drug_dose_PARPi*100), (int)no_circles, (int)(DR_fraction*100), (int)DR_type );

                    printf("Time: %d, t=%d, Cellcount: %d, AliveCount: %d.\n ", t, n_caption-100, no_cells, no_alive_cells );
                    n_caption++;
                    CreateOutputCellMaps(n_caption); // produce cell maps on paraview
                    CreateOutputCellMaps_resistant(n_caption); // produce drug-resistant-specific cell maps on paraview

                }
            }
        }
    }
    //end time loop
    Deallocate();
}
