/*********************************************************
 *
 * Handles multiple functionalities related to placing/removing cells on the lattice.
 *
 *********************************************************/

#ifndef H_CELL_PROLIFERATION_AND_DEATH
#define H_CELL_PROLIFERATION_AND_DEATH

/** functions*/
void SetInputOptions1(double mu_in, double sigma_in,   double mu_G1_in, double sigma_G1_in,   double mu_S_in, double sigma_S_in,   double mu_SD_in, double sigma_SD_in); //input parameters
int InitialCellCount(double mu_cells, double sigma_cells);
int CountCells(int t); //counts the total number of cells (alive and dead)
void CountCellsDetailed(int t); // counts the total number of cells in each state
int CountCyclingRepairableCells(int t); //count number of alive cells
int CountNCCells(int t); //count number of dead cells
void PlaceNewCell_Placeholder(int n); //places a new cell in no phase
void PlaceNewG1Cell(int n); //with added G1 time
void PlaceNewSCell(int n); // with added S time
void PlaceNewG2MCell(int n); // with added G2M time
void PlaceNewSDCell(int n); // with added SD time
void ResetMotherCell(int n);
void ScanForCellsToWakeUp();
void TryToWakeCell(int n);
void CellBecomesIrreparable(int n, int t);
void PutCellToSleep(int n_mother);
void CellDivision(int n_mother);
void ChooseNeighbourhood(int n_mother);
void MooreNeighbourhood(int n_mother);
int CheckMooreNbh1(int n_mother);
void VN_Neighbourhood(int n_mother);
int CheckVN_Nbh1(int n_mother);
int CheckLatticePoint(int n_check);
void CheckResistantPhase(int n);
int CircleRadius(int no_circles);
double mu_ds;
double sigma_ds;
double mu_dr;
double sigma_dr;
double mu_G1;
double sigma_G1;
double mu_S;
double sigma_S;
double mu_SD;
double sigma_SD;




int InitialCellCount(double mu_cells, double sigma_cells)
{
	int max_cells = 545;
	int min_cells = 419;
	
	double u = ((double)rand() / RAND_MAX);  // Random number between 0 and 1
    double uu = min_cells + u * (max_cells - min_cells);  // Transform to the desired range
    return uu;
} 


int CircleRadius(int no_circles)
{
		
int circle_rad;
if (no_circles == 1) 
{
    circle_rad = 11;
} 
else if (no_circles == 2) 
{
    circle_rad = 8;
} 
else if (no_circles == 4) 
{
    circle_rad = 5;
} 
else if (no_circles == 8) 
{
    circle_rad = 4;
} 
else if (no_circles == 16) 
{
    circle_rad = 3;
} 
else if (no_circles == 32) 
{
    circle_rad = 2;
} 
else if (no_circles == 64) 
{
    circle_rad = 1;
} 
else if (no_circles == 128) 
{
    circle_rad = 0;
} 
else if (no_circles == 256) 
{
    circle_rad = 0;
}
else
{
    circle_rad = 0;
}	
return circle_rad;
}


void SetInputOptions1(double mu_ds_in, double sigma_ds_in,  double mu_dr_in, double sigma_dr_in,   double mu_G1_in, double mu_S_in, double mu_SD_in)
{
    mu_ds=mu_ds_in;
    sigma_ds=sigma_ds_in;
	mu_dr=mu_dr_in;
    sigma_dr=sigma_dr_in;
    mu_G1 = mu_G1_in;
    mu_S = mu_S_in;
    mu_SD = mu_SD_in;
}


int CountCells(int t)
{
    /** Count number of cells on the lattice including dead cells*/
    int nocells=0;
    for(int n=0; n<NN; n++)
    {
        if(cell_cycle_phase[n]>0)
        {
            nocells++;
        }
    }
    return nocells;
}

void CountCellsDetailed(int t)
{
    /** Count (cell-cycle detailed) number of cells on the lattice */
    double no_g1_cells=0;
    double no_s_cells=0;
    double no_sd_cells=0;
    double no_g2m_cells=0;
    double no_cells=0;
    double no_all_cells=0;
    double no_nc_cells=0;

    for(int n=0; n<NN; n++)
    {
        if(cell_cycle_phase[n]>0)
        {
            no_all_cells++;

            if(cell_cycle_phase[n]==1)
            {
                no_g1_cells++;
                no_cells++;
            }
            if(cell_cycle_phase[n]==2)
            {
                no_s_cells++;
                no_cells++;
            }
            if(cell_cycle_phase[n]==3)
            {
                no_sd_cells++;
                no_cells++;
            }
            if(cell_cycle_phase[n]==4)
            {
                no_g2m_cells++;
                no_cells++;
            }
            if(cell_cycle_phase[n]==20)
            {
                no_nc_cells++;
				no_cells++;
            }
        }
    }
    printf("At time %d, NoAllCells=%lf, Nocells = %lf. G1=%lf. S=%lf. DS=%lf. G2M=%lf.\n ",t,no_all_cells, no_cells,no_g1_cells/no_cells*100,no_s_cells/no_cells*100,no_sd_cells/no_cells*100,no_g2m_cells/no_cells*
100);
}

int CountCyclingRepairableCells(int t)
{
    /** Count number of cells on the lattice */
    int nocyclingcells=0;
    for(int n=0; n<NN; n++)
    {
        if(cell_cycle_phase[n]>0 && cell_cycle_phase[n]!=20 && cell_cycle_phase[n]!=500)//(cell_lives_here[n]>0)
        {
            nocyclingcells++;
        }
    }
    return nocyclingcells;
}

int CountNCCells(int t)
{
    /** Count number of NC cells on the lattice */
    int nonccells=0;
    for(int n=0; n<NN; n++)
    {
        if(cell_cycle_phase[n]==20)
        {
            nonccells++;
        }
    }
    return nonccells;
}

void PlaceNewCell_Placeholder(int n)
{
    cell_cycle_phase[n] = 400;
    cell_lives_here[n]=5;
	cell_resistant_type[n]=-1;
}


void PlaceNewG1Cell(int n)
{
    {
        const double pi = 3.14159265358979323846;
        double u1, u2, uu,    u1_G1, u2_G1, uu_G1,    u1_S, u2_S, uu_S,    u1_SD, u2_SD, uu_SD;

        //Cell Cycle
        u1=( (double) rand() / (RAND_MAX) );
        u2=( (double) rand() / (RAND_MAX) );


        if(u1==0)
        {
            uu=0;
        }
        else
        {
			if(cell_resistant_type[n]==0)
			{
				uu = mu_ds + sigma_ds*sqrt(-2*log(u1))*cos(2*pi*u2);
			}
			else
			{
				uu = mu_dr + sigma_dr*sqrt(-2*log(u1))*cos(2*pi*u2);
			}
        }

        cell_cycle_length[n]=uu;
    }
	
	cell_cycle_length_G1[n] = mu_G1; 
	cell_cycle_length_S[n] = mu_S;
	cell_cycle_length_SD[n] = mu_SD;

    extra_G1_time[n] = (double)rand() / RAND_MAX * (cell_cycle_length_G1[n]);

    /** Set the initial variables for the newborn cell in lattice point n */
    cell_lives_here[n]=1;
    cell_cycle_phase[n]=1;
    cell_cycle_clock[n]= extra_G1_time[n]*cell_cycle_length[n];

    int valueofG1 = cell_cycle_length_G1[n]*cell_cycle_length[n];
}

void PlaceNewCell(int n) //this places a new cell in G1 on the grid with no extra time in G1 phase.
{
    {
        const double pi = 3.14159265358979323846;
        double u1, u2, uu,    u1_G1, u2_G1, uu_G1,    u1_S, u2_S, uu_S,    u1_SD, u2_SD, uu_SD;

        //Cell Cycle
        u1=( (double) rand() / (RAND_MAX) );
        u2=( (double) rand() / (RAND_MAX) );


        if(u1==0)
        {
            uu=0;
        }
        else
        {
			if(cell_resistant_type[n]==0)
			{
				uu = mu_ds + sigma_ds*sqrt(-2*log(u1))*cos(2*pi*u2);
			}
			else
			{
				uu = mu_dr + sigma_dr*sqrt(-2*log(u1))*cos(2*pi*u2);
			}
        }

        cell_cycle_length[n]=uu;
    }
	
	cell_cycle_length_G1[n] = mu_G1; 
	cell_cycle_length_S[n] = mu_S;
	cell_cycle_length_SD[n] = mu_SD;

    /** Set the initial variables for the newborn cell in lattice point n */
    cell_lives_here[n]=1;
    cell_cycle_phase[n]=1;
    cell_cycle_clock[n]= 1;
}

void PlaceNewSCell(int n)
{
    {
        const double pi = 3.14159265358979323846;
        double u1, u2, uu,    u1_G1, u2_G1, uu_G1,    u1_S, u2_S, uu_S,    u1_SD, u2_SD, uu_SD;

        //Cell Cycle
        u1=( (double) rand() / (RAND_MAX) );
        u2=( (double) rand() / (RAND_MAX) );


        if(u1==0)
        {
            uu=0;
        }
        else
        {
			if(cell_resistant_type[n]==0)
			{
				uu = mu_ds + sigma_ds*sqrt(-2*log(u1))*cos(2*pi*u2);
			}
			else
			{
				uu = mu_dr + sigma_dr*sqrt(-2*log(u1))*cos(2*pi*u2);
			}
        }

        cell_cycle_length[n]=uu;

    }
	cell_cycle_length_G1[n] = mu_G1; 
	cell_cycle_length_S[n] = mu_S;
	cell_cycle_length_SD[n] = mu_SD;
	
    extra_S_time[n] = (double)rand() / RAND_MAX * (cell_cycle_length_S[n]);

    /** Set the initial variables for the newborn cell in lattice point n */
    cell_lives_here[n]=1;
    cell_cycle_phase[n]=2;
    cell_cycle_clock[n]= (cell_cycle_length_G1[n]+extra_S_time[n])*cell_cycle_length[n];

    int valueofS = cell_cycle_length_G1[n]*cell_cycle_length[n] + cell_cycle_length_S[n]*cell_cycle_length[n];

}

void PlaceNewSDCell(int n)
{
    {
        const double pi = 3.14159265358979323846;
        double u1, u2, uu,    u1_G1, u2_G1, uu_G1,    u1_S, u2_S, uu_S,    u1_SD, u2_SD, uu_SD;

        //Cell Cycle
        u1=( (double) rand() / (RAND_MAX) );
        u2=( (double) rand() / (RAND_MAX) );


        if(u1==0)
        {
            uu=0;
        }
        else
        {
			if(cell_resistant_type[n]==0)
			{
				uu = mu_ds + sigma_ds*sqrt(-2*log(u1))*cos(2*pi*u2);
			}
			else
			{
				uu = mu_dr + sigma_dr*sqrt(-2*log(u1))*cos(2*pi*u2);
			}
        }

        cell_cycle_length[n]=uu;
    }
	
	cell_cycle_length_G1[n] = mu_G1; 
	cell_cycle_length_S[n] = mu_S;
	cell_cycle_length_SD[n] = mu_SD;

    extra_SD_time[n] = (double)rand() / RAND_MAX * (cell_cycle_length_SD[n]);


    /** Set the initial variables for the newborn cell in lattice point n */
    cell_lives_here[n]=1;
    cell_cycle_phase[n]=3;
    cell_cycle_clockSD[n]= (extra_SD_time[n])*cell_cycle_length[n];
	 cell_cycle_clock[n]= (cell_cycle_length_G1[n])*cell_cycle_length[n];



    int valueofSD = cell_cycle_length_SD[n]*cell_cycle_length[n];

}

void PlaceNewG2MCell(int n)
{
    {
        const double pi = 3.14159265358979323846;
        double u1, u2, uu,    u1_G1, u2_G1, uu_G1,    u1_S, u2_S, uu_S,    u1_SD, u2_SD, uu_SD;

        //Cell Cycle
        u1=( (double) rand() / (RAND_MAX) );
        u2=( (double) rand() / (RAND_MAX) );


        if(u1==0)
        {
            uu=0;
        }
        else
        {
			if(cell_resistant_type[n]==0)
			{
				uu = mu_ds + sigma_ds*sqrt(-2*log(u1))*cos(2*pi*u2);
			}
			else
			{
				uu = mu_dr + sigma_dr*sqrt(-2*log(u1))*cos(2*pi*u2);
			}
        }

        cell_cycle_length[n]=uu;
    }
    cell_cycle_length_G1[n] = mu_G1; 
	cell_cycle_length_S[n] = mu_S;
	cell_cycle_length_SD[n] = mu_SD;

	extra_G2_time[n] = (double)rand() / RAND_MAX * (1-cell_cycle_length_G1[n] - cell_cycle_length_S[n]);


    /** Set the initial variables for the newborn cell in lattice point n */
    cell_lives_here[n]=1;
    cell_cycle_phase[n]=4;
    cell_cycle_clock[n]= (cell_cycle_length_G1[n]+cell_cycle_length_S[n]+extra_G2_time[n])*cell_cycle_length[n];
	
}


void ResetMotherCell(int n)
{
    /** Reset mother cell after it has divided and produced a daughter cell */
    if(cell_cycle_phase[n]==4)
    {
        cell_cycle_phase[n]=1;
    }
    cell_cycle_clock[n]=1;
}

void ScanForCellsToWakeUp()
{
    /** Check if sleeping cell should wake up */
    for(int n=0; n<NN; n++)
    {
        if(cell_cycle_phase[n]==5)
        {
            TryToWakeCell(n);
        }
    }
}


void TryToWakeCell(int n)
{
    // Check (first order) neighbourhood to wake up sleeping cell
    bool wakeup = false;

    if(cell_lives_here[n-N-1]==0)
    {
        wakeup=true;
    }
    else if(cell_lives_here[n-N]==0)
    {
        wakeup=true;
    }
    else if(cell_lives_here[n-N+1]==0)
    {
        wakeup=true;
    }
    else if(cell_lives_here[n-1]==0)
    {
        wakeup=true;
    }
    else if(cell_lives_here[n+1]==0)
    {
        wakeup=true;
    }
    else if(cell_lives_here[n+N-1]==0)
    {
        wakeup=true;
    }
    else if(cell_lives_here[n+N]==0)
    {
        wakeup=true;
    }
    else if(cell_lives_here[n+N+1]==0)
    {
        wakeup=true;
    }

    if(wakeup==true)
    {
        ResetMotherCell(n);
    }

}

void CellBecomesIrreparable(int n, int t)
{
    /** Cell in lattice point n becomes irreparable. Note cell still remains on lattice */
    cell_cycle_clock[n]=1;
    cell_cycle_phase[n]=20;
    cell_lives_here[n]=1;
    scheduled_lattice_removal[n]=t;
}


void PutCellToSleep(int n_mother)
{
    if(cell_cycle_phase[n_mother] >0 && cell_cycle_phase[n_mother]<5)
    {
        cell_cycle_phase[n_mother]=5;
    }
}

void CellDivision(int n_mother)
{
    ChooseNeighbourhood(n_mother);
}

void ChooseNeighbourhood(int n_mother)
{
    int nbh_dice = rand() %2;
    if(nbh_dice==1)
    {
        MooreNeighbourhood(n_mother);
        //VN_Neighbourhood(n_mother);
    }
    else
    {
        //MooreNeighbourhood(n_mother);
         VN_Neighbourhood(n_mother);
    }
}

void MooreNeighbourhood(int n_mother)
{
    int n_daughter=-1;
    int mother_cell=0;
	

    // Check first order neighbourhood (NBH) //
    n_daughter = CheckMooreNbh1(n_mother);
	

    if(n_daughter==-1)
    {
        PutCellToSleep(n_mother);
    }


    if(n_daughter!=-1)
    {
        if(cell_cycle_phase[n_mother]==4)
        {
			cell_resistant_type[n_daughter] = cell_resistant_type[n_mother];
            PlaceNewCell(n_daughter);
        }
        ResetMotherCell(n_mother);
		
		
    }
}




int CheckMooreNbh1(int n_mother) 
{
    int daughter_dice = rand() % 9; // 8 directions
    int daughter_position = -1;
    int counter = 0;

    while (daughter_position == -1 && counter <= 8) 
	{
        int dx_new2 = (daughter_dice % 3) - 1;
        int dy_new2 = (daughter_dice / 3) - 1;

        int x = (n_mother % N + dx_new2 + N) % N;
        int y = (n_mother / N + dy_new2 + N) % N;

        int potential_daughter_position = x + N * y;

		daughter_position = CheckLatticePoint(potential_daughter_position);

        daughter_dice = (daughter_dice + 1) % 9;
        counter++;
    }

    return daughter_position;
}

int CheckLatticePoint(int check_n)
{
    if (cell_lives_here[check_n] == 0 )
    {
        
        return check_n;
    }
    else
    {
        return -1;
    }
}

void VN_Neighbourhood(int n_mother)
{
    int n_daughter=-1;
    int mother_cell=0;
	

    // Check first order neighbourhood (NBH) 
    n_daughter = CheckVN_Nbh1(n_mother);


    if(n_daughter==-1)
    {
        PutCellToSleep(n_mother);
    }
    if(n_daughter!=-1)
    {
         if(cell_cycle_phase[n_mother]==4)
        {
			cell_resistant_type[n_daughter] = cell_resistant_type[n_mother];
            PlaceNewCell(n_daughter);
        }
        ResetMotherCell(n_mother);
		
    }
}

int CheckVN_Nbh1(int n_mother) 
{
    int daughter_dice = rand() % 4;
    int daughter_position = -1;
    int counter = 0;

    while (daughter_position == -1 && counter <= 3) 
	{
        int dx_new2, dy_new2;
		
		if (daughter_dice == 0) 
		{
            dx_new2 = 0;
            dy_new2 = -1;
        } 
		else if (daughter_dice == 1) 
		{
            dx_new2 = 0;
            dy_new2 = 1;
        } 
		else if (daughter_dice == 2) 
		{
            dx_new2 = -1;
            dy_new2 = 0;
        } 
		else 
		{
            dx_new2 = 1;
            dy_new2 = 0;
        }

        int x = (n_mother % N + dx_new2 + N) % N;
        int y = (n_mother / N + dy_new2 + N) % N;

        int potential_daughter_position = x + N * y;

        daughter_position = CheckLatticePoint(potential_daughter_position);

        daughter_dice = (daughter_dice + 1) % 4;
        counter++;
    }

    return daughter_position;
}

void CheckResistantPhase(int n)
{
    for(int n=0; n<NN; n++)
    {
    if (cell_cycle_phase[n] == 1) 
	{
        cell_cycle_resistant_phase[n] = cell_resistant_type[n] * 20 + 1;
    } 
	else if (cell_cycle_phase[n] == 2) 
	{
        cell_cycle_resistant_phase[n] = cell_resistant_type[n] * 20 + 2;
    } 
	else if (cell_cycle_phase[n] == 3) 
	{
        cell_cycle_resistant_phase[n] = cell_resistant_type[n] * 20 + 3;
    } 
	else if (cell_cycle_phase[n] == 4) 
	{
        cell_cycle_resistant_phase[n] = cell_resistant_type[n] * 20 + 4;
    } 
	else if (cell_cycle_phase[n] == 5) 
	{
        cell_cycle_resistant_phase[n] = cell_resistant_type[n] * 20 + 5;
    } 
	else if (cell_cycle_phase[n] == 20) 
	{
        cell_cycle_resistant_phase[n] = cell_resistant_type[n] * 20 + 20;
    }
    }
}



#endif
