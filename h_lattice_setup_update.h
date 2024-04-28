/*********************************************************
*
* Initiate the in silico setup.
*
*********************************************************/

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>


#ifndef H_LATTICE_SETUP
#define H_LATTICE_SETUP


void place_cells(int no_circles, int circle_rad, double DR_type, double DR_fraction, double frac_G1_in, double frac_S_in, double frac_SD_in, double frac_G2M_in);


void place_cells(int no_circles, int circle_rad, double DR_type, double DR_fraction, double frac_G1_in, double frac_S_in, double frac_SD_in, double frac_G2M_in) 
{

std::vector<int> placedCells;

for (int circleIndex = 0; circleIndex < no_circles; circleIndex++) 
{
    bool canPlaceCircle = false;

    while (!canPlaceCircle)  
    {
        int middlecell = rand() % NN;
        canPlaceCircle = true;

        for (int i = -circle_rad; i <= circle_rad; i++) 
        {
            for (int j = -circle_rad; j <= circle_rad; j++) 
            {
                int x = (middlecell % N + i + N) % N;  // Apply periodic boundary conditions in x-direction
                int y = (middlecell / N + j + N) % N;  // Apply periodic boundary conditions in y-direction

                int cell = x + N * y;

                if (i * i + j * j <= circle_rad * circle_rad) 
                {
                    if (cell_lives_here[cell] != 0) 
                    {
                        canPlaceCircle = false;
                        break;
                    }
                }
            }
        }

        if (canPlaceCircle) 
        {
            for (int i = -circle_rad; i <= circle_rad; i++) 
            {
                for (int j = -circle_rad; j <= circle_rad; j++) 
                {
                    int x = (middlecell % N + i + N) % N;  // Apply periodic boundary conditions in x-direction
                    int y = (middlecell / N + j + N) % N;  // Apply periodic boundary conditions in y-direction

                    int cell = x + N * y;

                    if (i * i + j * j <= circle_rad * circle_rad) 
                    {
                        PlaceNewCell_Placeholder(cell);
                        placedCells.push_back(cell);
                    }
                }
            }
        }
    }    
}

    
 
int counter_G1 = 0; 
int counter_S = 0;
int counter_SD = 0;
int counter_G2M = 0;
int counter_res = 0;
int counter = 0;


	while (counter_res < placedCells.size() * DR_fraction)
	{
		int randomIndex = rand() % placedCells.size();
		
		int currentValue = placedCells[randomIndex];
		
		if(cell_resistant_type[currentValue]==-1)
		{
		// Check the condition before placing a new cell in each phase
			cell_resistant_type[currentValue] = DR_type;
			counter_res++;
		}
	}
	
for (int i = 0; i < placedCells.size(); ++i) 
{
    int currentValue = placedCells[i];
    
    if (cell_resistant_type[currentValue] == -1) 
	{
        // Assign a value of 0 to cell_resistant_type[currentValue]
        cell_resistant_type[currentValue] = 0;
    }
}

    
	int rounded_G1 = round(frac_G1_in * placedCells.size());
    int rounded_S = round(frac_S_in * placedCells.size());
    int rounded_SD = round(frac_SD_in * placedCells.size());
    int rounded_G2M = placedCells.size() - rounded_G1 - rounded_S - rounded_SD;
	
        while(counter_G1 < rounded_G1) 
		{
			int randomIndex = rand() % placedCells.size();
			int currentValue = placedCells[randomIndex];
			if (cell_lives_here[currentValue] == 5) 
			{
            PlaceNewG1Cell(currentValue);
            counter_G1++;
			}
        } 
		while(counter_S < rounded_S) 
		{
			int randomIndex = rand() % placedCells.size();
			int currentValue = placedCells[randomIndex];
			if (cell_lives_here[currentValue] == 5) 
			{
            PlaceNewSCell(currentValue);
            counter_S++;
			}

        } 
		while (counter_SD < rounded_SD) 
		{
			int randomIndex = rand() % placedCells.size();
			int currentValue = placedCells[randomIndex];
			if (cell_lives_here[currentValue] == 5) 
			{
            PlaceNewSDCell(currentValue);
            counter_SD++;
			}

        } 
		while (counter_G2M < rounded_G2M) 
		{
			int randomIndex = rand() % placedCells.size();
			int currentValue = placedCells[randomIndex];
			if (cell_lives_here[currentValue] == 5) 
			{
            PlaceNewG2MCell(currentValue);
            counter_G2M++;
			}
        }
    



	printf("number of cells %ld \n", placedCells.size());

}


#endif
