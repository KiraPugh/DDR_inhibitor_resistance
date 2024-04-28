/*********************************************************
 *
 * Handles cell progression through the cell cycle.
 *
 *********************************************************/
 
#ifndef H_CELL_CYCLE
#define H_CELL_CYCLE

void Compute_CellCycle(int t);
void Compute_CellCycle_CP(int t);
void EvaluateDrugs();
void SetInputOptions2(double hc_ATRi_in, double EC50_ATRi_in, double hc_PARPi_in, double EC50_PARPi_in, double prob_healthy_in, double prob_repair_in);


int dice_G1_fate, dice_SD_fate, protect_dice;
int attack_drugeffect_1, attack_drugeffect_2;

double Emax_ATRi;
double Emax_PARPi;
double hc_ATRi; 
double EC50_ATRi;
double hc_PARPi;
double EC50_PARPi;
double prob_healthy;
double prob_repair;


int GetDrugAttack_1(int n); // PARPi inhibits cells transitioning from G1 to S
int GetDrugAttack_2(int n); // ATRi and PARPi inhibits cells transitioning from SD to S 


/** Set input options to global variables */
void SetInputOptions2(double Emax_ATRi_in, double Emax_PARPi_in, double hc_ATRi_in, double EC50_ATRi_in, double hc_PARPi_in, double EC50_PARPi_in, double prob_healthy_in, double prob_repair_in)
{
	Emax_ATRi = Emax_ATRi_in;
	Emax_PARPi = Emax_PARPi_in;
    hc_ATRi = hc_ATRi_in; 
    EC50_ATRi = EC50_ATRi_in;
    hc_PARPi = hc_PARPi_in; 
    EC50_PARPi = EC50_PARPi_in;
    prob_healthy = prob_healthy_in;
    prob_repair = prob_repair_in;
}

/** Progress cell-cycle by time steps */
void Compute_CellCycle(int t)
{
    int length_cc;
    double length_G1;
    double length_S;
    double length_SD;
	
	int nocells = CountCells(t);

    for(int n=0; n<NN; n++)
    {
		
		
        length_cc=cell_cycle_length[n];		
        length_G1=cell_cycle_length_G1[n];
        length_S=cell_cycle_length_S[n];
        length_SD=cell_cycle_length_SD[n];

        /** Handle NC cells */
        if(cell_cycle_phase[n]==20)
        {
            cell_cycle_clock[n]++;
        }

        /** Handle cycling cell */
        if(cell_cycle_phase[n]>0  && cell_cycle_phase[n]<5 && cell_cycle_phase[n]!=3)
        {
            /** increment cell cyle clock +1 */
            cell_cycle_clock[n]++;

            if(cell_cycle_phase[n]==1) //If cell is in G1
            {
                if(cell_cycle_clock[n]>length_G1*length_cc) //If cell has reached end of G1
                {
                    /** Roll a dice to see if cell should enter S or SD state */
                    dice_G1_fate = rand() % 101;
					attack_drugeffect_1 = prob_healthy*(100-GetDrugAttack_1(n));
					
                    if(dice_G1_fate <= attack_drugeffect_1)
                    {
                        cell_cycle_phase[n]=2; //Go to S
                    }
                    else
                    {
                        cell_cycle_phase[n]=3;//Go to SD
                    }
                }
            }

            else if(cell_cycle_phase[n]==2) //If cell is in S
            {
                if(cell_cycle_clock[n]>(length_G1+length_S)*length_cc) //If cell has reached end of S
                {
                    cell_cycle_phase[n]=4; //Go to G2/M
                }
            }

            else if(cell_cycle_phase[n]==4) //If cell in G2/M
            {
                if(cell_cycle_clock[n] >= length_cc)
                {
                    {
                        CellDivision(n);
                    }

                }
            }
        }
		
		
        if(cell_cycle_phase[n]==3) //If cell is in SD
        {
            cell_cycle_clockSD[n]++;

            //If cell has reached end of SD state: Try to go to S. Otherwise go to SD.
            if( cell_cycle_clockSD[n] > length_SD*length_cc )
            {
                protect_dice = rand() % 101;
                attack_drugeffect_2 = prob_repair*(100-GetDrugAttack_2(n));

                if(protect_dice<=attack_drugeffect_2)
                {
                    /** Cell repairs */
                    cell_cycle_phase[n]=2; //Go to S
                }
                else
                {
                    /** Cell becomes irreparable */
                    CellBecomesIrreparable(n,t);
                }
            }
        }
    }
}






int GetDrugAttack_1(int n)
{
    double attack_E=0;
	
	if(cell_resistant_type[n] == 0 || cell_resistant_type[n] == 1 )
	{
        attack_E=100*Emax_PARPi * pow(drug_PARPi[n],hc_PARPi)/( pow(EC50_PARPi,hc_PARPi) + pow(drug_PARPi[n],hc_PARPi) );
	}
	else
	{
		attack_E = 0;
	}

    return (int)attack_E;

}
int GetDrugAttack_2(int n)

{
    double attack_PARPi=0;
    double attack_ATRi=0;
	double attack_E = 0;
	
	attack_PARPi=Emax_PARPi * pow(drug_PARPi[n],hc_PARPi)/( pow(EC50_PARPi,hc_PARPi) + pow(drug_PARPi[n],hc_PARPi) ) ;
    attack_ATRi=Emax_ATRi * pow(drug_ATRi[n],hc_ATRi)/( pow(EC50_ATRi,hc_ATRi) + pow(drug_ATRi[n],hc_ATRi) ) ;
	
	if(cell_resistant_type[n] == 0)
	{
        attack_E = 100 * ( attack_PARPi + attack_ATRi - attack_PARPi * attack_ATRi );
	}
	else if(cell_resistant_type[n] ==1 )
	{
		attack_E = 100 * ( attack_PARPi );
	}
	else if(cell_resistant_type[n] == 2)
	{
		attack_E =100 * ( attack_ATRi );
	}
	else if(cell_resistant_type[n] == 3)
	{
		attack_E = 0;
	}		
		//printf("attack_2 %f \n", attack_E);


    return (int)attack_E;

}



void EvaluateDrugs_c()
{
    double maxdrug =0;
    double mindrug =10000;
    double totaldrug =0;
    double nocells =0;

    for(int n=0; n<NN; n++)
    {
        if(cell_cycle_phase[n]>0 && cell_cycle_phase[n]<6)
        {
            nocells++;
            totaldrug=totaldrug+drug_ATRi[n];
            if(drug_ATRi[n]>maxdrug)
            {
                maxdrug=drug_ATRi[n];
            }
            if(drug_ATRi[n]<mindrug)
            {
                mindrug=drug_ATRi[n];
            }

        }

    }
    printf("MAX: %lf, MIN: %lf, AVG: %lf, TOT:%lf, nocells:%lf \n", maxdrug, mindrug, totaldrug/nocells, totaldrug, nocells);
}

void EvaluateDrugs_o()
{
    double maxdrug =0;
    double mindrug =10000;
    double totaldrug =0;
    double nocells =0;

    for(int n=0; n<NN; n++)
    {
        if(cell_cycle_phase[n]>0 && cell_cycle_phase[n]<6)
        {
            nocells++;
            totaldrug=totaldrug+drug_PARPi[n];
            if(drug_PARPi[n]>maxdrug)
            {
                maxdrug=drug_PARPi[n];
            }
            if(drug_PARPi[n]<mindrug)
            {
                mindrug=drug_PARPi[n];
            }

        }

    }
    printf("MAX: %lf, MIN: %lf, AVG: %lf, TOT:%lf, nocells:%lf \n", maxdrug, mindrug, totaldrug/nocells, totaldrug, nocells);
}
#endif
