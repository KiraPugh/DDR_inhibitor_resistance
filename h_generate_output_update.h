/*********************************************************
 *
 *  Create output. Scalar data (text/csv files) and cell/field maps (vtk).
 *
 *********************************************************/

#ifndef H_GENERATE_OUTPUT
#define H_GENERATE_OUTPUT
using namespace std;

/** Functions */
void CreateOutputCellMaps(int n_caption); //vtk file for cells
void CreateOutputCellMaps_resistant(int n_caption); //vtk file for cells
void CreateOutputScalarFieldMaps(int n_caption); //vtk file for fields (oxygen/drug)
void CreateOutputScalarFieldMaps2(int n_caption); //vtk file for fields (oxygen/drug)
void CreateOutputScalarData(int t, int mu_ds_in, int sigma_ds_in, int mu_dr_in, int sigma_dr_in, int drug_dose_ATRi, int drug_dose_PARPi, int no_circles, int DR_fraction, int DR_type );
void CreateOutputScalarData_LHC(int t, int prob_Damaged_S_in, int T_in_Sd_in, int hc_in, int EC50_in, int T_Death_Delay_in, int mu_ds_in, int sigma_ds_in,  int mu_dr_in, int sigma_dr_in, int drug_dose_ATRi, int drug_dose_PARPi, int LHC_sample );



// CREATE THE CELL MAPS ON PARAVIEW
void CreateOutputCellMaps(int n_caption)
{
  char fullname[256]="invitro_cellmap_3_"; //Choose file name
  sprintf(fullname, "invitro_cellmap_3_%d.vtk",n_caption); //Choose file name
  FILE *tempVTKfile;
  tempVTKfile = fopen(fullname,"w");

  /**Header, required format*/
  fprintf(tempVTKfile, "# vtk DataFile Version 2.0 \nTumourTREAT \nASCII \nDATASET POLYDATA \n");

  /**count how many cells exist */
  int nocells=0;
  for(int n = 0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0)
      {
          nocells++;
      }
    }

  /**Header, required format*/
  fprintf(tempVTKfile, "POINTS %d float \n", nocells);

  /**Write coordinate data*/
  int x_pos = 0;
  int y_pos = 0;
  for(int n=0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0)
	{
            x_pos=n/N;
            y_pos=n%N;
            fprintf(tempVTKfile, "%d %d %d \n",x_pos,y_pos,0);
	}
    }

  /** Header, required format*/
  fprintf(tempVTKfile, "POINT_DATA %d \nSCALARS sample_scalars float 1 \nLOOKUP_TABLE my_table \n", nocells);

  double temp=0;
  /** Write cell state data */
  for(int n=0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0)
	{
	  temp = cell_cycle_phase[n];
	  fprintf(tempVTKfile, "%lf \n", temp);
	}
    }
  fclose(tempVTKfile);
}

// CREATE THE CELL MAPS ON PARAVIEW
void CreateOutputCellMaps_resistant(int n_caption)
{
  char fullname[256]="invitro_cellmap_res_3_"; //Choose file name
  sprintf(fullname, "invitro_cellmap_res_3_%d.vtk",n_caption); //Choose file name
  FILE *tempVTKfile;
  tempVTKfile = fopen(fullname,"w");

  /**Header, required format*/
  fprintf(tempVTKfile, "# vtk DataFile Version 2.0 \nTumourTREAT \nASCII \nDATASET POLYDATA \n");

  /**count how many cells exist */
  int nocells=0;
  for(int n = 0; n<NN; n++)
    {
      if(cell_cycle_resistant_phase[n]>0)
      {
          nocells++;
      }
    }

  /**Header, required format*/
  fprintf(tempVTKfile, "POINTS %d float \n", nocells);

  /**Write coordinate data*/
  int x_pos = 0;
  int y_pos = 0;
  for(int n=0; n<NN; n++)
    {
      if(cell_cycle_resistant_phase[n]>0)
	{
            x_pos=n/N;
            y_pos=n%N;
            fprintf(tempVTKfile, "%d %d %d \n",x_pos,y_pos,0);
	}
    }

  /** Header, required format*/
  fprintf(tempVTKfile, "POINT_DATA %d \nSCALARS sample_scalars float 1 \nLOOKUP_TABLE my_table \n", nocells);

  double temp=0;
  /** Write cell state data */
  for(int n=0; n<NN; n++)
    {
      if(cell_cycle_resistant_phase[n]>0)
	{
	  temp = cell_cycle_resistant_phase[n];
	  fprintf(tempVTKfile, "%lf \n", temp);
	}
    }
  fclose(tempVTKfile);
}


void CreateOutputScalarFieldMaps(int n_caption)
{
    char fullname[256]="invitro3_oxy_";//Choose file name
    sprintf(fullname, "invitro3_oxy_%d.vtk",n_caption);//Choose file name
    FILE *tempVTKfile;
    tempVTKfile = fopen(fullname,"w");

    /**Header, required format*/
    fprintf(tempVTKfile, "# vtk DataFile Version 2.0 \nTumourTREAT \nASCII \nDATASET POLYDATA \n");

    /**count how many cells exist */
    int nolatticepoints=0;

    //   int nocells=0;
    for(int n = 0; n<NN; n++)
      {
	if(cell_cycle_phase[n]>0)
	  {
	    nolatticepoints++;
	    //	    nocells++;
	  }
      }

    /**Header, required format*/
    fprintf(tempVTKfile, "POINTS %d float \n", nolatticepoints);

    /**Write coordinate data*/
    int x_pos = 0;
    int y_pos = 0;
    for(int n=0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0){
            x_pos=n/N;
            y_pos=n%N;
            fprintf(tempVTKfile, "%d %d %d \n",x_pos,y_pos,0);
      }
    }

    /** Header, required format*/
    fprintf(tempVTKfile, "POINT_DATA %d \nSCALARS sample_scalars float 1 \nLOOKUP_TABLE my_table \n", nolatticepoints);

    double temp=0;
    /** Write cell state data */
    for(int n=0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0){
            temp = oxygen_scaled[n];
            fprintf(tempVTKfile, "%lf \n", temp);
      }
    }
    fclose(tempVTKfile);
}

void CreateOutputScalarFieldMaps2(int n_caption)
{
  char fullname[256]="invitro3_drug_ATRi_";//Choose file name
  sprintf(fullname, "invitro3_drug_ATRi_%d.vtk",n_caption);//Choose file name
  FILE *tempVTKfile;
  tempVTKfile = fopen(fullname,"w");

  /**Header, required format*/
  fprintf(tempVTKfile, "# vtk DataFile Version 2.0 \nTumourTREAT \nASCII \nDATASET POLYDATA \n");

  /**count how many cells exist */
  // int nolatticepoints=NN;


  int nolatticepoints=0;

  //   int nocells=0;
  for(int n = 0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0)
	{
	  nolatticepoints++;
	  //      nocells++;
	}
    }



  /**Header, required format*/
  fprintf(tempVTKfile, "POINTS %d float \n", nolatticepoints);

  /**Write coordinate data*/
  int x_pos = 0;
  int y_pos = 0;
  for(int n=0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0)
        {
      x_pos=n/N;
      y_pos=n%N;
      fprintf(tempVTKfile, "%d %d %d \n",x_pos,y_pos,0);
	}
    }

  /** Header, required format*/
  fprintf(tempVTKfile, "POINT_DATA %d \nSCALARS sample_scalars float 1 \nLOOKUP_TABLE my_table \n", nolatticepoints);

  double temp=0;
  /** Write cell state data */
  for(int n=0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0)
        {
      temp = drug_ATRi[n];
      fprintf(tempVTKfile, "%lf \n", temp);
	}
    }
  fclose(tempVTKfile);
}

void CreateOutputScalarFieldMaps3(int n_caption)
{
  char fullname[256]="invitro3_drug_PARPi_";//Choose file name
  sprintf(fullname, "invitro3_drug_PARPi_%d.vtk",n_caption);//Choose file name
  FILE *tempVTKfile;
  tempVTKfile = fopen(fullname,"w");

  /**Header, required format*/
  fprintf(tempVTKfile, "# vtk DataFile Version 2.0 \nTumourTREAT \nASCII \nDATASET POLYDATA \n");

  /**count how many cells exist */
  // int nolatticepoints=NN;


  int nolatticepoints=0;

  //   int nocells=0;
  for(int n = 0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0)
	{
	  nolatticepoints++;
	  //      nocells++;
	}
    }



  /**Header, required format*/
  fprintf(tempVTKfile, "POINTS %d float \n", nolatticepoints);

  /**Write coordinate data*/
  int x_pos = 0;
  int y_pos = 0;
  for(int n=0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0)
        {
      x_pos=n/N;
      y_pos=n%N;
      fprintf(tempVTKfile, "%d %d %d \n",x_pos,y_pos,0);
	}
    }

  /** Header, required format*/
  fprintf(tempVTKfile, "POINT_DATA %d \nSCALARS sample_scalars float 1 \nLOOKUP_TABLE my_table \n", nolatticepoints);

  double temp=0;
  /** Write cell state data */
  for(int n=0; n<NN; n++)
    {
      if(cell_cycle_phase[n]>0)
        {
      temp = drug_PARPi[n];
      fprintf(tempVTKfile, "%lf \n", temp);
	}
    }
  fclose(tempVTKfile);
}



void CreateOutputScalarData(int t, int mu_ds_in, int sigma_ds_in, int mu_dr_in, int sigma_dr_in, int drug_dose_ATRi, int drug_dose_PARPi, int no_circles, int DR_fraction, int DR_type)
{
    FILE *scalardata_file;

        string filename ="oct_B" ;//


    filename+="_mu_ds"+std::to_string(mu_ds_in);
    filename+="_sigma_ds"+std::to_string(sigma_ds_in);
	filename+="_mu_dr"+std::to_string(mu_dr_in);
    filename+="_sigma_dr"+std::to_string(sigma_dr_in);
    filename+="_DoseC"+std::to_string(drug_dose_ATRi);
    filename+="_DoseO"+std::to_string(drug_dose_PARPi);
	filename+="_NoCircles"+std::to_string(no_circles);
	filename+="_DRfrac"+std::to_string(DR_fraction);
	filename+="_DRtype"+std::to_string(DR_type);

    scalardata_file=fopen(filename.c_str(),"a");

    double g1cells=0;
    double scells=0;
    double sdcells=0;
    double g2mcells=0;
    double nccells=0;
    double g0cells=0;
    double allcyclingcells=0;
    double allhealthycells = 0;
    double allhealthycells_wg0 = 0;
    double allcells=0;
    double op1=0;
    double op2=0;
    double nocells_PARPi_res=0;
    double nocells_ATRi_res=0;
    double nocells_comb_res=0;
    double sensitivecells=0;
    double resistantcells=0;
	
	double nocells_PARPi_res_all=0;
    double nocells_ATRi_res_all=0;
    double nocells_comb_res_all=0;
    double sensitivecells_all=0;
    double resistantcells_incdead=0;




    for(int n=0; n<NN; n++)
        {
            if(cell_lives_here[n]>0)
            {
                if(cell_cycle_phase[n]==1)
                {
                    g1cells++;
                }
                else if(cell_cycle_phase[n]==2)
                {
                    scells++;
                }
                else if(cell_cycle_phase[n]==3)
                {
                    sdcells++;
                }
                else if(cell_cycle_phase[n]==4)
                {
                    g2mcells++;
                }
                else if(cell_cycle_phase[n]==20)
                {
                    nccells++;
                }
                else if(cell_cycle_phase[n]==5)
                {
                    g0cells++;
                }
            }
        }


    for(int n=0; n<NN; n++)
        {
            if(cell_lives_here[n]>0)
            {

              if(cell_resistant_type[n] == 0 && cell_cycle_phase[n] != 20)
                {
                	sensitivecells++;
                }
                else if(cell_resistant_type[n] == 1 && cell_cycle_phase[n] != 20)
                {
                	nocells_ATRi_res++;
                }
                else if(cell_resistant_type[n] == 2 && cell_cycle_phase[n] != 20)
                {
                	nocells_PARPi_res++;
                }
                else if(cell_resistant_type[n] == 3 && cell_cycle_phase[n] != 20)
                {
                	nocells_comb_res++;
                }
           }
        }

    for(int n=0; n<NN; n++)
        {
            if(cell_lives_here[n]>0)
            {

              if(cell_resistant_type[n] == 0)
                {
                	sensitivecells_all++;
                }
                else if(cell_resistant_type[n] == 1)
                {
                	nocells_ATRi_res_all++;
                }
                else if(cell_resistant_type[n] == 2)
                {
                	nocells_PARPi_res_all++;
                }
                else if(cell_resistant_type[n] == 3)
                {
                	nocells_comb_res_all++;
                }
           }
        }

        allcyclingcells=g1cells+scells+sdcells+g2mcells+nccells; //cycling cells with NC cells
        allhealthycells=g1cells+scells+sdcells+g2mcells; // just cycling cells
        allhealthycells_wg0=g1cells+scells+sdcells+g2mcells+g0cells; // cycling cells with G0

        allcells=g1cells+scells+sdcells+g2mcells+nccells+g0cells; //all cells on the lattice (cycling, G0, and NC)

        op1= 100*(nccells+sdcells)/allcyclingcells;
        op2= 100*(nccells+sdcells)/allcells;

        double g1cells_p = g1cells/(g1cells+g2mcells+scells+sdcells+g0cells)*100;
        double g2mcells_p = g2mcells/(g1cells+g2mcells+scells+sdcells+g0cells)*100;
        double scells_p = (scells+sdcells)/(g1cells+g2mcells+scells+sdcells+g0cells)*100;
        double sdcells_p = sdcells/(g1cells+g2mcells+scells+sdcells+g0cells)*100;

        resistantcells = nocells_PARPi_res + nocells_ATRi_res + nocells_comb_res;
		resistantcells_incdead = nocells_PARPi_res_all + nocells_ATRi_res_all + nocells_comb_res_all;
        double resistant_p = resistantcells/allcells*100;
        double sensitive_p = sensitivecells/allcells*100;
        double nocells_ATRi_res_p = nocells_ATRi_res/allcells*100;
        double nocells_PARPi_res_p = nocells_PARPi_res/allcells*100;
        double nocells_comb_res_p = nocells_comb_res/allcells*100;


        fprintf(scalardata_file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",t, op1, op2, allhealthycells, allhealthycells_wg0, allcells, g1cells_p, scells_p, sdcells_p, g2mcells_p, g0cells, nccells, sensitivecells, sensitive_p, resistantcells, resistant_p, nocells_ATRi_res, nocells_PARPi_res, nocells_comb_res, nocells_ATRi_res_all, nocells_PARPi_res_all, nocells_comb_res_all, sensitivecells_all,  resistantcells_incdead );

    fclose(scalardata_file);
}

void CreateOutputScalarData_LHC(int t, int prob_Damaged_S_in, int T_in_Sd_in, int hc_in, int EC50_in, int T_Death_Delay_in, int mu_ds_in, int sigma_ds_in, int mu_dr_in, int sigma_dr_in,  int drug_dose_ATRi, int drug_dose_PARPi, int LHC_sample )
{
    FILE *scalardata_file;

       string filename ="/Users/sara/Desktop/DDRinhibitor/SensitivityAnalysisDataLHC/lhc_";


 
    filename+="_sample_"+std::to_string(LHC_sample);

    scalardata_file=fopen(filename.c_str(),"a");


    double g1cells=0;
    double scells=0;
    double sdcells=0;
    double g2mcells=0;
    double nccells=0;
    double g0cells=0;
    double allcyclingcells=0;
    double allhealthycells = 0;
    double allhealthycells_wg0 = 0;
    double allcells=0;
    double op1=0;
    double op2=0;

    for(int n=0; n<NN; n++)
        {
            if(cell_lives_here[n]>0)
            {
                if(cell_cycle_phase[n]==1 || cell_cycle_phase[n]==21 || cell_cycle_phase[n]==41 || cell_cycle_phase[n]==61)
                {
                    g1cells++;
                }
                else if(cell_cycle_phase[n]==2 || cell_cycle_phase[n]==22 || cell_cycle_phase[n]==42 || cell_cycle_phase[n]==62)
                {
                    scells++;
                }
                else if(cell_cycle_phase[n]==3 || cell_cycle_phase[n]==23 || cell_cycle_phase[n]==43 || cell_cycle_phase[n]==63)
                {
                    sdcells++;
                }
                else if(cell_cycle_phase[n]==4 || cell_cycle_phase[n]==24 || cell_cycle_phase[n]==44 || cell_cycle_phase[n]==64)
                {
                    g2mcells++;
                }
                else if(cell_cycle_phase[n]==20 || cell_cycle_phase[n]==40 || cell_cycle_phase[n]==60 || cell_cycle_phase[n]==80)
                {
                    nccells++;
                }
                else if(cell_cycle_phase[n]==5 || cell_cycle_phase[n]==25 || cell_cycle_phase[n]==45 || cell_cycle_phase[n]==65)
                {
                    g0cells++;
                }
            }
        }
        allcyclingcells=g1cells+scells+sdcells+g2mcells+nccells;
        allhealthycells=g1cells+scells+sdcells+g2mcells;
        allhealthycells_wg0=g1cells+scells+sdcells+g2mcells+g0cells;

        allcells=g1cells+scells+sdcells+g2mcells+nccells+g0cells;

        op1= 100*(nccells+sdcells)/allcyclingcells;
        op2= 100*(nccells+sdcells)/allcells;

        fprintf(scalardata_file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",t, op1, op2, allhealthycells, allhealthycells_wg0, allcells, g1cells, scells, sdcells, g2mcells, g0cells, nccells);

    fclose(scalardata_file);
}

#endif
