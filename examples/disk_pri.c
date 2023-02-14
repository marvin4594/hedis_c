#include <math.h>
#include <stdio.h>
#include "../code/global.h"
#include "../code/nrutil.h"


//A galactic accretion disk with beta-viscosity using the Pringle-model


int main(int argc, char **argv)
{


//----------physical parameters-----------

  disk_inner_radius     =  1.E-3;
  disk_outer_radius     =  300.0;
  log_mass_central      =  4.0;
  log_mass_disk         =  10.0;
  beta                  =  0.01;
  cs                    =  9.28E-7;
  accr_eff              =  0.1;
  time_end              =  1.0E+9;

  SW_ag                 =  TRUE;                   //allow growth of central object
  SW_cem                =  TRUE;                   //consider escaped mass for gravitational acceleration

  mode_gravitation      =  3;                      //1: Kepler, 2: monopole (low order), 3: monopole (high order), 4: selfgravitation
  k                     =  1;                      //specify disk-model, 1: simple model (one equation), 3: extended model (three equations)
  mode_viscosity        =  1;                      //1: beta-viscosity, 2: Pringle-viscosity
  mode_initial_cond     =  1;

//----------numerical parameters----------

  N                     =  50;                     //number of grid points
  part                  =  0.8;                    //part of grid points between disk_inner_radius and disk_outer_radius
  c_cfl                 =  0.1;
  c_itw                 =  0.01;
  d_time_min            =  1.E-10;
  mu                    =  0.8;                    //Crank-Nicholson-Parameter

  mode_flux             =  1;                      //1: averages, 2: upwind-scheme
  mode_average          =  1;                      //1: geometric, 2: arithmetic
  mode_boundarys        =  4;                      //1: non-reflecting, 2: reflecting, 3: no-gradients, 4: simple model, 5: pringle-disk, 6: set by initial conditions
  mode_time_integr      =  2;                      //1: explicit, 2: implicit

  accuracy              =  1.E-4;
  epsil                 =  1.E-4;                  //for calculating H=dA/du in subroutine "henyey"
  iter_max              =  10;                     //maximum iterations for solving nonlinear system of equations
  num_visc              =  0.0;

//----------store options-----------------

  store_by_time         =  FALSE;                  //store after specified time or after specified number of timesteps
  time_store            =  1.E+6;                  //store after this time
  steps_store           =  10;                     //store after this number of timesteps
  time_store_general    =  1.E+9;                  //store general results after this time

  SW_store_distrib      =  TRUE;
  SW_store_timedev      =  TRUE;
  SW_store_general      =  FALSE;

  filename              =  "disk_pri";

//----------------------------------------------------------------------------------------------------



  char filename2[200];


  if (SW_store_general)
  {
    sprintf(filename2,"results/%s_general.dat",filename);
    file_general = fopen(filename2,"w");
  }


  disk();


  if (SW_store_general) fclose(file_general);

  return 0;

  
}






void initial_conditions(void)
{


  int i, j;
  double x, sum;
  double *sigma = dvector(1,N);


  
  switch (mode_initial_cond)
  {

    case 1:
    
      for(i = 1; i <= N_disk; i++)
        sigma[i] = sqrt(1 - pow((radius[i]-disk_inner_radius)/(disk_outer_radius-disk_inner_radius),2));
   
      for(i = 1, sum = 0; i <= N_disk; i++)
        sum += sigma[i]*ring_area[i];
     
      x = mass_disk*(1-0.001)/sum;

      for(i = 1; i <= N_disk; i++)
        sigma[i] *= x;


      for(i = N_disk+1; i <= N; i++)
        sigma[i] = 1/pow(radius[i],2);            //distribution for outer region
   
      for(i = N_disk+1, sum = 0; i <= N; i++)
        sum += sigma[i]*ring_area[i];
     
      x = mass_disk*0.001/sum;                    //put 1 promille of total mass into outer region

      for(i = N_disk+1; i <= N; i++)
        sigma[i] *= x;

      
      for(i = 1; i <= N; i++)
        u[1][i] = sigma[i]*radius[i];
     
      calc_gravacc();

      for(i = 1; i <= N; i++)
      {
        u[3][i] = u[1][i]*sqrt(c_grav*mass_central*radius[i] - a_g[i]*pow(radius[i],3));
        u[2][i] = 0;
      }

    break;
  
  
    case 2:

      for(i = 1; i <= N_disk; i++)
        sigma[i] = mass_disk*(1-0.001)/(c_pi*(pow(disk_outer_radius,2) - pow(disk_inner_radius,2)));     //constant distribution for the disk

      for(i = N_disk+1; i <= N; i++)
        sigma[i] = 1/pow(radius[i],2);                                                            //distribution for outer region


      for(i = N_disk+1, sum = 0; i <= N; i++)
        sum += sigma[i]*ring_area[i];
     
      x = mass_disk*0.001/sum;                                     //put 1 promille of total mass into outer region
     
      for(i = N_disk+1; i <= N; i++)
        sigma[i] *= x;

 
      for(i = 1; i <= N; i++)
        u[1][i] = sigma[i]*radius[i];
     
      calc_gravacc();

      for(i = 1; i <= N; i++)
      {
        u[3][i] = u[1][i]*sqrt(c_grav*mass_central*radius[i] - a_g[i]*pow(radius[i],3));
        u[2][i] = 0;
      }

    break;

  }


  for(j = 1; j <= 3; j++)
    for(i = 0; i <= N+1; i++)
      u_old[j][i] = u[j][i];

  free_dvector(sigma,1,N);


}
