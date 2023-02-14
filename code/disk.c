#include <stdio.h>
#include "global.h"


void disk(void)
{
  
  char filename2[200];

  alloc_arrays();
  init_parameters();
  init_grid();
  init_gravitation();
  initial_conditions();
  
  
  calc_variables();
  

  if (SW_store_timedev) {
    sprintf(filename2,"results/%s_timedev.dat",filename);
    file_timedev = fopen(filename2,"w");
  }
  
  if (SW_store_distrib) {
    sprintf(filename2,"results/%s_distrib.dat",filename);
    file_distrib = fopen(filename2,"w");
  }
  
  
  store_results();



  do
  {

    calc_timestep();

    timestep++;
    time += d_time;


    switch (mode_time_integr) {
      case 1:
        solve_explicit();
      break;
  
      case 2:
        solve_implicit();
      break;
    }
 

    calc_masses();
    store_results();

    
  } while ((!(time > time_end)) && !error);

  

  if (error) printf("ERROR, calculation aborted. Error-code: %i\n",error_code);              //error-code=0: timestep too small, else: solving sle failed


  free_arrays();


  if (SW_store_distrib) fclose(file_distrib);
  if (SW_store_timedev) fclose(file_timedev);

  
}
