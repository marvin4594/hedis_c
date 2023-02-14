#include <math.h>
#include <stdio.h>
#include "global.h"


void store_results(void)
{

  int i;

  

  if (!store_by_time)
  {
    index_store_results++;
    if (index_store_results == steps_store) store_flag = TRUE;
  }



  if (store_flag)
  {

    store_flag = FALSE;

    if (store_by_time) index_store_results++;
    if (!store_by_time) index_store_results = 0;


    if (SW_store_distrib) {

      for(i = 1; i <= N; i++)
        fprintf(file_distrib, "% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e\n", 
                radius[i], u[1][i], u[2][i], u[3][i], u[2][i]/(cs*u[1][i]), - c_grav*mass_bh/pow(radius[i],2), + a_g[i]);

      fprintf(file_distrib, "\n\n");
    }


    if (SW_store_timedev) fprintf(file_timedev, "%9i% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e\n", 
                                  timestep, time, d_time, d_time_cfl, mass_disk, mass_central, mass_bh, mass_lost, Mdot_disk, Mdot_edd, Mdot_bh);

      printf("%9i% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e\n",
             timestep, time, d_time, d_time_cfl, mass_disk, mass_central, mass_bh, ONE-(mass_lost+mass_disk)/mass_disk_init);
  }



  if ((time >= time_store_general) && (SW_store_general) && (!general_stored_yet))
  {
    general_stored_yet = TRUE;
    fprintf(file_general, "% 14.3e% 14.3e% 14.3e% 14.3e% 14.3e\n", disk_outer_radius, mass_disk_init, mass_bh, mass_central, mass_lost_edd);
  }


}
