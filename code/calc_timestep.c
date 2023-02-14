#include <math.h> 
#include "global.h"


void calc_timestep(void)
{


  double eps, t_1, t_2, t_3;
  int i;


  t_1 = d_radius[0]/( fabs(u[2][0]/u[1][0]) + cs );
  t_2 = sqrt(fabs(d_radius[0]/( - c_grav*mass_central/pow(radius[0],2) + a_g[0] )));
  t_3 = fabs( pow(d_radius[0],2)*u[1][0]/(beta*u[3][0]) );

  eps = (ONE - u_old[1][0]/u[1][0]);


  for(i = 1; i <= N+1; i++){ 
    t_1 = min(t_1, d_radius[i]/( fabs(u[2][i]/u[1][i]) + cs ));                                    //cfl-condition for advection
    t_2 = min(t_2, sqrt(fabs(d_radius[i]/( - c_grav*mass_central/pow(radius[i],2) + a_g[i] ))));   //cfl-condition for grav. acceleration
    t_3 = min(t_3, fabs( pow(d_radius[i],2)*u[1][i]/(beta*u[3][i]) ));                             //cfl-condition for beta-viscosity

    eps = max(eps, fabs(ONE - u_old[1][i]/u[1][i]));
  }



  d_time_cfl = c_cfl*min(t_1,min(t_2,t_3));

  d_time_itw = max(d_time_cfl, (c_itw/eps)*d_time);



  switch (mode_time_integr)
  {
    case 1:
      d_time = d_time_cfl;
    break;
  
    case 2:
      d_time = d_time_itw;
    break;
  }



  if (d_time < d_time_min)
    error = TRUE;


  if ((store_by_time) && (( time+d_time ) > index_store_results*time_store))                       //passed time to store?
  {
    d_time = index_store_results*time_store - time;                                                //cut timestep
    store_flag = TRUE;                                                                             //store results
  }


  if (SW_store_general && (!general_stored_yet) && ((time+d_time) > time_store_general))           //cut timestep for storing general results
     d_time = time_store_general - time;


}