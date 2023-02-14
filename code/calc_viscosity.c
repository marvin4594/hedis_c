#include "global.h"


void calc_viscosity(void)
{

  int i;


  switch (mode_viscosity) {

    case 1:     //beta-viscosity; bulk-viscosity=zero

        for(i = 0; i <= N+1; i++)
          eta[i] = beta*u[3][i]/radius[i];
      
        for(i = 0; i <= N; i++)
          eta_tilde[i] = beta*u_tilde[3][i]/radius_tilde[i];
      
    break;


    case 2:     //Pringle-viscosity; bulk-viscosity=zero

      for(i = 0; i <= N+1; i++)
        eta[i] = beta*u[1][i];
      
      for(i = 0; i <= N; i++)
        eta_tilde[i] = beta*u_tilde[1][i];

    break;
  
  }


}