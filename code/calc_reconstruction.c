#include <math.h>
#include "global.h"


void calc_reconstruction(void)
{


  int i, j;


  switch (mode_average)
  {

    case 1:     //geometric mean

      for (i = 0; i <= N; i++)
      {
        u_tilde[1][i] = sqrt( u[1][i]*u[1][i+1] );
        u_tilde[2][i] = sqrt( max(u[2][i]*u[2][i+1],ZERO) ) * sign(u[2][i+1]);
        u_tilde[3][i] = sqrt( u[3][i]*u[3][i+1] );
      }

    break;
    
  
    case 2:     //arithmetic mean

      for (j = 1; j <= 3; j++)
        for (i = 0; i <= N; i++)
          u_tilde[j][i] = HALF*( u[j][i]*u[j][i+1] );

    break;
  
  }




  switch (mode_flux)
  {

    case 1:     //use averages

      for (j = 1; j <= 3; j++)
        for (i = 0; i <= N; i++) {
          u_tilde_p[j][i] = u_tilde[j][i];
          u_tilde_m[j][i] = u_tilde[j][i];
        }
      
    break;
  

    case 2:     //use upwind scheme

      for (j = 1; j <= 3; j++)
        for (i = 1; i <= N; i++)
          u_dev[j][i] = minmod((u[j][i+1]-u[j][i])/d_radius_tilde[i],(u[j][i]-u[j][i-1])/d_radius_tilde[i-1]);          //u_dev = 0 at i=0 and i=N+1

      for (j = 1; j <= 3; j++)
        for (i = 0; i <= N; i++)
        {
          u_tilde_p[j][i] = u[j][i+1] + u_dev[j][i+1] *(radius_tilde[i]-radius[i+1]);
          u_tilde_m[j][i] = u[j][i]   + u_dev[j][i]   *(radius_tilde[i]-radius[i]);
        }

    break;
  
  }


}