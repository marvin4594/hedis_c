#include "global.h"


void solve_explicit(void)
{

  
  int i, j;
  

  for (j = 1; j <= 3; j++)
  {
    for (i = 0; i <= N+1; i++)
      u_old[j][i] = u[j][i];

    for (i = 1; i <= N; i++)
      D_old[j][i] = D[j][i];
  }


//-----Runge-Kutta 4th order----

  for (j = 1; j <= k; j++)
    for (i = 1; i <= N; i++)
    {
      k1[j][i] = - d_time*D[j][i];
      u[j][i]  = u_old[j][i] + HALF*k1[j][i];
    }

  calc_variables();

  for (j = 1; j <= k; j++)
    for (i = 1; i <= N; i++)
    {
      k2[j][i] = -d_time*D[j][i];
      u[j][i]   = u_old[j][i] + HALF*k2[j][i];
    }

  calc_variables();

  for (j = 1; j <= k; j++)
    for (i = 1; i <= N; i++)
    {
      k3[j][i] = -d_time*D[j][i];
      u[j][i]   = u_old[j][i] + k3[j][i];
    }

  calc_variables();

  for (j = 1; j <= k; j++)
    for (i = 1; i <= N; i++)
    {
      k4[j][i] = -d_time*D[j][i];
      u[j][i]   = u_old[j][i] + (ONE/SIX)*( k1[j][i] + TWO*k2[j][i] + TWO*k3[j][i] + k4[j][i] );
    }

  calc_variables();


}