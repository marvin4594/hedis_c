#include <math.h>
#include "global.h"


void solve_implicit(void)
{


  double maxval;
  int i, j;


  for(j = 1; j <= 3; j++)
    for(i = 0; i <= N+1; i++)
    {
      d_u[j][i] = u[j][i] - u_old[j][i];
      u_old[j][i] = u[j][i];
      u[j][i] = u_old[j][i] + d_u[j][i];
      d_u[j][i] = ZERO;
    }


  for(j = 1; j <= k; j++)
    for(i = 1; i <= N; i++)
      D_old[j][i] = D[j][i];


  iter = 0;


  
  do
  {

    iter++;

    calc_variables();

    henyey();

    solve_sle();

    for(j = 1; j <= 3; j++)
      for(i = 0; i <= N+1; i++)
        u[j][i] = u[j][i] + d_u[j][i];

    for(j = 1, maxval = 0; j <= k; j++)
      for(i = 0; i <= N+1; i++)
        maxval = max(maxval,fabs(d_u[j][i]/u[j][i]));
  
    
  } while (!(maxval < accuracy) && !(iter = iter_max));



  calc_variables();


}