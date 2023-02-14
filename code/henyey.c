#include <gsl/gsl_linalg.h>
#include <math.h>
#include "global.h"


void henyey(void)
{
//computes the Henyey-Matrix H


  int i, j, l, m, r, q, p;
  double delta;


  for(i = 1; i <= k*N; i++)
    for(j = 1; j <= k; j++)
      for(p = 1; p <= N; p++)
        gsl_matrix_set(H, (j-1)*N+p-1, i-1, A[j][p]);

    

  for(q = 1; q <= k*N; q++)
  {
    j = (q-1)/N + 1;
    i = q - (j-1)*N;
 
    delta = epsil*fabs(u[j][i]);
    if (u[j][i] == ZERO) delta = epsil;
    u[j][i] = u[j][i] + delta;
  
    calc_variables();

    for(r = 1; r <= k*N; r++)
    {
      l = (r-1)/N + 1;
      m = r - (l-1)*N;
      gsl_matrix_set(H, r-1, q-1, (A[l][m]-gsl_matrix_get(H, r-1, q-1))/delta);
    }

    u[j][i] = u[j][i] - delta;
 
  }


  calc_variables();

    
    
}
