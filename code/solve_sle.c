#include <gsl/gsl_linalg.h>
#include "global.h"
     

void solve_sle (void)
{
  
  int i, j, l, s;
  
  
  for(j = 1, l = 0; j <= k; j++)
    for(i = 1; i <= N; i ++, l++)
      gsl_vector_set (right, l, -A[j][i]);
  
  
  error_code = gsl_linalg_LU_decomp (H, permut, &s);
  
  if (error_code == 0) error_code = gsl_linalg_LU_solve (H, permut, right, x);
     
  if (error_code != 0) error = TRUE;


  for(j = 1, l = 0; j <= k; j++)
    for(i = 1; i <= N; i++, l++)
      d_u[j][i] = gsl_vector_get(x, l);
 

}
