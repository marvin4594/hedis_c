#include "global.h"


void calc_gravacc(void)
{
  
  int i, j;
  
  for(i = 0; i <= N+1; i++)
  {
    a_g[i] = 0;
    for(j = 0; j <= N+1; j++)
      a_g[i] += G[i][j]*u[1][j]/radius[j];
  }
  
  
}

