#include <math.h>
#include "global.h"


void init_grid(void)
{

  int i;
  
  
  for( i = -2; i <= N+1; i++)
    radius_tilde[i] = disk_inner_radius * pow(TEN, (i) * log10(disk_outer_radius/disk_inner_radius) / N_disk );
  
  
  for( i = -1; i <= N+1; i++)
    radius[i] = sqrt( radius_tilde[i]*radius_tilde[i-1] );
  
  
  for( i = 0; i <= N+1; i++)
    d_radius[i] = radius_tilde[i] - radius_tilde[i-1] ;
  
  
  for( i = -1; i <= N; i++)
    d_radius_tilde[i] = radius[i+1] - radius[i] ;
  
  
  for( i = 0; i <= N+1; i++)
    ring_area[i] = TWO*c_pi*radius[i]*d_radius[i];
  
  
  for( i = 0; i <= N; i++)
  {
    d2[1][i] =     ONE/(d_radius[i]*d_radius_tilde[i]);
    d2[2][i] = - ( ONE/(d_radius[i]*d_radius_tilde[i]) + ONE/(d_radius[i]*d_radius_tilde[i-1]) );
    d2[3][i] =     ONE/(d_radius[i]*d_radius_tilde[i-1]);
  }
  
  
  for( i = 1; i <= N; i++)
    f[i] = HALF*pow(d_radius_tilde[i-1],2)*radius_tilde[i-1];
  

}