#include <math.h> 
#include "global.h"


void init_gravitation(void)
{

  
  int i, j;

  
  
  for(i = 0; i <= N+1; i++)
    for(j = 0; j <= N+1; j++)
      G[i][j] = ZERO;                               //Kepler-gravitation



  switch (mode_gravitation)
  {

    case 2:     //monopole approximation, low order

      for(i = 0; i <= N+1; i++)
        for(j = 0; j <= i; j++)
          G[i][j] = - c_grav*ring_area[j]/pow(radius[i],2);

    break;


    case 3:     //monopole approximation, high order

      for(i = 1; i <= N; i++) {
        for(j = 1; j <= i; j++)
          G[i][j-1] = TWO*c_pi*(  - (ONE/THREE)*(pow(radius[j],3)-pow(radius[j-1],3))/(radius[j]-radius[j-1])
                  + (ONE/FOUR)*(pow(radius[j],2)-pow(radius[j-1],2))
                  + (ONE/FOUR)*(radius[j]+radius[j-1])*(pow(radius[j],2)-pow(radius[j-1],2))/(radius[j]-radius[j-1]) );
        for(j = 1; j <= i; j++)
          G[i][j] = G[i][j] + TWO*c_pi*(  + (ONE/THREE)*(pow(radius[j],3)-pow(radius[j-1],3))/(radius[j]-radius[j-1])
                     + (ONE/FOUR)*(pow(radius[j],2)-pow(radius[j-1],2))
                     - (ONE/FOUR)*(radius[j]+radius[j-1])*(pow(radius[j],2)-pow(radius[j-1],2))/(radius[j]-radius[j-1])  );
        for(j = 0; j <= N+1; j++)
          G[i][j] = -c_grav*G[i][j]/pow(radius[i],2);
      }
    
    break;
    

    case 4:     //full selfgravitation

      for(i = 0; i <= N+1; i++)
        for(j = 0; j <= N+1; j++)
          G[i][j] = - TWO*c_pi*c_grav*
                     (  (radius[i]*pow(radius_tilde[j],2)/(TWO*pow((radius[i]+radius_tilde[j]),3))) *
                             H2F1( FOUR*radius[i]*radius_tilde[j]/pow((radius[i]+radius_tilde[j]),2))
                      - (radius[i]*pow(radius_tilde[j-1],2)/(TWO*pow((radius[i]+radius_tilde[j-1]),3))) *
                             H2F1( FOUR*radius[i]*radius_tilde[j-1]/pow((radius[i]+radius_tilde[j-1]),2))  );

    break;
  
  }


}