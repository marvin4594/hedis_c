#include <math.h>
#include "global.h"


void calc_boundarys(void)
{


  int i, j;
  double v_in, v_out, xi, tau;



  switch (mode_boundarys) {

    case 1:     //non-reflecting boundary conditions

      for (j = 1; j <= 3; j++)
        for (i = 0; i <= N+1; i++)
          charak[j][i] = ZERO;
      
      i = 2;
      charak[1][i] = cs*log((u[1][i]*radius[i-1])/(u[1][i-1]*radius[i])) - (u[2][i]/u[1][i]-u[2][i-1]/u[1][i-1]);
      charak[2][i] = ( (u[3][i]/(u[1][i]*radius[i])) - (u[3][i-1]/(u[1][i-1]*radius[i-1])) );
      charak[3][i] = cs*log((u[1][i]*radius[i-1])/(u[1][i-1]*radius[i])) + (u[2][i]/u[1][i]-u[2][i-1]/u[1][i-1]);

      i = N;
      charak[1][i] =    cs*log((u[1][i]*radius[i-1])/(u[1][i-1]*radius[i])) - (u[2][i]/u[1][i]-u[2][i-1]/u[1][i-1]);
      charak[2][i] =    ( (u[3][i]/(u[1][i]*radius[i])) - (u[3][i-1]/(u[1][i-1]*radius[i-1])) );
      charak[3][i] =    cs*log((u[1][i]*radius[i-1])/(u[1][i-1]*radius[i])) + (u[2][i]/u[1][i]-u[2][i-1]/u[1][i-1]);


      v_in = u[2][2]/u[1][2] + (u[2][2]/u[1][2] - u[2][1]/u[1][1])*(radius_tilde[0]-radius[1])/d_radius_tilde[1];          //velocity at inner radius

      if ((v_in-cs) <= ZERO) charak[1][1] = charak[1][2];
      if ( v_in     <= ZERO) charak[2][1] = charak[2][2];
      if ((v_in+cs) <= ZERO) charak[3][1] = charak[3][2];


      v_out = u[2][N]/u[1][N] + (u[2][N]/u[1][N] - u[2][N-1]/u[1][N-1])*(radius_tilde[N]-radius[N])/d_radius_tilde[N-1];   //velocity at outer radius

      if (( v_out-cs) >= ZERO) charak[1][N+1] = charak[1][N];
      if (  v_out     >= ZERO) charak[2][N+1] = charak[2][N];
      if (( v_out+cs) >= ZERO) charak[3][N+1] = charak[3][N];


      u[1][0]   = ( radius[0]/radius[1]) * u[1][1] * exp( - (HALF/cs) * (charak[1][1] + charak[3][1]) );
      u[2][0]   = ( u[2][1]/u[1][1]             - HALF * (charak[3][1] - charak[1][1]) ) * u[1][0];
      u[3][0]   = ( u[3][1]/(u[1][1]*radius[1]) - charak[2][1] ) * radius[0] * u[1][0];

      u[1][N+1] = (radius[N+1]/radius[N]) * u[1][N] * exp ( + (HALF/cs) * (charak[3][N+1] + charak[1][N+1]) );
      u[2][N+1] = ( u[2][N]/u[1][N]               + HALF * (charak[3][N+1] - charak[1][N+1])) * u[1][N+1];
      u[3][N+1] = ( (u[3][N]/(radius[N]*u[1][N])) + charak[2][N+1] )*radius[N+1]*u[1][N+1];

    break;

  
    case 2:     //reflecting boundary conditions

      u[1][0] =   (radius[0]/radius[1]) * u[1][1];
      u[2][0] = - (u[2][1]/u[1][1]) * u[1][0];
      u[3][0] =   (radius[0]/radius[1]) * (u[1][0]/u[1][1]) * u[3][1];

      u[1][N+1] =   (radius[N+1]/radius[N]) * u[N][1];
      u[2][N+1] = - (u[2][N]/u[1][N]) * u[1][N+1];
      u[3][N+1] =   (radius[N+1]/radius[N]) * (u[1][N+1]/u[1][N]) * u[3][N];

    break;

  
    case 3:     //no-gradients boundary conditions

      u[1][0] = (radius[0]/radius[1]) * u[1][1];
      u[2][0] = (u[2][1]/u[1][1]) * u[1][0];
      u[3][0] = (radius[0]/radius[1]) * (u[1][0]/u[1][1]) * u[3][1];

      u[1][N+1] = (radius[N+1]/radius[N]) * u[N][1];
      u[2][N+1] = (u[2][N]/u[1][N]) * u[1][N+1];
      u[3][N+1] = (radius[N+1]/radius[N]) * (u[1][N+1]/u[1][N]) * u[3][N];

    break;


    case 4:     //boundarys for the pringle-model

      u[1][0]   = pow(10,-10) * u[1][1];
      u[1][N+1] = pow(10,-10) * u[1][N];

      u[2][0]   = pow(10,-10) * u[2][1];
      u[2][N+1] = pow(10,-10) * u[1][N];

    break;

  
    case 5:     //boundarys for the pringle-disk

      tau     = (THREE*beta/(FOUR*s_0))*time;

      u[1][0]   = radius[0]  *(sigma_0/(TWO*pow(sqrt(radius[0]/s_0),3)*sqrt(c_pi*tau)))
                   * (exp(-pow((sqrt(radius[0]/s_0)-1),2)/(FOUR*tau)) - exp(-pow((sqrt(radius[0]/s_0)+1),2)/(FOUR*tau)));
      u[1][N+1] = radius[N+1]  *(sigma_0/(TWO*pow(sqrt(radius[N+1]/s_0),3)*sqrt(c_pi*tau)))
                   * (exp(-pow((sqrt(radius[N+1]/s_0)-1),2)/(FOUR*tau)) - exp(-pow((sqrt(radius[N+1]/s_0)+1),2)/(FOUR*tau)));

      xi = sqrt(radius[0]/s_0);
      u[2][0]   = - ((THREE*beta*s_0*sigma_0)/(EIGHT*sqrt(c_pi*pow(tau,3))))
                    * (  (xi+ONE)*exp(- (pow((xi+ONE),2)/(FOUR*tau))) - (xi-ONE)*exp(- (pow((xi-ONE),2)/(FOUR*tau))));

      xi = sqrt(radius[N+1]/s_0);
      u[2][N+1] = - ((THREE*beta*s_0*sigma_0)/(EIGHT*sqrt(c_pi*pow(tau,3))))
                  * (  (xi+ONE)*exp(- (pow((xi+ONE),2)/(FOUR*tau))) - (xi-ONE)*exp(- (pow((xi-ONE),2)/(FOUR*tau))));

    break;
  
  }


}
