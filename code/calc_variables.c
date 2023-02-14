#include <math.h> 
#include "global.h"
#include "nrutil.h"


void calc_variables(void)
{

  
  int i, j;
  double a_p, a_m;


  calc_boundarys();
  calc_gravacc();



  switch (k)
  {

    case 3:

      calc_reconstruction();
      calc_viscosity();


      for (i = 0; i <= N; i++)
      {

        a_p = max( max(u_tilde_p[2][i]/u_tilde_p[1][i]+cs, u_tilde_m[2][i]/u_tilde_m[1][i]+cs), ZERO );
        a_m = min( min(u_tilde_p[2][i]/u_tilde_p[1][i]-cs, u_tilde_m[2][i]/u_tilde_m[1][i]-cs), ZERO );
 

        F_num[1][i] = ( a_p*u_tilde_m[2][i] - a_m*u_tilde_p[2][i] - a_p*a_m*num_visc*(u_tilde_m[1][i] - u_tilde_p[1][i]) )/(a_p - a_m);
        F_num[2][i] =   (    a_p*( pow(u_tilde_m[2][i],2)/u_tilde_m[1][i] + (cs*cs)*u_tilde_m[1][i] ) 
                           - a_m*( pow(u_tilde_p[2][i],2)/u_tilde_p[1][i] + (cs*cs)*u_tilde_p[1][i] )
                        - a_p*a_m*num_visc*(u_tilde_m[2][i] - u_tilde_p[2][i]) )/(a_p - a_m);
        F_num[3][i] =   (  a_p*( u_tilde_m[2][i]*u_tilde_m[3][i]/u_tilde_m[1][i] ) - a_m*( u_tilde_p[2][i]*u_tilde_p[3][i]/u_tilde_p[1][i] )
                        - a_p*a_m*num_visc*(u_tilde_m[3][i] - u_tilde_p[3][i]) )/(a_p - a_m);
      }


      for (i = 0; i <= N; i++)
      {

        V[2][i] =  + (zeta_tilde[i] + (4/3)*eta_tilde[i])*radius_tilde[i]
                            *(u_tilde[1][i]*(u[2][i+1]-u[2][i]) - u_tilde[2][i]*(u[1][i+1]-u[1][i]))/(pow(u_tilde[1][i],2)*d_radius_tilde[i])
                         + (zeta_tilde[i] - TWO_THIRD*eta_tilde[i])*u_tilde[2][i]/u_tilde[1][i];

        V[3][i] =  eta_tilde[i]*( radius_tilde[i]*(u_tilde[1][i]*(u[3][i+1]-u[3][i]) - u_tilde[3][i]*(u[1][i+1]-u[1][i]))
                                                /(d_radius_tilde[i]*pow(u_tilde[1][i],2)) - TWO*u_tilde[3][i]/u_tilde[1][i] );
 
      }
    

      for (i = 1; i <= N; i++)
      {

        W[2][i] =  - (zeta[i] + (4/3)*eta[i])*u[2][i]/(u[1][i]*radius[i]) - (zeta[i] - TWO_THIRD*eta[i])
                          *(u[1][i]*(u_tilde[2][i] - u_tilde[2][i-1]) - u[2][i]*(u_tilde[1][i] - u_tilde[1][i-1]))/(d_radius[i]*pow(u[1][i],2));

        S[2][i] =  u[1][i]*( + pow((u[3][i]/u[1][i]),2)/pow(radius[i],3) + (cs*cs)/radius[i] - c_grav*mass_central/pow(radius[i],2) + a_g[i] );

        for (j = 1; j <= 3; j++)
          D[j][i] = + ((F_num[j][i] - F_num[j][i-1]) - (V[j][i] - V[j][i-1])) /d_radius[i] - W[j][i] - S[j][i];

      }


      for (j = 1; j <= 3; j++)
        for (i = 1; i <= N; i++)
          A[j][i] = (u[j][i]-u_old[j][i])  +  d_time*( mu*D[j][i] + (1-mu)*D_old[j][i] );

    break;




    case 1:

      angmom[-1]    = sqrt(c_grav*mass_central*radius[-1]);
  
      for(i = 0; i <= N+1; i++)
      {
        angmom[i] = sqrt(c_grav*mass_central*radius[i] - a_g[i]*pow(radius[i],3));
        u[3][i] = angmom[i]*u[1][i];
      }
  
      calc_viscosity();

      for(i = 0; i <= N; i++)
        d2_angmom[i] = d2[1][i]*angmom[i+1] + d2[2][i]*angmom[i] + d2[3][i]*angmom[i-1];      //second derivative of angmom

      for(i = 1; i <= N; i++)
          u[2][i] = (   (eta[i]*radius[i] - eta[i-1]*radius[i-1])*(angmom[i] - angmom[i-1])
                      - 2*d_radius_tilde[i-1]*(eta[i]*angmom[i] - eta[i-1]*angmom[i-1])
                      + f[i]*( eta[i]*d2_angmom[i] + eta[i-1]*d2_angmom[i-1] ) )
                    /((angmom[i]-angmom[i-1])*d_radius_tilde[i-1]);

      for(i = 1; i <= N; i++)
      {
        D[1][i] = ( mu*(u[2][i+1]-u[2][i]) + (ONE-mu)*(u_old[2][i+1]-u_old[2][i]) )/d_radius[i];
        A[1][i] =   (u[1][i]-u_old[1][i]) + d_time*D[1][i];
      }


    break;
  
  }


}