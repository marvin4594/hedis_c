#include "global.h"


void calc_masses(void)
{

  int i;


  Mdot_i_old = Mdot_i;
  Mdot_o_old = Mdot_o;

  Mdot_i = TWO*c_pi*u[2][1];
  Mdot_o = TWO*c_pi*u[2][N+1];


  Mdot_disk = - HALF*(Mdot_i+Mdot_i_old);
  Mdot_edd = mass_bh/(accr_eff*tau_s);

  Mdot_bh = max(min(Mdot_disk, Mdot_edd),ZERO);


  mass_bh += d_time*Mdot_bh;


  if (SW_ag)
  {
    if (SW_cem)
    {
      mass_central += d_time*Mdot_disk;
    }
    else
    {
      mass_central = mass_bh;
    }
  }


  for(i = 1, mass_disk = 0; i <= N; i++)
    mass_disk += u[1][i]*d_radius[i];

  mass_disk *= TWO*c_pi;

  mass_lost += d_time*HALF*( Mdot_o - Mdot_i + Mdot_o_old - Mdot_i_old );       //total mass loss

  mass_lost_edd += d_time*max(Mdot_disk-Mdot_edd,ZERO);                         //mass loss due to Eddington-limit


}