#include <math.h>
#include <stdio.h> 
#include "global.h"


void init_parameters(void)
{

  int i, j;

  printf("\n");
  printf(" timestep     time/yr       d_t/yr        d_t(cfl)/yr   M_disk/M_sun  M_z/M_sun     M_BH/M_sun    Mass conserv.\n");
  printf(" --------------------------------------------------------------------------------------------------------------\n");


  N_disk = part*N;


  index_store_results = 0;
  store_flag = TRUE;
  general_stored_yet = FALSE;

  error = FALSE;

  timestep = 0;

  error_code = 0;

  d_time     = ZERO;
  d_time_cfl = ZERO;
  d_time_itw = ZERO;
  time       = ZERO;


  mass_disk    = pow(TEN,log_mass_disk);
  mass_central = pow(TEN,log_mass_central);

  mass_disk_init = mass_disk;

  mass_bh = mass_central;

  mass_lost     = ZERO;
  mass_lost_edd = ZERO;

  Mdot_i = ZERO;
  Mdot_o = ZERO;

  
  for(i = 0; i <= N+1; i++)
  {
    eta[i]        = ZERO;
    zeta[i]       = ZERO;
  }

  for(i = 0; i <= N; i++)
  {
    eta_tilde[i]  = ZERO;
    zeta_tilde[i] = ZERO;
  }


  for(j = 1; j <= 3; j++)
  {
    for(i = 1; i <= N; i++)
    {
      S[j][i]  = ZERO;
      W[j][i]  = ZERO;
    }

    for(i = 0; i <= N; i++)
    {
      V[j][i]  = ZERO;
      u_tilde[j][i] = ZERO;
    }
  
    for(i = 0; i <= N+1; i++)
    {
      d_u[j][i]     = ZERO;
      u[j][i]       = ZERO;
      u_dev[j][i]   = ZERO;
      u_old[j][i]   = ZERO;
    }
  }


  for(j = 1; j <= k; j++)
    for(i = 1; i <= N; i++)
    {
      A[j][i]     = ZERO;
      D[j][i]     = ZERO;
      D_old[j][i] = ZERO;
    }


}

