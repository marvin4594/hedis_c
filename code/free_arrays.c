#include "global.h"
#include "nrutil.h"


void free_arrays(void)
{
  

  free_dvector(a_g,0,N+1);

  free_dvector(angmom,   -1,N+1);
  free_dvector(d2_angmom, 0,N);

  free_dvector(eta,       0,N+1);
  free_dvector(eta_tilde, 0,N);
  free_dvector(zeta,      0,N+1);
  free_dvector(zeta_tilde,0,N);

  free_dvector(d_radius,       0,N+1);
  free_dvector(d_radius_tilde,-1,N); 
  free_dvector(radius,        -1,N+1);
  free_dvector(radius_tilde,  -2,N+1);
  free_dvector(ring_area,      0,N+1);

  free_dmatrix(d2,1,3,0,N);
  free_dvector(f,1,N);

  free_dmatrix(charak,1,3,0,N+1);
  
  free_dmatrix(k1,1,3,1,N);
  free_dmatrix(k2,1,3,1,N);
  free_dmatrix(k3,1,3,1,N);
  free_dmatrix(k4,1,3,1,N);
  
  free_dmatrix(d_u,  1,3,0,N+1);
  free_dmatrix(u,    1,3,0,N+1);
  free_dmatrix(u_dev,1,3,0,N+1);
  free_dmatrix(u_old,1,3,0,N+1);

  free_dmatrix(u_tilde,  1,3,0,N);
  free_dmatrix(u_tilde_m,1,3,0,N);
  free_dmatrix(u_tilde_p,1,3,0,N);

  free_dmatrix(F_num,1,3,0,N);
  free_dmatrix(V,    1,3,0,N);
  free_dmatrix(W,    1,3,1,N);
  free_dmatrix(S,    1,3,1,N);

  free_dmatrix(A,    1,k,1,N);
  free_dmatrix(D,    1,k,1,N);
  free_dmatrix(D_old,1,k,1,N);

  free_dmatrix(G,0,N+1,0,N+1);

  
  gsl_matrix_free(H);
  
  gsl_vector_free(right);  
  gsl_vector_free (x);
  
  gsl_permutation_free (permut);


}