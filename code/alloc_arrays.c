#include "global.h"
#include "nrutil.h"


void alloc_arrays(void)
{


  a_g            = dvector(0,N+1);

  angmom         = dvector(-1,N+1);
  d2_angmom      = dvector( 0,N);

  eta            = dvector(0,N+1);
  eta_tilde      = dvector(0,N);
  zeta           = dvector(0,N+1);
  zeta_tilde     = dvector(0,N);

  d_radius       = dvector( 0,N+1);
  d_radius_tilde = dvector(-1,N);
  radius         = dvector(-1,N+1);
  radius_tilde   = dvector(-2,N+1);
  ring_area      = dvector( 0,N+1);

  d2             = dmatrix(1,3,0,N);
  f              = dvector(1,N);

  charak         = dmatrix(1,3,0,N+1);

  k1             = dmatrix(1,3,1,N);
  k2             = dmatrix(1,3,1,N);
  k3             = dmatrix(1,3,1,N);
  k4             = dmatrix(1,3,1,N);

  d_u            = dmatrix(1,3,0,N+1);
  u              = dmatrix(1,3,0,N+1);
  u_dev          = dmatrix(1,3,0,N+1);
  u_old          = dmatrix(1,3,0,N+1);

  u_tilde        = dmatrix(1,3,0,N);
  u_tilde_m      = dmatrix(1,3,0,N);
  u_tilde_p      = dmatrix(1,3,0,N);

  F_num          = dmatrix(1,3,0,N);
  V              = dmatrix(1,3,0,N);
  W              = dmatrix(1,3,1,N);
  S              = dmatrix(1,3,1,N);

  A              = dmatrix(1,k,1,N);
  D              = dmatrix(1,k,1,N);
  D_old          = dmatrix(1,k,1,N);

  G              = dmatrix(0,N+1,0,N+1);
  
  
  H              = gsl_matrix_alloc (k*N, k*N);

  right          = gsl_vector_alloc(k*N);
  x              = gsl_vector_alloc(k*N);
  
  permut         = gsl_permutation_alloc(k*N);


}