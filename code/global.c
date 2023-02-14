#include <gsl/gsl_linalg.h>
#include "global.h"


FILE* file_general;
FILE* file_distrib;
FILE* file_timedev;

char *filename;


int mode_average, mode_boundarys, mode_flux, mode_gravitation, mode_initial_cond, mode_time_integr, mode_viscosity;
int N, N_disk, k;
int index_store_results, steps_store;
int error_code;
int iter_max, iter;
int timestep;

bool error;
bool general_stored_yet, store_by_time, store_flag;
bool SW_store_distrib, SW_store_general, SW_store_timedev;
bool SW_ag, SW_cem;

double accuracy, epsil, mu, num_visc, part;
double accr_eff, beta, cs;
double c_cfl, c_itw, d_time, d_time_cfl, d_time_itw, d_time_min, time, time_end, time_store, time_store_general;
double disk_inner_radius, disk_outer_radius;
double mass_bh, mass_central, mass_disk, mass_disk_init, mass_lost, mass_lost_edd;
double Mdot_i, Mdot_o, Mdot_i_old, Mdot_o_old, Mdot_bh, Mdot_disk, Mdot_edd;
double log_mass_central, log_mass_disk;
double s_0, sigma_0;


double *a_g;

double *angmom;
double *d2_angmom;

double *eta;
double *eta_tilde;
double *zeta;
double *zeta_tilde;

double *d_radius;
double *d_radius_tilde;
double *radius;
double *radius_tilde;
double *ring_area;

double **d2;
double *f;

double **charak;

double **k1;
double **k2;
double **k3;
double **k4;

double **d_u;
double **u;
double **u_dev;
double **u_old;

double **u_tilde;
double **u_tilde_m;
double **u_tilde_p;

double **F_num;
double **V;
double **W;
double **S;

double **A;
double **D;
double **D_old;

double **G;


gsl_matrix *H;

gsl_vector *right;
gsl_vector *x;

gsl_permutation * permut;
