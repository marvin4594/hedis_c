#include <gsl/gsl_linalg.h>


typedef int bool;
#define TRUE   (1)
#define FALSE  (0)


//-----constants-----

extern const double c_pi;

extern const double c_grav;
extern const double tau_s;

extern const double ZERO;
extern const double HALF;
extern const double ONE;
extern const double TWO;
extern const double THREE;
extern const double FOUR;
extern const double SIX;
extern const double EIGHT;
extern const double TEN;

extern const double TWO_THIRD;
extern const double FOUR_THIRD;


//-----function prototypes-----

double H2F1(double z);
double min(double a, double b);
double max(double a, double b);
double minmod(double a, double b);
double sign(double a);

void alloc_arrays(void);
void calc_boundarys(void);
void calc_gravacc(void);
void calc_masses(void);
void calc_reconstruction(void);
void calc_timestep(void);
void calc_variables(void);
void calc_viscosity(void);
void disk(void);
void free_arrays(void);
void henyey(void);
void init_gravitation(void);
void init_grid(void);
void init_parameters(void);
void initial_conditions(void);
void solve_explicit(void);
void solve_implicit(void);
void solve_sle (void);
void store_results(void);


//-----global variables-----

extern FILE* file_timedev;
extern FILE* file_distrib;
extern FILE* file_general;

extern char *filename;


extern int mode_average, mode_boundarys, mode_flux, mode_gravitation, mode_initial_cond, mode_time_integr, mode_viscosity;
extern int N, N_disk, k;
extern int index_store_results, steps_store;
extern int error_code;
extern int iter_max, iter;
extern int timestep;

extern bool error;
extern bool general_stored_yet, store_by_time, store_flag;
extern bool SW_store_distrib, SW_store_general, SW_store_timedev;
extern bool SW_ag, SW_cem;

extern double accuracy, epsil, mu, num_visc, part;
extern double accr_eff, beta, cs;
extern double c_cfl, c_itw, d_time, d_time_cfl, d_time_itw, d_time_min, time, time_end, time_store, time_store_general;
extern double disk_inner_radius, disk_outer_radius;
extern double mass_bh, mass_central, mass_disk, mass_disk_init, mass_lost, mass_lost_edd;
extern double Mdot_i, Mdot_o, Mdot_i_old, Mdot_o_old, Mdot_bh, Mdot_disk, Mdot_edd;
extern double log_mass_central, log_mass_disk;
extern double s_0, sigma_0;


extern double *a_g;

extern double *angmom;
extern double *d2_angmom;

extern double *eta;
extern double *eta_tilde;
extern double *zeta;
extern double *zeta_tilde;

extern double *d_radius;
extern double *d_radius_tilde;
extern double *radius;
extern double *radius_tilde;
extern double *ring_area;

extern double **d2;
extern double *f;

extern double **charak;

extern double **k1;
extern double **k2;
extern double **k3;
extern double **k4;

extern double **d_u;
extern double **u;
extern double **u_dev;
extern double **u_old;

extern double **u_tilde;
extern double **u_tilde_m;
extern double **u_tilde_p;

extern double **F_num;
extern double **V;
extern double **W;
extern double **S;

extern double **A;
extern double **D;
extern double **D_old;

extern double **G;


extern gsl_matrix *H;

extern gsl_vector *right;
extern gsl_vector *x;

extern gsl_permutation * permut;
