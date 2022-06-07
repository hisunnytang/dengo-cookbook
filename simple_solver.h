/*

The generalized rate data type holders.

*/


/* stdlib, hdf5, local includes */

#include "omp.h"

#include "time.h"
#include "sys/time.h"
#include "stdlib.h"
#include "math.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "stdio.h"
#include "string.h"

/* header files for CVODES/SUNDIALS */
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#ifdef CVKLU
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_klu.h>
#endif


/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


#ifndef MAX_NCELLS
#define MAX_NCELLS 1024
#endif

#ifndef NTHREADS
#define NTHREADS 8
#endif

#define NSPECIES 3
#define DMAX(A,B) ((A) > (B) ? (A) : (B))
#define DMIN(A,B) ((A) < (B) ? (A) : (B))

 

int simple_main(int argc, char **argv);



typedef struct simple_data {
    /* All of the network bins will be the same width */
    double dbin;
    double idbin;
    double bounds[2];
    int nbins;

    /* These will be for bins in redshift space */
    double d_zbin;
    double id_zbin;
    double z_bounds[2];
    int n_zbins;

    /* For storing and passing around
       redshift information */
    double current_z;
    double zdef;
    double dz;

    double Ts[NTHREADS][MAX_NCELLS];
    double Tdef[NTHREADS][MAX_NCELLS]; /* t1, t2, tdef */
    double dT[NTHREADS][MAX_NCELLS]; /* t1, t2, tdef */
    double logTs[NTHREADS][MAX_NCELLS];
    double invTs[NTHREADS][MAX_NCELLS];
    double dTs_ge[NTHREADS][MAX_NCELLS];

    /* Now we do all of our cooling and chemical tables */
    double r_k13[1024];
    double rs_k13[NTHREADS][MAX_NCELLS];
    double drs_k13[NTHREADS][MAX_NCELLS];
    
    double r_k22[1024];
    double rs_k22[NTHREADS][MAX_NCELLS];
    double drs_k22[NTHREADS][MAX_NCELLS];
    
    double c_h2formation_h2mcool[1024];
    double cs_h2formation_h2mcool[NTHREADS][MAX_NCELLS];
    double dcs_h2formation_h2mcool[NTHREADS][MAX_NCELLS];
    double c_h2formation_h2mheat[1024];
    double cs_h2formation_h2mheat[NTHREADS][MAX_NCELLS];
    double dcs_h2formation_h2mheat[NTHREADS][MAX_NCELLS];
    double c_h2formation_ncrd1[1024];
    double cs_h2formation_ncrd1[NTHREADS][MAX_NCELLS];
    double dcs_h2formation_ncrd1[NTHREADS][MAX_NCELLS];
    double c_h2formation_ncrd2[1024];
    double cs_h2formation_ncrd2[NTHREADS][MAX_NCELLS];
    double dcs_h2formation_ncrd2[NTHREADS][MAX_NCELLS];
    double c_h2formation_ncrn[1024];
    double cs_h2formation_ncrn[NTHREADS][MAX_NCELLS];
    double dcs_h2formation_ncrn[NTHREADS][MAX_NCELLS];
    
    int bin_id[NTHREADS][MAX_NCELLS];
    int ncells;
    


    // gamma as a function of temperature
    double g_gammaH2_1[1024];
    double g_dgammaH2_1_dT[1024];

    // store the gamma for that particular step
    double gammaH2_1[NTHREADS][MAX_NCELLS];
    double dgammaH2_1_dT[NTHREADS][MAX_NCELLS];
    
    // store 1 / (gamma - 1)
    double _gammaH2_1_m1[NTHREADS][MAX_NCELLS];
    
    double g_gammaH2_2[1024];
    double g_dgammaH2_2_dT[1024];

    // store the gamma for that particular step
    double gammaH2_2[NTHREADS][MAX_NCELLS];
    double dgammaH2_2_dT[NTHREADS][MAX_NCELLS];
    
    // store 1 / (gamma - 1)
    double _gammaH2_2_m1[NTHREADS][MAX_NCELLS];
    

    double scale[NTHREADS][3 * MAX_NCELLS ];
    double inv_scale[NTHREADS][3 * MAX_NCELLS];
    

    int nstrip;
    double mdensity[NTHREADS][MAX_NCELLS];
    double inv_mdensity[NTHREADS][MAX_NCELLS];

    
    
    
    const char *dengo_data_file;
    
    double reltol;
    double floor_value;
} simple_data;


/* Declare ctype RHS and Jacobian */
typedef int(*rhs_f)( realtype, N_Vector , N_Vector , void * );
#ifndef CVSPILS
typedef int(*jac_f)( realtype, N_Vector  , N_Vector , SUNMatrix , void *, N_Vector, N_Vector, N_Vector);
#endif
#ifdef CVSPILS
typedef int(*jac_f)(N_Vector , N_Vector , realtype,
             N_Vector, N_Vector,
             void *user_data, N_Vector);
#endif


void *setup_cvode_solver( rhs_f f, jac_f Jac,  int NEQ, 
        simple_data *data, SUNLinearSolver LS, SUNMatrix A, N_Vector y, double reltol, N_Vector abstol);

int cvode_solver( void *cvode_mem, double *output, int NEQ, double *dt, simple_data * data, N_Vector y, double reltol, N_Vector abstol );


simple_data *simple_setup_data( const char *FileLocation, int *NumberOfFields, char ***FieldNames);
void simple_read_rate_tables(simple_data*);
void simple_read_cooling_tables(simple_data*);
void simple_read_gamma(simple_data*);
void simple_interpolate_gamma(simple_data*, int );

void setting_up_extra_variables( simple_data * data, double * input, int nstrip );

int dengo_evolve_simple (double dtf, double &dt, double z,
                                     double *input, double *rtol,
                                     double *atol, unsigned long dims,
                                     simple_data *data, double *temp);

double evolve_in_batches( void * cvode_mem, N_Vector y_vec, N_Vector abstol,  
                          double reltol,double *input, int v_size, int d, int start_idx, 
                          int MAX_ITERATION, double dtf, simple_data *data );


 




#ifndef CVSPILS
#ifdef  CVKLU
int calculate_sparse_jacobian_simple( realtype t,
                                        N_Vector y, N_Vector fy,
                                        SUNMatrix J, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3);
#else
int calculate_jacobian_simple( realtype t,
               N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
#endif

#ifdef CVSPILS
int calculate_JacTimesVec_simple
            (N_Vector v, N_Vector Jv, realtype t,
             N_Vector y, N_Vector fy,
             void *user_data, N_Vector tmp);
#endif

int calculate_rhs_simple(realtype t, N_Vector y, N_Vector ydot, void *user_data);
void ensure_electron_consistency(double *input, double *equil_array, unsigned long nstrip, int nchem);
void temperature_from_mass_density(double *input, int nstrip, int nchem, 
                                   double *strip_temperature);

int simple_calculate_temperature(simple_data *data, double *input, int nstrip, int nchem);




typedef struct code_units
{

  int comoving_coordinates;
  double density_units;
  double length_units;
  double time_units;
  double velocity_units;
  double a_units;
  double a_value;

} code_units;


typedef struct dengo_field_data
{

  unsigned long int nstrip;
  unsigned long int ncells; 
  // let's just pass them passively through field_data
  double reltol;
  double floor_value;
  // This should be updated dynamically 
  // with dengo
  double *density;
  double *H2_1_density;
  double *H_1_density;
  double *ge_density;
    
  double *CoolingTime;
  double *MolecularWeight;
  double *temperature;
  double *Gamma;
  double *Pressure;

  int *grid_start;
  int *grid_end;
  int *grid_dimension;

  const char *dengo_data_file;
  code_units *units;
} dengo_field_data;


// we can call this and pass reltol and floor_value to solve_chemistry
int simple_solve_chemistry( code_units *units, dengo_field_data *field_data, double dt );
int simple_solve_chemistry_dt( code_units *units, dengo_field_data *field_data, double* rtol, double *atol, double dt );

int dengo_estimate_cooling_time( code_units* units, dengo_field_data * field_data );

int simple_calculate_cooling_timescale( double *cooling_time, double *input, int nstrip, simple_data *data);

int dengo_calculate_pressure_enzo( code_units*, dengo_field_data* );
int dengo_calculate_gamma_enzo( code_units*, dengo_field_data* );
int dengo_calculate_temperature_enzo( code_units*, dengo_field_data* );

int dengo_calculate_pressure( code_units*, dengo_field_data* );
int dengo_calculate_temperature( code_units*, dengo_field_data* );
int dengo_calculate_gamma( double* gamma_eff, simple_data*, double* input, int nstrip);
int dengo_calculate_mean_molecular_weight( code_units*, dengo_field_data * );

int calculate_equilibrium_abundance( simple_data *data, double *input, 
	int nstrip, unsigned long d, unsigned long dims, double *equil_array);
// Enzo interface
int dengo_estimate_cooling_time_enzo( code_units* units, dengo_field_data *field_data );

int simple_solve_chemistry_enzo( code_units *units, dengo_field_data *field_data, double dt );


int reshape_to_dengo_field_data( code_units* units, dengo_field_data *field_data, double* input );
int flatten_dengo_field_data( code_units* units, dengo_field_data *field_data, double *input );
int reshape_to_dengo_field_data_enzo( code_units* units, dengo_field_data *field_data, double* input, double *temp );
int flatten_dengo_field_data_enzo( code_units* units, dengo_field_data *field_data, double *input );


void ensure_species_conservation(double *input, double *mdensity, double *equil_array, simple_data *data, int nstrip, unsigned long d, unsigned long dims, int nchem);