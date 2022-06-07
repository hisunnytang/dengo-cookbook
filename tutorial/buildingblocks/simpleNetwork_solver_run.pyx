cimport numpy as np
import numpy as np
import time
from libc.stdlib cimport malloc, free
from libcpp cimport bool
from cpython cimport array

# NSPECIES here is N in the .C.template file
DEF NSPECIES = 4
DEF MAX_NCELLS=1

cdef extern from "simpleNetwork_solver.h":
    cdef int _MAX_NCELLS  "MAX_NCELLS"
    cdef int _NSPECIES "NSPECIES"

    ctypedef struct dengo_field_data:
        unsigned long int nstrip;
        unsigned long int ncells;
        double *H_1_density;
        double *H_2_density;
        double *de_density;
        double *ge_density;
        double *density;
        double *CoolingTime;
        double *MolecularWeight;
        double *temperature;
        double *Gamma;
        double reltol;
        double floor_value;
        const char *dengo_data_file;

    ctypedef struct code_units:
        int comoving_coordinates
        double density_units
        double length_units
        double time_units
        double velocity_units
        double a_units
        double a_value


    ctypedef struct simpleNetwork_data:
        double dbin
        double idbin
        double bounds[2]
        int nbins

        double d_zbin
        double id_zbin
        double z_bounds[2]
        int n_zbins

        double current_z
        double zdef
        double dz

        double Ts[MAX_NCELLS]
        double Tdef[MAX_NCELLS]
        double dT[MAX_NCELLS]
        double logTs[MAX_NCELLS]
        double dTs_ge[MAX_NCELLS]
        double r_k01[1024]
        double rs_k01[MAX_NCELLS]
        double drs_k01[MAX_NCELLS]
        double r_k02[1024]
        double rs_k02[MAX_NCELLS]
        double drs_k02[MAX_NCELLS]
        double c_reHII_reHII[1024]
        double cs_reHII_reHII[MAX_NCELLS]
        double dcs_reHII_reHII[MAX_NCELLS]
        
        int bin_id[MAX_NCELLS]
        int ncells

    double dengo_evolve_simpleNetwork (double dtf, double &dt, double z,
                                         double *input, double *rtol,
                                         double *atol, int dims,
                                         simpleNetwork_data *data, double *temp)
    int simpleNetwork_solve_chemistry_dt( code_units *units, dengo_field_data *field_data, double* rtol, double* atol, double dt)
    int simpleNetwork_solve_chemistry( code_units *units, dengo_field_data *field_data, double dt );
    int simpleNetwork_main(int argc, char **argv);



def main_run_simpleNetwork():
    t1 = time.time()
    simpleNetwork_main(0, NULL)
    t2 = time.time()
    print("Total elapsed time: %0.3e" % (t2-t1))

def run_simpleNetwork(ics, double tf, int niter = 10000,
                        int intermediate = 1, float z = -1.0, bool adaptive_step = True,
			float reltol = 1.0e-5, float floor_value = 1.0e-20):
    #assert(_MAX_NCELLS == MAX_NCELLS)
    assert(_NSPECIES == NSPECIES)
    cdef np.ndarray[np.float64_t, ndim=1] H_1_arr = ics["H_1"]
    if not H_1_arr.flags['C_CONTIGUOUS']:
        H_1_arr = np.ascontiguousarray( H_1_arr )
    cdef double[::1] H_1_memview = H_1_arr
    cdef np.ndarray[np.float64_t, ndim=2] H_1_int
    cdef np.ndarray[np.float64_t, ndim=1] H_2_arr = ics["H_2"]
    if not H_2_arr.flags['C_CONTIGUOUS']:
        H_2_arr = np.ascontiguousarray( H_2_arr )
    cdef double[::1] H_2_memview = H_2_arr
    cdef np.ndarray[np.float64_t, ndim=2] H_2_int
    cdef np.ndarray[np.float64_t, ndim=1] de_arr = ics["de"]
    if not de_arr.flags['C_CONTIGUOUS']:
        de_arr = np.ascontiguousarray( de_arr )
    cdef double[::1] de_memview = de_arr
    cdef np.ndarray[np.float64_t, ndim=2] de_int
    cdef np.ndarray[np.float64_t, ndim=1] ge_arr = ics["ge"]
    if not ge_arr.flags['C_CONTIGUOUS']:
        ge_arr = np.ascontiguousarray( ge_arr )
    cdef double[::1] ge_memview = ge_arr
    cdef np.ndarray[np.float64_t, ndim=2] ge_int
    cdef np.ndarray[np.float64_t, ndim=1] CoolingTime;
    cdef np.ndarray[np.float64_t, ndim=1] Gamma;
    cdef np.ndarray[np.float64_t, ndim=1] Temperature;
    cdef np.ndarray[np.float64_t, ndim=1] MolecularWeight;

    cdef np.ndarray[np.float64_t, ndim=2] CoolingTime_int;
    cdef np.ndarray[np.float64_t, ndim=2] Gamma_int;
    cdef np.ndarray[np.float64_t, ndim=2] Temperature_int;
    cdef np.ndarray[np.float64_t, ndim=2] MolecularWeight_int;

    cdef np.ndarray[np.float64_t, ndim=1] result_int
    cdef np.ndarray[np.float64_t, ndim=2] temp_int
    cdef np.ndarray[np.float64_t, ndim=1] t_int
    cdef np.ndarray[np.float64_t, ndim=1] dt_int

    cdef int i, j, k, iter
    cdef int dims = ge_arr.shape[0]
    cdef int NTOT = NSPECIES * dims
    cdef double *input = <double *> malloc(NTOT * sizeof(double))
    cdef double *prev = <double *> malloc(NTOT * sizeof(double))

    cdef dengo_field_data *field_data = <dengo_field_data *> malloc(sizeof(dengo_field_data));
    cdef code_units       * units     = <code_units *> malloc(sizeof(code_units));
    units.density_units = 1.67e-24;
    units.length_units  = 1.0;
    units.time_units    = 1.0;
    units.velocity_units = 1.0;
    CoolingTime = np.ascontiguousarray(np.zeros((dims), dtype=np.double))
    Gamma       = np.ascontiguousarray(np.zeros((dims), dtype=np.double))
    Temperature = np.ascontiguousarray(np.zeros((dims), dtype=np.double))
    MolecularWeight = np.ascontiguousarray(np.zeros((dims), dtype=np.double))
    cdef double[::1] ctime_memview = CoolingTime
    cdef double[::1] gamma_memview = Gamma
    cdef double[::1] temp_memview  = Temperature
    cdef double[::1] mweight_memview = MolecularWeight
    field_data.H_1_density = &H_1_memview[0];
    field_data.H_2_density = &H_2_memview[0];
    field_data.de_density = &de_memview[0];
    field_data.ge_density = &ge_memview[0];
    field_data.CoolingTime = &ctime_memview[0]
    field_data.Gamma       = &gamma_memview[0]
    field_data.temperature = &temp_memview[0]
    field_data.MolecularWeight = &mweight_memview[0]
    field_data.ncells  = dims;

    cdef np.ndarray[np.float64_t, ndim=1] atol = np.ones( (NTOT) )*floor_value*reltol
    cdef np.ndarray[np.float64_t, ndim=1] rtol = np.array([reltol])

    if intermediate == 1:
        # allocate memory for the intermediate results
        H_1_int = np.zeros((dims, niter), "float64")
        H_2_int = np.zeros((dims, niter), "float64")
        de_int = np.zeros((dims, niter), "float64")
        ge_int = np.zeros((dims, niter), "float64")
        CoolingTime_int = np.zeros((dims, niter), "float64")
        Gamma_int       = np.zeros((dims, niter), "float64")
        Temperature_int = np.zeros((dims, niter), "float64")
        MolecularWeight_int = np.zeros((dims, niter), "float64")

        temp_int = np.zeros((dims, niter), "float64")
        result_int = np.zeros(niter, "float64")
        t_int = np.zeros(niter, "float64")
        dt_int = np.zeros(niter, "float64")

    j = 0
    for i in range(dims):
        input[j] = prev[j] = H_1_arr[i]
        j += 1
        input[j] = prev[j] = H_2_arr[i]
        j += 1
        input[j] = prev[j] = de_arr[i]
        j += 1
        input[j] = prev[j] = ge_arr[i]
        j += 1

    cdef double ttot = 0.0
    cdef int status
    # Allocate some temporary data
    # Now we manually evolve
    #ttot = dengo_evolve_simpleNetwork(tf, dt, input, rtol, atol, dims, data)
    cdef double *t_now = <double *> malloc( sizeof(double) )
    cdef double *dt_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *success_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *ttot_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *temp = <double *> malloc(sizeof(double) * niter)

    cdef int niter_cvodes = niter
    cdef double dt  = tf / float(niter)
    field_data.dengo_data_file = "simpleNetwork_tables.h5"

    field_data.reltol = reltol;
    field_data.floor_value = floor_value;

    for iter in range(niter):
        status = simpleNetwork_solve_chemistry( units, field_data, dt);
        j = 0;

        for i in range(dims):
            H_1_int[i, iter] = field_data.H_1_density[i]
            j += 1
            H_2_int[i, iter] = field_data.H_2_density[i]
            j += 1
            de_int[i, iter] = field_data.de_density[i]
            j += 1
            ge_int[i, iter] = field_data.ge_density[i]
            j += 1
            temp_int[ i, iter ] = field_data.temperature[i]

        if status == 0:
            result_int[iter] = 1
            ttot += dt
        elif status == 1:
            result_int[iter] = 0

        t_int[iter] = ttot
        dt_int[iter] = dt

        if status == 0:
            if iter % 100 == 0:
                print("Successful iteration[% 5i]: (%0.3e) %0.3e / %0.3e" % (iter, dt, ttot, tf))

            if adaptive_step: 
                 dt *= 1.1

            if tf - ttot < dt:
                dt = tf - ttot
        elif status == 1:
            dt /= 2.0;
            # copy_array(prev, input, NTOT)
            # Reset the scaling array to match the new values
            # copy_array(input, scale, NTOT)
            if dt < 1e-20 * tf:
                print("dt too small (%0.3e / %0.3e) so breaking" % (dt, tf))
                break
            continue
        if ttot >= tf: break

    free(input)
    free(prev)
    free(t_now)
    free(temp)
    free(dt_arr)
    free(ttot_arr)
    free(success_arr)
    free(field_data)
    free(units)

    print("End in %s iterations: %0.5e / %0.5e (%0.5e)" % (iter + 1, ttot, tf, tf - ttot))

    rv, rv_t = {}, {}
    H_1_arr = rv["H_1"] = np.zeros(dims, "float64")
    H_2_arr = rv["H_2"] = np.zeros(dims, "float64")
    de_arr = rv["de"] = np.zeros(dims, "float64")
    ge_arr = rv["ge"] = np.zeros(dims, "float64")
    if intermediate:
        rv_t["H_1"] = H_1_int[:niter]
        rv_t["H_2"] = H_2_int[:niter]
        rv_t["de"] = de_int[:niter]
        rv_t["ge"] = ge_int[:niter]
        rv_t["successful"] = result_int.astype("bool")
        rv_t['T'] = temp_int
        rv_t['t'] = t_int
        rv_t['dt'] = dt_int

    j = 0
    for i in range(dims):
        H_1_arr[i] = input[j] #* 1.00794
        j += 1
        H_2_arr[i] = input[j] #* 1.00794
        j += 1
        de_arr[i] = input[j] #* 1.0
        j += 1
        ge_arr[i] = input[j] #* 1.0
        j += 1
    return rv, rv_t

cdef copy_array(double *input, double *output, int dims):
    cdef int i
    for i in range(dims):
        output[i] = input[i]

cdef floor_values(double *input, int dims, double floor):
    cdef int i
    for i in range(dims):
        if input[i] < floor:
            input[i] = floor