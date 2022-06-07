
/* THIS FILE HAS BEEN AUTO-GENERATED.  DO NOT EDIT. */

/* This is C++ code to read HDF5 files for
   reaction rates, cooling rates, and initial
   conditions for the chemical network defined
   by the user.  In addition, this contains
   code for calculating temperature from the
   gas energy and computing the RHS and the
   Jacobian of the system of equations which
   will be fed into the solver.
*/


#include "simpleNetwork_solver.h"

///////////////////////////////////////////////////////////////////////////////
/////////// Setup the reaction, cooling rate data table ///////////////////////
///////////////////////////////////////////////////////////////////////////////
simpleNetwork_data *simpleNetwork_setup_data( const char *FileLocation, int *NumberOfFields, char ***FieldNames)
{

    //-----------------------------------------------------
    // Function : simpleNetwork_setup_data
    // Description: Initialize a data object that stores the reaction/ cooling rate data 
    //-----------------------------------------------------

    int i, n;
    
    simpleNetwork_data *data = (simpleNetwork_data *) malloc(sizeof(simpleNetwork_data));
    
    // point the module to look for simpleNetwork_tables.h5
    data->dengo_data_file = FileLocation;

    /* allocate space for the scale related pieces */

    // Number of cells to be solved in a batch 
    data->nstrip = MAX_NCELLS;
    /*initialize temperature so it wont crash*/
    for ( i = 0; i < MAX_NCELLS; i++ ){
        for( n = 0; n < NTHREADS; n++ ){
            data->Ts[n][i]    = 1000.0;
            data->logTs[n][i] = log(1000.0);
        }
    }

    /* Temperature-related pieces */
    data->bounds[0] = 1.0;
    data->bounds[1] = 100000000.0;
    data->nbins = 1024 - 1;
    data->dbin = (log(data->bounds[1]) - log(data->bounds[0])) / data->nbins;
    data->idbin = 1.0L / data->dbin;

    /* Redshift-related pieces */
    data->z_bounds[0] = 0.0;
    data->z_bounds[1] = 10.0;
    data->n_zbins = 0 - 1;
    data->d_zbin = (log(data->z_bounds[1] + 1.0) - log(data->z_bounds[0] + 1.0)) / data->n_zbins;
    data->id_zbin = 1.0L / data->d_zbin;
    
    simpleNetwork_read_rate_tables(data);
    //fprintf(stderr, "Successfully read in rate tables.\n");

    simpleNetwork_read_cooling_tables(data);
    //fprintf(stderr, "Successfully read in cooling rate tables.\n");
    
    simpleNetwork_read_gamma(data);
    //fprintf(stderr, "Successfully read in gamma tables. \n");

    if (FieldNames != NULL && NumberOfFields != NULL) {
        NumberOfFields[0] = 4;
        FieldNames[0] = new char*[4];
        i = 0;
        FieldNames[0][i++] = strdup("H_1");
        
        FieldNames[0][i++] = strdup("H_2");
        
        FieldNames[0][i++] = strdup("de");
        
        FieldNames[0][i++] = strdup("ge");
        
    }

    data->dengo_data_file = NULL;

    return data;

}


void simpleNetwork_read_rate_tables(simpleNetwork_data *data)
{
    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "simpleNetwork_tables.h5";   
    }

    hid_t file_id = H5Fopen( filedir , H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_double(file_id, "/k01", data->r_k01);
    H5LTread_dataset_double(file_id, "/k02", data->r_k02);
    
    H5Fclose(file_id);
}


void simpleNetwork_read_cooling_tables(simpleNetwork_data *data)
{

    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "simpleNetwork_tables.h5";   
    }
    hid_t file_id = H5Fopen( filedir , H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_double(file_id, "/reHII_reHII",
                            data->c_reHII_reHII);

    H5Fclose(file_id);
}

void simpleNetwork_read_gamma(simpleNetwork_data *data)
{

    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "simpleNetwork_tables.h5";   
    }
    
    hid_t file_id = H5Fopen( filedir , H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */

    H5Fclose(file_id);

}
 


/*
   This setup may be different than the user may anticipate, as a result
   of the lockstep timestep we use for a pencil beam through the grid.
   As such, it accepts the number of things to interpolate and makes
   assumptions about the sizes of the rates.
*/

/* This also requires no templating other than for the solver name...*/
void simpleNetwork_interpolate_rates(simpleNetwork_data *data,
                    int nstrip)
{
    int i, bin_id, zbin_id;
    double lb, t1, t2;
    double lbz, z1, z2, Tdef, zdef;
    int no_photo = 0;
    lb = log(data->bounds[0]);
    lbz = log(data->z_bounds[0] + 1.0);


    i = 0;
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    for ( i = 0; i < nstrip; i++ ){
        data->bin_id[threadID][i] = bin_id = (int) (data->idbin * (data->logTs[threadID][i] - lb));
        if (data->bin_id[threadID][i] <= 0) {
            data->bin_id[threadID][i] = 0;
        } else if (data->bin_id[threadID][i] >= data->nbins) {
            data->bin_id[threadID][i] = data->nbins - 1;
        }
        t1 = (lb + (bin_id    ) * data->dbin);
        t2 = (lb + (bin_id + 1) * data->dbin);
        data->Tdef[threadID][i] = (data->logTs[threadID][i] - t1)/(t2 - t1);
        data->dT[threadID][i] = (t2 - t1);
        /*fprintf(stderr, "INTERP: %d, bin_id = %d, dT = % 0.16g, T = % 0.16g, logT = % 0.16g\n",
                i, data->bin_id[i], data->dT[i], data->Ts[i],
                data->logTs[i]);*/
    
    if ((data->current_z >= data->z_bounds[0]) && (data->current_z < data->z_bounds[1])) {
        zbin_id = (int) (data->id_zbin * (log(data->current_z + 1.0) - lbz));
        if (zbin_id <= 0) {
            zbin_id = 0;
        } else if (zbin_id >= data->n_zbins) {
            zbin_id = data->n_zbins - 1;
        }
        z1 = (lbz + (zbin_id    ) * data->d_zbin);
        z2 = (lbz + (zbin_id + 1) * data->d_zbin);
        data->zdef = (log(data->current_z + 1.0) - z1)/(z2 - z1);
        data->dz = (exp(z2) - exp(z1)); //note: given this, we don't have to divide rate of change by z
    } else {
        no_photo = 1;
    }
    }

    zdef   = data->zdef;
    
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k01[threadID][i] = data->r_k01[bin_id] +
            Tdef * (data->r_k01[bin_id+1] - data->r_k01[bin_id]);
        data->drs_k01[threadID][i] = (data->r_k01[bin_id+1] - data->r_k01[bin_id]);
        data->drs_k01[threadID][i] /= data->dT[threadID][i];
        data->drs_k01[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k02[threadID][i] = data->r_k02[bin_id] +
            Tdef * (data->r_k02[bin_id+1] - data->r_k02[bin_id]);
        data->drs_k02[threadID][i] = (data->r_k02[bin_id+1] - data->r_k02[bin_id]);
        data->drs_k02[threadID][i] /= data->dT[threadID][i];
        data->drs_k02[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_reHII_reHII[threadID][i] = data->c_reHII_reHII[bin_id] +
            Tdef * (data->c_reHII_reHII[bin_id+1] - data->c_reHII_reHII[bin_id]);
        data->dcs_reHII_reHII[threadID][i] = (data->c_reHII_reHII[bin_id+1] - data->c_reHII_reHII[bin_id]);
        data->dcs_reHII_reHII[threadID][i] /= data->dT[threadID][i];
        data->dcs_reHII_reHII[threadID][i] *= data->invTs[threadID][i];
    }
    
    
}
 


void simpleNetwork_interpolate_gamma(simpleNetwork_data *data,
                    int i)
{   

    /*
     * find the bin_id for the given temperature 
     * update dT for i_th strip
     */

    int bin_id, zbin_id;
    double lb, t1, t2;
    double lbz, z1, z2;
    int no_photo = 0;
    lb = log(data->bounds[0]);
    lbz = log(data->z_bounds[0] + 1.0);
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    data->bin_id[threadID][i] = bin_id = (int) (data->idbin * (data->logTs[threadID][i] - lb));
    if (data->bin_id[threadID][i] <= 0) {
        data->bin_id[threadID][i] = 0;
    } else if (data->bin_id[threadID][i] >= data->nbins) {
        data->bin_id[threadID][i] = data->nbins - 1;
    }
    t1 = (lb + (bin_id    ) * data->dbin);
    t2 = (lb + (bin_id + 1) * data->dbin);
    data->Tdef[threadID][i] = (data->logTs[threadID][i] - t1)/(t2 - t1);
    data->dT[threadID][i] = (t2 - t1);

    
       
    }




///////////////////////////////////////////////////////////////////////////////
/////////////////// Main Evolution Routines            ////////////////////////
///////////////////////////////////////////////////////////////////////////////

int simpleNetwork_main(int argc, char** argv )
{
    //-----------------------------------------------------
    // Function : simpleNetwork_main
    // Description: this will look for initial condition files from the CLI, 
    //              evolve the ODE system to dtf specified from the CLI, (if not it's set to freefall time),
    //              and write the result to simpleNetwork_solution.h5 if output file name is not specified
    // Parameter:   argv[1]   : initial condition file name (in hdf5 format)
    //              argv[2]   : output file name (in hdf5 format)
    //              argv[3]   : desired final time reached by solver (in seconds)
    //
    // Note:        Units of initial conditions/ output 
    //              Baryons: mass density in a.m.u * cm^-3 ( 1mH = 1.00794 amu )
    //              de     : number density of electrons in cm^-3 
    //              ge     : internal energy per mass density of the cell (erg / g )
    //-----------------------------------------------------
    simpleNetwork_data *data = simpleNetwork_setup_data(NULL, NULL, NULL);

    /* Initial conditions */
    hid_t file_id;
    if (argc < 2){
    file_id = H5Fopen("simpleNetwork_initial_conditions.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {fprintf(stderr, "Failed to open "
        "simpleNetwork_initial_conditions.h5 so dying.\n");
        return(1);}
    } else {
        file_id = H5Fopen( argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {fprintf(stderr, "Failed to open  your initial_conditions file so dying.\n");
        return(1);}
       
            
    }

    /* Allocate the correct number of cells */
    hsize_t dims; /* We have flat versus number of species */

    /* Check gas energy to get the number of cells */
    fprintf(stderr, "Getting dimensionality from ge:\n");
    herr_t status = H5LTget_dataset_info(file_id, "/ge", &dims, NULL, NULL);
    if(status == -1) {
        fprintf(stderr, "Error opening initial conditions file.\n");
        return 1;
    }
    fprintf(stderr, "  ncells = % 3i\n", (int) dims);
    data->ncells = dims;

    int N = 4;

    double *atol, *rtol;
    atol = (double *) malloc(N * dims * sizeof(double));
    rtol = (double *) malloc(N * dims * sizeof(double));

    double *tics = (double *) malloc(dims * sizeof(double));
    double *ics = (double *) malloc(dims * N * sizeof(double));
    double *input = (double *) malloc(dims * N * sizeof(double));
    double *temp  = (double *) malloc(dims * sizeof(double) );

    unsigned int i = 0, j;
    
    fprintf(stderr, "Reading I.C. for /H_1\n");
    H5LTread_dataset_double(file_id, "/H_1", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "H_1[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /H_2\n");
    H5LTread_dataset_double(file_id, "/H_2", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "H_2[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /de\n");
    H5LTread_dataset_double(file_id, "/de", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "de[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /ge\n");
    H5LTread_dataset_double(file_id, "/ge", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "ge[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    
    double *density = (double *) malloc(dims *sizeof(double) );
    H5LTread_dataset_double(file_id, "/density", density);

    H5Fclose(file_id);
    
    // double dtf = 31557000000000.0;
    double dtf, t0;
    t0 = 2.992e15;
    dtf = t0 / sqrt(density[0]);
    
    // if the output time is specified,
    // it overrides the freefall time
    if (argc > 3){
        dtf = atof( argv[3] ); 
    }

    double dt = -1.0;
    double z = -1.0;
    for (i = 0; i < dims * N; i++) input[i] = ics[i];
    double ttot;
    int flag = dengo_evolve_simpleNetwork(dtf, dt, z, input, rtol, atol, dims, data, temp);
    if (flag > 0){
        fprintf(stderr, "solver failed, Time reached by the solver: %0.5g \n", dt);    
    }

    /* Write results to HDF5 file */

    if (argc < 3){
        file_id = H5Fcreate("simpleNetwork_solution.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    } else{
        file_id = H5Fcreate( argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    hsize_t dimsarr[1];
    dimsarr[0] = dims;
    i = 0;
    
    double H_1[dims];
    for (j = 0; j < dims; j++) {
        H_1[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /H_1\n");
    H5LTmake_dataset_double(file_id, "/H_1", 1, dimsarr, H_1);
    i++;
    
    double H_2[dims];
    for (j = 0; j < dims; j++) {
        H_2[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /H_2\n");
    H5LTmake_dataset_double(file_id, "/H_2", 1, dimsarr, H_2);
    i++;
    
    double de[dims];
    for (j = 0; j < dims; j++) {
        de[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /de\n");
    H5LTmake_dataset_double(file_id, "/de", 1, dimsarr, de);
    i++;
    
    double ge[dims];
    for (j = 0; j < dims; j++) {
        ge[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /ge\n");
    H5LTmake_dataset_double(file_id, "/ge", 1, dimsarr, ge);
    i++;
    


    H5LTmake_dataset_double(file_id, "/T", 1, dimsarr, temp);

    double time[1];
    time[0] = ttot;
    double timestep[1];
    timestep[0] = dt;
    H5LTset_attribute_double(file_id, "/", "time", time, 1); 
    H5LTset_attribute_double(file_id, "/", "timestep", timestep, 1);
    H5Fclose(file_id);
    
    free(temp);
    free(tics);
    free(ics);
    free(data);
    free(rtol);
    free(atol);
    free(input);
    free(density);

    return 0;
}
 



int dengo_evolve_simpleNetwork (double dtf, double &dt, double z, double *input,
            double *rtol, double *atol, unsigned long dims, simpleNetwork_data *data, double *temp_array ){

    //-----------------------------------------------------
    // Function     : dengo_evolve_simpleNetwork
    // Description  : Main ODE solver function in dengo
    
    // Parameter    :   dtf     : Desired time to be reached by the solver
    //                  dt      : Pointer to the actual time reached by the solver
    //                  z       : Current redshift
    //                  input   : Array to store the initial value of each species, 
    //                            it will be updated with the value at time dt
    //                  rtol    : relative tolerance required for convergenece in each internal CVODE timesteps
    //                  atol    : absolute tolerance required for convergence in each interanl CVODE timesteps
    //                  dims    : dimension of the input array, i.e. no. of species * no. of cells
    //                  data    : simpleNetwork_data object that relay the reaction/cooling rates, and normalizations 
    //                  temp_array: temperature of each cell by the end of the evolution
    //                           
    //-----------------------------------------------------
    unsigned long i, j;
    hid_t file_id;
    /* fprintf(stderr, "  ncells = % 3i\n", (int) dims); */
    int N = 4;

    data->reltol = rtol[0];
    double H_1;
    double H_2;
    double de;
    double ge;
    for (i = 0; i < dims; i++) {
        j = i * N;
        input[j] /= 1.00794; // H_1
        j++;
        input[j] /= 1.00794; // H_2
        j++;
        input[j] /= 1.0; // de
        j++;
        j++;
    }

    double floor_value = data->floor_value;
    for (i = 0; i < dims*N; i++) input[i] = fmax(floor_value, input[i]);
    
    // when partial equilibrium solver is used
    // the equilibrium abundance is passed through the temp_array
    double *equil_array = &temp_array[dims];

    // TODO:
    // Need to consider the cases where 
    // the incoming array assumes some equilibrium species,
    // - should first calculate temperature, then rates
    // - then equilibrium species
    // - then electron conservation
    // calculate_equilibrium_abundance(data, input, dims, 0, dims, equil_array);
    // ensure_electron_consistency(input, equil_array, dims, N);

    rhs_f f = calculate_rhs_simpleNetwork;

    #ifndef CVSPILS
    #ifdef  CVKLU
    jac_f jf = calculate_sparse_jacobian_simpleNetwork;
    #else
    jac_f jf = calculate_jacobian_simpleNetwork;
    #endif
    #endif
    
    #ifdef CVSPILS
    jac_f jf = calculate_JacTimesVec_simpleNetwork;
    #endif

    if (dt < 0) dt = dtf / 1e0;
    data->current_z = z;
    int niter = 0;
    int siter = 0;


    // Initialize a CVODE object, memory spaces
    // and attach rhs, jac to them
    int flag;
    double reltol = rtol[0];
    void *cvode_mem;
    int MAX_ITERATION = 10;
    double mh = 1.66054e-24;
    int nstrip = data->nstrip;

    // inputs are now grouped into a batch of nstrip
    // and send into the CVode solver
    //
    // ntimes: the number of times the CVode solver will be called
    // v_size: the size of the input vector for the solver
    //       : v_size = N * nstrip
    // v_size_res: leftover strip that doesn't fit into a v_size batch
    // N     : Number of species
    
    int v_size      = N * nstrip;
    int nstrip_res  = dims % nstrip;
    int v_size_res  = N * nstrip_res;
    unsigned long ntimes      = dims / nstrip;
    
    SUNLinearSolver LS;
    SUNMatrix A;
    N_Vector y_vec, abstol;

    y_vec = NULL;   
    LS = NULL;
    A  = NULL;
    
    double *yvec_ptr;
    double *atol_ptr;

    // these objects should be initialize once !!
    // in each separate thread!!
    double *ttot = (double *) malloc( (ntimes + 1)* sizeof(double));
    int d;
    int sum;
    int NSPARSE = 15; // no. of sparse jacobian compoents
    
    int threadID;

    // this routine is also called from the cythonized routine too
    // where MAX_NCELLS > input strip length
    // to avoid error message, we will catch it here
    if (ntimes > 0){

    #pragma omp parallel private (A, LS, cvode_mem, threadID, y_vec, abstol) num_threads(NTHREADS)
    {

    y_vec  = N_VNew_Serial(v_size);
    abstol = N_VNew_Serial(v_size); 
    yvec_ptr = N_VGetArrayPointer(y_vec);
    atol_ptr = N_VGetArrayPointer(abstol);
    
    // Need to be initialized before feeding into A, LS
    for (i = 0; i < v_size; i++) {
            yvec_ptr[i]   = 1.0;
            atol_ptr[i]   = reltol;
    }
    
    #ifdef CVKLU
    A         = SUNSparseMatrix      ( v_size, v_size, nstrip * NSPARSE, CSR_MAT );
    LS        = SUNKLU( y_vec, A);
    #else
    A         = SUNDenseMatrix      ( v_size, v_size );
    LS        = SUNDenseLinearSolver( y_vec, A);
    #endif
    
    cvode_mem = setup_cvode_solver  ( f, jf, v_size , data, LS, A, y_vec, reltol, abstol);
    
    // d: d th path going into the solver
    #pragma omp for private (sum, i, d, siter, threadID) schedule(static, 1)
    for (int d = 0; d < ntimes; d++){
        #ifdef _OPENMP
        threadID = omp_get_thread_num();
        #else
        threadID = 0;
        #endif
        ttot[d] = evolve_in_batches( cvode_mem, y_vec, abstol, reltol, input, v_size, d, d*v_size, MAX_ITERATION, dtf, data );

        // re-calculate temperature at the final output
        // updated in the data->Ts[threadID]
        simpleNetwork_calculate_temperature(data, &input[d*v_size], nstrip, N);
        
        // fprintf(stderr, "%d th strip from thread %d = %0.5g\n", d, threadID, ttot[d]);
        if (ttot[d] < dtf){
            fprintf(stderr, "FAILED FROM thread %d at ttot[%d] = %0.5g \n", threadID, d, ttot[d]);    
        } 

        for ( i = 0; i < nstrip; i ++){
            temp_array[ d * nstrip + i ] = data->Ts[threadID][i];
        }
    calculate_equilibrium_abundance(data, &input[d*v_size], nstrip, d, dims, equil_array);

    } // for d dims loop
    
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    #ifdef _OPENMP
    N_VDestroy(y_vec);
    N_VDestroy(abstol);
    #endif
    }
    } // if (ntimes > 0)

    if ( v_size_res > 0 ){
        #ifdef _OPENMP
        threadID = omp_get_thread_num();
        #else
        threadID = 0;
        #endif

        data->nstrip = nstrip_res;
        d = ntimes;
    
        y_vec  = N_VNew_Serial( v_size_res );
        abstol = N_VNew_Serial( v_size_res ); 

        yvec_ptr = N_VGetArrayPointer(y_vec);
        atol_ptr = N_VGetArrayPointer(abstol);

        for (i = 0; i < v_size_res; i++) {
            yvec_ptr[i]   = 1.0;
            atol_ptr[i]   = reltol;
        }
       
        #ifdef CVKLU
        A  = SUNSparseMatrix(v_size_res, v_size_res, nstrip_res * NSPARSE, CSR_MAT);
        LS = SUNKLU(y_vec,A); 
        #else
        A  = SUNDenseMatrix(v_size_res, v_size_res);
        LS = SUNDenseLinearSolver(y_vec, A);
        #endif
        
        cvode_mem = setup_cvode_solver( f, jf, v_size_res, data, LS, A, y_vec, reltol, abstol );
        ttot[d] = evolve_in_batches(cvode_mem, y_vec, abstol, reltol, input, v_size_res, d, d * v_size,  MAX_ITERATION, dtf, data);
        simpleNetwork_calculate_temperature(data, &input[d*v_size], nstrip_res, N);
        for ( i = 0; i < nstrip_res; i ++){
            temp_array[ d * nstrip + i ] = data->Ts[threadID][i];
        }
        
        
        calculate_equilibrium_abundance(data, &input[d*v_size], nstrip_res, d, dims, equil_array);
        //fprintf(stderr, "%d th strip = %0.5g\n", d, ttot[d]);
        
        // free all unused memory;
        CVodeFree(&cvode_mem);
        SUNLinSolFree(LS);
        SUNMatDestroy(A);
        N_VDestroy(y_vec);
        N_VDestroy(abstol);    
    } else{
      ttot[ntimes] = dtf;      
    }

    // inputs are in `number density`
    for (i = 0; i < dims; i++) {
      j = i * N;
      H_1 = input[j];
      input[j] *= 1.00794; // H_1
      j++;
      H_2 = input[j];
      input[j] *= 1.00794; // H_2
      j++;
      de = input[j];
      input[j] *= 1.0; // de
      j++;
      j++;
    }

    double dt_final = dtf;
    
    for (int d = 0; d < (ntimes + 1); d++){
        if (ttot[d] < dt_final) dt_final = ttot[d];    
    }

    // fprintf(stderr, "Fraction of completion (dt (%0.3g) / dtf (%0.3g)): %0.3g\n", dt, dt_final, dt_final/dtf);
    free(ttot);

    dt = dt_final;
    if (dt_final < dtf) return 1;
    return 0;

}
 



double evolve_in_batches( void * cvode_mem, N_Vector y_vec, N_Vector abstol,  
                          double reltol, double *input, int v_size, int d, int start_idx, 
                          int MAX_ITERATION, double dtf, simpleNetwork_data *data ){ 
    // Function   :     evolve_in_batches
    // Description:     this would evolve the ODE system in bataches of size v_size
    //                  and return the final time reached by the solver
    //
    // Parameter  :     cvode_mem   : CVODE memory object
    //                  y_vec           : Array to store/relay the data to the CVODE solver                  
    //                  abstol          : Array of absolute tolerance for each internal timestep
    //                  reltol          : Relative tolerance requested for each internal timestep
    //                  input           : Array to store the input number density/ energy of the species
    //                  v_size          : Size of the ODE system, i.e. No. of Species * No. of Strips/cells
    //                  d               : Batch count
    //                  start_idx       : Index to the first element from "input" for the current batch
    //                  MAX_ITERATION   : Maximum Call/retry to the CVODE solver
    //                  dtf             : Desired final output time
    //                  data            : simpleNetwork_data object that stores rate/ temperature/ scale data
    // Return    :      ttot            : final time reached by the solver

    int i, siter, flag;
    double dt, ttot;
    double y[v_size];
    int nstrip = data->nstrip;

    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    // access the array pointer of sundials vector object
    double *yvec_ptr = N_VGetArrayPointer(y_vec);
    double *atol_ptr = N_VGetArrayPointer(abstol);

    // by default: scale the input array
    #ifdef SCALE_INPUT
    double *scale = data->scale[threadID];
    double *inv_scale = data->inv_scale[threadID];
    for (i = 0; i < v_size; i++){ 
        scale[i]      = input[ start_idx + i];
        inv_scale[i]  = 1.0 / scale[i];
        yvec_ptr[i]   = 1.0;
        atol_ptr[i]   = reltol*reltol; //atol_ptr[i]*inv_scale[i];
    }
    setting_up_extra_variables(data, data->scale[threadID], nstrip );
    #else
    for (i = 0; i < v_size; i++){ 
        yvec_ptr[i]   = input[ start_idx + i];
        atol_ptr[i]   = reltol*reltol*yvec_ptr[i]; //atol_ptr[i]*inv_scale[i];
    }
    setting_up_extra_variables(data, yvec_ptr, nstrip );
    #endif

    // initialize a dt for the solver  
    dt = dtf;
    ttot = 0.0;
    siter = 0;
            
    while (ttot < dtf) { 
        // fprintf(stderr, "%d th strip: %d iterations, time: %0.5g\n", d, siter, ttot );    
        flag = cvode_solver( cvode_mem, y, v_size , &dt, data, y_vec, reltol, abstol);

        for (i = 0; i < v_size; i++) {
            if (y[i] < 0) {
            // this catches negatives values smaller than the abstol
            // and replaces them
            // need to check the convergence of this approach
            if (y[i] + atol_ptr[i] > 0 ){
                y[i] = atol_ptr[i]; 
            } else{
                fprintf(stderr, "negative \n");
            flag = 1;
                    break;
	    }
            }
        }
            
        if (flag < 1){
            // flag = 0 => success
            // we reset the scale of each component 
            // with the solution returned by the solver
            #ifdef SCALE_INPUT
            for (i = 0; i < v_size ; i++){
                yvec_ptr[i]   = 1.0;
        	//atol_ptr[i]   = fmin(reltol, atol_ptr[i]*inv_scale[i]);
                scale[i]     = y[i] * scale[i];
                inv_scale[i] = 1.0 /  scale[i];
            }
            #else
            for (i = 0; i < v_size ; i++){
                yvec_ptr[i]   = y[i];
            }
            #endif
            ttot += dt;
            dt = DMIN(dt * 2.0, dtf - ttot);
        } else{
            dt /= 2.0;
            dt = DMIN(dt , dtf - ttot);
        }
            if (siter == MAX_ITERATION) break;
            siter++;
    } // while loop for each strip
   
    // copy the results back to the input array
    // regardless the SCALE_INPUT
    // the returned array is still in number density
    // i.e. not scaled
    #ifdef SCALE_INPUT
    for (i = 0; i < v_size; i++){ 
        input[ start_idx + i] = scale[i];
    }
    #else
    for (i = 0; i < v_size; i++){ 
        input[ start_idx + i] = y[i];
    }
    #endif


    return ttot;
}


///////////////////////////////////////////////////////////////////////////////
//////////// Evaluate Temperature /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


int simpleNetwork_calculate_temperature(simpleNetwork_data *data,
                        double *input, int nstrip, int nchem)
{
    int i, j;
    double density, T, Tnew;
    double kb = 1.3806504e-16; // Boltzmann constant [erg/K]
    double mh = 1.66054e-24;
    double gamma = 5.e0/3.e0;
    double _gamma_m1 = 1.0 / (gamma - 1);

    double dge_dT;

    /* Calculate total density */
    double H_1;
    double H_2;
    double de;
    double ge;
    
    i = 0;
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else 
    int threadID = 0;
    #endif

    for ( i = 0; i < nstrip; i++ ){
        j = i * nchem;
        H_1 = input[j];
        j++;
        H_2 = input[j];
        j++;
        de = input[j];
        j++;
        ge = input[j];
        j++;
    
        /*
        */
	
        // TODO: pull the rates from simpleNetwork_data
        // these species usually contribute negligbly to the number density (?)
	// it is a little tricky here,
	// since these species depends on the temperature 
	// and the abundance of the rest of the species
	// BUT, without their abundance, temperature CANNOT be evaluated....
	// FOR NOW, a not entirely correct physically, 
	// BUT a not-too-bad surrogate is:
	// assume these species has negligible abundance....

        density = 1.0079400000000001*H_1 + 1.0079400000000001*H_2;
        
        // Requires iteration on the convergence of temperature
        // since the gammaH2 is not fixed 
        data->Ts[threadID][i] = density*ge*mh/(kb*(H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + de/(gamma - 1.0)));
        

        if (data->Ts[threadID][i] < data->bounds[0]) {
            data->Ts[threadID][i] = data->bounds[0];
        } else if (data->Ts[threadID][i] > data->bounds[1]) {
            data->Ts[threadID][i] = data->bounds[1];
        }
        data->logTs[threadID][i] = log(data->Ts[threadID][i]);
        data->invTs[threadID][i] = 1.0 / data->Ts[threadID][i];

        dge_dT = kb*(H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + de/(gamma - 1.0))/(density*mh);
        data->dTs_ge[threadID][i] = 1.0 / dge_dT;
    } // for i in nstrip loop
    return 0;
         
}
 



///////////////////////////////////////////////////////////////////////////////
////////// Ensure Conservation of mass, charge, species ///////////////////////
///////////////////////////////////////////////////////////////////////////////


void ensure_electron_consistency(double *input, double *equil_array, unsigned long nstrip, int nchem)
{
    // inputs are assumed to be in number density
    unsigned long i, j;

    /* Now we set up some temporaries */
    double H_1;
    double H_2;
    double de;
    double ge;

    double total_e = 0.0;
    int e_indx;

    for (i = 0; i<nstrip; i++) {
	    total_e = 0.0;
        j = i * nchem;
        H_1 = input[j];
        
        total_e += H_1 * 0.0;
        j++;
        H_2 = input[j];
        
        total_e += H_2 * 1.0;
        j++;
        de = input[j];
        e_indx = j;
        j++;
        ge = input[j];
        
        j++;
        input[e_indx] = total_e;
    }
}



///////////////////////////////////////////////////////////////////////////////
//////////////////// Auxillary Functions that estimate state variables ////////
///////////////////////////////////////////////////////////////////////////////
// Again these are exposed to the user
// and can be called basically with units, field_data objects
///////////////////////////////////////////////////////////////////////////////

int dengo_calculate_pressure_enzo( code_units* units, dengo_field_data *field_data){

    unsigned long int i, j, k, d, dims;
    int nchem = 4;
    
    // calculate total number of cells 
    int is, ie, js, je, ks, ke;
    // avoid ghost zones
    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    // number of cells that actually required calculations
    dims = ie - is + 1;
    dims*= je - js + 1;
    dims*= ke - ks + 1;

    const char * FileLocation = field_data->dengo_data_file;
    simpleNetwork_data *data = simpleNetwork_setup_data( FileLocation, NULL, NULL);

    // fills the `input` array with field data
    // note that input are still some kind of mass density, rho/ m_amu  
    double *input = (double *) malloc(dims*nchem*sizeof(double));
    flatten_dengo_field_data_enzo(units, field_data, input);

    // now we should call internal dengo function to calculate gamma
    // once, gamma is updated, we can use P = (gamma -1)* rho * u
    double *gamma_eff = (double *) malloc(dims*sizeof(double));
    double *pressure;
    double *density;
    double *thermal_energy;
    double *input_batch;
    double *gamma_eff_batch;

    int nstrip           = data->nstrip;
    unsigned long ntimes = dims / nstrip;
    int nstrip_res       = dims % nstrip;
    int threadID; 

    #pragma omp parallel for private (i, j ,d, input_batch, gamma_eff_batch, pressure, thermal_energy, density, threadID) num_threads(NTHREADS) schedule(static, 1)
    for ( d = 0; d < ntimes; d++ ){
        #ifdef _OPENMP
        threadID = omp_get_thread_num();
        #else
        threadID = 0;
        #endif

        input_batch        = &input[d*nstrip*nchem];
        gamma_eff_batch    = &gamma_eff[d*nstrip];
        pressure           = &field_data->Pressure[d*nstrip];
        thermal_energy     = &field_data->ge_density[d*nstrip];
        density            = &field_data->density[d*nstrip];

        setting_up_extra_variables( data, input_batch, nstrip );
        simpleNetwork_calculate_temperature(data,  input_batch , nstrip, nchem );
        dengo_calculate_gamma( gamma_eff_batch, data, input_batch, nstrip);
        for ( i = 0; i < nstrip; i++ ){
            pressure[i] = (gamma_eff_batch[i] - 1.0)* density[i]* thermal_energy[i];   
        }
    }

    if (nstrip_res > 0){
        input_batch        = &input[ntimes*nstrip*nchem];
        gamma_eff_batch    = &gamma_eff[ntimes*nstrip];
        pressure           = &field_data->Pressure[ntimes*nstrip];
        thermal_energy     = &field_data->ge_density[ntimes*nstrip];
        density            = &field_data->density[ntimes*nstrip];
        
        setting_up_extra_variables( data, input_batch, nstrip_res );
        simpleNetwork_calculate_temperature(data,  input_batch , nstrip_res , nchem );
        dengo_calculate_gamma( gamma_eff_batch, data, input_batch, nstrip_res);
        for ( i = 0; i < nstrip_res; i++ ){
            pressure[i] = (gamma_eff_batch[i] - 1.0)* density[i]* thermal_energy[i];   
        }
    }
    
    double pressure_unit = pow(units->length_units, 4)* units->density_units / pow(units->time_units,2); 
    for (unsigned long i = 0; i < dims; i++ ){
        field_data->Pressure[i] /= pressure_unit; 
    }

    free(input);
    free(gamma_eff);
    return 1;

}



int dengo_calculate_gamma_enzo(code_units* units, dengo_field_data *field_data){
    unsigned long i, j, k, d, dims;
    int nchem = 4;
    
    // calculate total number of cells 
    int is, ie, js, je, ks, ke;
    // avoid ghost zones
    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    // number of cells that actually required calculations
    dims = ie - is + 1;
    dims*= je - js + 1;
    dims*= ke - ks + 1;

    const char * FileLocation = field_data->dengo_data_file;
    simpleNetwork_data *data = simpleNetwork_setup_data( FileLocation, NULL, NULL);

    // fills the `input` array with field data
    // note that input are still some kind of mass density, rho/ m_amu  
    double *input = (double *) malloc(dims*nchem*sizeof(double));
    flatten_dengo_field_data_enzo(units, field_data, input);

    // now we should call internal dengo function to calculate gamma
    // once, gamma is updated, we can use P = (gamma -1)* rho * u
    double *input_batch;
    double *gamma_eff_batch;

    int nstrip     = data->nstrip;
    unsigned long ntimes = dims / nstrip;
    int nstrip_res = dims % nstrip;
    int threadID; 
    
    #pragma omp parallel for private (i, j ,d, input_batch, gamma_eff_batch, threadID) num_threads(NTHREADS) schedule(static, 1)
    for ( d = 0; d < ntimes; d++ ){
        #ifdef _OPENMP
        threadID = omp_get_thread_num();
        #else
        threadID = 0;
        #endif
        input_batch        = &input[d*nstrip*nchem];
        gamma_eff_batch    = &field_data->Gamma[d*nstrip];

        setting_up_extra_variables( data, input_batch, nstrip );
        simpleNetwork_calculate_temperature(data,  input_batch , nstrip, nchem );
        dengo_calculate_gamma( gamma_eff_batch, data, input_batch, nstrip);
    }

    if (nstrip_res > 0){
        input_batch        = &input[ntimes*nstrip*nchem];
        gamma_eff_batch    = &field_data->Gamma[ntimes*nstrip];
        
        setting_up_extra_variables( data, input_batch, nstrip_res );
        simpleNetwork_calculate_temperature(data,  input_batch , nstrip_res , nchem );
        dengo_calculate_gamma( gamma_eff_batch, data, input_batch, nstrip_res);
    }
    
    free(input);
    return 1;

}



int dengo_calculate_temperature_enzo(code_units* units, dengo_field_data *field_data){
    unsigned long int i, j, k, d, dims;
    int nchem = 4;
    
    // calculate total number of cells 
    int is, ie, js, je, ks, ke;
    // avoid ghost zones
    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    // number of cells that actually required calculations
    dims = ie - is + 1;
    dims*= je - js + 1;
    dims*= ke - ks + 1;

    const char * FileLocation = field_data->dengo_data_file;
    simpleNetwork_data *data = simpleNetwork_setup_data( FileLocation, NULL, NULL);

    // fills the `input` array with field data
    // note that input are still some kind of mass density, rho/ m_amu  
    double *input = (double *) malloc(dims*nchem*sizeof(double));
    flatten_dengo_field_data_enzo(units, field_data, input);
    
    // since we dont know if memory are allocated to the gamma pointer
    // lets just keep it local here
    double *gamma_eff = (double*) malloc(sizeof(double)*dims);
    
    double *input_batch;
    double *cooling_time_batch;
    double *gamma_eff_batch;
    int nstrip      = data->nstrip;
    unsigned long ntimes      = dims / nstrip;
    int nstrip_res            = dims % nstrip;

    int threadID; 
    /* Now we set up some temporaries */

    #pragma omp parallel for private (i, j ,d, input_batch, gamma_eff_batch, threadID) num_threads(NTHREADS) schedule(static, 1)
    for ( d = 0; d < ntimes; d++ ){
        #ifdef _OPENMP
        threadID = omp_get_thread_num();
        #else
        threadID = 0;
        #endif
        input_batch        = &input[d* nchem];
        gamma_eff_batch    = &gamma_eff[d*nstrip];

        setting_up_extra_variables( data, input_batch, nstrip );
        simpleNetwork_calculate_temperature(data,  input_batch , nstrip, nchem );
        dengo_calculate_gamma( gamma_eff_batch, data, input_batch, nstrip);
        for ( i = 0; i < nstrip; i++ ){
            field_data->temperature[d*nstrip + i] = data->Ts[threadID][i]; 
        }
    }

    if (nstrip_res > 0){
        input_batch        = &input[ntimes * nchem];
        gamma_eff_batch    = &gamma_eff[ntimes*nstrip];
        
        setting_up_extra_variables( data, input_batch, nstrip_res );
        simpleNetwork_calculate_temperature(data,  input_batch , nstrip_res , nchem );
        dengo_calculate_gamma( gamma_eff_batch, data, input_batch, nstrip_res);
       
        for ( i = 0; i < nstrip_res; i++ ){
            field_data->temperature[ntimes*nstrip + i] = data->Ts[threadID][i]; 
        }
    }
    
    free(input);
    free(data);

    return 1;

}





int dengo_estimate_cooling_time_enzo( code_units* units, dengo_field_data *field_data ){
    

    unsigned long int i, j, k, d, dims;
    int nchem = 4;
    
    // calculate total number of cells 
    int is, ie, js, je, ks, ke;
    // avoid ghost zones
    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    // number of cells that actually required calculations
    dims = ie - is + 1;
    dims*= je - js + 1;
    dims*= ke - ks + 1;

    const char * FileLocation = field_data->dengo_data_file;
    simpleNetwork_data *data = simpleNetwork_setup_data( FileLocation, NULL, NULL);

    // fills the `input` array with field data
    // note that input are still some kind of mass density, rho/ m_amu  
    double *input = (double *) malloc(dims*nchem*sizeof(double));
    flatten_dengo_field_data_enzo(units, field_data, input);

    int nstrip      = data->nstrip;
    unsigned long ntimes      = dims / nstrip;
    int nstrip_res            = dims % nstrip;
    
    double *input_batch;
    double *cooling_time_batch;

    // update the redshift
    double a = units->a_value * units->a_units;
    double z = 1./a - 1.;
    data->current_z = z;

    #pragma omp parallel for private (i, j ,d, input_batch, cooling_time_batch) num_threads(NTHREADS) schedule(static, 1) 
    for ( d = 0; d < ntimes; d++ ){
        input_batch        = &input[d* nstrip* nchem];
        cooling_time_batch = &field_data->CoolingTime[d * nstrip];
        simpleNetwork_calculate_cooling_timescale( cooling_time_batch, input_batch, nstrip, data);
    }

    if (nstrip_res > 0){
        input_batch        = &input[ntimes* nstrip * nchem];
        cooling_time_batch = &field_data->CoolingTime[ntimes*nstrip]; 
        simpleNetwork_calculate_cooling_timescale( cooling_time_batch, input_batch, nstrip_res, data );    
    }

    for (i = 0; i < dims; i++ ){
        field_data->CoolingTime[i] /= units->time_units; 
    }
    
    free(input);
    free(data);

    // in the world of Enzo
    //FAIL = 0
    return 1;
}



int dengo_estimate_cooling_time( code_units* units, dengo_field_data *field_data ){
    
    unsigned long int i, j, d, dims;
    int nchem     = 4;
    dims = field_data->ncells;
    
    const char * FileLocation = field_data->dengo_data_file;
    simpleNetwork_data *data = simpleNetwork_setup_data( FileLocation, NULL, NULL);
    
    int nstrip      = data->nstrip;

    unsigned long ntimes      = dims / nstrip;
    int nstrip_res            = dims % nstrip;
    
    // flatten field_data entry to a 1d input array
    double *input = (double *) malloc( sizeof(double) * nchem * dims );
    flatten_dengo_field_data( units, field_data, input );
    
    double *input_batch;
    double *cooling_time_batch;
    
    #pragma omp parallel for private (i, j ,d, input_batch, cooling_time_batch) num_threads(NTHREADS) schedule(static, 1) 
    for ( d = 0; d < ntimes; d++ ){
        input_batch        = &input[d* nstrip *nchem];
        cooling_time_batch = &field_data->CoolingTime[d * nstrip];
        simpleNetwork_calculate_cooling_timescale( cooling_time_batch, input_batch, nstrip, data);
    }

    if (nstrip_res > 0){
        input_batch        = &input[ntimes* nstrip* nchem];
        cooling_time_batch = &field_data->CoolingTime[ntimes*nstrip]; 
        simpleNetwork_calculate_cooling_timescale( cooling_time_batch, input_batch, nstrip_res, data );    
    }

    for (i = 0; i < dims; i++ ){
        field_data->CoolingTime[i] /= units->time_units; 
    }
    
    free(input);
    free(data);
}




int simpleNetwork_calculate_cooling_timescale( double *cooling_time, double *input, int nstrip, simpleNetwork_data *data){
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num(); 
    #else
    int threadID = 0;
    #endif

    unsigned long i, j, dims;
    int flag;
    int nchem = 4;
    /* Now we set up some temporaries */
    // make sure the input are in number density
    for (i = 0; i < nstrip; i++) {
        j = i * nchem;
        input[j] /= 1.00794; // H_1
        j++;
        input[j] /= 1.00794; // H_2
        j++;
        input[j] /= 1.0; // de
        j++;
        j++;
    }
    
    // calculate temperature and cooling rate
    setting_up_extra_variables(data, input, nstrip );
    flag = simpleNetwork_calculate_temperature(data,  input , nstrip, nchem );
    if (flag > 0){
        // check if the temperature failed to converged
        return -1;    
    }
    simpleNetwork_interpolate_rates(data, nstrip);
    double H_2;
    double H_1;
    double ge;
    double de;

    // this might be redundant
    // should only select relavant rates
    double *k01 = data->rs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *reHII_reHII = data->cs_reHII_reHII[threadID];
    
    
    
    
    double z;
    double T;
    double mdensity, inv_mdensity, dge_dt;

    for ( i = 0; i < nstrip; i++ ){
        
        T            = data->Ts[threadID][i];
        z            = data->current_z;
        mdensity     = data->mdensity[threadID][i];
        inv_mdensity = data->inv_mdensity[threadID][i];
        
        

        j = i * nchem;
        H_1 = input[j];
        j++;
        H_2 = input[j];
        j++;
        de = input[j];
        j++;
        ge = input[j];
        j++;
   	
        // obtain a quasis-equilibrium estimate

        //
        // Species: ge
        //
        dge_dt = -reHII_reHII[i]*H_2*de;
        dge_dt *= inv_mdensity;
        cooling_time[i] = fabs( ge / dge_dt);
    
    //fprintf(stderr, "----------------\n");
    }
}







int dengo_calculate_temperature( code_units *units, dengo_field_data *field_data){
    
    unsigned long int i, j, d, dims;
    int nchem     = 4;
    dims = field_data->ncells;
    
    const char * FileLocation = field_data->dengo_data_file;
    simpleNetwork_data *data = simpleNetwork_setup_data( FileLocation, NULL, NULL);
    
    int nstrip      = data->nstrip;
    unsigned long ntimes      = dims / nstrip;
    int nstrip_res            = dims % nstrip;
    
    // flatten field_data entry to a 1d input array
    double *input = (double *) malloc( sizeof(double) * nchem * dims );
    flatten_dengo_field_data( units, field_data, input );

    double *input_batch;
    double *cooling_time_batch;
    double *gamma_eff_batch;

    int threadID; 
    /* Now we set up some temporaries */

    #pragma omp parallel for private (i, j ,d, input_batch, gamma_eff_batch, threadID) num_threads(NTHREADS) schedule(static, 1)
    for ( d = 0; d < ntimes; d++ ){
        #ifdef _OPENMP
        threadID = omp_get_thread_num();
        #else
        threadID = 0;
        #endif
        input_batch        = &input[d* nchem];
        gamma_eff_batch    = &field_data->Gamma[d*nstrip];

        setting_up_extra_variables( data, input_batch, nstrip );
        simpleNetwork_calculate_temperature(data,  input_batch , nstrip, nchem );
        dengo_calculate_gamma( gamma_eff_batch, data, input_batch, nstrip);
        for ( i = 0; i < nstrip; i++ ){
            field_data->temperature[d*nstrip + i] = data->Ts[threadID][i]; 
        }
    }

    if (nstrip_res > 0){
        input_batch        = &input[ntimes * nchem];
        gamma_eff_batch    = &field_data->Gamma[ntimes*nstrip];
        
        setting_up_extra_variables( data, input_batch, nstrip_res );
        simpleNetwork_calculate_temperature(data,  input_batch , nstrip_res , nchem );
        dengo_calculate_gamma( gamma_eff_batch, data, input_batch, nstrip_res);
       
        for ( i = 0; i < nstrip_res; i++ ){
            field_data->temperature[ntimes*nstrip + i] = data->Ts[threadID][i]; 
        }
    }
    
    free(input);
    free(data);
}





int dengo_calculate_gamma( double* gamma_eff, simpleNetwork_data *data, double *input, int nstrip  ){
    unsigned long int i, j, d, dims;
    int nchem     = 4;
    double H_1;
    double H_2;
    double de;
    double ge;
    
    double gamma = 5.0/3.0;
    
    double n_gamma_m1, ndensity;
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    // both temperature and gamma are updated before calling this
    for (i = 0; i < nstrip; i++ ){
        ndensity = 0;
        j = i * nchem;
        H_1 = input[j] / 1.00794;
        ndensity += H_1;
        H_2 = input[j] / 1.00794;
        ndensity += H_2;
        de = input[j] / 1.0;
        ndensity += de;
        ge = input[j] / 1.0;
      
        n_gamma_m1 = H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + de/(gamma - 1.0);
        gamma_eff[i] = ndensity / n_gamma_m1 + 1.0; 
    }

}



int dengo_calculate_mean_molecular_weight( code_units *units, dengo_field_data *field_data ){
    unsigned long int i, j, d, dims;
    int nchem     = 4;
    dims = field_data->ncells;
    
    
    // flatten field_data entry to a 1d input array
    double *input = (double *) malloc( sizeof(double) * nchem * dims );
    flatten_dengo_field_data( units, field_data, input );
    double H_2;
    double H_1;
    double ge;
    double de;
    double ndensity, mdensity;

    for (i = 0; i < dims; i++ ){
        ndensity = 0.0;
        mdensity = 0.0;
        j = i * nchem;
        // Species: H_1 
        mdensity += input[j];
        ndensity += input[j] /1.00794;
        j++;
        // Species: H_2 
        mdensity += input[j];
        ndensity += input[j] /1.00794;
        j++;
        // Species: de
        ndensity += input[j] /1.0;
        j++;
        // Species: ge
        j++;
        field_data->MolecularWeight[i] = mdensity / ndensity;   
    }
    
    free(input);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////Solve Chemistry   ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// basically all the routines that are exposed to the user
// that can be called from dengo


int simpleNetwork_solve_chemistry( code_units *units, dengo_field_data *field_data, double dt ){
    // to use this, reltol and floor_value must be specified through 
    // dengo_field_data 
    unsigned long int i, j, d, dims;
    int N = 4;
    dims = field_data->ncells; // total number of strips to be evaluated
    
    const char * FileLocation = field_data->dengo_data_file;
    simpleNetwork_data *data = simpleNetwork_setup_data( FileLocation, NULL, NULL);

    double *input = (double *) malloc(dims * N * sizeof(double));
    double *temp  = (double *) malloc(dims * (0+1) * sizeof(double) );
    double *atol  = (double *) malloc(dims * N * sizeof(double));
    double *rtol  = (double *) malloc(sizeof(double));
    
    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;

    rtol[0]            = field_data->reltol;
    data->reltol       = field_data->reltol;
    double floor_value = field_data->floor_value;

    flatten_dengo_field_data(units, field_data, input);
    /*
    for (int i = 0; i < N; i++){
        fprintf(stderr, "input[%d] = %0.5g\n", i, input[i] );
    }
    */
    
    int flag;
    double z;
    double dtf;
    
    // convert code time to seconds 
    dt *= units->time_units;
    dtf = dt;
    
    // update the rate table location
    data->dengo_data_file = field_data->dengo_data_file; 

    flag = dengo_evolve_simpleNetwork(dtf, dt, z, input, rtol, atol, dims, data, temp);
    
    reshape_to_dengo_field_data(units, field_data, input);  
    
    for ( d = 0; d< dims; d++  ){
        field_data->temperature[d] = temp[d];
    }

    for ( d = 0; d< dims; d++  ){
    }

    //dengo_estimate_cooling_time( units, field_data );
    //dengo_calculate_temperature(units, field_data);
    //dengo_calculate_mean_molecular_weight( units, field_data );

    free(input);
    free(temp);
    free(data);
    free(atol);
    free(rtol);

    if (flag > 0) return 1;

    return 0;
}




int simpleNetwork_solve_chemistry_enzo( code_units *units, dengo_field_data *field_data, double dt ){
    //-------------------------------------------------------------------------
    // Function: simpleNetwork_solve_chemistry_enzo
    // Description: takes the same input as simpleNetwork_solve_chemistry
    //              BUT field_data needs to be specified with things like
    //              grid_start, grid_end, grid_dimension
    //              such that dengo can avoid evaluating ghost zones
    //-------------------------------------------------------------------------

    // to use this, reltol and floor_value must be specified through 
    // dengo_field_data 
    unsigned long int i, j, k, d, dims;
    int N = 4;
    
    // calculate total number of cells 
    int is, ie, js, je, ks, ke;
    // avoid ghost zones
    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    // number of cells that actually required calculations
    dims = ie - is + 1;
    dims*= je - js + 1;
    dims*= ke - ks + 1;

    const char * FileLocation = field_data->dengo_data_file;
    simpleNetwork_data *data = simpleNetwork_setup_data( FileLocation, NULL, NULL);

    double *input = (double *) malloc(dims * N * sizeof(double));
    double *temp  = (double *) malloc(dims * (0+1) * sizeof(double) );
    double *atol  = (double *) malloc(dims * N * sizeof(double));
    double *rtol  = (double *) malloc(sizeof(double));
    
// fills the `input` array with field data
    // note that input are still some kind of mass density, rho/ m_amu  
    flatten_dengo_field_data_enzo(units, field_data, input);

    // specify the relative tolerance and floor value
    rtol[0]            = field_data->reltol;
    data->reltol       = field_data->reltol;
    double floor_value = field_data->floor_value;
    
    // specifiy the redshift
    int flag;
    double a = units->a_value * units->a_units;
    double z = 1./a - 1.;
    double dtf;
    
    // convert code time to seconds 
    dt *= units->time_units;
    dtf = dt;
    
    // update the rate table location
    data->dengo_data_file = field_data->dengo_data_file; 
    data->floor_value = floor_value;
    
    // evolve chemistry in dengo
    flag = dengo_evolve_simpleNetwork(dtf, dt, z, input, rtol, atol, dims, data, temp);

    // fill results in `input` back to field data
    // in appropriate units
    reshape_to_dengo_field_data_enzo(units, field_data, input, temp);

    // free all pointers
    free(input);
    free(temp);
    free(data);
    free(atol);
    free(rtol);

    // in Enzo; 0 == FAIL
    if (flag > 0) return 0;
    // 1 == SUCCESS
    return 1;
}



int simpleNetwork_solve_chemistry_dt( code_units *units, dengo_field_data *field_data, double* reltol, double* abstol, double dt ){
    // TODO:
    // this is not only called from the python modules 
    // but this is dumb...
    // this is almost the same as simpleNetwork_solve_chemistry!!!
    // should be replaced later....

    unsigned long int i, j, d, dims;
    int N = 4;
    dims = field_data->ncells; // total number of strips to be evaluated

    // turned the field data into a long chain of 
    // 1-D array
    //
    // N: number of species
    // d: d th number of strip to be evalulated
    // i: index for each species
    //    should be in correct order and handled 
    //    by dengo templates
    // i.e.
    // input[d*N + i] = field_data->HI_density[]
    //
    
    const char * FileLocation = field_data->dengo_data_file;
    simpleNetwork_data *data = simpleNetwork_setup_data( FileLocation, NULL, NULL);

    double *input = (double *) malloc(dims * N * sizeof(double));
    double *temp  = (double *) malloc(dims * sizeof(double) );
    double *atol  = (double *) malloc(dims * N * sizeof(double));
    double *rtol  = (double *) malloc(sizeof(double));
    
    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;


    #pragma omp parallel for private (i, j ,d) num_threads(NTHREADS) schedule(static,1)
    for ( d = 0; d< dims; d++  ){
        j = d*N;
        // this should be the normalized 
        // by the input units later
        // atol = input * rtol;
        // which again should be set by dengo
        // input in *mass density per amu* 
        // and energy in the units of (erg / g)
        input[j]  = field_data->H_1_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->H_2_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->de_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->ge_density[d] ;
        input[j] *= UNIT_E_per_M;
        j++;
    }

    /*
    for (int i = 0; i < N; i++){
        fprintf(stderr, "input[%d] = %0.5g\n", i, input[i] );
    }
    */
    
    int flag;
    double z;
    double dtf;
    
    // convert code time to seconds 
    dt *= units->time_units;
    dtf = dt;
    
    // update the rate table location
    data->dengo_data_file = field_data->dengo_data_file; 

    flag = dengo_evolve_simpleNetwork(dtf, dt, z, input, reltol, abstol, dims, data, temp);
    
    #pragma omp parallel for private (i, j ,d) num_threads(NTHREADS) schedule(static, 1)
    for ( d = 0; d< dims; d++  ){
        j = d*N;
        // this should be the normalized 
        // by the input units later
        // atol = input * rtol;
        // which again should be set by dengo
        field_data->H_1_density[d] = input[j];
        field_data->H_1_density[d] /= units->density_units / m_amu;
        j++;
        field_data->H_2_density[d] = input[j];
        field_data->H_2_density[d] /= units->density_units / m_amu;
        j++;
        field_data->de_density[d] = input[j];
        field_data->de_density[d] /= units->density_units / m_amu;
        j++;
        field_data->ge_density[d] = input[j];
        field_data->ge_density[d] /= UNIT_E_per_M;
        j++;
    field_data->temperature[d] = temp[d];
    }
    
    dengo_estimate_cooling_time( units, field_data );
    //dengo_calculate_temperature(units, field_data);
    dengo_calculate_mean_molecular_weight( units, field_data );


    free(input);
    free(temp);
    free(data);
    free(atol);
    free(rtol);

    if (flag > 0) return 1;

    return 0;
}



int calculate_equilibrium_abundance( simpleNetwork_data *data, double *input, 
    int nstrip, unsigned long d, unsigned long dims, double *equil_array){

    //-----------------------------------------------------
    // Function     : calculate_equilibrium_abundance
    // Parameter    :   
    //                  input   : Array to store the initial value of each species, 
    //                            it will be updated with the value at time dt
    //                  data    : simpleNetwork_data object that relay the reaction/cooling rates, and normalizations 
    //                  equil_array: temperature of each cell by the end of the evolution
    //                  d          : batch count  
    //                  dims       : total number of cells to be solved
    //-----------------------------------------------------
    
    // update the rates given the updated temperature from simpleNetwork_data
    simpleNetwork_interpolate_rates(data, nstrip);
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    int nchem = 4;
    unsigned long i, j;
    double H_1;
    double H_2;
    double de;
    double ge;
    double *k01 = data->rs_k01[threadID];
    double *k02 = data->rs_k02[threadID];

    for ( i = 0; i < nstrip; i++ ){
        j = i * nchem;
        H_1 = input[j];
        j++;
        H_2 = input[j];
        j++;
        de = input[j];
        j++;
        ge = input[j];
        j++;
    }
    return 0;
}




/*

int calculate_JacTimesVec_simpleNetwork
            (N_Vector v, N_Vector Jv, realtype t,
             N_Vector y, N_Vector fy,
             void *user_data, N_Vector tmp)
{
    // TODO:
    // as of now this is utterly useless, 
    // cos it runs even slower than the actual dense linear solver ...
    // BUT! kept in mind that autodiff could easily replace this 
    // but some linear call to rhs/ f evauluations O(n) time
    // but not O(n^2) i think...
    
    // We iterate over all of the rates
    // Calcuate temperature first
    int nstrip = 1;
    int nchem = 4;
    simpleNetwork_data *data = (simpleNetwork_data*)user_data; 
    

    int i, j;
    j = 0;

    // change N_Vector back to an array 
    double y_arr[ 4 ];
    y_arr[0] = Ith(y , 1);
    y_arr[1] = Ith(y , 2);
    y_arr[2] = Ith(y , 3);
    y_arr[3] = Ith(y , 4);
    // Abundances are scaled in the calculate temperature module
    int flag;
    flag = simpleNetwork_calculate_temperature(data, y_arr, nstrip, nchem);
    if (flag > 0){
        return 1;    
    }
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    simpleNetwork_interpolate_rates(data, nstrip);
    
    // Now We set up some temporaries
    double *Tge = data->dTs_ge[threadID];
    
    // Define the reaction rates
    double *k01 = data->rs_k01[threadID];
    double *rk01 = data->drs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *rk02 = data->drs_k02[threadID];
    double *reHII_reHII = data->cs_reHII_reHII[threadID];
    double *rreHII_reHII = data->dcs_reHII_reHII[threadID];
    
    

    double scale;
    // Define the species
    double H_1, v0;
    double H_2, v1;
    double de, v2;
    double ge, v3;
    scale = data->scale[threadID][0];
    v0 = Ith( v, 1 );
    v0 *= scale;
    scale = data->scale[threadID][1];
    v1 = Ith( v, 2 );
    v1 *= scale;
    scale = data->scale[threadID][2];
    v2 = Ith( v, 3 );
    v2 *= scale;
    scale = data->scale[threadID][3];
    v3 = Ith( v, 4 );
    v3 *= scale;

    double z;
    double T;

    double mh = 1.67e-24;
    double mdensity;
    
    int jj;
    jj = 0;

    j = i*nchem;
    mdensity = 0.0;
    z = data->current_z;
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    H_1 = Ith( y, 1  )*scale;
    
    mdensity += H_1;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    H_2 = Ith( y, 2  )*scale;
    
    mdensity += H_2;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    de = Ith( y, 3  )*scale;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    ge = Ith( y, 4  )*scale;
    
    j++;
    


    mdensity *= mh;
        
    j = 0;
    //
    // Species: H_1
    //
    Ith(Jv, 1 ) = -k01[i]*de*v0 + k02[i]*de*v1 + Tge[i]*v3*(-H_1*de*rk01[i] + H_2*de*rk02[i]) + v2*(-k01[i]*H_1 + k02[i]*H_2);

    scale = data->scale[threadID][0];
    Ith(Jv, 1) /= scale;

    
    //
    // Species: H_2
    //
    Ith(Jv, 2 ) = k01[i]*de*v0 - k02[i]*de*v1 + Tge[i]*v3*(H_1*de*rk01[i] - H_2*de*rk02[i]) + v2*(k01[i]*H_1 - k02[i]*H_2);

    scale = data->scale[threadID][1];
    Ith(Jv, 2) /= scale;

    
    //
    // Species: de
    //
    Ith(Jv, 3 ) = k01[i]*de*v0 - k02[i]*de*v1 + Tge[i]*v3*(H_1*de*rk01[i] - H_2*de*rk02[i]) + v2*(k01[i]*H_1 - k02[i]*H_2);

    scale = data->scale[threadID][2];
    Ith(Jv, 3) /= scale;

    
    //
    // Species: ge
    //
    Ith(Jv, 4 ) = -reHII_reHII[i]*H_2*v2/mdensity - reHII_reHII[i]*de*v1/mdensity - H_2*Tge[i]*de*rreHII_reHII[i]*v3/mdensity;

    scale = data->scale[threadID][3];
    Ith(Jv, 4) /= scale;

    
    Ith(Jv, 4) /= mdensity;
    
    return 0;
}

*/


///////////////////////////////////////////////////////////////////////////////
////////////////// RHS functions and Jacobian Functions ///////////////////////
///////////////////////////////////////////////////////////////////////////////

int calculate_rhs_simpleNetwork(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    simpleNetwork_data *data = (simpleNetwork_data* ) user_data;
    int i, j;

    int nchem = 4;
    int nstrip = data->nstrip;
    
    /* change N_Vector back to an array */
    double y_arr[ nchem * nstrip ];
    double H_1;
    double H_2;
    double de;
    double ge;
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif
   
    double *yvec_ptr = N_VGetArrayPointer(y);
    double *ydot_ptr = N_VGetArrayPointer(ydot);

    #ifdef SCALE_INPUT
    double *scale     = data->scale[threadID];
    double *inv_scale = data->inv_scale[threadID];   
    for ( i = 0; i < nstrip*nchem; i++ ){
        y_arr[i] = yvec_ptr[i]*scale[i];
    }
    #else
    for ( i = 0; i < nstrip*nchem; i++ ){
        y_arr[i] = yvec_ptr[i];
    }
    #endif

    int flag;
    flag = simpleNetwork_calculate_temperature(data, y_arr , nstrip, nchem );
    if (flag > 0){
        // check if the temperature failed to converged
        return -1;    
    }
    simpleNetwork_interpolate_rates(data, nstrip);

    /* Now we set up some temporaries */
    double *k01 = data->rs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *reHII_reHII = data->cs_reHII_reHII[threadID];

    double z;
    double T;

    double mh = 1.66054e-24;
    double mdensity, inv_mdensity;

   
    for ( i = 0; i < nstrip; i++ ){
        
        T            = data->Ts[threadID][i];
        z            = data->current_z;
        mdensity     = data->mdensity[threadID][i];
        inv_mdensity = data->inv_mdensity[threadID][i];

        j = i * nchem;
        H_1 = y_arr[j];
        j++;
        H_2 = y_arr[j];
        j++;
        de = y_arr[j];
        j++;
        ge = y_arr[j];
        j++;
    
        j = i * nchem;
        // Species: H_1
        ydot_ptr[j] = -k01[i]*H_1*de + k02[i]*H_2*de;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "H_1: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: H_2
        ydot_ptr[j] = k01[i]*H_1*de - k02[i]*H_2*de;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "H_2: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: de
        ydot_ptr[j] = k01[i]*H_1*de - k02[i]*H_2*de;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "de: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: ge
        ydot_ptr[j] = -reHII_reHII[i]*H_2*de;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        ydot_ptr[j] *= inv_mdensity;
        
        //fprintf(stderr, "ge: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
    
    //fprintf(stderr, "----------------\n");
    }
    return 0;
    }



int calculate_jacobian_simpleNetwork( realtype t,
                                        N_Vector y, N_Vector fy,
                                        SUNMatrix J, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3)
{
    /* We iterate over all of the rates */
    /* Calcuate temperature first */
    

    simpleNetwork_data *data = (simpleNetwork_data*)user_data; 
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    int nchem = 4;
    int nstrip = data->nstrip;
    int i, j;

    /* change N_Vector back to an array */
    double *yvec_ptr = N_VGetArrayPointer(y);
    double y_arr[ 4 * nstrip ];
 
    #ifdef SCALE_INPUT
    double *scale     = data->scale[threadID];
    double *inv_scale = data->inv_scale[threadID];   
    double scale1, inv_scale2;
    for ( i = 0; i < nstrip*nchem; i++ ){
        y_arr[i] = yvec_ptr[i]*scale[i];
    }
    #else
    for ( i = 0; i < nstrip*nchem; i++ ){
        y_arr[i] = yvec_ptr[i];
    }
    #endif


    /*
    int flag;
    flag = simpleNetwork_calculate_temperature(data, y_arr, nstrip, nchem );
    if (flag > 0){
        // check if the temperature failed to converged
        return -1;    
    }
    simpleNetwork_interpolate_rates(data, nstrip);
    */

    // simpleNetwork_calculate_temperature(data, y_arr, nstrip, nchem);
    // simpleNetwork_interpolate_rates(data, nstrip);

    /* Now We set up some temporaries */
    double *Tge = data->dTs_ge[threadID];
    double *k01 = data->rs_k01[threadID];
    double *rk01= data->drs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *rk02= data->drs_k02[threadID];
    double *reHII_reHII = data->cs_reHII_reHII[threadID];
    double *rreHII_reHII = data->dcs_reHII_reHII[threadID];
    double H_1;
    double H_2;
    double de;
    double ge;
    double z;
    double T;

    double mh = 1.66054e-24;
    double mdensity, inv_mdensity;
    
    
  

    j = 0;
    mdensity = 0.0;
    z = data->current_z;

    for ( i = 0; i < nstrip; i++ ){
        j = i * nchem;
        H_1 = y_arr[j];
        j++;
        H_2 = y_arr[j];
        j++;
        de = y_arr[j];
        j++;
        ge = y_arr[j];
        j++;
	
        mdensity = data->mdensity[threadID][i];
        inv_mdensity = 1.0 / mdensity;
        
        

        j = i * nchem;
        //
        // Species: H_1
        //
        
        
        // H_1 by H_1
        
        SM_ELEMENT_D( J, j + 0, j + 0 ) = -k01[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 0, j + 0) *= inv_scale2*scale1;
        #endif
        
        // H_1 by H_2
        
        SM_ELEMENT_D( J, j + 0, j + 1 ) = k02[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 0, j + 1) *= inv_scale2*scale1;
        #endif
        
        // H_1 by de
        
        SM_ELEMENT_D( J, j + 0, j + 2 ) = -k01[i]*H_1 + k02[i]*H_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 0, j + 2) *= inv_scale2*scale1;
        #endif
        
        // H_1 by ge
        
        SM_ELEMENT_D( J, j + 0, j + 3 ) = -H_1*de*rk01[i] + H_2*de*rk02[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 0, j + 3) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 0, j + 3 ) *= Tge[i];
        //
        // Species: H_2
        //
        
        
        // H_2 by H_1
        
        SM_ELEMENT_D( J, j + 1, j + 0 ) = k01[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 1, j + 0) *= inv_scale2*scale1;
        #endif
        
        // H_2 by H_2
        
        SM_ELEMENT_D( J, j + 1, j + 1 ) = -k02[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 1, j + 1) *= inv_scale2*scale1;
        #endif
        
        // H_2 by de
        
        SM_ELEMENT_D( J, j + 1, j + 2 ) = k01[i]*H_1 - k02[i]*H_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 1, j + 2) *= inv_scale2*scale1;
        #endif
        
        // H_2 by ge
        
        SM_ELEMENT_D( J, j + 1, j + 3 ) = H_1*de*rk01[i] - H_2*de*rk02[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 1, j + 3) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 1, j + 3 ) *= Tge[i];
        //
        // Species: de
        //
        
        
        // de by H_1
        
        SM_ELEMENT_D( J, j + 2, j + 0 ) = k01[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 2, j + 0) *= inv_scale2*scale1;
        #endif
        
        // de by H_2
        
        SM_ELEMENT_D( J, j + 2, j + 1 ) = -k02[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 2, j + 1) *= inv_scale2*scale1;
        #endif
        
        // de by de
        
        SM_ELEMENT_D( J, j + 2, j + 2 ) = k01[i]*H_1 - k02[i]*H_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 2, j + 2) *= inv_scale2*scale1;
        #endif
        
        // de by ge
        
        SM_ELEMENT_D( J, j + 2, j + 3 ) = H_1*de*rk01[i] - H_2*de*rk02[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 2, j + 3) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 2, j + 3 ) *= Tge[i];
        //
        // Species: ge
        //
        
        
        // ge by H_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 3, j + 0 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 3, j + 0) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 3, j + 0 ) *= inv_mdensity;
        
        // ge by H_2
        
        SM_ELEMENT_D( J, j + 3, j + 1 ) = -reHII_reHII[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 3, j + 1) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 3, j + 1 ) *= inv_mdensity;
        
        // ge by de
        
        SM_ELEMENT_D( J, j + 3, j + 2 ) = -reHII_reHII[i]*H_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 3, j + 2) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 3, j + 2 ) *= inv_mdensity;
        
        // ge by ge
        
        SM_ELEMENT_D( J, j + 3, j + 3 ) = -H_2*de*rreHII_reHII[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 3, j + 3) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 3, j + 3 ) *= inv_mdensity;
        SM_ELEMENT_D( J, j + 3, j + 3 ) *= Tge[i];
    }
    return 0;
}



#ifdef CVKLU
int calculate_sparse_jacobian_simpleNetwork( realtype t,
                                        N_Vector y, N_Vector fy,
                                        SUNMatrix J, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3)
{
    /* We iterate over all of the rates */
    /* Calcuate temperature first */
    

    simpleNetwork_data *data = (simpleNetwork_data*)user_data; 
    
    int nchem = 4;
    int nstrip = data->nstrip;
    int i, j;
    int NSPARSE = 15;
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif
    
    /* change N_Vector back to an array */
    double y_arr[ 4 * nstrip ];
    double *scale     = data->scale[threadID];
    double *inv_scale = data->inv_scale[threadID];

    //TODO: Here we assumed during the evaluation of jacobian
    // temperature is approximately constant, 
    // i.e. close enough to the point evaluation of f(y)
    // such that the rates and temperature need not be interpolated or evalulated 
    // again during the jacobian evaluation.
    // We havent really fully explored the effect of this `assumption`...
    // But it definitely boost the performance 
    
    /*
    int flag;
    flag = simpleNetwork_calculate_temperature(data, y_arr , nstrip, nchem );
    if (flag > 0){
        // check if the temperature failed to converged
        return -1;    
    }
    simpleNetwork_interpolate_rates(data, nstrip);
    */

    // simpleNetwork_calculate_temperature(data, y_arr, nstrip, nchem);
    // simpleNetwork_interpolate_rates(data, nstrip);

    /* Now We set up some temporaries */
    // CSR is what we choose
    sunindextype *rowptrs = SUNSparseMatrix_IndexPointers(J);
    sunindextype *colvals = SUNSparseMatrix_IndexValues(J);
    realtype *matrix_data = SUNSparseMatrix_Data(J);
    
    SUNMatZero(J);
   
    double *Tge = data->dTs_ge[threadID];
    double *k01 = data->rs_k01[threadID];
    double *rk01= data->drs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *rk02= data->drs_k02[threadID];
    double *reHII_reHII = data->cs_reHII_reHII[threadID];
    double *rreHII_reHII = data->dcs_reHII_reHII[threadID];
    double H_1;
    double H_2;
    double de;
    double ge;
    double z;
    double T;

    double mh = 1.66054e-24;
    double mdensity, inv_mdensity;
    
    double scale2, inv_scale1;

    j = 0;
    mdensity = 0.0;
    z = data->current_z;
   
    int k = 0;
    
    double *yvec_ptr = N_VGetArrayPointer(y);

    for ( i = 0; i < nstrip; i++ ){

        #ifdef SCALE_INPUT
        j = i * nchem;
        H_1 = yvec_ptr[j]*scale[j];
        j++;
        
        H_2 = yvec_ptr[j]*scale[j];
        j++;
        
        de = yvec_ptr[j]*scale[j];
        j++;
        
        ge = yvec_ptr[j]*scale[j];
        j++;
        
        #else
        j = i * nchem;
        H_1 = yvec_ptr[j];
        j++;
        
        H_2 = yvec_ptr[j];
        j++;
        
        de = yvec_ptr[j];
        j++;
        
        ge = yvec_ptr[j];
        j++;
        
        #endif

        mdensity = data->mdensity[threadID][i];
        inv_mdensity = 1.0 / mdensity; 
        
        
        

        j = i * NSPARSE;
        // H_1 by H_1
        colvals[j + 0] = i * nchem + 0 ;
        matrix_data[ j + 0 ] = -k01[i]*de;

        
        // H_1 by H_2
        colvals[j + 1] = i * nchem + 1 ;
        matrix_data[ j + 1 ] = k02[i]*de;

        
        // H_1 by de
        colvals[j + 2] = i * nchem + 2 ;
        matrix_data[ j + 2 ] = -k01[i]*H_1 + k02[i]*H_2;

        
        // H_1 by ge
        colvals[j + 3] = i * nchem + 3 ;
        matrix_data[ j + 3 ] = -H_1*de*rk01[i] + H_2*de*rk02[i];

        
        matrix_data[ j + 3] *= Tge[i];
        // H_2 by H_1
        colvals[j + 4] = i * nchem + 0 ;
        matrix_data[ j + 4 ] = k01[i]*de;

        
        // H_2 by H_2
        colvals[j + 5] = i * nchem + 1 ;
        matrix_data[ j + 5 ] = -k02[i]*de;

        
        // H_2 by de
        colvals[j + 6] = i * nchem + 2 ;
        matrix_data[ j + 6 ] = k01[i]*H_1 - k02[i]*H_2;

        
        // H_2 by ge
        colvals[j + 7] = i * nchem + 3 ;
        matrix_data[ j + 7 ] = H_1*de*rk01[i] - H_2*de*rk02[i];

        
        matrix_data[ j + 7] *= Tge[i];
        // de by H_1
        colvals[j + 8] = i * nchem + 0 ;
        matrix_data[ j + 8 ] = k01[i]*de;

        
        // de by H_2
        colvals[j + 9] = i * nchem + 1 ;
        matrix_data[ j + 9 ] = -k02[i]*de;

        
        // de by de
        colvals[j + 10] = i * nchem + 2 ;
        matrix_data[ j + 10 ] = k01[i]*H_1 - k02[i]*H_2;

        
        // de by ge
        colvals[j + 11] = i * nchem + 3 ;
        matrix_data[ j + 11 ] = H_1*de*rk01[i] - H_2*de*rk02[i];

        
        matrix_data[ j + 11] *= Tge[i];
        // ge by H_2
        colvals[j + 12] = i * nchem + 1 ;
        matrix_data[ j + 12 ] = -reHII_reHII[i]*de;

        
        matrix_data[j + 12] *= inv_mdensity;
        // ge by de
        colvals[j + 13] = i * nchem + 2 ;
        matrix_data[ j + 13 ] = -reHII_reHII[i]*H_2;

        
        matrix_data[j + 13] *= inv_mdensity;
        // ge by ge
        colvals[j + 14] = i * nchem + 3 ;
        matrix_data[ j + 14 ] = -H_2*de*rreHII_reHII[i];

        
        matrix_data[j + 14] *= inv_mdensity;
        matrix_data[ j + 14] *= Tge[i];
        rowptrs[ i * nchem +  0] = i * NSPARSE + 0;
        rowptrs[ i * nchem +  1] = i * NSPARSE + 4;
        rowptrs[ i * nchem +  2] = i * NSPARSE + 8;
        rowptrs[ i * nchem +  3] = i * NSPARSE + 12;
       
        #ifdef SCALE_INPUT
        j = i * nchem;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 0 ];
        matrix_data[ i * NSPARSE + 0]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 1 ];
        matrix_data[ i * NSPARSE + 1]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 2]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 3]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 0 ];
        matrix_data[ i * NSPARSE + 4]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 1 ];
        matrix_data[ i * NSPARSE + 5]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 6]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 7]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 0 ];
        matrix_data[ i * NSPARSE + 8]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 1 ];
        matrix_data[ i * NSPARSE + 9]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 10]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 11]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 1 ];
        matrix_data[ i * NSPARSE + 12]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 13]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 14]  *= inv_scale1*scale2;
        #endif

    }

    rowptrs[ i * nchem ] = i * NSPARSE ;
    return 0;
}

#endif




void setting_up_extra_variables( simpleNetwork_data * data, double * input, int nstrip ){
    //-------------------------------------------------------------------------    
    // Function: setting_up_extra_variables
    // Desciption: calculating variables that are independent on the state with time
    //             to avoid repeated evaluation. Examples here are h2 optical depth
    //             and cie_optical depth. Functions that depends only on density
    //-------------------------------------------------------------------------    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    int i, j;
    double mh = 1.66054e-24;
    double mdensity;
    // TODO: maybe this should be filled out by Dengo as well
    for ( i = 0; i < nstrip; i++){
        data->mdensity[threadID][i] = 0;
        j = i * 4;
        // species: H_1
        data->mdensity[threadID][i] += input[j] * 1.00794;
        j++;
        // species: H_2
        data->mdensity[threadID][i] += input[j] * 1.00794;
        j++;
        j++;
        j++;
        // TODO: update temperature and rates to obtain abundances of equilibrium states
	// for more detail: go to the _calculate_temperature 
        data->mdensity[threadID][i] *= mh;
        data->inv_mdensity[threadID][i] = 1.0 / data->mdensity[threadID][i];
    }
}



///////////////////////////////////////////////////////////////////////////////
/////////////////// Sturcture Incoming Data ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Reshape/ Flatten data to match the data shape/ type 
// required by the solver in dengo


int flatten_dengo_field_data(code_units *units, dengo_field_data *field_data, double *input){

    //-----------------------------------------------------
    // Function     :   flatten_dengo_field_data 
    // Parameter    :   
    //                  code_units: units from the incoming field_data 
    //                  field_data: dengo_field_data class that contains pointer to species array 
    //                  input     : 1D array that flattens the field_data, 
    //                              i.e. 
    //                              s = species, s0 = 0th species, with d dimensions
    //                              [s0, s1, s2..., , s0_1, s1_1, ... s0_d, s1_d, ...] 
    //-----------------------------------------------------
    //
    unsigned long d, dims, i, j;
    dims = field_data->ncells;

    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;
    double dom = units->density_units/m_amu;
    int N = 4;

    #pragma omp parallel for private (i, j ,d) num_threads(NTHREADS) schedule(static,1)
    for ( d = 0; d< dims; d++  ){
        j = d*N;
        // this should be the normalized 
        // by the input units later
        // atol = input * rtol;
        // which again should be set by dengo
        // input in *mass density per amu* 
        // and energy in the units of (erg / g)
        input[j]  = field_data->H_1_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->H_2_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->de_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->ge_density[d] ;
        input[j] *= UNIT_E_per_M;
        j++;
    }
}



int reshape_to_dengo_field_data( code_units* units, dengo_field_data *field_data, double* input ){
    //------------------------------------------------------------------------------------
    // Function   :     reshape_to_dengo_field_data
    // Description:     reshape the 1d output array from solver to a dengo_field_data object  
    //                  and covert them to code units 
    //                  i.e. ge_density in erg /g
    //                       H_1_density in g / cm^-3 / amu (mass density per amu)
    // Parameter  :     code_units
    //                  dengo_field_data
    //                  input
    //------------------------------------------------------------------------------------

    unsigned long int i, j, d, dims;
    int N = 4;
    dims = field_data->ncells; // total number of strips to be evaluated
    

    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;
    double dom = units->density_units/ m_amu;
   
    #pragma omp parallel for private (i, j ,d) num_threads(NTHREADS) schedule(static, 1)
    for ( d = 0; d< dims; d++  ){
        j = d*N;
        field_data->H_1_density[d] = input[j];
        field_data->H_1_density[d] /= dom;
        j++;
        field_data->H_2_density[d] = input[j];
        field_data->H_2_density[d] /= dom;
        j++;
        field_data->de_density[d] = input[j];
        field_data->de_density[d] /= dom;
        j++;
        field_data->ge_density[d] = input[j];
        field_data->ge_density[d] /= UNIT_E_per_M;
        j++;
    }
   
    return 0;
}



// and a enzo version
//

int flatten_dengo_field_data_enzo(code_units *units, dengo_field_data *field_data, double *input){

    //-----------------------------------------------------
    // Function     :   flatten_dengo_field_data_enzo
    // Description  :   To read in data from Enzo pointers 
    // Parameter    :   
    //                  code_units: units from the incoming field_data 
    //                  field_data: dengo_field_data class that contains pointer to species array 
    //                  input     : 1D array that flattens the field_data, 
    //                              i.e. 
    //                              s = species, s0 = 0th species, with d dimensions
    //                              [s0, s1, s2..., , s0_1, s1_1, ... s0_d, s1_d, ...] 
    //                              abundances in units of mass density / m_amu
    //                              m_amu is in atomic mass units                           
    //-----------------------------------------------------
    //
    
    int is, ie, js, je, ks, ke;
    int i, j, k, N;
    int ni, nj, nk, idim, jdim, kdim;
    unsigned long dims, ccount, c, idx;

    N = 4;
    // avoid ghost zones
    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    idim = field_data->grid_dimension[0];
    jdim = field_data->grid_dimension[1];
    kdim = field_data->grid_dimension[2];

    // number of cells that actually required calculations
    ni = ie - is + 1;
    nj = je - js + 1;
    nk = ke - ks + 1;
    dims = ni*nj*nk;
    field_data->ncells = dims;

    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;
    double dom = units->density_units/m_amu;


    ccount = 0;
    for (k = ks; k <= ke; k++){
    for (j = js; j <= je; j++){
    for (i = is; i <= ie; i++){
        c = ccount * N;
	idx = ((k* jdim + j)*idim + i);
        input[c]  = field_data->H_1_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->H_2_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->de_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->ge_density[idx];
        input[c] *= UNIT_E_per_M;
    c++;
    ccount += 1;

    }}}
}



int reshape_to_dengo_field_data_enzo( code_units* units, dengo_field_data *field_data, double* input, double *temp ){
    //------------------------------------------------------------------------------------
    // Function   :     reshape_to_dengo_field_data
    // Description:     reshape the 1d output array from solver to a dengo_field_data object  
    //                  and covert them to code units 
    //                  i.e. ge_density in erg /g
    //                       H_1_density in g / cm^-3 / amu (mass density per amu)
    // Parameter  :     code_units
    //                  dengo_field_data
    //                  input
    //------------------------------------------------------------------------------------

    unsigned long i, j, k, d; 
    unsigned long idx, ccount, c,dims;
    int is, ie, js, je, ks, ke;
    int ni, nj,nk;
    int N = 4;

    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    int idim = field_data->grid_dimension[0];
    int jdim = field_data->grid_dimension[1];
    int kdim = field_data->grid_dimension[2];

    ni = ie - is + 1;
    nj = je - js + 1;
    nk = ke - ks + 1;
    dims = ni*nj*nk;
    field_data->ncells = dims;
    dims = field_data->ncells; // total number of strips to be evaluated

    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;
    double dom = units->density_units/ m_amu;
   
    // extra bits for calculating  conservation
    ////////////////////////////////////////////
    ccount = 0;
    for (k = ks; k <= ke; k++){
    for (j = js; j <= je; j++){
    for (i = is; i <= ie; i++){
        c = ccount * N;
	idx = ((k* jdim + j)*idim + i);
        field_data->H_1_density[idx] = input[c];
        field_data->H_1_density[idx] /= dom;
	c++;
        field_data->H_2_density[idx] = input[c];
        field_data->H_2_density[idx] /= dom;
	c++;
        field_data->de_density[idx] = input[c];
        field_data->de_density[idx] /= dom;
	c++;
        field_data->ge_density[idx] = input[c];
        field_data->ge_density[idx] /= UNIT_E_per_M;
	c++;
	ccount += 1;
    }}}

    return 0;
}




int read_init_data_to_dengo( dengo_field_data *field_data, char const *filename){

    // this reads the initial abundances of the data from
    // a hdf5 file, and initialize a field_data object

    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    if (file_id < 0){
        fprintf(stderr, "failed to open %s so dying. \n", filename);
        return (1);
        }

    hsize_t dims;
    /* Check gas energy to get the number of cells */
    fprintf(stderr, "Getting dimensionality from ge:\n");
    herr_t status = H5LTget_dataset_info(file_id, "/ge", &dims, NULL, NULL);
    if(status == -1) {
        fprintf(stderr, "Error opening initial conditions file.\n");
        return 1;
    }
    fprintf(stderr, "  ncells = % 3i\n", (int) dims);
    
    field_data->ncells = (int) dims;
    int N = 4;
    double *atol, *rtol;
    atol = (double *) malloc(N * dims * sizeof(double));
    rtol = (double *) malloc(N * dims * sizeof(double));

    double *tics = (double *) malloc(dims * sizeof(double));
    double *ics = (double *) malloc(dims * N * sizeof(double));
    double *input = (double *) malloc(dims * N * sizeof(double));
    
    unsigned int i = 0, j;
    double *H_1 = (double *) malloc(dims * sizeof(double));    
    field_data->H_1_density = H_1;
    
    double *H_2 = (double *) malloc(dims * sizeof(double));    
    field_data->H_2_density = H_2;
    
    double *de = (double *) malloc(dims * sizeof(double));    
    field_data->de_density = de;
    
    double *ge = (double *) malloc(dims * sizeof(double));    
    field_data->ge_density = ge;
    
    fprintf(stderr, "Reading I.C. for /H_1\n");
    H5LTread_dataset_double(file_id, "/H_1", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->H_1_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "H_1[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /H_2\n");
    H5LTread_dataset_double(file_id, "/H_2", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->H_2_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "H_2[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /de\n");
    H5LTread_dataset_double(file_id, "/de", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->de_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "de[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /ge\n");
    H5LTread_dataset_double(file_id, "/ge", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->ge_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "ge[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    
    H5Fclose(file_id);
    free(input);
    free(ics);    
    return 1;
}

