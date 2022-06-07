#include "simpleNetwork_solver.h"

int solve_chemistry(int argc, char **argv);
int solve_chemistry_enzo(int argc, char **argv);

// Handy Function to compare identical array
int compare_array(double *array, double rtol, unsigned long N){
    /* Check if the array element are the same within tolerance level */
    double ref, diff;
    ref = array[0];
    for (unsigned long i = 0; i < N; i++){
        diff = fabs(ref - array[i])/ref;
        if (diff > rtol){
            printf("Exceeded tolerance level (%0.5g); ref = %0.5g; array[%lu] = %0.5g\n", rtol, ref, i, array[i]);
            return 1;
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc > 1){
       simpleNetwork_main(argc, argv);
       return 0;
    }
    return 0;
    //}else{
    //    solve_chemistry_enzo(argc, argv);
    //}
}

// Sample on how to use dengo with primordial chemistry
// with simpleNetwork_solve_chemistry_enzo
/*
int solve_chemistry_enzo(int argc, char **argv) {
    dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));

    int Nx = 16;
    int Ny = 16;
    int Nz = 16;
    int N = Nx*Ny*Nz;
    field_data->ncells = N; 
    double density = 1.0e0; // in cm^-3
    double T = 1000.0; // in K
    double mH, k, tiny;

    mH = 1.67e-24;
    k  = 1.380e-16;
    tiny = 1.0e-20; 
    density *= mH;
    double G = 6.67259e-8;
    code_units *units = (code_units *) malloc(sizeof(code_units));
    units->density_units = 1.0;
    units->length_units = 1.0;
    units->time_units = 1.0;
    units->velocity_units = 1.0;
    double *H_1_density = (double*) malloc(N * sizeof(double));
    double *H_2_density = (double*) malloc(N * sizeof(double));
    double *de_density = (double*) malloc(N * sizeof(double));
    double *ge_density = (double*) malloc(N * sizeof(double));
    double *cooling_time = (double *) malloc( N * sizeof(double) );
    double *gamma = (double * ) malloc( N * sizeof(double) );
    double *temperature = (double *) malloc( N * sizeof(double) );
    double *mean_molecular_weight = (double *) malloc( N * sizeof(double) );
    double *density_arr = (double *) malloc( N * sizeof(double) );

    double *abstol = (double*) malloc(sizeof(double)* N * 4);
    double reltol = 1.0e-5;

    for ( int i = 0; i < N; i++){
        H2_1_density[i] = 1.0e-5 * density;
        H_1_density[i]  = 0.76   * density;
        He_1_density[i] = 0.24   * density;
	    de_density[i]    = 1.0e-5 * density;
	    H_2_density[i]   = 1.0e-5 * density;

	    // ge ~ nkT / (gamma - 1)/ rho; gamaa ~ 5/3
        ge_density[i]   = 3.0/2.0 * k * T / mH;
	    density_arr[i] = (1+2.0*1e-5)*density;
    }

    for ( int i = 0; i < N * 4; i++ ){
    	abstol[i] = tiny * reltol;
    }
    field_data->H_2_density = H_2_density;
    field_data->H_1_density = H_1_density;
    field_data->ge_density = ge_density;
    field_data->de_density = de_density;

    field_data->density         = density_arr;
    field_data->CoolingTime     = cooling_time;
    field_data->Gamma           = gamma;
    field_data->temperature     = temperature;
    field_data->MolecularWeight = mean_molecular_weight;

    int gstart[3];
    int gend[3];
    int gd[3];

    gstart[0] = 0;
    gstart[1] = 0;
    gstart[2] = 0;

    gend[0] = Nx-1;
    gend[1] = Ny-1;
    gend[2] = Nz-1;

    gd[0] = Nx;
    gd[1] = Ny;
    gd[2] = Nz;

    field_data->grid_start = &gstart[0];
    field_data->grid_end   = &gend[0];
    field_data->grid_dimension = &gd[0];


    const char *fileloc = "/home/kwoksun2/dengo_install/simpleNetwork_tables.h5";
    field_data->dengo_data_file = fileloc;
    field_data->reltol = reltol;

    units->a_value = 1.0;
    units->a_units = 1.0;

    double dt = 1.0 / sqrt(G * density) ;
    fprintf(stderr, "MAX_NCELLS = %d \n", MAX_NCELLS);
    simpleNetwork_solve_chemistry_enzo( units, field_data, dt );
    dengo_estimate_cooling_time_enzo( units, field_data);
    fprintf(stderr, "H_1 = %0.5g\n", field_data->H_1_density[0] / mH );
    fprintf(stderr, "H_2 = %0.5g\n", field_data->H_2_density[0] / mH );
    fprintf(stderr, "de = %0.5g\n", field_data->de_density[0] / mH );
    fprintf(stderr, "ge = %0.5g\n", field_data->ge_density[0] );
    fprintf(stderr, "CoolingTime = %0.5g\n", field_data->CoolingTime[0]);

    double *Pressure = (double *) malloc(sizeof(double)*N);
    double *Temperature = (double *) malloc(sizeof(double)*N);
    double *Gamma  = (double *) malloc(sizeof(double)*N);
    field_data->Gamma       = Gamma;
    field_data->Pressure    = Pressure;
    field_data->temperature = Temperature;

    dengo_calculate_pressure_enzo   (units, field_data);
    dengo_calculate_temperature_enzo(units, field_data);
    dengo_calculate_gamma_enzo      (units, field_data);


    for (int i = 0; i < 1; i++){
    fprintf(stderr, "Gamma    = %0.5g ", field_data->Gamma[i]);
    fprintf(stderr, "Pressure = %0.5g ", field_data->Pressure[i]);
    fprintf(stderr, "Temperature = %0.5g\n", field_data->temperature[i]);
    }

    fprintf(stderr, "\nGamma    = %0.5g ", field_data->Gamma[0]);
    compare_array(field_data->Gamma, 1.0e-4, N);
    fprintf(stderr, "\nPressure = %0.5g ", field_data->Pressure[0]);
    compare_array(field_data->Pressure, 1.0e-4, N);
    fprintf(stderr, "\n Temperature = %0.5g\n", field_data->temperature[0]);
    compare_array(field_data->temperature, 1.0e-4, N);
    unsigned long d;
    // lets just compare everything!!!!!!
    double ref0, frac;
    if (compare_array(field_data->H_2_density, reltol, N) == 1)
        printf("H_2 is not consistent\n");
    if (compare_array(field_data->H_1_density, reltol, N) == 1)
        printf("H_1 is not consistent\n");
    if (compare_array(field_data->ge_density, reltol, N) == 1)
        printf("ge is not consistent\n");
    if (compare_array(field_data->de_density, reltol, N) == 1)
        printf("de is not consistent\n");if (compare_array(field_data->CoolingTime, reltol, N) == 1)
        printf("CoolingTime is not consistent\n");

    free(field_data);
    free(H_1_density);
    free(H_2_density);
    free(de_density);
    free(ge_density);
    free(cooling_time);
    free(gamma);
    free(temperature);
    free(mean_molecular_weight);
    free(abstol);
}

*/