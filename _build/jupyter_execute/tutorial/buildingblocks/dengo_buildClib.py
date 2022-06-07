#!/usr/bin/env python
# coding: utf-8

# # Build C-library
# This template comes with a `Makefile` that lets the user to compile the library. Apart from specifying the `CVODE_PATH`, `SUITESPARSE_PATH`, `DENGO_INSTALL_PATH`, `HDF5_PATH`, there are additional options the user can choose. One can choose either the `CVDLS` which refers to the CVODE dense linear solver, or `CVKLU` which refers to the CVODE KLU Sparse solver. `CVKLU` is oftentimes preferred as it leverages the emptyness in the Jacobian matrix, and utilizes the SuiteSparse library to perform the newton iterations. `DENGO_INSTALL_PATH` is where the library will be installed. `MAX_NCELLS` determines how many cells would be batched together when the `CVODE` solver is invoked.
# ```
# # Use the Dense Linear Solver interface (CVDLS)
# #OPTIONS+= -DMAX_NCELLS=1
# # Use the Sparse KLU direct solver
# OPTIONS = -DCVKLU -DMAX_NCELLS=256
# ```
# The `Make` script uses `libtool` to compile the files and build the libraries. The command `Make` would build and install the libraries. 
# Now you have built a `libdengo.so` in your specified installed path, and you are ready to link it to your favorite program with the sample below, or checkout the `Makefile:run_dengo`.
# 
# ```
# $(CC) -o $@ $^ $(OPTIONS) -I$(DENGO_INSTALL_PATH)/include -I$(CVODE_PATH)/include -I$(INCLUDES_KLU) -I$(HDF5_PATH)/include -L$(HDF5_PATH)/lib -lhdf5_hl -lhdf5 -L$(DENGO_INSTALL_PATH)/lib -ldengo -lm
# ```

# ## Main components that drives `Dengo` C-solver
# 
# In brief, similar to the example we shown before with the python RHS function, the C templates generate primarily solvers
# 
# - <font color='blue'>{{solver_name}}_solver.h </font>
# - <font color='blue'>{{solver_name}}_solver.C </font>(major modules in Dengo)
# - <font color='blue'>{{solver_name}}_solver_main.h </font>
# - <font color='blue'>{{solver_name}}_solver_main.C </font>(example script to use the C library) 
# - <font color='blue'>initialize_cvode_solver.C </font>(wrapper function for the CVode library)
# - <font color='blue'>Makefile </font>(to compile the dengo library `libdengo.a`)
# 
# 
# An example usage of the built dengo library can be found in the `{{solver_name}}_solver_main.C`. The solver is accessed through 
# ```C
# int simple_solve_chemistry_enzo( code_units *units, dengo_field_data *field_data, double dt );
# ```
# This API follows closely the `Grackle` solver `solver_chemistry` which is tailored to work with `Enzo` simulations. `dengo_field_data` is a data structure that points to the respective abundance pointer `{species.name}_density`. It also informs the setup of the incoming grid. `code_units` informs the dengo solver how to convert from code_units into physical cgs units. `dt` is the delta time the user want to advance the system for. The structure definition taken from `{{solver_name}}_solver.h` are shown below. 
# `dengo_field_data.CoolingTime`, `dengo_field_data.MolecularWeight`, `dengo_field_data.temperature`, `dengo_field_data.Gamma`, `dengo_field_data.Pressure` are **not** required inputs. User however needs to allocate space to these pointers if one want Dengo to compute these quantities alongside the chemical evolution. `dengo_data_file` points to the reaction rates `hdf5` data for on-the-fly interpolations.
# 
# ```C
# typedef struct dengo_field_data
# {
# 
#   unsigned long int nstrip;
#   unsigned long int ncells; 
#   // let's just pass them passively through field_data
#   double reltol;
#   double floor_value;
#   // This should be updated dynamically 
#   // with dengo
#   double *density;
#   double *H_1_density;
#   double *H_2_density;
#   double *de_density;
#   double *ge_density;
#     
#   double *CoolingTime;
#   double *MolecularWeight;
#   double *temperature;
#   double *Gamma;
#   double *Pressure;
# 
#   int *grid_start;
#   int *grid_end;
#   int *grid_dimension;
# 
#   const char *dengo_data_file;
#   code_units *units;
# } dengo_field_data;
# 
# typedef struct code_units
# {
# 
#   int comoving_coordinates;
#   double density_units;
#   double length_units;
#   double time_units;
#   double velocity_units;
#   double a_units;
#   double a_value;
# 
# } code_units;
# ```

# ## Conclusion
# 1. We outlie how `ChemicalNetwork.write_solver` works.
# 2. How to build the cython module from the templates generated.
# 3. We also outline how the C-library can be built.
# 
# In the next chapter, we will outline the procedures needed to link hydrodynamical simulations like `enzo` to `dengo`.
