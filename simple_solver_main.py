import numpy
import pyximport
import os

# write the solver for various network

# test the cython intallation
os.environ["HDF5_DIR"] = /home/kwoksun2/anaconda3
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import simple_solver_run

# this runs the `simple_main` function in the simple_solver
# reading the simple_initial_conditions.h5
simple_solver_run.main_run_simple()

# test the run_simple(init, )