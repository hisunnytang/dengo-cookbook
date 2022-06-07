import numpy
import pyximport
import os

# write the solver for various network

# test the cython intallation
os.environ["HDF5_DIR"] = /home/kwoksun2/anaconda3
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import simpleNetwork_solver_run

# this runs the `simpleNetwork_main` function in the simpleNetwork_solver
# reading the simpleNetwork_initial_conditions.h5
simpleNetwork_solver_run.main_run_simpleNetwork()

# test the run_simpleNetwork(init, )