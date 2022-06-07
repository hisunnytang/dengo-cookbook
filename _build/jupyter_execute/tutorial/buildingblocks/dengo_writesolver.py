#!/usr/bin/env python
# coding: utf-8

# # Making Use of Dengo templates
# In this chapter, we will walk through how to build the C library with the template provided in `dengo.templates`. This is particularly useful when one wants to extend the solver to solve the chemistry alongside hydrodynamical simulations like `enzo`. 
# 
# In this example, we will stick with the 2-species primordial chemistry model we used in the last two chapters.

# ## Generate templates with `ChemicalNetwork.write_solver`
# With the `ChemicalNetwork` object, `Dengo` can write the C solver.
# The corresponding auxillary library paths are needed to be set as the **environomental variables**, in order to compile our python modules.
# 
# - <font color='green'>HDF5_DIR </font>          (HDF5 installation path)
# - <font color='green'>CVODE_PATH </font>        (CVode installation path)
# - <font color='green'>SUITESPARSE_PATH </font>  (SuiteSparse library which is optional unless we use `KLU` option)
# - <font color='green'> DENGO_INSTALL_PATH</font> (Installation Path of Dengo)
# 
# `solver_template` are the existing templates under `dengo.templates`. Currently there are two major ones `be_chem_solve`

# In[1]:


import dengo
from dengo.chemical_network import \
 ChemicalNetwork, \
 reaction_registry, \
 cooling_registry, species_registry
import dengo.primordial_rates
import dengo.primordial_cooling

dengo.primordial_rates.setup_primordial()

simpleNetwork = ChemicalNetwork()
simpleNetwork.add_reaction("k01")
simpleNetwork.add_reaction("k02")
simpleNetwork.add_cooling("reHII")
simpleNetwork.init_temperature((1e0, 1e8))


# In[2]:


solver_name = "simpleNetwork"
simpleNetwork.write_solver(
    solver_name, 
    output_dir = ".", 
    solver_template = "cv_omp/sundials_CVDls", 
    ode_solver_source = "initialize_cvode_solver.C"
)


# In[3]:


get_ipython().system('ls simpleNetwork_*')


# ## `ChemicalNetwork.write_solver` 
# Look in your directory there are 9 extra files, from `C` to `python` codes!
# 
# ### Main components that drives `Dengo` C-solver
# 
# - <font color='blue'>{{solver_name}}_solver.h </font>
# - <font color='blue'>{{solver_name}}_solver.C </font>(major modules in Dengo)
# - <font color='blue'>{{solver_name}}_solver_main.h </font>
# - <font color='blue'>{{solver_name}}_solver_main.C </font>(example script to use the C library) 
# - <font color='blue'>initialize_cvode_solver.C </font>(wrapper function for the CVode library)
# - <font color='blue'>Makefile </font>(to compile the dengo library `libdengo.a`)
# 
# ### Helper function to compile `Dengo` C files for `Python` wrapper
# - <font color='blue'>{{solver_name}}_solver_run.pyxbld </font>
# - <font color='blue'>{{solver_name}}_solver_run.pyxdep </font>
# - <font color='blue'>{{solver_name}}_solver_run.pxd </font>
# - <font color='blue'>{{solver_name}}_solver_run.pyx </font> (major Python wrapper)
