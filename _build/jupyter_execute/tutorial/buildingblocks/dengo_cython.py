#!/usr/bin/env python
# coding: utf-8

# # Build Cython Module from Templates
# In the last chapter, we used our `simpleNework` to generate 9 files from `dengo`. We outline how can we build the library from here and use cython to invoke the function written in `C`.
# 
# ```
# simpleNetwork_solver.C	      simpleNetwork_solver_run.pyx
# simpleNetwork_solver.h	      simpleNetwork_solver_run.pyxbld
# simpleNetwork_solver_main.C   simpleNetwork_solver_run.pyxdep
# simpleNetwork_solver_main.py  simpleNetwork_tables.h5
# simpleNetwork_solver_run.pxd
# ```

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


# ## Compile the Cython Module 
# The solver templates come with a `.pyx` files that lets you build a python interface to the C-library. A more well-rounded `Cython` tutorial can be found here [Cython Tutorial].(https://cython.readthedocs.io/en/latest/src/tutorial/cython_tutorial.html).

# In[2]:


import pyximport
import numpy as np
solver_name = "simpleNetwork"

pyximport.install(
    setup_args={"include_dirs":np.get_include()},
    reload_support=True, 
    inplace=True, 
    language_level=3
)

simple_solver_run = pyximport.load_module(
    f"{solver_name}_solver_run",
    f"{solver_name}_solver_run.pyx",
    build_inplace = True, 
    pyxbuild_dir = "_dengo_temp", 
    )


# ## Invoking the Cython Solver
# In the below we show how the cython module can be built and used.
# `{solver_name}_solver_run.run_{solver_name}(init_values, dt, niter)` is the entry point for the built solver.
# It takes an dictionary that contains the abundances, and thermal energy, an $\rm dt$, time to advance the fluid parcel, and `niter` the maximum number of iterations as arguments.
# It assumes that abundances are in number density with the units of $\rm cm^{-3}$, and thermal energy in $\rm erg/g$. `niter` implicitly sets the initial timestep for the solver, i.e. $dt_{solver} = \rm dt / \rm niter$.
# Our `cython` module gives the same trajectories as what we had done from scratch with the previous chapters! 

# In[16]:


NCELLS = 1
density = 1e-2

init_array = np.ones(NCELLS) * density
init_values = dict()
init_values['H_1']     = init_array 
init_values['H_2']     = init_array 
init_values['de']      = init_array 
init_values['ge']      = np.ones(NCELLS)*1e13

total_density = simpleNetwork.calculate_total_density(init_values)
init_values = simpleNetwork.convert_to_mass_density(init_values)
init_values['de'] = simpleNetwork.calculate_free_electrons(init_values)
init_values['density'] = simpleNetwork.calculate_total_density(init_values)
number_density = simpleNetwork.calculate_number_density(init_values)

rv, rv_int = simple_solver_run.run_simpleNetwork(init_values, 1e16, niter = 1e5)


# ## Visualize the Trajectories
# They match with what we have built from scratch before!

# In[17]:


import matplotlib.pyplot as plt
flags = rv_int['t'] > 0
sp_names = [s.name for s in sorted(simpleNetwork.required_species)]
output = np.zeros((len(sp_names),sum(flags)))
                  
for i, s in enumerate(sp_names):
    output[i] = rv_int[s][:,flags]

H_1_traj, H_2_traj, de_traj, ge_traj = output
T_traj    = rv_int['T'][0,flags]
timesteps = rv_int['t'][flags]


f, ax = plt.subplots()
l1= ax.loglog(timesteps, H_1_traj)
l2= ax.loglog(timesteps, H_2_traj)
l3= ax.loglog(timesteps, de_traj)
ax2 = ax.twinx()
l4= ax2.loglog(timesteps, T_traj, color='C3')


ax.set_ylabel(r"Abundance ($\rm cm^{-3}$)")
ax2.set_ylabel("Temperature (K)")

ax.set_xlabel("Time (s)")

f.legend(
    [l1,l2,l3,l4], 
    labels=[r'$\rm H$',r'$\rm H^+$',r'$\rm e^-$',r'$\rm T$'],
    loc=[0.15,0.4]
)
plt.subplots_adjust(right=0.85)


# In[ ]:




