#!/usr/bin/env python
# coding: utf-8

# In[1]:


# importing necessary libraries
import pyximport
import numpy as np
import os

from   dengo.chemistry_constants import tiny, kboltz, mh, G
import dengo.primordial_rates
import dengo.primordial_cooling

from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry, \
    species_registry
import dengo.primordial_rates
import dengo.primordial_cooling

import sympy
from sympy import lambdify
import matplotlib.pyplot as plt


# In[2]:


def setup_primordial_network():
    """Initial a ChemicalNetwork object
    for primordial network 9-species model
    Return:
     primordial: ChemicalNetwork with primordial reactions and cooling
    """
    # this register all the rates specified in `primordial_rates.py`
    dengo.primordial_rates.setup_primordial()

    # initialize the chmical network object
    primordial = ChemicalNetwork()

    # add all the reactions
    primordial.add_reaction("k01")
    primordial.add_reaction("k02")
    primordial.add_reaction("k03")
    primordial.add_reaction("k04")
    primordial.add_reaction("k05")
    primordial.add_reaction("k06")
    primordial.add_reaction("k07")
    primordial.add_reaction("k08")
    primordial.add_reaction("k09")
    primordial.add_reaction("k10")
    primordial.add_reaction("k11")
    primordial.add_reaction("k12")
    primordial.add_reaction("k13")
    primordial.add_reaction("k14")
    primordial.add_reaction("k15")
    primordial.add_reaction("k16")
    primordial.add_reaction("k17")
    primordial.add_reaction("k18")
    primordial.add_reaction("k19")
    primordial.add_reaction("k21")
    primordial.add_reaction("k22")
    primordial.add_reaction("k23")

    primordial.add_cooling("brem")
    primordial.add_cooling("reHII")
    primordial.add_cooling("reHeIII")
    primordial.add_cooling("gloverabel08")
    primordial.add_cooling("ceHI")
    primordial.add_cooling("h2formation")
    primordial.add_cooling("reHeII2")
    primordial.add_cooling("reHeII1")
    primordial.add_cooling("ciHeIS")
    primordial.add_cooling("ceHeII")
    primordial.add_cooling("ciHI")
    primordial.add_cooling("ceHeI")
    primordial.add_cooling("gammah")
    primordial.add_cooling("ciHeI")
    primordial.add_cooling("ciHeII")
    primordial.add_cooling("cie_cooling")
    primordial.add_cooling("compton")

    # This defines the temperature range for the rate tables
    primordial.init_temperature((1e0, 1e8))
    
    return primordial
#     primordial.enforce_conservation = True
#     primordial.set_equilibrium_species("H2_2")

cn_simple = setup_primordial_network()


# In[3]:


output_dir = "."
solver_name = "simple"
use_omp = True
use_cvode = True
use_suitesparse = True

# specify the library path
os.environ["HDF5_DIR"] = "/home/kwoksun2/anaconda3"
os.environ["CVODE_PATH"] = "/home/kwoksun2/dengo-merge/cvode-3.1.0/instdir"
os.environ["HDF5_PATH"]  = "/home/kwoksun2/anaconda3"
os.environ["SUITESPARSE_PATH"] = "/home/kwoksun2/dengo-merge/suitesparse"
os.environ["DENGO_INSTALL_PATH"] = "/home/kwoksun2/dengo_install"

# write the solver
# cn_simple.write_solver(solver_name, output_dir=output_dir,
#                 solver_template="be_chem_solve/rates_and_rate_tables",
#                 ode_solver_source="BE_chem_solve.C")


# In[4]:


# install the library
pyximport.install(setup_args={"include_dirs": np.get_include()},
                  reload_support=True, inplace=True)


# In[5]:


# compile a BE chem solve one
solver_nameBE = "simpleBE"
cn_simple.write_solver(solver_nameBE, output_dir=output_dir,
                solver_template="be_chem_solve/rates_and_rate_tables",
                ode_solver_source="BE_chem_solve.C")
simple_solver_run_beChem = pyximport.load_module(
    "{}_solver_run".format(solver_nameBE),
    "{}_solver_run.pyx".format(solver_nameBE),
    build_inplace=True, pyxbuild_dir="_dengo_temp")


# In[6]:


solver_nameCV = "simpleCV"
cn_simple.write_solver(solver_nameCV, output_dir=output_dir,
                solver_template="cv_omp/sundials_CVDls",
                ode_solver_source="initialize_cvode_solver.C")
simple_solver_run_CV = pyximport.load_module(
    "{}_solver_run".format(solver_nameCV),
    "{}_solver_run.pyx".format(solver_nameCV),
    build_inplace=True, pyxbuild_dir="_dengo_temp")


# In[9]:


def setup_initial_conditions(network, density, temperature, h2frac, NCELLS):
    # setting initial conditions
    temperature = np.ones((NCELLS))*temperature
    init_array = np.ones(NCELLS) * density
    init_values = dict()
    init_values["H_1"] = init_array * 0.76 * ( 1-h2frac )
    init_values['H_2'] = init_array * tiny
    init_values['H_m0'] = init_array * tiny
    init_values['He_1'] = init_array * 0.24
    init_values['He_2'] = init_array * tiny
    init_values['He_3'] = init_array * tiny
    init_values['H2_1'] = init_array * 0.76 * h2frac
    init_values['H2_2'] = init_array * tiny
    init_values['de'] = init_array * tiny

    # update and calculate electron density and etc with the handy functions
    # init_values = primordial.convert_to_mass_density(init_values)
    init_values['de'] = network.calculate_free_electrons(init_values)
    # init_values['density'] = primordial.calculate_total_density(init_values)
    init_values['density'] = np.ones((NCELLS))*density
    number_density = network.calculate_number_density(init_values)

    # set up initial temperatures values used to define ge
    init_values['T'] = temperature

    # calculate ge (very crudely, no H2 help here)
    gamma = 5.0/3.0
    mH = 1.67e-24
    init_values["ge"] = 3.0 / 2.0 * temperature * kboltz / mH

    return init_values


# In[72]:


# Initial Conditions
density = 1e10 # number density
temperature = 1000.0 # K
H2Fraction  = 0.001 # molecular mass fraction
ncells = 1 # number of cells

# freefall timescale
dtf = 10.0/ np.sqrt(G*mh*density)

states =  setup_initial_conditions(cn_simple, density, temperature, H2Fraction, ncells)
rv, rv_intBE = simple_solver_run_beChem.run_simpleBE(states, dtf, niter=1e6, reltol = 1.0e-5, z = 0.0);


# In[73]:


# Initial Conditions
density = 1e10 # number density
temperature = 1000.0 # K
H2Fraction  = 0.001 # molecular mass fraction
ncells = 1 # number of cells

# freefall timescale
dtf = 10.0/ np.sqrt(G*mh*density)

states =  setup_initial_conditions(cn_simple, density, temperature, H2Fraction, ncells)
rv, rv_intCV = simple_solver_run_CV.run_simpleCV(states, dtf, niter=1e6, reltol = 1.0e-5, z = 0.0)


# In[74]:


f, ax = plt.subplots(figsize=(6,4))
ax2 = ax.twinx()

def plot_results(ax, ax2, rv_int, label, ls='-'):
    flag = rv_int["successful"]
    t    = rv_int["t"][flag]
    H2I  = rv_int["H2_1"][0][flag]
    HI   = rv_int["H_1"][0][flag]
    T    = rv_int['ge'][0][flag]


    ax.loglog(t, HI, label='HI '+label, alpha=0.5, ls=ls)
    ax.loglog(t, H2I, label="H2I " +label, alpha=0.5, ls =ls)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("number density $(\mathrm{cm^{-3}})$")
    
    ax2.loglog(t, T, label='T', color='k', ls =ls)
    ax2.set_ylabel("temperature (K)")
    return ax, ax2 

ax, ax2 = plot_results(ax, ax2, rv_intCV, label='CV')
ax, ax2 = plot_results(ax, ax2, rv_intBE, label='BE', ls ='--')

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2)


# In[34]:


# extract_final_step(rv_intBE), 
extract_final_step(rv_intCV)


# In[33]:


def extract_final_step(result):
    flag = result['successful']
    skip = ['successful', 't', 'dt']
    out = {}
    for k, v in result.items():
        if k in skip: continue
        out[k] = v[0][flag][-1]
    out['t'] = result['t'][flag][-1]
    return out

# extract_final_step(rv_intBE)


# In[70]:


def run_grid(solver, density, temperature, h2frac):
    results = []
    for d in density:
        for T in temperature:
            for f in h2frac:
                init_values = setup_initial_conditions(cn_simple, d, T, f, 1)
                dtf = float(1.0/ np.sqrt(G*mh*d))
                rv, rv_int = solver(init_values, dtf, niter=1e4, reltol = 1.0e-5)
                results.append(extract_final_step(rv_int))
    return results


# In[142]:


darray = np.logspace(0,10,11)
tarray = np.logspace(2, 3.5,11)
h2frac = [1.0e-6, 1.0e-5, 1.0e-4]
cv_grid = run_grid(simple_solver_run_CV.run_simpleCV, darray, tarray, h2frac)


# In[143]:


darray = np.logspace(0,10,11)
tarray = np.logspace(2, 3.5,11)
h2frac = [1.0e-6, 1.0e-5, 1.0e-4]
be_grid = run_grid(simple_solver_run_beChem.run_simpleBE, darray, tarray, h2frac)


# In[148]:


def reshape_grid(grid):
    keys = grid[0].keys()
    
    out = {}
    for k in keys: out[k] = []
    for g in grid:
        for k in keys:
            out[k].append(g[k])
    return {k: np.array(v) for k, v in out.items()}


# In[149]:


be_grid_ =reshape_grid(be_grid)
cv_grid_ =reshape_grid(cv_grid)


# In[150]:


ratio = {}
for k in cv_grid_.keys():
    ratio[k] = be_grid_[k]/ cv_grid_[k]


# In[152]:


for k in ratio.keys():
    print(k)
    err = np.abs(ratio[k].reshape(11,11,3)[:,:,0] -1 + 1e-12)
    plt.pcolormesh(tarray, darray, np.log10(err), vmin=-3, vmax=0, cmap='jet')
    plt.xscale('log')
    plt.yscale('log')
    plt.colorbar()
    plt.show()


# In[160]:


err = np.abs(ratio['T'].reshape(11,11,3)[:,:,1] -1 + 1e-12)


# In[161]:


err.max()


# In[7]:


import h5py
import matplotlib.pyplot as plt


# In[2]:


f= h5py.File("simpleBE_tables.h5", 'r')


# In[3]:


f.keys()


# In[9]:


plt.plot(f['dgammaH2_1_dT'][:])


# In[ ]:




