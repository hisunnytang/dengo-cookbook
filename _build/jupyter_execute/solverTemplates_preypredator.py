#!/usr/bin/env python
# coding: utf-8

# # How Does Dengo Generate Solvers?
# 
# This tutorial tells you how the solver are generated with the help of `dengo.ChemicalNetwork`, and `Jinja2`. In short, `dengo.ChemicalNetwork` carries the full information of the chemical reactions and cooling actions of interest. It internally generates the symbolic representation of the dynamics of each chemical species. This can be exported as `C++` or `python` code with a pre-written templates which can be found under `dengo/templates`. In this example we will be demonstrating how to generate rhs and solve the initial value problem with`scipy.odeint`

# # Prey-Predator Model
# We took the Prey-Predator model as our motivating example. It is also known as [Lotka-Volterra Equations (Wikipedia)](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations).
# It describes the dynamics of the two species (predator and prey) with 2 first order differential equations.
# $$
# \begin{align*}
# \frac{dx}{dt} &= \alpha x - \beta xy \\
# \frac{dy}{dt} &= \delta xy - \gamma y, \\
# \end{align*}
# $$
# where $x$ is the number of prey, $y$ is the number of predator, $t$ represents time.
# 
# Taken from wikipedia:
# ![](https://upload.wikimedia.org/wikipedia/commons/thumb/1/16/Lotka_Volterra_dynamics.svg/676px-Lotka_Volterra_dynamics.svg.png)

# # Importing Libraries
# 

# In[1]:


import numpy as np
from dengo.reaction_classes import \
 reaction, \
 ChemicalSpecies, \
 registry_setup, \
 species_registry
from dengo.chemical_network import ChemicalNetwork
import os
import pyximport
import matplotlib
import matplotlib.pyplot as plt
import copy
import pytest
import h5py


# ## Parameters Setup

# In[2]:


# parameters for the prey-predator model
α = 2./3.
β = 4./3.
γ = 1.0/2
δ = 2.0
# initial conditions
predator0 = 2.0
prey0 = 2.0


# ## Initialize `ChemicalNetwork`
# Our **network** consists of 4 species and 4 reactions.
# 
# ### Register the species
# ```python
# predator = ChemicalSpecies("predator", 1.0)
# dead_predator = ChemicalSpecies("dead_predator", 1.0)
# prey = ChemicalSpecies("prey", 1.0)
# dead_prey = ChemicalSpecies("dead_prey", 1.0)
# ```
# 
# ### Register the reactions
# 
# $$
# \begin{align*}
# \rm prey &\rightarrow \rm prey + prey\\
# \rm predator &\rightarrow \rm dead \rm ~ predator \\
# \rm prey + \rm predator &\rightarrow \rm dead ~ prey + \rm ~ predator + \rm ~   \frac{\gamma}{\beta} predator\\
# \end{align*}
# $$
# 
# ```python
# @reaction("exp_growth_prey", [(1, prey), ], [(2, prey), ])
# def rxn(state):
#     return α * np.ones_like(state.T)
# @reaction("predation", [(1, predator), (1, prey)], [
#        (1, dead_prey), (γ/β + 1, predator)])
# def rxn(state):
#     return β * np.ones_like(state.T)
# @reaction(
#  "natural_death_predator", [(1, predator), ],
#  [(1, dead_predator), ])
# def rxn(state):
#     return γ * np.ones_like(state.T)
# ```
# 
# ### Adding the reactions to the `ChemicalNetwork`
# ```python
# cN = ChemicalNetwork()
# cN.add_reaction("exp_growth_prey")
# cN.add_reaction("predation")
# cN.add_reaction("natural_death_predator")
# ```

# In[3]:


@registry_setup
def setup_predator_prey_rates():
    predator = ChemicalSpecies("predator", 1.0)
    dead_predator = ChemicalSpecies("dead_predator", 1.0)
    prey = ChemicalSpecies("prey", 1.0)
    dead_prey = ChemicalSpecies("dead_prey", 1.0)

    # predator-prey model
    @reaction("exp_growth_prey", [(1, prey), ], [(2, prey), ])
    def rxn(state):
        return α * np.ones_like(state.T)
    @reaction("predation", [(1, predator), (1, prey)], [
            (1, dead_prey), ((δ/β + 1), predator)])
    def rxn(state):
        return β*np.ones_like(state.T)
    @reaction(
     "natural_death_predator", [(1, predator), ],
     [(1, dead_predator), ])
    def rxn(state):
        return γ * np.ones_like(state.T)

def predator_prey_network():
    setup_predator_prey_rates()
    cN = ChemicalNetwork()
    cN.add_reaction("exp_growth_prey")
    cN.add_reaction("predation")
    cN.add_reaction("natural_death_predator")

    # this shouldnt be compulsory...
    cN.init_temperature((1e0, 1e8))
    return cN


# In[4]:


cn = predator_prey_network()


# ## Building the solver
# In this example, we will walk you through how to write a template from scratch that can be fed into a scipy solver. This can be done with `ChemicalNetwork.write_solver` and the combination of templates available under `dengo/templates`.
# 
# ### Evaluate the reaction rates
# Reaction rates usually have a temperature dependence. For example, for reactions following the (Arrhenius equation)[https://en.wikipedia.org/wiki/Arrhenius_equation] usually have the forms of $$k(T) = A e^{-\frac{E_a}{RT}}$$, where $k$ is the reaction rate, $E_a$ is the activation energy of the reaction, $T$ is the temperature, $A$, $R$ are the pre-exponential factor, and the universal gas constant respectively. $A$ is sometimes dependent further on temperature in (Modified Arrhenius equation)https://en.wikipedia.org/wiki/Arrhenius_equation#Modified_Arrhenius_equation].
# 
# Evaluating these rates on the fly would be computationally expensive. One possible way of reducing the computational time is to interpolate from a pre-calculated reaction rates table. The rates are specified when the reactions `rxn` are first created with the `@reaction` decorator. They can be evaluated handily with `rxn.coeff_fn(chemicalnetwork)`. The range of temperature of interest for example $T = \rm (1, 10^8) K$ can be first specified with `ChemicalNetwork.init_temperature(T_bounds=(1e0, 1e8), n_bins=1024)`. The added reaction objects can be accessed with `ChemicalNetwork.reactions`. For example, the reaction rates of `exp_growth_prey` can the accessed with the snippet below
# ```python
# rxn_rate = cn.reactions['exp_growth_prey'].coeff_fn(ChemicalNetwork)
# ```
# The output `rxn_rate` is an numpy array with a length of `[n_bins]`.
# 
# A reaction rate table is generated and exported to a `hdf5` file below.

# In[5]:


solver_name = 'prey_predator_solver'
output_dir  = "."


# In[6]:


ofn = os.path.join(output_dir, f"{solver_name}_tables.h5")
f = h5py.File(ofn, "w")

for rxn in sorted(cn.reactions.values()):
    f.create_dataset(
        f"/{rxn.name}", data=rxn.coeff_fn(cn).astype("float64")
    )
if hasattr(rxn, "tables"):
    for tab in rxn.tables:
        print(rxn.name, tab, rxn)
        f.create_dataset(
            f"/{rxn.name}_{tab}",
            data=rxn.tables[tab](self).astype("float64"),
        )
f.close()


# ### Evaluate the temperature
# 
# The temperature $T$ as we have seen above is critical to the rate at which the reaction proceeds. The temperature can be evaluated from the internal energy term `ge`. 
# Internal energy of an ideal gas is:
# $$ E = c_V T = \frac{nkT}{\gamma -1}$$
# For monoatomic gas $\gamma$ is $5/3$, and diatomic gas $\gamma$ is $7/5$. $\gamma$ refers to the adiabatic constant, and it is directly related to the degree of freedom available to the species $f = \frac{2}{\gamma -1}$. 
# 
# The total internal energy in the mixture of ideal gas is:
# $$E = \sum_s \frac{n_s kT}{\gamma_s -1}$$.
# $T$ can be thus be calculated from $E$ and the abundance of all the avaialble species $n_s$.

# ### The RHS function
# 
# The dynamics is specified by the set of ODE equations.
# $$ \frac{d \bf y}{dt} = f(\bf y) $$
# where $\bf y$ corresponds to the abundance vector for the species of interest, and $f(\bf y)$ describes the dynamics.
# 
# `Dengo` aggreates the reactions specific to each species $s$ with `ChemicalNetwork.species_total(s)` with `sympy` internally. These sympy expression can be exported to various different code style with `sympy.printing` to `C`, `python` for example.

# ### Ordinary Differential Equation
# 
# Here we outline the steps needed for a first order backward-euler integration. This is under the umbrella of a wider class of integration methods called implicit methods. 
# 
# > Implicit methods require an extra computation (solving the above equation), and they can be much harder to implement. Implicit methods are used because many problems arising in practice are stiff, for which the use of an explicit method requires impractically small time steps $\Delta t$ to keep the error in the result bounded (see numerical stability). That said, whether one should use an explicit or implicit method depends upon the problem to be solved.
# >
# > [Explicit and implicit methods: Computation](https://en.wikipedia.org/wiki/Explicit_and_implicit_methods) 
# 
# 
# #### Backward Euler Method
# $$ 
# \begin{align*}
# \frac{d \bf y}{dt} &= f(\bf y) \\
# y(t_{i+1}) &\approx y(t_i) + h f(y_{i+1}) \\
# F(x) &= x - y(t_i) - h f(x)
# \end{align*}
# $$
# where h is step-size of the integration. The solution of  $F(x) = 0$ gives straightforwardly $y_{i+1}$. The solution to $F(x)$ can be found iteratively by the newtons method
# $$
# x_{k+1} = x_k - \frac{F(x_k)}{F'(x_k)} \\
# F'(x) = 1 - h \frac{\partial f}{\partial x}
# $$
# Here $k$ corresponds to the step taken. The iteration is stopped when the difference $|x_{k+1} - x_{k}|$ is less than the given tolerance level. $F'(x)$ is the derivative of the function $F$ and requires the Jacobian.

# ### The Jacobian Function
# In cases where the set of reactions are stiff to evolve, the backward differentiation formulas and newton's method are often employed in conjunction with to integrate the system. The availability of an exact jacobian is beneficial to solving the stiff system efficiently. Note that this is also optional, as modern solver package could also approximate the jacobian numerically by finite difference. 
# 
# $$J = \frac{\partial \bf f}{ \partial \bf y} $$
# 
# In `Dengo`, the reactions are handled internally through the `sympy` engine. The derivative can be easily obtained from the sympy analytical derivatives `sympy.diff`. 
# 
# ```python
# import sympy
# x = sympy.symbols('x')
# xcube = x**3
# sympy.diff(xcube,x) == 3*x*x
# ```

# In[7]:


def f(state):
    """RHS function of each chemical species

    Parameters
    ----------
    state : ndarray with shape [NSPECIES + 1, N]
    Abundances sorted by name, and the last dimension corresponds to the current check_time

    Returns
    -------
    dy/dt: rate of change of each species
    """

# retreive the species 
    dead_predator,dead_prey,ge,predator,prey,current_time= state

# calculate temperature
    T = calculate_temperature(state)
    
# calculate mass density
    mdensity = 1.0*dead_predator + 1.0*dead_prey + 1.0*predator + 1.0*prey*mh;
    inv_mdensity = 1/mdensity;
        
# calculate the h2 optical depth approximation        
    h2_optical_depth_approx  = min( 1.0, pow( (mdensity / (1.34e-14) )  , -0.45) );
    
    tau      = pow( (mdensity / 3.3e-8 ), 2.8);
    tau      = max( tau, 1.0e-5 );
    cie_optical_depth_approx = min( 1.0, (1.0 - exp(-tau) ) / tau );

# interpolate the rates
    exp_growth_prey,predation,natural_death_predator, = interpolate_rates(T)

    
#     = interpolate_cooling_rates(T)

# rhs function
    
    ddead_predator = natural_death_predator[i]*predator*np.ones_like(ge)
    ddead_prey = predation[i]*predator*prey*np.ones_like(ge)
    dge = 0*np.ones_like(ge)
    dpredator = -natural_death_predator[i]*predator + 1.5*predation[i]*predator*prey*np.ones_like(ge)
    dprey = exp_growth_prey[i]*prey - predation[i]*predator*prey*np.ones_like(ge) 

    return np.array([ddead_predator,ddead_prey,dge,dpredator,dprey,0.0*current_time
    ])


# In[8]:


import numpy as np
import h5py

mh      = 1.67e-24
kb      = 1.38e-16
gamma   = 5./3.
gammaH2_1 = 7./5.
gammaH2_2 = 7./5.
_gamma_m1 = 1./ (gamma-1.)

# read rates in as global variables
rates_table = 'reaction_rates.h5'
ratef = h5py.File(rates_table, 'r')

# Reaction Rates
{% for k in network.reactions.keys()%}
out{{k}}dev = ratef['{{k}}'][:]
{%- endfor %} 

# Cooling Rates
{%- for name, rate in network.cooling_actions | dictsort %}
{%- for name2 in rate.tables | sort %}
out_{{name}}_{{name2}} = ratef["{{name}}_{{name2}}"][:]
{%- endfor %}
{%- endfor %}
tdev = ratef['T'][:]
ratef.close()

def interpolate_rates(T):
    """Interpolate all the reaction rates based on temperature
    """
    {% for k in network.reactions.keys()%}
    {{k}} = np.interp(T, tdev, out{{k}}dev)
    {%- endfor %} 
    return (
    {%- for k in network.reactions.keys() -%}
    {{k}}, 
    {%- endfor -%}
    )
def interpolate_cooling_rates(T):
    """Interpolate all the cooling rates based on temperature
    """
    {%- for name, rate in network.cooling_actions | dictsort %}
    {%- for name2 in rate.tables | sort %}
    {{name}}_{{name2}} = np.interp(T, tdev, out_{{name}}_{{name2}})
    {%- endfor %}
    {%- endfor %}
    return (
    {%- for name, rate in network.cooling_actions | dictsort -%}
    {%- for name2 in rate.tables | sort -%}
    {{name}}_{{name2}}, 
    {%- endfor -%}
    {%- endfor -%}
    )

def calculate_temperature(state):
    """calculate temperature based on the N different input state

    Parameters
    ----------
    state : ndarray with shape [NSPECIES + 1, N]
    Abundances sorted by name, and the last dimension corresponds to the current check_time

    Returns
    -------
    Temperature: ndarray

    """
# retreive the species 
    {% for s in network.required_species | sort -%}
    {{s.name}}, 
    {%- endfor -%}
    _= state


    density = {{network.print_mass_density()}}

    return {{network.temperature_calculation()}}



def f(state):
    """RHS function of each chemical species

    Parameters
    ----------
    state : ndarray with shape [NSPECIES + 1, N]
    Abundances sorted by name, and the last dimension corresponds to the current check_time

    Returns
    -------
    dy/dt: rate of change of each species
    """

# retreive the species 
    {% for s in network.required_species | sort -%}
    {{s.name}}, 
    {%- endfor -%}
    current_time= state

# calculate temperature
    T = calculate_temperature(state)
    
# calculate mass density
    mdensity = {{network.print_mass_density()}}*mh;
    inv_mdensity = 1/mdensity;
        
# calculate the h2 optical depth approximation        
    h2_optical_depth_approx  = min( 1.0, pow( (mdensity / (1.34e-14) )  , -0.45) );
    
    tau      = pow( (mdensity / 3.3e-8 ), 2.8);
    tau      = max( tau, 1.0e-5 );
    cie_optical_depth_approx = min( 1.0, (1.0 - exp(-tau) ) / tau );

# interpolate the rates
    {% for k in network.reactions.keys() -%}
    {{k}}, 
    {%- endfor %} = interpolate_rates(T)

    
    {% for name, rate in network.cooling_actions | dictsort -%}
    {%- for name2 in rate.tables | sort -%}
    {{name}}_{{name2}}, 
    {%- endfor -%}
    {%- endfor -%} = interpolate_cooling_rates(T)

# rhs function
    {% for s in network.required_species | sort %}
    d{{s.name}} = {{rhs_dict[s]}}*np.ones_like(ge)
    {%- endfor %} 

    return np.array([
    {%- for s in network.required_species | sort -%}
    d{{s.name}}, 
    {%- endfor -%}
    0.0*current_time
    ])


# ### Jinja2 Template Writer
# `Jinja2` is a popular templating engine. For example, if we have a "template" file as below. 
# 
# ```
# %%writefile solver.py
# {{name}} had a little {{animal}}
# ```
# 
# 

# In[150]:


get_ipython().run_cell_magic('writefile', 'jinja_example.txt', '{{name}} had a little {{animal}}\n')


# In[148]:


from jinja2 import Environment, FileSystemLoader
file_loader = FileSystemLoader(".")
env = Environment(loader=file_loader)
template = env.get_template("jinja_example.txt")


# In[94]:


cn.species_total('predator')


# In[80]:


cn.species_total('dead_prey')


# In[121]:


from sympy.printing.pycode import pycode


# In[122]:


rhs_dict = {}
for s in cn.required_species:
    rhs_dict[s] = pycode(cn.species_total(s))


# In[123]:





# In[124]:


import jinja2


# In[125]:


get_ipython().system('ls ~/data/dengo-merge/dengo/templates/scipy')


# In[126]:


get_ipython().system('cat /mnt/gv0/homes/kwoksun2/dengo-merge/dengo/templates/scipy/dengo_scipy.py.template')


# In[104]:


import jinja2

templateLoader = jinja2.FileSystemLoader(searchpath="/mnt/gv0/homes/kwoksun2/dengo-merge/dengo/templates/scipy")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "dengo_scipy.py.template"
template = templateEnv.get_template(TEMPLATE_FILE)


# In[127]:


template_vars = dict(
             network=cn, solver_name="solver_name", init_values={}, rhs_dict=rhs_dict
         )


# In[129]:


with open('solver.py', 'w') as f:
    f.write(template.render(template_vars))


# In[132]:


import solver


# In[95]:


import pyximport
import os
import numpy as np

solver_name = "test_grackle"
output_dir  = '.'
network = dengo_network

# specify the library path
os.environ["HDF5_DIR"] = "/home/kwoksun2/anaconda3"
os.environ["CVODE_PATH"] = "/home/kwoksun2/dengo-merge/cvode-3.1.0/instdir"
os.environ["HDF5_PATH"]  = "/home/kwoksun2/anaconda3"
os.environ["SUITESPARSE_PATH"] = "/home/kwoksun2/dengo-merge/suitesparse"
os.environ["DENGO_INSTALL_PATH"] = "/home/kwoksun2/dengo_install"

# install the library
pyximport.install(setup_args={"include_dirs": np.get_include()},
                  reload_support=True, inplace=True)

network.write_solver(solver_name, output_dir=output_dir,
                solver_template="cv_omp/sundials_CVDls",
                ode_solver_source="initialize_cvode_solver.C")


# In[ ]:




