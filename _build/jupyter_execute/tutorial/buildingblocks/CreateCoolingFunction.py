#!/usr/bin/env python
# coding: utf-8

# # Create Cooling Functions
# 
# The release/ absorption of energy via various physical process would modify the thermal state of the medium, which might in turn accelerates/ terminates further reactions. Keeping track of the thermal energy through the course of evolution is quintessential to the accurate modeling of the process. 
# 
# ```{note}
# While it is termed cooling actions, or cooling function in `dengo`, it is **not** limited only to cooling functions. This module can account for general modifications to the thermal energy of the fluid parcel $\frac{d\epsilon}{dt}$ whether it is injecting/ releasing energy from the fluid parcel. For example, when molecular hydrogen is formed from 3-body reactions, $4.48 eV$ is released. We will demonstrate that how to account for the chemical heating in one of the chapters to follow. 
# ```
# 
# 
# $k01$ is the reaction that describes the recombination of an ionized hydrogen atom with an electron. Energy is released as a result of such process.
#  
#  
# ## Register new cooling functions with the `cooling_action` decorator
# Similar to reactions, the cooling functions can be specified with a decorator `dengo.reaction_classes.cooling_action`.
# a `cooling_action` to `dengo` is primarily composed of:
# - **name** of the cooling function, i.e. `reHII`
# - **Cooling Equation**, `-reHII * HII * de`
# - **Cooling Rate**, that often times dependent on the temperature of the gas parcel, and sometimes on external factors as well.
# 
# In the example below, `HII` and `de` are defined and registered in `dengo` species registry. `reHII` the cooling rate associated to this particular cooling function. The cooling rate takes `eq` as an input. `@eq.table` register the cooling rates that are dependent on temperature. From then on `dengo` would associate the symbol `reHII` appeared in the cooling equation to the user-supplied cooling rate `def reHII(state)` internally.

# In[1]:


from dengo.reaction_classes import cooling_action
from dengo.chemical_network import cooling_registry

# -- reHII --
@cooling_action("reHII", "-reHII * HII * de")
def cool(eq):
    @eq.table
    def reHII(state):
        # (from Hui and Gnedin 1997)
        lambdaHI = 2.0 * 157807e0 / state.T
        vals = 1.778e-29 * state.T * lambdaHI**1.965 / \
         (1.0e0 + (lambdaHI/0.541)**0.502)**2.697
        return vals


# In[2]:


cooling_registry['reHII'], type(cooling_registry['reHII'])

