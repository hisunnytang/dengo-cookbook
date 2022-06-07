#!/usr/bin/env python
# coding: utf-8

# Creating a Network
# ------------------
# We can now assemble them together in this `ChemicalNetwork`. This object helps us do all the neccessary computations to arrive at the symbolic **rhs** and **jacobian** functions, which ultimately eases us of the process of deriving them by hand. 
# 
# In our example with only two reactions, that consisted of $\rm H, H^+, e^-$. Note that in the below they represent the **number density $(\rm cm^{-3})$** of different species.
# 
# $$
#     \mathrm{ H^+ + e^- \rightarrow H}
# $$
# 
# $$
#     \mathrm{ H + e^- \rightarrow H^+ + 2 e^-}
# $$
# 
# $$
#     \begin{align*}
#     \mathrm{
#     \frac{d H}{dt} = k_{02}(T) \,  H^+ \, e^-  - k_{01}(T) \,  H \, e^- \\
#     \frac{d  H^+}{dt} = - k_{02}(T) \,  H^+ \, e^- + k_{01}(T) \,  H \, e^-  \\
#     \frac{d e^-}{dt} = - k_{02}(T) \,  H^+ \, e^-  + k_{01}(T) \,  H \, e^- \\
#     }
#     \end{align*}
# $$
# 
# We can do a quick sanity check on the conservation of the species $H$ and the charge $e^-$ and the equation above apparently satisfies our conservation law.
# 
# $$
#     \frac{d}{dt} (\rm H + H^+) = 0
# $$ 
# 
# $$
#     \frac{d}{dt} (\rm H^+ -e^-) = 0 
# $$

# ## Import Libraries and Create the Network
# 
# Primordial rates and cooling for the 9-species network are included in the default dengo library in `dengo.primordial_rates` and `dengo.primordial_cooling`. The reactions and cooling are added automatically to the `reaction_registry`, `cooling_registry` and `species_registry` with the call to `dengo.primordial_rates.setup_primordial`.  Here we setup the same sample network we demonstrated in the last chapter with `k01`, `k02` and `reHII`.

# In[3]:


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


# `simpleNetwork` is fed with the reactions and cooling actions. We are now ready to use `dengo` to build a solver.
