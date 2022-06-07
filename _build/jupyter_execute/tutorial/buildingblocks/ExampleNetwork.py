#!/usr/bin/env python
# coding: utf-8

# # An example of 2-species chemical network
# 
# In our example with only two reactions, that consisted of $\rm H, H^+, e^-$. The first reaction corresponds to the recombination of ionized hydrogen and an electron. The second reaction refers to the collisional ionization of neural hydrogen atom to positively charged atoms by an electron. These are 
# 
# $$
#     \begin{align*}
#     \rm H^+ + e^- &\rightarrow \rm H \quad &(k01)  \\
#     \rm H + e^- &\rightarrow \rm H^+ + 2 e^- \quad &(k02)
#     \end{align*}
# $$
# 
# These two reactions can be written as the set of ODE equations below. 
# Here $k_{01}(T)$ and $k_{02}(T)$ are the temperature dependent reaction rates. 
# Note that in the below they represent the **number density $(\rm cm^{-3})$** of different species.
# $$
#     \begin{align*}
#     \rm \frac{d H}{dt} &=  \rm k_{02}(T) \,  H^+ \, e^-  - k_{01}(T) \,  H \, e^- \\
#     \rm \frac{d  H^+}{dt} &= \rm - k_{02}(T) \,  H^+ \, e^- + k_{01}(T) \,  H \, e^-  \\
#     \rm  \frac{d e^-}{dt} &=  \rm - k_{02}(T) \,  H^+ \, e^-  + k_{01}(T) \,  H \, e^- \\
#     \end{align*}
# $$ 

# In[ ]:




