#!/usr/bin/env python
# coding: utf-8

# # Interfacing with Hydro Simulations
# We have come up with a sample script for the combining `dengo` and [enzo](https://github.com/enzo-project/enzo-dev) together. 
# 
# 
# https://github.com/hisunnytang/dengo-templates
# It is a helper function that helps generate the tmeplates for enzo simulations.
# 
# 
# ## Outline of what `generateEnzoTemplates.py` does
# - it takes `ChemicalNetwork` as input
# - Based on the reactions and chemical cooling/ heating, the required scripts are generated from pre-written enzo templates files. 
# - `enzo-templates/templates` contains the `Jinja2` templates files for `enzo`
# 
# ## Compiling Dengo-Enabled Enzo
# 
# - these can be placed directly to the `enzo/src/enzo` directory
# - Dengo-enabled enzo can be built by make dengo-yes grackle-no, below shows the example snippit from make show-config
# ```Make
# CONFIG_GRACKLE  [grackle-{yes,no}]                        : no
# CONFIG_DENGO    [dengo-{yes,no}]                          : yes
# ```
# 
# - Specify paths to various libraries
# ```
# LOCAL_DENGO_INSTALL = {{network._dengo_install_path}}
# LOCAL_CVODE_INSTALL = {{network._cvode_path}}
# LOCAL_SUITESPARSE_INSTALL = {{network._suitesparse_path}}
# ```
# 
# - Compile `Enzo`!
# 
# ```{note}
# A example workflow is documented under https://github.com/hisunnytang/test-enzo-dengo.
# ```
# 
# ```{note}
# Currently it works with the particular commit: `git checkout d84c2415b7914c7bab729c2d818b5af8c85c1918`. Checkout this particular commit for a working integration. We will work to incorporate it to more recent versions of Enzo
# ```
