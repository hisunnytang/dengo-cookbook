#!/usr/bin/env python
# coding: utf-8

# # Creating a New Chemical Species
# 
# We start by defining individual species.  This can be done inside a python
# module of your choosing, which we will run at the end to output our new
# network.  Species are defined by a small number of attributes:
# 
#  * Name (which will be used in cooling functions and internally to the solver)
#  * Number: Mostly unused except when handling ions.
#  * Atomic weight (in integer AMU)
#  * Number of free electrons that is contributes
# 
# This information is used when calculating things like the contribution of a
# species to the number density.
# We now have three symbolic "species" objects for hydrogen, ionized hydrogen,
# and electrons.  Note that Dengo will happily create ions for species defined in
# the CHIANTI database.

# In[1]:


from dengo.reaction_classes import AtomicSpecies, MolecularSpecies
from dengo.chemical_network import species_registry

HI = AtomicSpecies('H', free_electrons=0)
HII = AtomicSpecies("H", free_electrons=1)
H2I = MolecularSpecies("H2", weight=2.01588, free_electrons=0)
de  = species_registry['de']


# Currently in `dengo`, `AtomicSpecies` and `MolecularSpecies` are the two primary ways to initialize a chemical species. `Dengo.reaction_classes.AtomicSpecies` will try to look for the constituent elements based on the expression from the periodic table, and assign them with appropriate weights. These elements are automatically register to `dengo.reaction_classes.species_registry`. 
# 
# Yeah you have successfully created your `ChemicalSpecies` in `Dengo`!
