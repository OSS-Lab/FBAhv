from copy import deepcopy

# Function Definition
# genHVM.py takes a user-supplied model file and creates an intergrated host-
# virus model, given a VBOF

# Inputs:
# Model             User-supplied model (cobra.io.core.model.Model instance)
# VBOF              Virus biomass objective function created by genVBOF.py

# Outputs:
# hvm               Integrated host-virus model

def genHVM(Model,VBOF):
    "Generate_HVM"

    # Integrate the VBOF Reaction
    HVM = deepcopy(Model)
    HVM.add_reaction(VBOF)

    # Outputs
    return HVM

