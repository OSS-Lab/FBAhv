# Description

This project contains scripts allowing to add a "virus biomass function" to a
cell metabolic model (a SBML model), and then to perform an analysis of this
"Host-Virus Model" (HVM).

The scripts allowing to create the HVM and to analyse it are stored in the
modelAnalysis module. The `analyse_model.py` script in the root directory of
this project consists in a command-line interface to this module. An example
of invokation (using the Recon2.2 model) would be

> `python analyse_model.py SARS_CoV_2.txt recon2.xml biomass_reaction`

where "SARS_CoV_2.txt" is the genbank file for SARS_CoV_2, "recon2.xml" is
the SBML file for the Recon2.2 model and "biomass_reaction" is the id of the
biomass reaction in the Recon2.2 model. For informations about the other
options of the command-line interface of `analyse_model.py`, refer to its
documentation.

The `minimalMediaCreation` directory contains a script to determine the
composition of the minimal medium allowing growth for a given model.

The scripts `create_reaction.py` and `get_tissue-specific_stoichiometry.sh` in
the root directory are used to produce modified versions of the Recon2.2 model
where the `biomass_protein` reaction is modified such that its stoichiometry
reflects the RNA expression data from a specific tissue, from the human
protein atlas project (https://www.proteinatlas.org/). In order to use it,
use the `get_tissue-specific_stoichiometry.sh` script.

# License informations

The analysis codes provided here are “AS IS" and covered under a BSD2 Licence.
They are intended solely for non-commercial, academic use in the hope that they
will be useful for analysing other host-virus pairs in the context of infection
or viral ecology. If you would like to explore application of the provided code
and analysis approach in a commercial setting, please get in touch with the
project developers or University Warwick tech transfer office at
“ventures@warwick.ac.uk” quoting reference ‘host-virus metabolic modelling
