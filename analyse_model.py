# Dependencies: in-built python
import  os
from copy import deepcopy
# Dependencies: required from user
import  cobra
from    cobra   import Model
from Bio import SeqIO
import pandas as pd
#import optlang.gurobi_interface as grb
# Dependencies: modelAnalysis module
from modelAnalysis.info import metDict
from modelAnalysis.genVBOF2 import genVBOF2, get_virus_names
from modelAnalysis.genHVM import genHVM
from modelAnalysis.hvmCompare import hvmCompare
from modelAnalysis.hvmKnockout import hvmKnockout
from modelAnalysis.hvmEnforce import hvmEnforce
from modelAnalysis.hvmDoubleKnockout import hvmDoubleKnockout
from modelAnalysis.hvmDoubleEnforcement import hvmDoubleEnforcement
from modelAnalysis.hvmFluxVar import hvmFluxVar

#######################################################################################
# Function Definition
# analyse_model.py uses the modelAnalysis module to create a host-virus integrated model (HVM),
# compare the host and virus optimised states, and then predict potential anti-viral
# targets

# Inputs:
# VirusGB           User-supplied GenBank file (NCBI) for desired virus
# Model             User-supplied model file (.mat, .sbml)
# HostRxn           Name of the host objective reaction

# Outputs:
# hvm               Integrated host-virus model
# hvmComp           Comparison of HOS and VOS flux values (as % total flux)
# koVirus           Vector of virus optima, for HVM with knockout constraints
# enfVirus          Vector of virus optima, for HVM with additional host-constraint

def cobra_model(file_path):
    """Open a metabolic model.
    This function is used by argparse for type validation."""
    if ".mat".lower() in file_path.lower():
        return cobra.io.load_matlab_model(file_path)
    elif ".xml".lower() in file_path.lower():
        return cobra.io.read_sbml_model(file_path)
    elif ".json".lower() in file_path.lower():
        return cobra.io.read_json_model(file_path)
    elif ".yaml".lower() in file_path.lower():
        return cobra.io.read_yaml_model(file_path)
    else:
        raise ValueError('Unable to load model: unsupported file type. Consult _README')

def genbank_file(file_path):
    """Open a genbank file.
    This function is used by argparse for type validation."""
    return next(SeqIO.parse(file_path, "genbank"))

def medium_file(file_path):
    """Open a medium file (headerless tsv, 2 columns (str, float))
    Return a generator of (str, float) tuples
    This function is used by argparse for type validation."""
    def row_generator(file_path):
        with open(file_path, "r") as fh:
            for line in fh:
                reaction_id, lower_bound, upper_bound = line.rstrip().split("\t")
                yield (reaction_id, float(lower_bound), float(upper_bound))
    return row_generator(file_path)

def get_hvm(VirusGB, Model, model_name, HostRxn):
    """Create a new metabolic network model (HVM, Host Virus Model) by adding a
    virus biomass function to the metabolic network model of a host cell.
    The objective of the HVM is set to be the host biomass synthesis reaction.
    
    parameters:
    - VirusGB: genbank file for the virus (generator of Bio.SeqRecord.SeqRecord
      as returned by Bio.SeqIO.parse on a genbank file)
    - Model: metabolic network model for the host (cobra.core.model.Model)
    - model_name: name of the model (str). This argument should be an empty
      string if Model has a name attribute.
    - hostRxn: id of the host's biomass synthesis reaction in Model (str)
    
    return:
    a 2-tuple containing:
    - the virus reaction (cobra.core.reaction.Reaction)
    - the Host Virus Model (cobra.core.model.Model)"""
    print("creating virus biomass function....")
    # Genereate the virus biomass objective function    (genVBOF.py)
    virusRxn = genVBOF2(VirusGB, Model, model_name=model_name)

    print("creating host-virus integrated model....\n")
    # Create the host-virus integrated model            (genHVM.py)
    HVM = genHVM(Model,virusRxn)
    hvm_host_biomass_reaction = HVM.reactions.get_by_id(HostRxn)
    HVM.objective = hvm_host_biomass_reaction
    return (virusRxn, HVM)

def analyse_model(VirusGB, Model, HostRxn, solver, model_name=None, prefix="", overwrite=False, defined_medium=None):
    """Create a HVM model from a host model and a virus genbank file (it works
    only with SARS-CoV-2 for the moment) and runs FBAhv model analysis on it.
    
    parameters:
    - VirusGB: genbank file for the virus (generator of Bio.SeqRecord.SeqRecord
      as returned by Bio.SeqIO.parse on a genbank file)
    - Model: metabolic network model for the host (cobra.core.model.Model)
    - hostRxn: id of the host's biomass synthesis reaction in Model (str)
    - solver: name of the solver to use (supported: glpk, gurobi) (str)
    - model_name: name of the model (str). This argument should be an empty
      string if Model has a name attribute. It is only necessary when the
      model is imported from a matlab file which has no name attribute.
    - prefix: prefix for the name of the output files to write (str)
    - overwrite: whether or not existing files should be overwritten (bool)
    - defined_medium: iterable of 3-tuples being the reaction id, the lower
      and the upper boundary for that reaction. If this parameter is not
      specified (ie left as None) the exchange reactions boundaries are set to
      (-1000, 1000). If this parameter is set however, the bounds of every
      exchange reactions not defined in the medium are set to (0, 1000) (ie no
      absorbtion) then the bounds of the reactions defined by this
      parameter are set to what they are defined in defined_medium.
      (Iterable[(str, float, float)]) """

    short_name, full_name = get_virus_names(VirusGB)
    print("Virus identified as \"{}\"".format(full_name))
    print("Abbreviated as \"{}\"".format(short_name))

    # importing the solver
    if solver == "gurobi":
        import optlang.gurobi_interface as grb

    # make sure the objective function is the model's biomass function
    host_biomass_reaction = next(reaction for reaction in Model.reactions
                                 if reaction.id == HostRxn)
    Model.objective = host_biomass_reaction

    # enforce virus medium
    if defined_medium:
        for reaction in Model.exchanges:
            Model.reactions.get_by_id(reaction.id).bounds = (0, 1000)
        for reaction_id, lower_bound, upper_bound in defined_medium:
            Model.reactions.get_by_id(reaction_id).bounds = (lower_bound, upper_bound)
    else:
        for reaction in Model.exchanges:
            Model.reactions.get_by_id(reaction.id).bounds = (-1000, 1000)

    print("creating host-virus integrated model....\n")
    # Create the host-virus integrated model            (genHVM.py)
    virusRxn, HVM = get_hvm(VirusGB, Model, model_name, HostRxn)
    HVM.objective = virusRxn
    HVM_file = ''.join([str(Model), "_virus.xml"])
    cobra.io.write_sbml_model(HVM, HVM_file)

    #output handling
    virus_rxn_file = ''.join([prefix, "results_virusRxn_", short_name])
    f = open(virus_rxn_file,"w+")
    f.write("Virus Reaction is:\n%s" % (virusRxn.reaction))
    f.close()
    print("Virus reaction saved in {}".format(virus_rxn_file))

    print("performing comparison among host- vs. virus-optimal model....\n")
    # Perform the comparison analysis                   (hvmCompare.py)
    (outputDf_flux,outputDf_summary) = hvmCompare(deepcopy(HVM),HostRxn,solver)
    flux_comp_file = ''.join([prefix,"results_hvmFluxComp_",solver,"_",short_name])
    outputDf_flux.to_csv(flux_comp_file, index=False)
    flux_comp_file = ''.join([prefix,"results_hvmFluxSummary_",solver,"_",short_name])
    outputDf_summary.to_csv(flux_comp_file, index=False)
    print("...just finished flux comparison and have put out results on file\n")

    # Perform the knockout analysis                     (hvmKnockout.py)
    ko_file = ''.join([prefix,"results_hvmKnockOut_", solver, "_", short_name])
    if os.path.isfile(ko_file) and not overwrite:
        print("The knockout file ({}) already exists: reading from it\n".format(ko_file))
        koVirus = pd.read_csv(ko_file)
    else:
        print("performing knockout analysis....\n")
        koVirus  = hvmKnockout(deepcopy(HVM), HostRxn, solver)
        koVirus.to_csv(ko_file, index=False)
    
    print("...just finished knockout analysis and have put out results on file\n")

    # Perform the double knockout analysis              (hvmDoubleKnockout.py)
    host_fva_file = ''.join([prefix,"results_hvmFluxVarHost_",solver,"_",short_name])
    virus_fva_file = ''.join([prefix,"results_hvmFluxVarVirus_",solver,"_",short_name])
    if os.path.isfile(host_fva_file) and os.path.isfile(host_fva_file) and not overwrite:
        print("The FVA files ({}; {}) already exist: reading from them\n".format(host_fva_file, virus_fva_file))
        fluxVarHost = pd.read_csv(host_fva_file)
        fluxVarVirus = pd.read_csv(virus_fva_file)
    else:
        print("performing FVA analysis....\n")
        (fluxVarHost,fluxVarVirus) = hvmFluxVar(deepcopy(HVM),HostRxn,solver)
        fluxVarHost.to_csv(host_fva_file, index=False)
        fluxVarVirus.to_csv(virus_fva_file, index=False)
        print("...just finished FVA analysis and have put out results on file\n")

    # Perform the host-derived enforcement analysis     (hvmEnforce.py)
    enforce_file = ''.join([prefix,"results_hvmEnforce_",solver,"_",short_name])
    if os.path.isfile(enforce_file) and not overwrite:
        print("The enforcement analysis file ({}) already exists".format(enforce_file))
    else:
        print("performing flux enforcement analysis....\n")
        enfVirus = hvmEnforce(deepcopy(HVM),HostRxn,fluxVarHost,fluxVarVirus,solver)
        enfVirus.to_csv(enforce_file, index=False)
        print("...just finished flux enforcement analysis and have put out results on file\n")

    #Perform the double knockout analysis              (hvmDoubleKnockout.py)
    double_ko_file = ''.join([prefix,"results_hvmDoubleKO_",solver,"_",short_name])
    if os.path.isfile(double_ko_file) and not overwrite:
        print("The double KO file ({}) already exists")
    else:
        print("performing double knockout analysis....\n")
        doubleKOvirus = hvmDoubleKnockout(deepcopy(HVM),HostRxn,koVirus,solver)
        doubleKOvirus.to_csv(double_ko_file, index=False)
        print("...just finished double knockout analysis and have put out results on file\n")

    # Perform the double knockout analysis              (hvmDoubleEnforcement.py)
    double_enforce_file = ''.join([prefix,"results_hvmDoubleEnf_",solver,"_",short_name])
    if os.path.isfile(double_enforce_file) and not overwrite:
        print("The double enforcement file ({}) already exists".format(double_enforce_file))
    else:
        print("performing double enforcement analysis....\n")
        doubleEnfvirus = hvmDoubleEnforcement(deepcopy(HVM),HostRxn,koVirus,fluxVarHost,fluxVarVirus,solver)
        doubleEnfvirus.to_csv(double_enforce_file, index=False)
        print("...just finished double enforcement analysis and have put out results on file\n")
    
    print("\n\nall done Sir!")
    # Create outputs
    return (virusRxn,HVM,hvmComp,koVirus,enfVirus,doubleKOvirus)

#main function call
#virus genbank file, metabolic model, name of biomass function for host in the model, linear optimisation solver
#  'glpk'                               #default
#  'gurobi'                             #default
#  'glpk_exact'                         #not tested
#   'cplex'                             #not tested

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Command line interface for FBAhv")
    parser.add_argument("genbank",
                        help="genbank file for a virus",
                        type=genbank_file)
    parser.add_argument("model",
                        help="model for a host (.mat, .xml, .yaml, .json)",
                        type=cobra_model)
    parser.add_argument("hostRxn",
                        help="name of the biomass production reaction in the model")
    parser.add_argument("--solver",
                        help="solver to use for computation",
                        choices=["glpk", "gurobi"],
                        default="glpk")
    parser.add_argument("--model-name",
                        help="""provide this parameter only if the model being
inputed has no `name` attribute.""",
                        default=None)
    parser.add_argument("--prefix",
                        help="prefix for the output files",
                        default="")
    parser.add_argument("--overwrite",
                        action="store_true",
                        help="if result files already exists, perform the "
                             "analysis anyway and ovewrite them")
    parser.add_argument("--medium",
                        type=medium_file,
                        help="path to a tsv medium file which specifies which "
                             "exchange reaction should be kept in the model "
                             "and with which boundaries "
                             "(lower bound: intake, upper bound: excretion). "
                             "NOTE: if this parameter is not set, the "
                             "boundaries of exchange reactions will be set "
                             "to (-1000, 1000). If this "
                             "parameter is set however, the upper bound of "
                             "each exchange reaction not defined through this "
                             "parameter will be set to zero.")
    args = parser.parse_args()
    analyse_model(args.genbank,
                  args.model,
                  args.hostRxn,
                  args.solver,
                  model_name=args.model_name,
                  prefix=args.prefix,
                  overwrite=args.overwrite,
                  defined_medium=args.medium)
    # test: `python FBAhv.py SARS_CoV2.txt recon2_lung.xml biomass_reaction`
