import collections
from Bio import SeqIO
import pickle
import cobra
import pickle
from modelAnalysis.info import aaDict, k_atp_protein

def get_tuple(line, sep="\t"):
    return line.rstrip().split(sep)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("tissue",
                        help="name of the tissue (in the protein atlas rna expression file)")
    args = parser.parse_args()

    rna_expression_file = "{}_rna_expression.tsv".format(args.tissue)
    uniprot_ensembl_mapping_file = "uniprot_to_ensembl.dat"
    human_proteome_swissprot_file = "UP000005640.xml"
    human_proteome_trembl_file = "UP000005640_unreviewed.xml"
    human_proteome_seq_dict = "human_proteome_sequences.pickle"
    complete_uniprot_ensembl_file = "complete_uniprot_to_ensembl.pickle"
    weighted_aminoacids_count_file = "weighted_aminoacid_count_{}.pickle".format(args.tissue)
    met_dict_file = "recon2_met_dict.txt"
    model_file = "recon2_fixed.xml"
    protein_reaction_id = "biomass_protein"
    protein_reaction_name = "protein component of biomass (based on RNA expression in {})".format(args.tissue)
    new_biomass_protein_pickle = "biomass_protein_{}.pickle".format(args.tissue)
    tissue_specific_model = "recon2_{}.xml".format(args.tissue)

# Store tissue rna expression in a dictionary mapping ensembl id to normalized
# expression level
    with open(rna_expression_file, "r") as fh:
        tuples = (get_tuple(line) for line in fh)
        ensembl_expression = dict()
        for ensembl_id, gene_name, expression_level in tuples:
            ensembl_expression[ensembl_id] = float(expression_level)

# Read the human proteome file
# NOTE: at the time this script is written, the "swiss"-formatted file
# reading through SeqIO.parse is broken.
# The dictionary is pickled because it takes quite some time to generate
    try:
        with open(human_proteome_seq_dict, "rb") as fh:
            sequences = pickle.load(fh)
    except:
        print("Generating the sequences dict")
        sequences = dict() # uniprot id to aminoacid sequence
        print("Reading Swissprot data")
        with open(human_proteome_swissprot_file, "r") as fh:
            human_proteome_swissprot = SeqIO.UniprotIO.UniprotIterator(fh)
            for record in human_proteome_swissprot:
                for accession_id in record.annotations["accessions"]:
                    sequences[accession_id] = record.seq
        print("Reading TrEMBL data")
        with open(human_proteome_trembl_file, "r") as fh:
            human_proteome_trembl = SeqIO.UniprotIO.UniprotIterator(fh)
            for record in human_proteome_trembl:
                for accession_id in record.annotations["accessions"]:
                    if accession_id not in sequences:
                        sequences[accession_id] = record.seq
        with open(human_proteome_seq_dict, "wb") as fh:
            print("saving it to {}".format(human_proteome_seq_dict))
            pickle.dump(sequences, fh)
        print("done")

    try:
        with open(complete_uniprot_ensembl_file, "rb") as fh:
            ensembl_to_uniprot = pickle.load(fh)
    except:
        print("Generating ensemble-to-uniprot id dictionary")
        # Store ensembl id to uniprot ids (one ensembl id actually points to
        # potentially multiple uniprot ids)
        ensembl_to_uniprot = collections.defaultdict(list)
        with open(uniprot_ensembl_mapping_file, "r") as fh:
            tuples = (get_tuple(line) for line in fh)
            for uniprot_id, ensembl_id in tuples:
                ensembl_to_uniprot[ensembl_id].append(uniprot_id)

        # Store supplementary ensembl_id -> uniprot id by looking at OpenTargets and
        # Bgee crossrefs
        print("Adding supplementary crossreferences")
        print("Reading from Swissprot")
        with open(human_proteome_swissprot_file, "r") as fh:
            human_proteome = SeqIO.UniprotIO.UniprotIterator(fh)
            for record in human_proteome:
                opentargets_crossref = [crossref for crossref in record.dbxrefs
                                        if crossref.startswith("OpenTargets:")]
                bgee_crossref = [crossref for crossref in record.dbxrefs
                                 if crossref.startswith("Bgee:")]
                if opentargets_crossref:
                    dbname, ensembl_id = opentargets_crossref[0].split(":")
                    ensembl_to_uniprot[ensembl_id].extend(record.annotations["accessions"])
                if bgee_crossref:
                    dbname, ensembl_id = bgee_crossref[0].split(":")
                    ensembl_to_uniprot[ensembl_id].extend(record.annotations["accessions"])
        print("Reading from TrEMBL")
        with open(human_proteome_trembl_file, "r") as fh:
            human_proteome = SeqIO.UniprotIO.UniprotIterator(fh)
            for record in human_proteome:
                opentargets_crossref = [crossref for crossref in record.dbxrefs
                                        if crossref.startswith("OpenTargets:")]
                bgee_crossref = [crossref for crossref in record.dbxrefs
                                 if crossref.startswith("Bgee:")]
                if opentargets_crossref:
                    dbname, ensembl_id = opentargets_crossref[0].split(":")
                    ensembl_to_uniprot[ensembl_id].extend(record.annotations["accessions"])
                if bgee_crossref:
                    dbname, ensembl_id = bgee_crossref[0].split(":")
                    ensembl_to_uniprot[ensembl_id].extend(record.annotations["accessions"])
        with open(complete_uniprot_ensembl_file, "wb") as fh:
            print("saving it to {}".format(complete_uniprot_ensembl_file))
            pickle.dump(ensembl_to_uniprot, fh)
        print("done")

# For every ensembl id in the tissue expression levels data,
# get it's protein sequence from uniprot
    print("Compute sum of aminoacids in proteome weighted by RNA expression")
    total_aminoacid_counter = collections.Counter()
    peptide_bonds = 0 # length of proteins weighted their by RNA expression
    nbof_rna = len(ensembl_expression)
    for i, expression_data in enumerate(ensembl_expression.items()):
        ensembl_id, expression = expression_data
        uniprot_ids = ensembl_to_uniprot[ensembl_id]
        # ensembl_id not associated to any uniprot id may exist in uniprot but are
        # not included in the proteome data probably because they are evidenced
        # at transcript level but not at protein level. Those putative proteins
        # are skipped from the aminoacid count
        if uniprot_ids:
            print("{} -> {}: {}/{}".format(ensembl_id, "; ".join(uniprot_ids), i, nbof_rna))
            valid_uniprot_id = next(uniprot_id for uniprot_id in uniprot_ids
                                    if uniprot_id in sequences)
            sequence = sequences[valid_uniprot_id]
            aminoacids_count = collections.Counter(sequence)
            weighted_aminoacid_count = {aa: count * expression for aa, count in aminoacids_count.items()}
            total_aminoacid_counter += weighted_aminoacid_count
            peptide_bonds += len(sequence) * expression
        else:
            print("{}: no protein information: skipped".format(ensembl_id))
    print("done")

    print("saving it to {}".format(weighted_aminoacids_count_file))
    with open(weighted_aminoacids_count_file, "wb") as fh:
        pickle.dump({"aminoacid count": total_aminoacid_counter,
                     "peptide bonds": peptide_bonds}, fh)

    print("Generate the biomass reaction")

# Read the metabolites dictionary
    with open(met_dict_file, "r") as fh:
        met_dict = dict()
        next(fh)
        next(fh)
        for line in fh:
            key, ID, name, chebi_ID = line.split(",")
            met_dict[key] = ID.lstrip()

# Normalize the aminoacid count.
# The reaction's stoichiometry is multiplied by a scaling factor such that the
# mass of aminoacid being consummed is the same as in the original reaction.
    print("Reading the model")
    model = cobra.io.read_sbml_model(model_file)
    print("Normalizing the stoichiometry")
    original_protein_reaction = model.reactions.get_by_id("biomass_protein")
# Create a reverse metabolite dictionary which maps the model's ids to their
# generic names (A, G, H, P...)
    reverse_met_dict = {value: key for key, value in met_dict.items()}
# Compute the molecular weight of the aminoacids in the default protein biomass
# reaction in the model (molecular weight of individual aminoacids weighted by
# their stoichiometric coefficient)
# NOTE: the conditional of this list comprehension excludes reagents of the
# biomass reaction which are not aminoacids (biomass, atp, h2o...)
    original_protein_reaction_weight = sum([-coef * aaDict[reverse_met_dict[met.id]]
                                            for met, coef
                                            in original_protein_reaction.metabolites.items()
                                            if met.id in reverse_met_dict
                                            and reverse_met_dict[met.id] in aaDict])
    original_total_aa_stoich = sum(-coef for met, coef in original_protein_reaction.metabolites.items()
                                   if met.id in reverse_met_dict
                                   and reverse_met_dict[met.id] in aaDict)
    original_atp_stoich = next(coef for met, coef in original_protein_reaction.metabolites.items() if met.id == met_dict["atp"])
    print("cost per aminoacid in original reaction: {} atp/aa".format(original_atp_stoich / original_total_aa_stoich))
# Compute the molecular weight of the aminoacids in the total aminoacid count
# from RNA
# NOTE: "X" corresponds to unidentified aminoacids in uniprot.
# NOTE: "U" corresponds to pyrrolidone carboxylic acid, which does not exist in
# the recon2 model.
    weight_from_count = sum(count * aaDict[aa_id]
                            for aa_id, count
                            in total_aminoacid_counter.items()
                            if aa_id != "X" and aa_id != "U")
    scaling_factor = original_protein_reaction_weight / weight_from_count

    reaction_stoichiometry = {aa: -count for aa, count in total_aminoacid_counter.items()}
# Add ATP, ADP and water to the stoichiometry
    virus_protein_atp_cost = peptide_bonds * k_atp_protein
    reaction_stoichiometry["atp"] = -virus_protein_atp_cost 
    reaction_stoichiometry["h2o"] = -virus_protein_atp_cost 
    reaction_stoichiometry["adp"] = virus_protein_atp_cost 
    reaction_stoichiometry["Pi"] = virus_protein_atp_cost 
    reaction_stoichiometry["h"] = virus_protein_atp_cost 

    normalized_count = {key: value * scaling_factor 
                        for key, value in reaction_stoichiometry.items()
                        if key != "X" and key != "U"}
    # DEBUG
    print("new protein reaction:")
    print({met: coef for met, coef in normalized_count.items()})
    print("cost per aminoacid in new reaction:")
    print(virus_protein_atp_cost / -sum(coef for met, coef in reaction_stoichiometry.items()
                                        if met in aaDict))

# Turn the count dictionary into an input suitable for a
# cobra.core.reaction.Reaction instance.
    print("Instanciating the reaction")
    metabolites_dict = {model.metabolites.get_by_id(met_dict[key]): value
                        for key, value in normalized_count.items()}
    protein_metabolite = model.metabolites.get_by_id("biomass_protein_c")
    metabolites_dict[protein_metabolite] = 1

# Create the new biomass_protein function
    biomass_protein = cobra.core.reaction.Reaction(id=protein_reaction_id,
                                                   name=protein_reaction_name)
    biomass_protein.add_metabolites(metabolites_dict)
    print(biomass_protein)
    print("Saving the new biomass protein reaction as {}".format(new_biomass_protein_pickle))
    with open(new_biomass_protein_pickle, "wb") as fh:
        pickle.dump(biomass_protein, fh)

# Create the updated model
    original_protein_reaction_index = next(i for i, reaction in enumerate(model.reactions)
                                           if reaction.id == original_protein_reaction.id)
    model.reactions[original_protein_reaction_index] = biomass_protein
    print("Saving the new model as {}".format(tissue_specific_model))
    cobra.io.write_sbml_model(model, tissue_specific_model)
