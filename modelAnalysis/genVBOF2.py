from cobra import Reaction, Metabolite
from collections import Counter
from modelAnalysis.info import ntpsDict, aaDict, N_A, k_atp_protein, k_atp_rna, k_ppi

#######################################################################################
# Function Setup

nbof_spikes = 90
# NOTE: the keys in this dictionary correspond to the name of the product in the genbank file
virus_composition = {
    "GU28": { # COVID19
# Sources for all numbers below are: 
# Source: Neuman, B.W., Adair, B. D., et.al. J. of Virology (2006, 80:16) -> 8M2:4N:1S3 stoichiometry per spike and 75 spike particles
# Source: Neuman, B.W., Joseph, J. S., et.al. J. of Virology (2008, 82:11).
# Source: Neuman, B.W., Kiss G., et.al. J. of Structural Biology (2011, 174) -> # 8M2:4N:1S3 stoichiometry per spike and 90 spike particles
# Source: Escors D., Camafeita E, et.al. J. of Virology (2001, 75:24) -> M:N from whole virion is 3:1
# Source: Barcena, M., Oostergetel, G. T., et.al. PNAS (2008, 106:2) 
        # Copy number for viral genome                      [Cg]                
        "Cg": 1,
        "proteins": {
        # Copy numbers for structural proteins N, M, S, and E
        # Copy number for viral M protein
        "membrane glycoprotein": 16*nbof_spikes+100,
        # Copy number for viral Spike protein - based on 8M2:4N:1S3 stoichiometry per spike and 75 spike particles       
        "surface glycoprotein": 3*nbof_spikes,
        # Copy number for viral N protein - based on 8M2:4N:1S3 stoichiometry per spike and 75 spike particles               
        "nucleocapsid phosphoprotein": 4*nbof_spikes+130,
        # Copy number for viral Envelope protein - no good estimate, but said to be low         
        "envelope protein": 5,
        # Copy numbers for nonstructural proteins pp1a and pp1ab
        # Copy number for viral nonstructural polyprotein pp1a [Cnp_1a] - no good estimate, but said to be low       
        "orf1a polyprotein": 5,
        # Copy number for viral nonstructural polyprotein pp1ab [Cnp_1ab] - no good estimate, but said to be low  
        "orf1ab polyprotein": 5,
        # Copy numbers for accessory proteins 
        # Copy number for viral accessory polyprotein ORF3a  [Cnp] - no good estimate, but said to be low 
        "ORF3a protein": 5,
        # Copy number for viral accessory polyprotein ORF6  [Cnp] - no good estimate, but said to be low    
        "ORF6 protein": 5,
        # Copy number for viral accessory polyprotein ORF7a  [Cnp] - no good estimate, but said to be low    
        "ORF7a protein": 5,
        # Copy number for viral accessory polyprotein ORF7b  [Cnp] - no good estimate, but said to be low    
        "ORF7b": 5,
         # Copy number for viral accessory polyprotein ORF8  [Cnp] - no good estimate, but said to be low    
        "ORF8 protein": 5,
         # Copy number for viral accessory polyprotein ORF10  [Cnp] - no good estimate, but said to be low    
        "ORF10 protein": 5}}
}

def multiply_counter(counter, n):
    """return a collection.Counter whose values have been multiplied by a
    scalar"""
    # the multiplicand does not need to be a collections.Counter, but at least
    # it needs to be a dict
    assert isinstance(counter, dict)
    assert float(n)
    return Counter({key: value * n for key, value in counter.items()})

def get_virus_names(virus_record):
    """Return a tuple containing a short name and a
    full name for the virus. This is an alternative to how it is done in Sean's
    genVBOF."""
    long_name = virus_record.description.rstrip(", complete genome")
    # try to get a short name for the virus by taking the prefix of the first CDS
    short_name = next(feature for feature in virus_record.features if feature.type == "CDS").qualifiers["locus_tag"][0].rstrip("_gp01")
    return (short_name, long_name)

model_name_to_met_dict_file = {
    "Swainston2016 - Reconstruction of human metabolic network (Recon 2.2)": "met_dicts/recon2_met_dict.txt",
    "Genome-scale metabolic model for hepatocytes, iHepatocytes2322": "met_dicts/iHepatocytes2322_met_dict.txt",
    "macrophage_model": "met_dicts/mac_met_dict.txt",
    "recon 2_LarsNielsen model": "met_dicts/recon2_LarsNielsen_met_dict.txt",
    "GSM_human model": "met_dicts/GSM_human_met_dict.txt", "lung model": "met_dicts/lung_met_dict.txt"
}

def load_metabolite_id_dict(model, model_name=None):
    """Provide a dictionary mapping generic metabolic name (a string like "A",
    "atp", "h2o", "h" or "PPi") to the metabolite object in the model.
    
    The model_name parameter allows to indicate the name of the model in case
    the model object has no name attribute. Do not set this parameter if the
    model object already has a name."""

    name = model_name or model.name
    if name in model_name_to_met_dict_file:
        met_dict_file = model_name_to_met_dict_file[name]
    else:
        raise NotImplementedError("This model is not covered: \"{}\"".format(name))

    with open(met_dict_file, "r") as fh:
        met_tuples = (line.split(", ")[:2] for line in fh.readlines()[2:])
    met_dict = {key: model.metabolites.get_by_id(met_id) for key, met_id in met_tuples}
    return met_dict

#######################################################################################
# Function Definition
# genVOBF.py takes user-supplied viral genome file (.gb) and creates a
# biomass objective function that is characterises the genomic, proteomic and
# energy requirements for the production of virus particles

# Inputs:
# VirusGB           User-supplied GenBank file (NCBI) for desired virus

# Outputs:
# VBOF              Virus biomass objective function for desired virus

def genVBOF2(virus_record, model, model_name=None):
    """New version of the genVBOF function by Hadrien.
    Builds a Virus Biomass Objective Function (basically a virus biomass
    production reaction, from aminoacids and nucleotides) from a genbank
    file.
    
    Params:
    - virus_record: genbank record of a virus (output from Bio.SeqIO.parse)
    - model: a cobra metabolic model (cobra.core.model.Model)

    Returns:
    - virus biomass objective function (cobra.core.reaction.Reaction)
    """
    met_dict = load_metabolite_id_dict(model, model_name=model_name)

    # VIRUS IDENTIFICATION
    taxonomy = " ".join([taxon.lower() for taxon in virus_record.annotations["taxonomy"]])
    if "betacoronavirus" not in taxonomy:
        raise NotImplementedError('Virus family is not supported: Unable to create VBOF. Consult _README')
    short_name, full_name = get_virus_names(virus_record)

    # AMINOACID COUNT
    all_cds = {feature for feature in virus_record.features if feature.type == "CDS"}
    # Check that our own virus_composition dict contain exactly the
    # proteins defined in the genbank file, no more, no less.
    protein_names_in_gb_file = {cds.qualifiers["product"][0] for cds in all_cds}
    protein_names_in_our_data = {protein_name for protein_name in virus_composition[short_name]["proteins"]}
    assert protein_names_in_gb_file == protein_names_in_our_data

    virus_aa_composition = Counter()
    # protein name -> number of atp involved in its peptide bonds formations
    # (accounting for the number of copies of protein)
    peptide_bond_formation = dict()
    for cds in all_cds:
        protein_name = cds.qualifiers["product"][0]
        aa_sequence = cds.qualifiers["translation"][0]
        aa_count = Counter(aa_sequence)
        copies_per_virus = virus_composition[short_name]["proteins"][protein_name]
        virus_aa_composition += multiply_counter(aa_count, copies_per_virus)
        peptide_bond_formation[protein_name] = (len(aa_sequence) * k_atp_protein - k_atp_protein) * copies_per_virus

    # [3] Precursor frequency
    # Genome                            [Nucleotides]
    Cg = virus_composition[short_name]["Cg"] # number of genome copies per virus
    virus_nucl_count = Counter(str(virus_record.seq))
    countA  = virus_nucl_count["A"]
    countC  = virus_nucl_count["C"]
    countG  = virus_nucl_count["G"]
    countU  = virus_nucl_count["T"]    # Base 'T' is pseudo for base 'U'
    antiA   = countU
    antiC   = countG
    antiG   = countC
    antiU   = countA
    # Count summation
    totNTPS     = (Cg * (countA + countC + countG + countU + antiA + antiC + antiG + antiU))
    totAA       = sum(count for count in virus_aa_composition.values())

    # [4] VBOF Calculations
    # Nucleotides
    # mol.ntps/mol.virus
    V_a = (Cg*(countA + antiA))
    V_c = (Cg*(countC + antiC))
    V_g = (Cg*(countG + antiG))
    V_u = (Cg*(countU + antiU))
    # g.ntps/mol.virus
    G_a = V_a * ntpsDict["atp"]
    G_c = V_c * ntpsDict["ctp"]
    G_g = V_g * ntpsDict["gtp"]
    G_u = V_u * ntpsDict["ttp"]

    # Amino Acids
    # g.a/mol.virus
    G_aa = {aa: count * aaDict[aa] for aa, count in virus_aa_composition.items()}
    # Total genomic and proteomic molar mass
    M_v     = (G_a + G_c + G_g + G_u) + sum(G_aa.values())

    # Stoichiometric coefficients
    # Nucleotides [mmol.ntps/g.virus] (for the genome)
    S_atp = 1000 * (V_a/M_v)
    S_ctp = 1000 * (V_c/M_v)
    S_gtp = 1000 * (V_g/M_v)
    S_utp = 1000 * (V_u/M_v)

    # Amino acids [mmol.aa/g.virus]
    S_aa = {aa: 1000 * V_aa / M_v for aa, V_aa in virus_aa_composition.items()}

    # Energy requirements
    # Genome: Phosphodiester bond formation products [Pyrophosphate]
    # SARS Cov 2 is a single stranded RNA virus: it has to first do an
    # intermediary reverse copy of itself and then replicate itself from
    # that intermediary strand.
    genTemp = (((countA + countC + countG + countU) * k_ppi) - k_ppi)
    genRep  = (((antiA + antiC + antiG + antiU) * k_ppi) - k_ppi)
    genTot  = genTemp + genRep
    V_ppi   = genTot
    S_ppi   = 1000 * (V_ppi / M_v)

    # Proteome: Peptide bond formation [ATP + H2O]
    # Note: ATP used in this process is denoated as ATPe/Ae [e = energy version]
    V_Ae = sum(peptide_bond_formation.values())
    S_Ae = 1000 * (V_Ae / M_v)
 
    # [5] VBOF Reaction formatting and output
    # Left-hand terms: Nucleotides
    # Note: ATP term is a summation of genome and energy requirements
    S_ATP   = (S_atp + S_Ae) * -1
    S_CTP   = S_ctp * -1
    S_GTP   = S_gtp * -1
    S_UTP   = S_utp * -1
 
    # Left-hand terms: Amino Acids
    S_AAf = {aa: -coef for aa, coef in S_aa.items()}
    # Left-hand terms: Energy Requirements
    S_H2O   = S_Ae * -1
    # Right-hand terms: Energy Requirements
    S_ADP   = S_Ae
    S_Pi    = S_Ae
    S_H     = S_Ae
    S_PPi   = S_ppi

    reaction_name       = short_name + '_prodrxn_VN'
    virus_reaction      = Reaction(reaction_name)
    virus_reaction.name = full_name + ' production reaction'
    virus_reaction.subsystem                = 'Virus Production'
    virus_reaction.lower_bound              = 0
    virus_reaction.upper_bound              = 1000

    virus_reaction.add_metabolites(({
        met_dict['atp']: S_ATP,
        met_dict['ctp']: S_CTP,
        met_dict['gtp']: S_GTP,
        met_dict['utp']: S_UTP,
        met_dict['A']: S_AAf['A'],
        met_dict['R']: S_AAf['R'],
        met_dict['N']: S_AAf['N'],
        met_dict['D']: S_AAf['D'],
        met_dict['C']: S_AAf['C'],
        met_dict['Q']: S_AAf['Q'],
        met_dict['E']: S_AAf['E'],
        met_dict['G']: S_AAf['G'],
        met_dict['H']: S_AAf['H'],
        met_dict['I']: S_AAf['I'],
        met_dict['L']: S_AAf['L'],
        met_dict['K']: S_AAf['K'],
        met_dict['M']: S_AAf['M'],
        met_dict['F']: S_AAf['F'],
        met_dict['P']: S_AAf['P'],
        met_dict['S']: S_AAf['S'],
        met_dict['T']: S_AAf['T'],
        met_dict['W']: S_AAf['W'],
        met_dict['Y']: S_AAf['Y'],
        met_dict['V']: S_AAf['V'],
        met_dict['h2o']: S_H2O,
        met_dict['adp']: S_ADP,
        met_dict['Pi']:  S_Pi,
        met_dict['h']:   S_H,
        met_dict['PPi']: S_PPi}))
    return virus_reaction
