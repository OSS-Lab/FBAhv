#!/bin/bash

E_BADARG=65
usage="$0 TISSUE_NAME"

# This script performs every operations necessary to obtain a new protein
# biomass function from the RNA expression data from the human protein atlas.

# Download the RNA expression data from the human protein atlas.
# Extract the lines corresponding to its expression in a specific tissue.

if [ -z $1 ]; then
    echo "ERROR: no tissue name provided"
    echo "$usage"
    exit $E_BADARG
fi

export tissue=$1

RNA_expression_file="rna_tissue_consensus.tsv"
RNA_expression_url="https://www.proteinatlas.org/download/${RNA_expression_file}.zip"
specific_rna_expression_file="${tissue}_rna_expression.tsv"
if [ -s $specific_rna_expression_file ]; then
    echo "$specific_rna_expression_file already exists"
else
    wget "$RNA_expression_url"
    unzip "${RNA_expression_file}.zip"
    awk -F $'\t' -e '{if ($3 == ENVIRON["tissue"]) {printf("%s\t%s\t%.1f\n", $1, $2, $4)}}' "$RNA_expression_file" \
       > "$specific_rna_expression_file"
    rm "${RNA_expression_file}.zip"
    if [ -s "$specific_rna_expression_file" ]; then
        echo "$specific_rna_expression_file"
        #rm "$RNA_expression_file"
    else
        echo "ERROR: no \"$tissue\" data in $RNA_expression_file"
        echo "Execution will stop here"
        echo "$RNA_expression_file will be kept so you can check inside"
        exit $E_BADARG
    fi
fi

# Download the uniprot mapping to other databases for the human genome.
# Extract the lines corresponding to the Ensembl database then remove the
# original file for the sake of space.

data_file="HUMAN_9606_idmapping.dat"
uniprot_id_mapping_url="ftp://ftp.uniprot.org/pub/databases/uniprot/current%5Frelease/knowledgebase/idmapping/by_organism/${data_file}.gz"
uniprot_to_ensembl_file="uniprot_to_ensembl.dat"

if [ -s $uniprot_to_ensembl_file ]; then
    echo "$uniprot_to_ensembl_file already exists"
else
    wget "$uniprot_id_mapping_url"
    gzip -d "${datafile}.gz"
    awk -F $'\t' -e '{if ($2 == "Ensembl") {printf("%s\t%s\n", $1, $3)}}' "$data_file" \
       > "$uniprot_to_ensembl_file"
    rm "$data_file"
fi

# Download the whole human proteome uniprot sequences

human_proteome_file_swissprot="UP000005640.xml"
human_proteome_url_swissprot="http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=proteome:UP000005640%20reviewed:yes&fil=&force=yes&format=xml"
if [ -s $human_proteome_file_swissprot ]; then
    echo "$human_proteome_file_swissprot already exists"
else
    wget -O "${human_proteome_file_swissprot}.gz" "$human_proteome_url_swissprot"
    gzip -d "${human_proteome_file_swissprot}.gz"
fi

human_proteome_file_trembl="UP000005640_unreviewed.xml"
human_proteome_url_trembl="http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=proteome:UP000005640%20reviewed:no&fil=&force=yes&format=xml"
if [ -s $human_proteome_file_trembl ]; then
    echo "$human_proteome_file_trembl already exists"
else
    wget -O "${human_proteome_file_trembl}.gz" "$human_proteome_url_trembl"
    gzip -d "${human_proteome_file_trembl}.gz"
fi

# Create the new protein biomass reaction
# NOTE: a lot of file names are defined at the beginning of this script
python create_reaction.py "$tissue"
