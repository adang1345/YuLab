"""Given cancer mutation data from COSMIC and TCGA, create a data file that organizes the data.

COSMIC mutation data were obtained from http://cancer.sanger.ac.uk/cosmic/download and are in
    ../Mutation Data/CosmicMutantExport.tsv
The mapping between gene name and UniProt ID was obtained from http://www.genenames.org/cgi-bin/download and is in
    ../Mutation Data/CosmicUniProt.txt

Data from COSMIC contain many types of mutations. I extract only the following types of mutations:
"Deletion - Frameshift"
"Deletion - In frame"
"Insertion - Frameshift"
"Insertion - In frame"
"Substitution - Missense"
"Substitution - coding silent"

TCGA raw data files were obtained from https://gdc-portal.nci.nih.gov/search/f?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22%3C%3D%22,%22content%22:%7B%22field%22:%22cases.diagnoses.age_at_diagnosis%22,%22value%22:%5B7305%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_category%22,%22value%22:%5B%22Simple%20Nucleotide%20Variation%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D%5D%7D&pagination=%7B%22files%22:%7B%22from%22:0,%22size%22:20,%22sort%22:%22data_category:desc,%22%7D%7D&facetTab=files
Data include all publicly-available simple nucleotide variation information on the website. Descriptions of raw data
files are available at https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification

Only the following mutation types are extracted: "Missense_Mutation", "Silent", "Frame_Shift_Del", "In_Frame_Del",
"Frame_Shift_Ins", "In_Frame_Ins".

Exported data are non-repetitive in the following sense: If two entries from the raw data refer to the same amino acid
change, then they are considered the same mutation. Also, any ambiguous mutations (for which there is a question mark in
the mutation specified) are not extracted."""

import os

# extract COSMIC data

# create mapping between gene name and UniProt ID for COSMIC
uniprot_file = open("../Mutation Data/COSMIC/CosmicUniProt.txt")
next(uniprot_file)
cosmic_to_uniprot = {}
for line in uniprot_file:
    line = line.split()
    if len(line) == 2:
        cosmic_to_uniprot[line[0]] = line[1]
    else:
        cosmic_to_uniprot[line[0]] = "-"
uniprot_file.close()

cosmic_data = open("../Mutation Data/COSMIC/CosmicMutantExport.tsv")
next(cosmic_data)
fwrite = open("CancerMutationData.txt", "w")
fwrite.write("GeneName\tUniProtID\tMutationID\tProteinMutation\tMutationType\tSource\n")

write_data = []
prev_mutations = set()
relevant_mutations = ("Deletion - Frameshift", "Deletion - In frame", "Insertion - Frameshift", "Insertion - In frame",
                      "Substitution - Missense", "Substitution - coding silent", "Missense_Mutation", "Silent",
                      "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins", "In_Frame_Ins")
mutation_type_convert = {"Deletion - Frameshift":           "Del_Frameshift",
                         "Deletion - In frame":             "Del_Inframe",
                         "Insertion - Frameshift":          "Ins_Frameshift",
                         "Insertion - In frame":            "Ins_Inframe",
                         "Substitution - Missense":         "Missense",
                         "Substitution - coding silent":    "Silent",
                         "Missense_Mutation":               "Missense",
                         "Silent":                          "Silent",
                         "Frame_Shift_Del":                 "Del_Frameshift",
                         "In_Frame_Del":                    "Del_Inframe",
                         "Frame_Shift_Ins":                 "Ins_Frameshift",
                         "In_Frame_Ins":                    "Ins_Inframe"}
source = "COSMIC"
for line in cosmic_data:
    line = line.rstrip("\n").split("\t")
    gene_name = line[0]
    if "_" in gene_name:
        gene_name = gene_name[:gene_name.index("_")]
    mutation_type = line[19]
    protein_mutation = line[18][2:]
    if mutation_type not in relevant_mutations or (gene_name, protein_mutation) in prev_mutations or "?" in protein_mutation:
        # If this mutation has an unwanted type or if it has been encountered before, then ignore it.
        continue
    # invariant: this mutation has not been previously encountered
    prev_mutations.add((gene_name, protein_mutation))
    mutation_type = mutation_type_convert[mutation_type]
    mutation_id = line[16]
    try:
        uniprot_id = cosmic_to_uniprot[gene_name]
    except KeyError:
        uniprot_id = "-"
    new_data = [gene_name, uniprot_id, mutation_id, protein_mutation, mutation_type, source]
    for x in range(len(new_data)):
        if not new_data[x]:
            new_data[x] = "-"
    write_data.append(new_data)


# extract TCGA data
raw_files = os.listdir("../Mutation Data/TCGA")
mutation_id = "-"
source = "TCGA"

for filename in raw_files:
    tcga_data = open("../Mutation Data/TCGA/" + filename)
    next(tcga_data)
    next(tcga_data)  # skip first 2 lines

    for line in tcga_data:
        line = line.rstrip("\n").split("\t")
        mutation_type = line[8]
        gene_name = line[60]
        protein_mutation = line[36][2:]
        if mutation_type not in relevant_mutations or (gene_name, protein_mutation) in prev_mutations or "?" in protein_mutation:
            # go to next mutation if we're not considering this one
            continue
        prev_mutations.add((gene_name, protein_mutation))
        mutation_type = mutation_type_convert[mutation_type]
        uniprot_id = line[67]
        new_data = [gene_name, uniprot_id, mutation_id, protein_mutation, mutation_type, source]
        for x in range(len(new_data)):
            if not new_data[x]:
                new_data[x] = "-"
        write_data.append(new_data)

write_data.sort(key=lambda x: x[4])  # sort data by mutation type
for line in write_data:
    fwrite.write("\t".join(line) + "\n")
fwrite.close()
