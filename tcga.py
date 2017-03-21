"""Obtain mutation data from TCGA

Raw data files were obtained from https://gdc-portal.nci.nih.gov/search/f?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22%3C%3D%22,%22content%22:%7B%22field%22:%22cases.diagnoses.age_at_diagnosis%22,%22value%22:%5B7305%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_category%22,%22value%22:%5B%22Simple%20Nucleotide%20Variation%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D%5D%7D&pagination=%7B%22files%22:%7B%22from%22:0,%22size%22:20,%22sort%22:%22data_category:desc,%22%7D%7D&facetTab=files
Data include all publicly-available simple nucleotide variation information on the website. Descriptions of raw data
files are available at https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification

Only the following mutation types are extracted: "Missense_Mutation", "Silent", "Frame_Shift_Del", "In_Frame_Del",
"Frame_Shift_Ins", "In_Frame_Ins". Extracted data are non-repetitive. Two mutations are collapsed to a single entry if
they refer to the same change in amino acid.
"""

import os

# make tuple of mutation types that we want to extract
mutation_types = ("Missense_Mutation", "Silent", "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins",
                  "In_Frame_Ins")
raw_files = os.listdir("../Mutation Data/TCGA")
fwrite = open("TCGAMutationData.txt", "w")
fwrite.write("GeneName\tUniProtID\tProteinMutation\tMutationType\n")

prev_mutations = set()
write_data = []
for filename in raw_files:
    tcga_data = open("../Mutation Data/TCGA/" + filename)
    next(tcga_data)
    next(tcga_data)  # skip first 2 lines

    for line in tcga_data:
        line = line.rstrip("\n").split("\t")
        mutation_type = line[8]
        gene_name = line[60]
        mutation = line[36][2:]
        if mutation_type not in mutation_types or (gene_name, mutation) in prev_mutations:
            # go to next mutation if we're not considering this one
            continue
        prev_mutations.add((gene_name, mutation))
        uniprot_id = line[67]
        new_data = [gene_name, uniprot_id, mutation, mutation_type]
        for x in range(len(new_data)):
            if not new_data[x]:
                new_data[x] = "None"
        write_data.append(new_data)

write_data.sort(key=lambda x: x[3])  # sort data by mutation type
for line in write_data:
    fwrite.write("\t".join(line) + "\n")
fwrite.close()
