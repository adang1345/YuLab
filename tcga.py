"""Obtain mutation data from TCGA

Raw data files were obtained from https://gdc-portal.nci.nih.gov/search/f?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22%3C%3D%22,%22content%22:%7B%22field%22:%22cases.diagnoses.age_at_diagnosis%22,%22value%22:%5B7305%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_category%22,%22value%22:%5B%22Simple%20Nucleotide%20Variation%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D%5D%7D&pagination=%7B%22files%22:%7B%22from%22:0,%22size%22:20,%22sort%22:%22data_category:desc,%22%7D%7D&facetTab=files
Data include all publicly-available simple nucleotide variation information on the website. Descriptions of raw data
files are available at https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification

Only the following mutation types are extracted: "Missense_Mutation", "Silent", "Frame_Shift_Del", "In_Frame_Del",
"Frame_Shift_Ins", "In_Frame_Ins".
"""

import os

# make tuple of mutation types that we want to extract
mutation_types = ("Missense_Mutation", "Silent", "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins",
                  "In_Frame_Ins")
raw_files = os.listdir("../Mutation Data/TCGA")
fwrite = open("TCGAMutationData.txt", "w")
fwrite.write("GeneName\tUniProtID\tMutation\tMutationType\n")

for curr_mutation in mutation_types:
    for filename in raw_files:
        tcga_data = open("../Mutation Data/TCGA/" + filename)
        next(tcga_data)
        next(tcga_data)  # skip first 2 lines

        for line in tcga_data:
            line = line.rstrip("\n").split("\t")
            mutation_type = line[8]
            if mutation_type != curr_mutation:  # go to next mutation if we're not looking for this one currently
                continue
            gene_name = line[60]
            uniprot_id = line[67]
            mutation = line[36][2:]
            new_data = [gene_name, uniprot_id, mutation, mutation_type]
            for x in range(len(new_data)):
                if not new_data[x]:
                    new_data[x] = "None"
            fwrite.write("\t".join(new_data) + "\n")
        print("Finished " + filename + " for " + curr_mutation)

fwrite.close()
