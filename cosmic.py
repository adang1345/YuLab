"""Given cancer mutation data from COSMIC, create a data file that organizes the data.

COSMIC mutation data were obtained from http://cancer.sanger.ac.uk/cosmic/download and are in
    ../Mutation Data/CosmicMutantExport.tsv
The mapping between gene name and UniProt ID was obtained from http://www.genenames.org/cgi-bin/download and is in
    ../Mutation Data/CosmicUniProt.txt"""

# create mapping between gene name and UniProt ID
uniprot_file = open("../Mutation Data/CosmicUniProt.txt")
next(uniprot_file)
cosmic_to_uniprot = {}
for line in uniprot_file:
    line = line.split()
    if len(line) == 2:
        cosmic_to_uniprot[line[0]] = line[1]
    else:
        cosmic_to_uniprot[line[0]] = "None"
uniprot_file.close()

cosmic_data = open("../Mutation Data/CosmicMutantExport.tsv")
next(cosmic_data)
fwrite = open("CosmicMutationData.txt", "w")
fwrite.write("GeneName\tUniProtID\tMutationID\tMutation\tMutationType\tPMID\n")

for line in cosmic_data:
    line = line.rstrip("\n").split("\t")
    gene_name = line[0]
    try:
        uniprot_id = cosmic_to_uniprot[gene_name]
    except KeyError:
        uniprot_id = "None"
    mutation_id = line[16]
    mutation = line[18][2:]
    mutation_type = line[19]
    pmid = line[30]
    new_data = [gene_name, uniprot_id, mutation_id, mutation, mutation_type, pmid]
    for x in range(len(new_data)):
        if not new_data[x]:
            new_data[x] = "None"
    fwrite.write("\t".join(new_data) + "\n")

fwrite.close()
