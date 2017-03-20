"""Given cancer mutation data from COSMIC, create a data file that organizes the data.

COSMIC mutation data were obtained from http://cancer.sanger.ac.uk/cosmic/download and are in
    ../Mutation Data/CosmicMutantExport.tsv
The mapping between gene name and UniProt ID was obtained from http://www.genenames.org/cgi-bin/download and is in
    ../Mutation Data/CosmicUniProt.txt

Data from COSMIC contain many types of mutations. I extract only

"Deletion - Frameshift"
"Deletion - In frame"
"Insertion - Frameshift"
"Insertion - In frame"
"Substitution - Missense"
"Substitution - coding silent"

mutations. I ignore

"Complex"
"Complex - compound substitution"
"Complex - deletion inframe"
"Complex - frameshift"
"Complex - insertion inframe"
"Frameshift"
"Nonstop extension"
"Substitution - Nonsense"
"Whole gene deletion"
"Unknown"."""

# create mapping between gene name and UniProt ID
uniprot_file = open("../Mutation Data/COSMIC/CosmicUniProt.txt")
next(uniprot_file)
cosmic_to_uniprot = {}
for line in uniprot_file:
    line = line.split()
    if len(line) == 2:
        cosmic_to_uniprot[line[0]] = line[1]
    else:
        cosmic_to_uniprot[line[0]] = "None"
uniprot_file.close()

cosmic_data = open("../Mutation Data/COSMIC/CosmicMutantExport.tsv")
next(cosmic_data)
fwrite = open("CosmicMutationData.txt", "w")
fwrite.write("GeneName\tUniProtID\tMutationID\tMutation\tMutationType\tPMID\n")

write_data = []
prev_mutation_id = ""
for line in cosmic_data:
    line = line.rstrip("\n").split("\t")
    mutation_id = line[16]
    pmid = line[30]
    if mutation_id == prev_mutation_id:  # if this line is a repeat, then add information to previous line and continue
        prev_pmids = write_data[-1][-1]
        if pmid and prev_pmids == "None":
            write_data[-1][-1] = pmid
        elif pmid and pmid not in prev_pmids:
            write_data[-1][-1] += ";" + pmid
        continue
    # invariant: this line is not a repeat of the previous
    gene_name = line[0]
    try:
        uniprot_id = cosmic_to_uniprot[gene_name]
    except KeyError:
        uniprot_id = "None"
    mutation = line[18][2:]
    mutation_type = line[19]
    new_data = [gene_name, uniprot_id, mutation_id, mutation, mutation_type, pmid]
    for x in range(len(new_data)):
        if not new_data[x]:
            new_data[x] = "None"
    if new_data[4] in ("Deletion - Frameshift", "Deletion - In frame", "Insertion - Frameshift", "Insertion - In frame",
                       "Substitution - Missense", "Substitution - coding silent"):
        write_data.append(new_data)
    prev_mutation_id = mutation_id

write_data.sort(key=lambda x: x[4])  # sort data by mutation type
for line in write_data:
    fwrite.write("\t".join(line) + "\n")
fwrite.close()
