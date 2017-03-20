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
"Unknown".

Exported data are non-repetitive in the following sense:
If two entries from the raw data refer to the same nucleotide variation, then they are considered the same mutation."""
# todo add check for same nucleotide variation when removing entries

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
fwrite.write("GeneName\tUniProtID\tMutationID\tMutation\tMutationType\n")

write_data = []
prev_mutation_ids = set()
relevant_mutations = ("Deletion - Frameshift", "Deletion - In frame", "Insertion - Frameshift", "Insertion - In frame",
                      "Substitution - Missense", "Substitution - coding silent")
for line in cosmic_data:
    line = line.rstrip("\n").split("\t")
    mutation_id = line[16]
    mutation_type = line[19]
    if mutation_type not in relevant_mutations or mutation_id in prev_mutation_ids:
        continue

    # invariant: this mutation has not been previously encountered
    prev_mutation_ids.add(mutation_id)
    gene_name = line[0]
    try:
        uniprot_id = cosmic_to_uniprot[gene_name]
    except KeyError:
        uniprot_id = "None"
    mutation = line[18][2:]
    new_data = [gene_name, uniprot_id, mutation_id, mutation, mutation_type]
    for x in range(len(new_data)):
        if not new_data[x]:
            new_data[x] = "None"
    write_data.append(new_data)

write_data.sort(key=lambda x: x[4])  # sort data by mutation type
for line in write_data:
    fwrite.write("\t".join(line) + "\n")
fwrite.close()
