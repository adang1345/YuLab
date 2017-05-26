"""Extract substitution disease mutations from HGMD. Include only those that do not have "R" in the Variant_class
column. Include only those for which there exists a mapping from HGMD to UniProt ID without any errors or warnings.

Data extracted are non-repetitive at the protein level."""

import io

# construct mapping from HGMD to UniProt
with open("../Mutation Data/HGMD/hgmd2uniprot.txt") as map_f:
    hgmd2uniprot = {}
    for x in map_f.readlines():
        xs = x.strip().split()
        if len(xs) != 1 and xs[1] != "BISQUE" and xs[1] != "WARNING:" and xs[1] != "WARNING":
            hgmd2uniprot[xs[0]] = xs[1]

# read mutations from HGMD and write file
fwrite = open("../Mutation Data/HGMDMutationData.txt", "w")
fwrite.write("GeneName\tUniProtID\tMutationID\tProteinMutation\n")

already_seen = set()
with io.open("../Mutation Data/HGMD/HGMD_substitutions_2015.3.tsv", encoding="utf8") as mut_file:
    for x in mut_file.readlines():
        xs = x.split("\t")
        gene_name = xs[8]
        protein_mutation = xs[7]
        if protein_mutation not in hgmd2uniprot:  # ignore mutations for which no mapping to UniProt can be found
            continue
        uniprotid = hgmd2uniprot[protein_mutation]
        mutationid = xs[2]
        if xs[1] != "R" and protein_mutation != "\\N" and protein_mutation not in already_seen:
            protein_mutation2 = protein_mutation.split(":")[1][2:]
            fwrite.write("\t".join((gene_name, uniprotid, mutationid, protein_mutation2)) + "\n")
        already_seen.add(protein_mutation)

fwrite.close()
