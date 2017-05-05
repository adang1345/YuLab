"""It seems that mutation_control1.py downloaded nearly all the data but is missing a few files. This script determines
which files are missing."""

import os

# obtain IDs of proteins
ids_file = open("uniprot_human_ids.txt")
ids = ids_file.read().split()
ids_file.close()

files = set(os.listdir("../Mutation Control/UniProt/"))
for x in ids:
    if (x + ".fasta") not in files:
        print(x)
