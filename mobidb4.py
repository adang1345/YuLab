"""It seems that mobidb3.py downloaded nearly all the data but is missing a few files. This script determines which
files are missing."""

import os

# obtain IDs of proteins
ids_file = open("mobidb_human_ids.txt")
ids = ids_file.read().split()
ids_file.close()

files = set(os.listdir("../Disordered Region Data/MobiDB"))
for x in ids:
    if (x + "_mobidb.json" not in files) or (x + "_pfam.json" not in files) or (x + "_disorder.json" not in files):
        print(x)
