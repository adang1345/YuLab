"""Generate random protein mutations. Create the file ../Mutation Control/Random.txt"""

import common_tools
import random
import os

NUM_MUTATIONS = 20_000
OUTPUT = "../Mutation Control/RandomInterface.txt"


def range_to_set(r):
    """Given a string r similar to the form "a;b-c", return the set {a, b, b+1, ... , c-1, c}"""
    r = [x.split("-") for x in r.split(";")]
    s = set()
    for x in r:
        if len(x) == 1:
            s.add(int(x[0]))
        else:
            s.update(list(range(int(x[0]), int(x[1])+1)))
    return s


with open("../Interface Data/hSIN_organized.txt") as idataf:
    idata = [x.split() for x in idataf.readlines() if x[:8] != "UniProtA"]
uniprot_interface = {}
for x in idata:
    uniprotid_a = x[0]
    uniprotid_b = x[1]
    iregion_a = x[3]
    iregion_b = x[5]
    if uniprotid_a in uniprot_interface:
        uniprot_interface[uniprotid_a].update(range_to_set(iregion_a))
    else:
        uniprot_interface[uniprotid_a] = range_to_set(iregion_a)
    if uniprotid_b in uniprot_interface:
        uniprot_interface[uniprotid_b].update(range_to_set(iregion_b))
    else:
        uniprot_interface[uniprotid_b] = range_to_set(iregion_b)

disorder_type = "DisProt"

with open("../Disordered Region Data/DisorderData.txt") as ddataf:
    ddata = []
    for x in ddataf.readlines():
        xl = x.rstrip().split("\t")
        if disorder_type == "DisProt":
            if xl[0] != "Source" and xl[0] == "DisProt":
                ddata.append(xl)
        elif disorder_type == "Exp":
            if xl[0] != "Source" and xl[7] == "Exp":
                ddata.append(xl)
        elif disorder_type == "All":
            if xl[0] != "Source":
                ddata.append(xl)
        else:
            assert False, "Incorrect disorder type indicated"
uniprot_disorder = {}
for x in ddata:
    uniprot = x[1]
    dregion = x[5]
    uniprot_disorder[uniprot] = range_to_set(dregion)

uniprot_all = list(x[:x.index(".")] for x in os.listdir("../UniProt Sequences"))  # list of all UniProt IDs
amino_acids = list(common_tools.amino_acid_abbrev.values())  # list of all 1-letter amino acid abbreviations
uniprot2seq = common_tools.uniprot2seq(uniprot_all)  # map from UniProt ID to protein sequence prepended with >

print("begin generating random mutations")

rand_mutations = set()  # set of (UniProt ID, mutation) tuples
while len(rand_mutations) < NUM_MUTATIONS:
    uniprotid = random.choice(uniprot_all)
    if uniprotid not in uniprot_interface or uniprotid not in uniprot_disorder:
        continue
    # generate mutation from the UniProt ID
    try:
        seq = uniprot2seq[uniprotid]
    except KeyError:
        continue
    mut_position = random.choice(list(uniprot_interface[uniprotid]))
    mut_from = seq[mut_position]
    mut_to = mut_from
    while mut_to == mut_from:
        mut_to = random.choice(amino_acids)
    mutation = mut_from + str(mut_position) + mut_to

    rand_mutations.add((uniprotid, mutation))

# write to file
output_file = open(OUTPUT, "w")
output_file.write("GeneName\tUniProtID\tMutationID\tProteinMutation\n")
i = 1
for rand_mutation in rand_mutations:
    output_file.write(f"-\t{rand_mutation[0]}\tRandom{i}\t{rand_mutation[1]}\n")
    i += 1
output_file.close()
