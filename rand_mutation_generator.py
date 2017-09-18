"""Generate random protein mutations. Create the file ../Mutation Control/Random.txt"""

import common_tools
import random
import os

NUM_MUTATIONS = 10_000_000
OUTPUT = "../Mutation Control/Random.txt"

uniprot_all = list(x[:x.index(".")] for x in os.listdir("../UniProt Sequences"))  # list of all UniProt IDs
amino_acids = list(common_tools.amino_acid_abbrev.values())  # list of all 1-letter amino acid abbreviations
uniprot2seq = common_tools.uniprot2seq(uniprot_all)  # map from UniProt ID to protein sequence prepended with >

print("begin generating random mutations")

# generate sample set of all (uniprot id, mutation position) pairs
uniprot_position_sampleset = []
for uniprotid in uniprot2seq:
    positions = range(1, len(uniprot2seq[uniprotid]))
    uniprot_position_sampleset.extend((uniprotid, position) for position in positions)

rand_mutations = set()  # set of (UniProt ID, mutation) tuples
while len(rand_mutations) < NUM_MUTATIONS:
    uniprotid, mut_position = random.choice(uniprot_position_sampleset)

    # generate mutation from the UniProt ID
    try:
        seq = uniprot2seq[uniprotid]
    except KeyError:
        continue
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
