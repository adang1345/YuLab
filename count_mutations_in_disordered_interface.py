"""Count the number of substitution mutations that are in disordered interface regions or in
structured interface regions. Also count the number of residues in disordered and structured regions for all
proteins for which there is at least one mutation, disorder data are available, and interface data are available.

This script takes the following arguments in the order given:

path the to file containing the mutations
    Must be tab-delimited with columns GeneName, UniProtID, MutationID, ProteinMutation in that order
type of disordered region to consider, must be "All", "Exp", or "DisProt"
"""

# TODO create file containing list of each thing mentioned in results, to debug

import sys
import os


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


# check that arguments are correct
assert len(sys.argv) == 3
assert os.path.isfile(sys.argv[1])
assert sys.argv[2] in ("All", "Exp", "DisProt")

mutation_filepath = sys.argv[1]
disorder_type = sys.argv[2]

# read interface data file and construct mapping from UniProt ID to set containing locations of interface regions
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

# Read disorder region data file and construct mapping from UniProt ID to set containing locations of disordered
# regions. Type of data included depends on user input.
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

# read mutation data and consider only mutations for proteins that show up in the interface data and
# the disorder data
with open(mutation_filepath) as mutf:
    mutdata = []
    for x in mutf.readlines():
        xl = x.split()
        if xl[1] in uniprot_interface and xl[1] in uniprot_disorder:
            mutdata.append(xl)

# count number of mutations in each category
in_disordered_interface = 0
in_structured_interface = 0
mut_uniprot = set()  #
for x in mutdata:
    uniprot = x[1]
    # try:
    mut_location = int(x[3][1:-1])
    # except ValueError:
    #     continue
    if mut_location in uniprot_interface[uniprot]:
        if mut_location in uniprot_disorder[uniprot]:
            print(x, "is in disordered interface")
            in_disordered_interface += 1
        else:
            print(x, "is in structured interface")
            in_structured_interface += 1
        mut_uniprot.add(uniprot)

# count total number of residues in each category
disordered_interface_size = 0
structured_interface_size = 0
total_protein_count = 0
for x in uniprot_interface:
    if x in uniprot_disorder and x in mut_uniprot:
        print(x, "has interface", uniprot_interface[x])
        print(x, "has disorder", uniprot_disorder[x])
        print(x, "has disordered interface", uniprot_interface[x].intersection(uniprot_disorder[x]))
        print(x, "has structured interface", uniprot_interface[x].difference(uniprot_disorder[x]))
        disordered_interface_residues = uniprot_interface[x].intersection(uniprot_disorder[x])
        disordered_interface_size += len(disordered_interface_residues)
        structured_interface_residues = uniprot_interface[x].difference(uniprot_disorder[x])
        structured_interface_size += len(structured_interface_residues)
        assert disordered_interface_residues.intersection(structured_interface_residues) == set()
        total_protein_count += 1

print("Mutations in disordered interface\t" + str(in_disordered_interface))
print("Mutations in structured interface\t" + str(in_structured_interface))
print("Disordered interface total # of residues\t" + str(disordered_interface_size))
print("Structured interface total # of residues\t" + str(structured_interface_size))
print("Total number of proteins considered\t" + str(total_protein_count))
