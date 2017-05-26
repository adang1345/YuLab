"""Count the number of dbSNP missense mutations that are in disordered interface regions or in
structured interface regions. Also count the number of residues in disordered and structured regions for all
proteins for which there is at least one mutation, disorder data are available, and interface data are available."""


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


# read interface data file and construct mapping from UniProt ID to set containing locations of interface domains
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
# regions. Include experimental and predicted data.
with open("../Disordered Region Data/DisorderData.txt") as ddataf:
    ddata = []
    for x in ddataf.readlines():
        xl = x.split()
        if xl[0] != "Source" and xl[0] == "DisProt":
            ddata.append(xl)
uniprot_disorder = {}
for x in ddata:
    uniprot = x[1]
    dregion = x[5]
    uniprot_disorder[uniprot] = range_to_set(dregion)

# read dbSNP mutation data and consider only missense mutations for proteins that show up in the interface data and
# the disorder data
with open("../Mutation Control/MutationControlNeutral.txt") as mutf:
    mutdata = []
    for x in mutf.readlines():
        xl = x.split()
        if xl[4] == "Missense" and xl[1] in uniprot_interface and xl[1] in uniprot_disorder:
            mutdata.append(xl)

# count number of mutations in each category
in_disordered_interface = 0
in_structured_interface = 0
mut_uniprot = set()
for x in mutdata:
    uniprot = x[1]
    try:
        mut_location = int(x[3][1:-1])
    except ValueError:
        continue
    if mut_location in uniprot_interface[uniprot]:
        if mut_location in uniprot_disorder[uniprot]:
            in_disordered_interface += 1
        else:
            in_structured_interface += 1
    mut_uniprot.add(uniprot)

# count total number of residues in each category
disordered_interface_size = 0
structured_interface_size = 0
total_interface_size = 0
total_protein_count = 0
for x in uniprot_interface:
    if x in uniprot_disorder and x in mut_uniprot:
        disordered_interface_size += len(uniprot_interface[x].intersection(uniprot_disorder[x]))
        structured_interface_size += len(uniprot_interface[x].difference(uniprot_disorder[x]))
        total_interface_size += len(uniprot_interface[x])
        total_protein_count += 1

print("Control mutations in disordered interface\t" + str(in_disordered_interface))
print("Control mutations in structured interface\t" + str(in_structured_interface))
print("Disordered interface total # of residues\t" + str(disordered_interface_size))
print("Structured interface total # of residues\t" + str(structured_interface_size))
print("Total number of proteins considered\t" + str(total_protein_count))
