"""Read nondisease mutations from dbSNP and construct a file containing a list of RefSeq protein IDs."""

# map from 3-letter abbreviation to 1-letter abbreviation
amino_acid_abbrev = {'Cys':'C', 'Asp':'D', 'Ser':'S', 'Gln':'Q', 'Lys':'K', 'Ile':'I', 'Pro':'P', 'Thr':'T', 'Phe':'F',
                     'Asn':'N', 'Gly':'G', 'His':'H', 'Leu':'L', 'Arg':'R', 'Trp':'W', 'Ala':'A', 'Val':'V', 'Glu':'E',
                     'Tyr':'Y', 'Met':'M'}


def mutation_abbrev(m):
    """Given mutation m as a string that has 3-letter amino acid abbreviations, return the mutation using 1-letter
    amino acid abbreviations. Example: "Ser180Asn" becomes "S180N"."""
    start_position = -1
    end_position = -1
    for x in range(len(m)):
        if m[x].isdigit() and start_position == -1:
            start_position = x
        elif start_position != -1 and end_position == -1 and not m[x].isdigit():
            end_position = x
            break
    return amino_acid_abbrev[m[:start_position]] + m[start_position: end_position] + amino_acid_abbrev[m[end_position:]]


# make set of RefSeqIDs extracted from dbSNP
nondisease_mutation_file = open("../Mutation Control/dbSNP/nondisease.xml")
refseqids = set()
for line in nondisease_mutation_file:
    mutation_start = line.index(":p.")
    refseqid_start = mutation_start
    while line[refseqid_start] != ">":  # RefSeq ID starts at the closest > sign left of mutation
        refseqid_start -= 1
        if refseqid_start == 0:
            raise IndexError("RefSeqID not found")
    refseqid = line[refseqid_start + 1: mutation_start]
    if refseqid[2] == "_":  # check that RefSeq ID is valid
        refseqids.add(refseqid)
    while True:
        try:
            mutation_start = line.index(":p.", mutation_start+4)
            refseqid_start = mutation_start
            while line[refseqid_start] != ">":
                refseqid_start -= 1
                if refseqid_start == 0:
                    raise IndexError("RefSeqID not found")
            refseqid = line[refseqid_start + 1: mutation_start]
            if refseqid[2] == "_":
                refseqids.add(refseqid)
        except (ValueError, KeyError):
            break
nondisease_mutation_file.close()

print("Number of RefSeq IDs: " + str(len(refseqids)))

f = open("../Mutation Control/dbSNP/refseq_list.txt", "w")
for x in refseqids:
    f.write(x + "\n")
f.close()
