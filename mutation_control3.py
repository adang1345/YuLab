"""nondisease.xml contains mutations obtained from dbSNP that are labeled as "Benign" or "Likely Benign". Extract these
mutations and select only the ones that can be mapped to UniProt protein positions using refseq2uniprot.txt.

Currently, my solution is to take the RefSeq ID of a mutation and map it to a UniProt ID. Due to the complication that
the RefSeq ID might pertain to a protein that is not the canonical isoform in UniProt, I check whether the 
"""

#  TODO finish control mutation extraction and comments

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


# make set of (RefSeq, protein mutation) tuples extracted from dbSNP
nondisease_mutation_file = open("../Mutation Control/dbSNP/nondisease.xml")
refseqid_mutation = set()
for line in nondisease_mutation_file:
    mutation_start = line.index(":p.") + len(":p.")
    refseqid_start = mutation_start
    while line[refseqid_start] != ">":  # RefSeq ID starts at the closest > sign left of mutation
        refseqid_start -= 1
        if refseqid_start == 0:
            raise IndexError("Entrez ID not found")
    refseqid = line[refseqid_start + 1: mutation_start - len(":p.")]
    mutation_end = line.index("</hgvs>", mutation_start)
    try:
        refseqid_mutation.add((refseqid, mutation_abbrev(line[mutation_start: mutation_end])))
    except KeyError:  # skip this mutation if there's unknown information
        pass
    while True:
        try:
            mutation_start = line.index(":p.", mutation_end) + len(":p.")
            refseqid_start = mutation_start
            while line[refseqid_start] != ">":
                refseqid_start -= 1
                if refseqid_start == 0:
                    raise IndexError("Entrez ID not found")
            refseqid = line[refseqid_start + 1: mutation_start - len(":p.")]
            mutation_end = line.index("</hgvs>", mutation_start)
            refseqid_mutation.add((refseqid, mutation_abbrev(line[mutation_start: mutation_end])))
        except (ValueError, KeyError):
            break
nondisease_mutation_file.close()


print("Total number of mutations: " + str(len(refseqid_mutation)))


mapping_file = open("../Mutation Control/dbSNP/uniprot_human2refP1may17.txt")
entrezids = set()
for line in mapping_file:
    entrezids.update(line.split()[1].split(";"))

c = 0
for x in refseqid_mutation:
    if x[0] in entrezids:
        c += 1
print(c)
