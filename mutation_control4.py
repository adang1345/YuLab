"""nonpathogenic.xml contains mutations obtained from dbSNP that are labeled as "Benign", "Likely Benign", "Other",
"Uncertain Significance", and "Untested". Extract these
mutations and select only the ones that can be mapped to UniProt protein positions using refseq2uniprot.txt.

Currently, my solution is to take the RefSeq ID of a mutation and map it to a UniProt ID. If a single RefSeq ID maps to
more than one UniProt ID, then I pick a UniProt ID arbitrarily. Due to the complication that the RefSeq ID might pertain
to a protein that is not the canonical isoform in UniProt, I then check whether the original amino acid specified in the
mutation matches that in the canonical sequence from UniProt. If not, I discard this mutation. This filtering method
has an error rate of about 2.5%.
"""

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


# make set of (rsID, RefSeq, protein mutation) tuples extracted from dbSNP
nonpathogenic_mutation_file = open("../Mutation Control/dbSNP/nonpathogenic.xml")
rsid_refseqid_mutation = set()
for line in nonpathogenic_mutation_file:
    rsid_start = line.index('<Rs rsId="') + len('<Rs rsId="')
    rsid_end = line.index('"', rsid_start)
    rsid = "rs" + line[rsid_start: rsid_end]
    mutation_start = line.index(":p.") + len(":p.")
    refseqid_start = mutation_start
    while line[refseqid_start] != ">":  # RefSeq ID starts at the closest > sign left of mutation
        refseqid_start -= 1
        if refseqid_start == 0:
            raise IndexError("Entrez ID not found")
    refseqid = line[refseqid_start + 1: mutation_start - len(":p.")]
    mutation_end = line.index("</hgvs>", mutation_start)
    try:
        rsid_refseqid_mutation.add((rsid, refseqid, mutation_abbrev(line[mutation_start: mutation_end])))
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
            rsid_refseqid_mutation.add((rsid, refseqid, mutation_abbrev(line[mutation_start: mutation_end])))
        except (ValueError, KeyError):
            break
    # break
nonpathogenic_mutation_file.close()


# Construct map from RefSeq ID to UniProt
refseq_uniprot = {}
with open("../Mutation Control/dbSNP/refseq2uniprot.txt") as refseq2uniprot_file:
    for line in refseq2uniprot_file.readlines():
        line_split = line.split()
        refseqids = line_split[0].split(",")
        uniprot = line_split[1]
        for x in refseqids:
            refseq_uniprot[x] = uniprot


# Of the UniProt IDs that exist in the RefSeq-to-UniProt mapping, construct a mapping from UniProt ID to protein
# sequence in FASTA format without the header, preceded by a ">". The reason for having a character inserted in front of
# the protein sequence is so that the sequence can be 1-indexed instead of 0-indexed.
uniprot_seq = {}
for uniprot in refseq_uniprot.values():
    try:
        with open("../UniProt Sequences/" + uniprot + ".fasta") as seq_file:
            seq = seq_file.read()
            seq = ">" + seq[seq.index("\n")+1:].replace("\n", "")
            uniprot_seq[uniprot] = seq
    except FileNotFoundError:
        pass


# Add UniProt information to each mutation and filter out the mutations that do not map to canonical sequences.
uniprot_rsid_refseqid_mutation = set()
uniprot_mutation = set()  # keep set of (UniProt ID, mutation) tuples to protect against repeats
for m in rsid_refseqid_mutation:
    refseqid = m[1]
    mutation = m[2]
    try:  # ignore this mutation if there's no RefSeq-UniProt mapping or if we don't have the UniProt sequence
        uniprot = refseq_uniprot[refseqid]
        seq = uniprot_seq[uniprot]
    except KeyError:
        continue
    mutation_original = mutation[0]
    mutation_position = int(mutation[1:-1])
    if (mutation_position < len(seq) and seq[mutation_position] == mutation_original
        and (uniprot, mutation) not in uniprot_mutation):
        uniprot_rsid_refseqid_mutation.add((uniprot, m[0], refseqid, mutation))
    uniprot_mutation.add((uniprot, mutation))

# convert data set into a list and sort by UniProt ID
uniprot_rsid_refseqid_mutation = list(uniprot_rsid_refseqid_mutation)
uniprot_rsid_refseqid_mutation.sort(key=lambda x: x[0])

# Write data file
gene_name = "-"
mutation_type = "Missense"
source = "dbSNP"
fwrite = open("../Mutation Control/MutationControlNonpathogenic.txt", "w")
fwrite.write("GeneName\tUniProtID\tMutationID\tProteinMutation\tMutationType\tSource\n")
for m in uniprot_rsid_refseqid_mutation:
    uniprot_id = m[0]
    mutation_id = m[1]
    protein_mutation = m[3]
    fwrite.write("\t".join((gene_name, uniprot_id, mutation_id, protein_mutation, mutation_type, source)) + "\n")
fwrite.close()
