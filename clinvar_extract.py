"""Extract the missense mutations of interest from ClinVar (date file at
ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz).

Bisque is used to map the RefSeq IDs to UniProt. In the case that a RefSeq ID can be mapped to multiple UniProt IDs, an
arbitrary one is chosen.

Argument to this script must be a comma-separated list of Clinical significance values to extract, surrounded by
quotes."""

import sys
import common_tools

clinvar_filepath = "../Mutation Control/ClinVar/variant_summary.txt"
output_filepath = "../Mutation Control/ClinVar.txt"

assert len(sys.argv) > 1, "Need Clinical Significance argument"
clinsig_filter = set(sys.argv[1].split(","))  # set of clinical significance to accept
print("Keeping only mutations {}".format(clinsig_filter))

fwrite = open(output_filepath, "w")
fwrite.write("GeneName\tUniProtID\tMutationID\tProteinMutation\n")
clinvar_file = open(clinvar_filepath)
next(clinvar_file)  # skip first line
genename_mutationid_proteinmutation = []  # list of tuple (GeneName, UniProtID, MutationID, ProteinMutation)
refseq_ids = set()  # list of RefSeq IDs that need to be mapped to UniProt

# extract relevant data from ClinVar file
for line in clinvar_file:
    line_l = line.rstrip().split("\t")
    type = line_l[1]
    name = line_l[2]
    gene_symbol = line_l[4]
    clinical_significance = line_l[6]
    assembly = line_l[16]
    if (type != "single nucleotide variant" or "(p." not in name or "=)" in name or "Ter)" in name or "p.Ter" in name or
                clinical_significance not in clinsig_filter or assembly != "GRCh38"):
        continue  # filter out mutations we don't want
    rs_num = line_l[9]
    protein_mutation = name[name.index("(p.")+3: name.rindex(")")]
    protein_mutation = common_tools.mutation_abbrev(protein_mutation)  # convert to 1-letter abbreviations
    refseq_id = name[:name.index("(")]
    # temporarily hold RefSeq ID instead of UniProt ID
    genename_mutationid_proteinmutation.append([gene_symbol, refseq_id, rs_num, protein_mutation])
    refseq_ids.add(refseq_id)


# get help from user to do mapping with Bisque
with open("to_map.txt", "w") as to_map:
    to_map.write("\n".join(refseq_ids))
    to_map.write("\n")
input("Use Bisque and store result in text format in from_map.txt. Make sure there are no headers.\n"
      "Press Enter when you are done.\n")

while True:
    try:
        from_map = open("from_map.txt")
        break
    except FileNotFoundError:
        input("from_map.txt not found. Press Enter to try again.\n")

refseq2uniprot = {}
for line in from_map:
    line_l = line.rstrip().split(",")
    refseq = line_l[0]
    uniprot = line_l[3]
    refseq2uniprot[refseq] = uniprot

# filter out unmappable ones, mismatches, and repeats
uniprot2seq = common_tools.uniprot2seq(refseq2uniprot.values())
filtered_mutations = []
map_count = 0  # number of successful maps
mismatch_count = 0  # number of mismatches
already_seen = set()  # set of (UniProtID, ProteinMutation)
for entry in genename_mutationid_proteinmutation:
    mutation_position = int(entry[3][1:-1])
    original = entry[3][0]
    try:
        entry[1] = refseq2uniprot[entry[1]]
    except KeyError:
        continue
    if entry[1] in uniprot2seq:
        seq = uniprot2seq[entry[1]]
        if (entry[1], entry[3]) not in already_seen:
            already_seen.add((entry[1], entry[3]))
            if mutation_position < len(seq) and seq[mutation_position] == original:
                filtered_mutations.append(entry)
                map_count += 1
            else:
                mismatch_count += 1

# write mutations to file
for entry in filtered_mutations:
    fwrite.write("\t".join(entry) + "\n")

print("Successfully mapped: {}".format(map_count))
print("Mismatches: {}".format(mismatch_count))

clinvar_file.close()
