"""Generate a document organizing the UniProt control missense mutations."""

import os


def find_missense_mutations(xml):
    """Given an xml string containing UniProt data, return a list of tuples consisting of a mutation ID and a protein
    missense mutation similar to the form T342Y."""
    seq_variant_start = xml.find('<feature type="sequence variant"')
    seq_variant_end = xml.find("</feature>", seq_variant_start)
    mutations = []
    mutations_set = set()
    while seq_variant_start != -1:
        try:
            # get the mutation id
            id_start = xml.index('id="', seq_variant_start, seq_variant_end) + len('id="')
            id_end = xml.index('"', id_start, seq_variant_end)
            mutation_id = xml[id_start: id_end]

            # get the original amino acid
            original_start = xml.index("<original>", id_end, seq_variant_end) + len("<original>")
            original_end = xml.index("</original>", original_start, seq_variant_end)
            original_amino_acid = xml[original_start: original_end]

            # get the new amino acid
            new_start = xml.index("<variation>", original_end, seq_variant_end) + len("<variation>")
            new_end = xml.index("</variation>", new_start, seq_variant_end)
            new_amino_acid = xml[new_start: new_end]

            if len(original_amino_acid) == len(new_amino_acid) == 1:  # consider only single amino acid changes
                # get the position
                position_start = xml.index('<position position="', new_end, seq_variant_end) + len('<position position="')
                position_end = xml.index('"/>', position_start, seq_variant_end)
                position = xml[position_start: position_end]
                mutation = original_amino_acid + position + new_amino_acid
                if mutation not in mutations_set:
                    mutations.append((mutation_id, mutation))
                    mutations_set.add(mutation)
        except ValueError:  # skip this mutation if data are incomplete
            pass
        seq_variant_start = xml.find('<feature type="sequence variant"', seq_variant_end)
        seq_variant_end = xml.find("</feature>", seq_variant_start)
    return mutations


fwrite = open("MutationControl.txt", "w")
fwrite.write("GeneName\tUniProtID\tMutationID\tProteinMutation\tMutationType\tSource\n")
gene_name = "-"
mutation_type = "Missense"
source = "UniProt"
for file in os.listdir("../Mutation Control/UniProt/"):
    with open("../Mutation Control/UniProt/" + file) as f:
        xml = f.read()
    uniprot = file[:file.index(".")]
    for m in find_missense_mutations(xml):
        mutation_id = m[0]
        protein_mutation = m[1]
        fwrite.write("\t".join((gene_name, uniprot, mutation_id, protein_mutation, mutation_type, source)) + "\n")
fwrite.close()
