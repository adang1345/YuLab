"""ExAC mutation data were obtained from ftp://ftp.broadinstitute.org/pub/ExAC_release/current/ExAC.r0.3.1.sites.vep.vcf.gz
and unzipped. Extract the missense mutation data from this file into the standard format."""

import common_tools


def hgvsp2standard(s):
    """Given a string s representing a mutation of the form 'ENSP00000334393.3:p.Gly9Asp', return the mutation using
    1-letter amino acid abbreviations.
    Example: 'ENSP00000334393.3:p.Gly9Asp' becomes 'G9D'."""
    i = s.rindex(".")
    return common_tools.mutation_abbrev(s[i+1:])


mutation_file = open("../Mutation Control/ExAC/ExAC.r0.3.1.sites.vep.vcf")
fwrite = open("../Mutation Control/ExAC.txt", "w")

header = ("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|"
         "CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|"
         "MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|DOMAINS|"
         "HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|"
         "MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF|context|ancestral")\
            .split("|")
num_columns = len(header)
consequence_index = header.index("Consequence")
genename_index = header.index("Gene")
uniprot_index = header.index("SWISSPROT")
protmutation_index = header.index("HGVSp")

fwrite.write("\t".join(header))

# uniprot2seq = common_tools.uniprot2seq_all()

# skip headers
for line in mutation_file:
    if line == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n":
        break

fwrite.write("GeneName\tUniProtID\tMutationID\tProteinMutation\n")
# i = 0
consequences = set()
for line in mutation_file:
    start = line.index(";CSQ=") + len(";CSQ=")
    end = line.index(";AC_POPMAX")

    # List of data for relevant mutations in this line. Each element of mutations is a list of data bits corresponding
    # to the header.
    mutations = [m.split("|") for m in line[start: end].split(",")]

    for m in mutations:
        assert len(m) == num_columns, m
        consequence = m[consequence_index]
        genename = m[genename_index]
        uniprot = m[uniprot_index]
        protmutation = m[protmutation_index]

        # skip mutations that are not missense or contain missing information
        if "missense_variant" not in consequence or uniprot == "" or protmutation == "":
            continue

        if genename == "":
            genename = "-"
        fwrite.write("\t".join((genename, uniprot, "-", hgvsp2standard(protmutation))))
        fwrite.write("\n")


    # print(mutations)
    # if i == 100:
    #     break
    # i += 1


mutation_file.close()
fwrite.close()
print(consequences)

# todo test what consequences we have in whole file
