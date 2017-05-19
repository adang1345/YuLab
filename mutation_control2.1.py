"""Create a file containing RefSeq IDs for all human missense SNP mutations that have been mapped to chromosomes or
mitochondria. I previously downloaded raw files containing these IDs separated by chromosome, and this script extracts
and synthesizes the RefSeq IDs."""

import os

rsid_list = open("../Mutation Control/dbSNP/refseq_list.txt", "w")
rsids = set()
for x in os.listdir("../Mutation Control/dbSNP/RefSeqID Lists"):
    if x == "downloader.sh":
        continue
    with open("../Mutation Control/dbSNP/RefSeqID Lists/" + x) as fread:
        for line in fread:
            rsid = line.split()[0]
            if rsid[0] == "B":
                rsid = rsid[rsid.index("rs"):]
            if rsid[:2] != "rs":
                print(line)
                input()
            rsids.add(rsid)
rsids = list(rsids)
rsids.sort(key=lambda x: int(x[2:]))
for x in rsids:
    rsid_list.write(x + "\n")
rsid_list.close()
