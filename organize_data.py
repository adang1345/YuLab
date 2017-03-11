"""Organize the data from DisProt and MobiDB and put into a standardized table. Columns are Source, UniProtID, LocalID,
PfamID, Seq, DStart, DEnd, DSeq, DMethodType, DMethod, DNotes. Each line represents a single disordered region. There
are potential data overlaps between two experimental regions if their evidence come from different papers.

Source - database of origin"""

import json
import io
import os


def disprot_pfamid(j):
    """Given DisProt json data j, return a comma-separated string of the Pfam accession ID(s) for the protein. If the
    protein has multiple domains, where each domain has its own Pfam ID, then the string will contain multiple IDs.
    If a protein has more than one copy of a domain, then the string may contain the same Pfam ID more than once. If
    the protein has no Pfam IDs, then "None" is returned."""
    pfam_info = j["protein"]["pfam"]
    if not pfam_info:
        return "None"
    return ",".join(pfam_info[x]["accession"] for x in range(len(pfam_info)))


def disprot_notes(r):
    """Given DisProt json region data r, return a comma-separated string of the supplemental notes regarding this
    region. The notes indicate whether experimental data for r might be dubious. If there are no notes, then "None" is
    returned."""
    if r["tags"][0] is None:
        return "None"
    return ",".join(r["tags"])


def mobidb_pfamid(j):
    """Given MobiDB Pfam json data j, return a comma-separated string of the Pfam accession ID(s) for the protein. If
    the protein has multiple domains, where each domain has its own Pfam ID, then the string will contain multiple IDs.
    If a protein has more than one copy of a domain, then the string may contain the same Pfam ID more than once. If
    the protein has no Pfam IDs, then "None" is returned."""
    if not j["graphic"] or not j["graphic"]["regions"]:
        return "None"
    pfam_info = j["graphic"]["regions"]
    return ",".join(pfam_info[x]["metadata"]["accession"] for x in range(len(pfam_info)))


def mobidb_experimental_regions(j):
    """Given MobiDB disorder json data j, return a list of the experimental disordered regions determined by PDB-NMR
    and PDB-XRay consensus data. Although MobiDB offers data obtained directly from DisProt, I will ignore these regions
    because I already have information from DisProt. Each element of the list returned is a dictionary with the
    following keys: 'start' (the start position), 'end' (the end position), and 'type' (the experimental method)."""
    # pdb_nmr and pdb_xray are lists, where each element is a dictionary that contains information about one consensus
    # disorder sequence. Keys to this dictionary are 'start' (the start position), 'end' (the end position), and 'ann'
    # (the type of region, which always 'D', meaning disorder)
    pdb_nmr = [x for x in j["consensus"]["pdb_nmr"] if x["ann"]=="D"]
    pdb_xray = [x for x in j["consensus"]["pdb_xray"] if x["ann"]=="D"]
    exp_regions = []
    for x in pdb_nmr:
        x["type"] = "PDB-NMR"
        del x["ann"]
        exp_regions.append(x)
    for x in pdb_xray:
        x["type"] = "PDB-XRay"
        del x["ann"]
        exp_regions.append(x)
    return exp_regions


def mobidb_predicted_regions(j):
    """Given MobiDB mobidb json data j, return a list of the long disorder regions. The MobiDB database creators used
    an amalgam of experimental and predicted data to generate the long disorder regions, as described in their paper.
    The database and I both consider these predicted regions, not experimental regions. Each element of the list
    returned is a dictionary with the following keys: 'start' (the start position) and 'end' (the end position)."""
    long = [x for x in j["consensus"]["long"] if x["ann"] == "d"]
    return long


# create file for writing organized data
fwrite = open("DisorderData.txt", "w")
fwrite.write("Source\tUniProtID\tLocalID\tPfamID\tSeq\tDStart\tDEnd\tDSeq\tDMethodType\tDMethod\tDNotes\n")

# organize DisProt data
for prot in os.listdir("../Disordered Region Data/DisProt"):
    with io.open("../Disordered Region Data/DisProt/"+prot, "r", encoding='utf8') as f:
        text = f.read()
    j = json.loads(text)
    source = "DisProt"
    uniprotid = j["protein"]["uniprot_accession"]
    localid = j["protein"]["disprot_id"]
    pfamid = disprot_pfamid(j)
    seq = j["protein"]["sequence"]
    regions = j["regions"]  # extract data about disordered regions
    for r in regions:
        dstart = str(r["start"])
        dend = str(r["end"])
        dseq = r["sequence"]
        dmethodtype = "Exp"
        dmethod = r["method"]["id"]
        dnotes = disprot_notes(r)
        fwrite.write(source + "\t" + uniprotid + "\t" + localid + "\t" + pfamid + "\t" + seq + "\t" + dstart + "\t" +
                     dend + "\t" + dseq + "\t" + dmethodtype + "\t" + dmethod + "\t" + dnotes + "\n")


# Organize MobiDB data. Not all predicted regions are used for each predictor; only the long disorder predicted
# consensus and experimental consensus data are considered. Additionally, some of the data I downloaded include entries
# that are not yet available in the current version (2.0) of MobiDB. These entries are ignored.
with open("mobidb_human_ids.txt") as ids_file:
    ids = ids_file.read().split()

for prot in ids:
    with io.open("../Disordered Region Data/MobiDB/" + prot + "_disorder.json", "r", encoding='utf8') as f:
        text = f.read()
    j_disorder = json.loads(text)
    if "error" in j_disorder:
        # If file containing disorder information has an error, then go onto next protein. This means that the
        # database has not been updated to contain detailed information about this protein.
        continue
    with io.open("../Disordered Region Data/MobiDB/"+prot+"_mobidb.json", "r", encoding='utf8') as f:
        text = f.read()
    j_mobidb = json.loads(text)
    with io.open("../Disordered Region Data/MobiDB/"+prot+"_pfam.json", "r", encoding='utf8') as f:
        text = f.read()
    j_pfam = json.loads(text)
    source = "MobiDB"
    uniprotid = j_mobidb["acc"]
    localid = j_mobidb["ename"]
    pfamid = mobidb_pfamid(j_pfam)
    seq = j_mobidb["seq"]
    for r in mobidb_experimental_regions(j_disorder):
        dstart = r["start"]
        dend = r["end"]
        dseq = seq[dstart-1: dend]
        dmethodtype = "Exp"
        dmethod = r["type"]
        dnotes = "None"
        fwrite.write(source + "\t" + uniprotid + "\t" + localid + "\t" + pfamid + "\t" + seq + "\t" + str(dstart) +
                       "\t" + str(dend) + "\t" + dseq + "\t" + dmethodtype + "\t" + dmethod + "\t" + dnotes + "\n")
    for r in mobidb_predicted_regions(j_disorder):
        dstart = r["start"]
        dend = r["end"]
        dseq = seq[dstart-1: dend]
        dmethodtype = "Pred"
        dmethod = "LongConsensus"
        dnotes = "None"
        fwrite.write(source + "\t" + uniprotid + "\t" + localid + "\t" + pfamid + "\t" + seq + "\t" + str(dstart) +
                       "\t" + str(dend) + "\t" + dseq + "\t" + dmethodtype + "\t" + dmethod + "\t" + dnotes + "\n")
