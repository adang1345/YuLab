"""Organize the data from DisProt and put into a standardized table. Columns are Source, UniProtID, LocalID, PfamID,
Seq, DStart, DEnd, DSeq, DMethodType, DMethod, DNotes. Each line represents a single disordered region. There are
potential data overlaps where two experimental regions

Source - database of origin

"""


import json
import io
import os
import time


def disprot_pfamid(j):
    """Given DisProt json data j, return a semicolon-separated string of the Pfam accession ID(s) for disordered regions
    of the protein. Pfam domains are included if they overlap with at least one amino acid of a disordered region of the
    protein. If a protein has multiple copies of a domain, then only one instance of the Pfam ID will be returned. If
    the protein has no Pfam IDs provided in the database, then "None" is returned."""
    pfam_info = j["protein"]["pfam"]
    if not pfam_info:
        return "None"
    pfam_ids = ""
    disordered_regions = [(j["regions"][r]["start"], j["regions"][r]["end"]) for r in range(len(j["regions"]))]
        # make list of (start,end) tuples for all disordered regions in this protein
    for x in range(len(pfam_info)):
        pfam_candidate = pfam_info[x]
        pfam_candidate_region = (pfam_candidate["start"], pfam_candidate["end"])
        pfam_candidate_id = pfam_candidate["accession"]
        for d in disordered_regions:
            if are_overlapping(pfam_candidate_region, d) and pfam_candidate_id not in pfam_ids:
                pfam_ids += pfam_candidate_id + ";"
    if not pfam_ids:  # if there are no Pfam ids that overlap with any disordered regions
        return "None"
    return pfam_ids[:-1]  # don't want last semicolon


def are_overlapping(a, b):
    """a and b are two tuples or lists of two elements each, where the first element is the start position of some
    sequence, and the second element is the end position of some sequence. The end position must always be greater than
    or equal to the start position. Return whether the two sequences overlap based on these endpoints."""
    assert type(a) == type(b) == tuple or type(a) == type(b) == list
    return a[0] <= b[1] and a[1] >= b[0]


def disprot_regions(j):
    """Given DisProt json data j, extract data about the disordered regions for this protein. Return a tuple (r,m,p).
    r is a string containing the region boundaries in the form "start-end;start-end;start-end".
    m is a semicolon-separated string containing the experimental methods for all regions. No repetition is present.
    p is a semicolon-separated string containing the PubMed IDs for the papers used as evidence for the regions."""
    regions = j["regions"]
    region_boundaries = []
    methods = ""
    pmid = ""
    for r in regions:
        dstart = r["start"]
        dend = r["end"]
        region_boundary = (dstart, dend)
        region_boundaries.append(region_boundary)
        dmethod = r["method"]["name"]
        if dmethod not in methods:
            methods += dmethod + ";"
        pmid += r["pmid"] + ";"
    for _ in range(100):  # 100 iterations should be much more than enough to remove all repeats
        region_boundaries = remove_repeats(region_boundaries)
    str_region_boundaries = [(str(x),str(y)) for x,y in region_boundaries]
    str_region_boundaries = ";".join("-".join(x) for x in str_region_boundaries)
    return str_region_boundaries, methods[:-1], pmid[:-1]


def remove_repeats(l):
    """l is a non-empty list of tuples (a0,a1) where a0,a1 are the endpoints of some subsequence. Return a list where
    most repeats are removed. (a0,a1) is considered a repeat of (b0,b1) if a0 >= b0 and a1 <= b1, i.e., if (a0,a1) is fully
    encompassed by (b0,b1). Regions where one does not fully encompass the other and the regions overlap are considered
    distinct. Note that due to this definition of a repeat, longer disordered regions tend to be kept while shorter ones
    tend to be removed when the data are repetitive.

    Currently, this function fails to handle the case where a sequence c fully encompasses two non-repetitive sequences
    a and b, and c is added to the return list after a and b. To remove all repeats, this function should be iterated n
    times, where n is the maximum number of non-repetitive sequences such that there exists a long sequence that
    encompasses none of these."""
    no_repeats = [l[0]]
    for a in l:
        for b in range(len(no_repeats)):
            if a[0] >= no_repeats[b][0] and a[1] <= no_repeats[b][1]:
                # new element is a repeat, so skip it
                break
            elif no_repeats[b][0] >= a[0] and no_repeats[b][1] <= a[1]:
                # old element is a repeat, so replace with new element and move on
                no_repeats[b] = a
                break
            elif b == len(no_repeats) - 1:
                # new element does not repeat with any other element in the list, so append it and move on
                no_repeats.append(a)
                break
            # new element and old element do not repeat, and we have not yet reached end of no_repeats list, so keep
            # checking next value of no_repeats and comparing it with new value
    return no_repeats


def disprot_notes(r):
    """Given DisProt json region data r, return a semicolon-separated string of the supplemental notes regarding this
    region. The notes indicate whether experimental data for r might be dubious. If there are no notes, then "None" is
    returned."""
    if r["tags"][0] is None:
        return "None"
    return ",".join(r["tags"])


def mobidb_pfamid_experimental(j_pfam, j_disorder):
    """Given MobiDB Pfam json data j, return a semicolon-separated string of the Pfam accession ID(s) for the
    experimental disordered regions of the protein. There are no repeats if the protein contains the same ID more than
    once. If the protein has no Pfam IDs, then "None" is returned."""
    if not j_pfam["graphic"] or not j_pfam["graphic"]["regions"]:
        return "None"
    pfam_info = j_pfam["graphic"]["regions"]
    pfam_ids = ""
    for x in range(len(pfam_info)):
        pfam_candidate = pfam_info[x]["metadata"]
        pfam_candidate_region = (int(pfam_candidate["start"]), int(pfam_candidate["end"]))
        pfam_candidate_id = pfam_candidate["accession"]
        for d in mobidb_experimental_regions(j_disorder)[2]:
            if are_overlapping(pfam_candidate_region, d) and pfam_candidate_id not in pfam_ids:
                pfam_ids += pfam_candidate_id + ";"
    if not pfam_ids:
        return "None"
    return pfam_ids[:-1]


def mobidb_experimental_regions(j):
    """Given MobiDB disorder json data j, return a data tuple for the experimental disordered regions determined by
    PDB-NMR and PDB-XRay consensus data. Although MobiDB offers data obtained directly from DisProt, I will ignore these
    regions because I already have information from DisProt. The first element of the tuple is a string containing all
    the regions, and the second element is a string containing all the experimental methods. The third element is a list
    of region boundaries."""
    pdb_nmr = [(x["start"],x["end"]) for x in j["consensus"]["pdb_nmr"] if x["ann"]=="D"]
    pdb_xray = [(x["start"],x["end"]) for x in j["consensus"]["pdb_xray"] if x["ann"]=="D"]
    region_boundaries = pdb_nmr + pdb_xray
    if not region_boundaries:
        return "", "", []
    for _ in range(100):
        region_boundaries = remove_repeats(region_boundaries)
    str_region_boundaries = [(str(x), str(y)) for x, y in region_boundaries]
    str_region_boundaries = ";".join("-".join(x) for x in str_region_boundaries)
    methods = ""
    if pdb_nmr:
        methods += "PDB-NMR;"
    if pdb_xray:
        methods += "PDB-XRay;"
    return str_region_boundaries, methods[:-1], region_boundaries


def mobidb_predicted_regions(j):
    """Given MobiDB mobidb json data j, return a list of the long disorder regions. The MobiDB database creators used
    an amalgam of experimental and predicted data to generate the long disorder regions, as described in their paper.
    The database and I both consider these predicted regions, not experimental regions. Each element of the list
    returned is a dictionary with the following keys: 'start' (the start position) and 'end' (the end position)."""
    long = [(str(x["start"]),str(x["end"])) for x in j["consensus"]["long"] if x["ann"] == "d"]
    str_long = ";".join("-".join(x) for x in long)
    return str_long


start_time = time.time()

# create file for writing organized data
fwrite = open("DisorderData2.txt", "w")
fwrite.write("Source\tUniProtID\tLocalID\tPfamID\tSeq\tDRegions\tDMethods\tDMethodType\tPMID\n")
used_uniprot_ids = set()  # track UniProt IDs of already-used UniProt IDs so that repeat entries can be ignored

# organize DisProt data
for prot in os.listdir("../Disordered Region Data/DisProt"):
    with io.open("../Disordered Region Data/DisProt/"+prot, "r", encoding='utf8') as f:
        text = f.read()
    j = json.loads(text)
    source = "DisProt"
    uniprotid = j["protein"]["uniprot_accession"]
    used_uniprot_ids.add(uniprotid)
    localid = j["protein"]["disprot_id"]
    pfamid = disprot_pfamid(j)
    seq = j["protein"]["sequence"]
    dregion_info = disprot_regions(j)
    dregions = dregion_info[0]
    dmethods = dregion_info[1]
    dmethodtype = "Exp"
    dnotes = dregion_info[2]
    fwrite.write(source + "\t" + uniprotid + "\t" + localid + "\t" + pfamid + "\t" + seq + "\t" + dregions + "\t" +
                 dmethods + "\t" + dmethodtype + "\t" + dnotes + "\n")


# Organize MobiDB data. Not all predicted regions are used for each predictor; only the long disorder predicted
# consensus and experimental consensus data are considered. Additionally, some of the data I downloaded include entries
# that are not yet available in the current version (2.0) of MobiDB. These entries are ignored.
with open("mobidb_human_ids.txt") as ids_file:
    ids = ids_file.read().split()

for prot in ids:
    # read and parse files
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

    # obtain experimental data
    source = "MobiDB"
    uniprotid = j_mobidb["acc"]
    if uniprotid in used_uniprot_ids:  # ignore this protein if DisProt has already given information
        continue
    localid = j_mobidb["ename"]
    pfamid = mobidb_pfamid_experimental(j_pfam, j_disorder)
    seq = j_mobidb["seq"]
    mobidb_experimental = mobidb_experimental_regions(j_disorder)
    if mobidb_experimental[0]:  # if any experimental regions exist, record their data
        dregions = mobidb_experimental[0]
        dmethods = mobidb_experimental[1]
        dmethodtype = "Exp"
        dnotes = "None"
        fwrite.write(source + "\t" + uniprotid + "\t" + localid + "\t" + pfamid + "\t" + seq + "\t" + dregions + "\t" +
                 dmethods + "\t" + dmethodtype + "\t" + dnotes + "\n")
    else:  # otherwise, search for predicted regions
        dregions = mobidb_predicted_regions(j_disorder)
        if dregions:  # if there are any predicted regions
            dmethods = "LongConsensus"
            dmethodtype = "Pred"
            dnotes = "None"
            fwrite.write(source + "\t" + uniprotid + "\t" + localid + "\t" + pfamid + "\t" + seq + "\t" + dregions +
                         "\t" + dmethods + "\t" + dmethodtype + "\t" + dnotes + "\n")

fwrite.close()

print("Total runtime (s): " + str(time.time()-start_time))
