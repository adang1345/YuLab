import io
import json

def mobidb_experimental_regions(j):
    """Given MobiDB disorder json data j, return a list of the experimental disordered regions determined by PDB-NMR
    and PDB-XRay consensus data. Although MobiDB offers data obtained directly from DisProt, I will ignore these regions
    because I already have information from DisProt. Each element of the list returned is a tuple (a,b), where a is the
    start position of a region and b is the end position of a region."""
    pdb_nmr = [(x["start"],x["end"]) for x in j["consensus"]["pdb_nmr"] if x["ann"]=="D"]
    pdb_xray = [(x["start"],x["end"]) for x in j["consensus"]["pdb_xray"] if x["ann"]=="D"]
    experimental_regions = pdb_nmr + pdb_xray
    if not experimental_regions:
        return []
    for _ in range(100):
        experimental_regions = remove_repeats(experimental_regions)
    return experimental_regions

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

# Organize MobiDB data. Not all predicted regions are used for each predictor; only the long disorder predicted
# consensus and experimental consensus data are considered. Additionally, some of the data I downloaded include entries
# that are not yet available in the current version (2.0) of MobiDB. These entries are ignored.
with io.open("../Disordered Region Data/MobiDB/P21333_disorder.json", "r", encoding='utf8') as f:
    text = f.read()
j_disorder = json.loads(text)
with io.open("../Disordered Region Data/MobiDB/P21333_mobidb.json", "r", encoding='utf8') as f:
    text = f.read()
j_mobidb = json.loads(text)
with io.open("../Disordered Region Data/MobiDB/P21333_pfam.json", "r", encoding='utf8') as f:
    text = f.read()
j_pfam = json.loads(text)
uniprotid = j_mobidb["acc"]
print(uniprotid)
print(mobidb_experimental_regions(j_disorder))