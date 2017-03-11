"""Generate refined file of MobiDB entries to query"""

f = open("mobidb_human_entries.txt")
data = f.read()
f.close()

g = open("mobidb_human_ids.txt", "w")
data2 = data.split(":")

for x in range(len(data2)):
    if '"acc"' in data2[x]:
        prot_id = data2[x+1]
        comma = prot_id.index(",")
        prot_id = prot_id[1:comma-1]
        g.write(prot_id + "\n")
g.close()
