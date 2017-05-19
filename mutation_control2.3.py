"""mutation_control2.2.py downloaded all human missense SNPs from dbSNP, split into many files. Combine all the data
into 1 file and remove repeat entries. Repeat entries can arise when two refSNP IDs have been merged but were downloaded
in separate queries. Save new large file in neutral.xml."""

import os

fwrite = open("../Mutation Control/dbSNP/neutral.xml", "w")
done_ids = set()
for f in os.listdir("../Mutation Control/dbSNP/XML"):
    with open("../Mutation Control/dbSNP/XML/" + f) as f2:
        for line in f2:
            if line == "\n":  # move to next line if this one is blank
                continue
            try:
                rsid_start = line.index('<Rs rsId="') + len('<Rs rsId="')
                rsid_end = line.index('"', rsid_start)
                rsid = "rs" + line[rsid_start: rsid_end]
                if rsid not in done_ids:
                    fwrite.write(line)
                done_ids.add(rsid)
            except ValueError:
                print("Error occurred at file " + f + " at line " + line)
fwrite.close()
