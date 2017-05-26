"""Organize the hSIN data such that all interface domains for a binary interaction between two proteins are in one line.
Ignore titin. Data are organized in the columns UniProtA, UniProtB, PfamA, ResiduesA, PfamB, ResiduesB. A single binary
interaction that contains more than 1 pair of interacting domains will have PfamA, ResiduesA, PfamB, and ResiduesB data
separated by semicolon for each domain and in the correct corresponding order. Overlapping residues are not
consolidated."""

data = {}

f = open("../Interface Data/hSIN_20170524.txt")
next(f)
for line in f:
    line_split = line.split()
    if "Q8WZ42" in line_split:  # ignore titin
        continue
    uniprot_a = line_split[2]
    uniprot_b = line_split[3]
    uab = (uniprot_a, uniprot_b)
    if uab in data:
        data[uab]["pfamA"] += ";" + line_split[4]
        data[uab]["residuesA"] += ";" + line_split[5] + "-" + line_split[6]
        data[uab]["pfamB"] += ";" + line_split[7]
        data[uab]["residuesB"] += ";" + line_split[8] + "-" + line_split[9]
    else:
        data[uab] = {}
        data[uab]["pfamA"] = line_split[4]
        data[uab]["residuesA"] = line_split[5] + "-" + line_split[6]
        data[uab]["pfamB"] = line_split[7]
        data[uab]["residuesB"] = line_split[8] + "-" + line_split[9]
f.close()

fwrite = open("../Interface Data/hSIN_organized.txt", "w")
fwrite.write("UniProtA\tUniProtB\tPfamA\tResiduesA\tPfamB\tResiduesB\n")
for uab in data:
    fwrite.write("\t".join((uab[0], uab[1], data[uab]["pfamA"], data[uab]["residuesA"], data[uab]["pfamB"],
                            data[uab]["residuesB"])) + "\n")
fwrite.close()
