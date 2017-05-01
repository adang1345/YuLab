import io

with io.open("../Mutation Data/HGMD/HGMD_substitutions_2015.3.tsv",encoding="utf8") as mut_file:
    mut_data = set()
    for x in mut_file.readlines():
        xs = x.split("\t")
        if xs[1] == "DM" and xs[7] != "\\N":
            if xs[7] in mut_data:
                print(xs[7])
            mut_data.add(xs[7])
