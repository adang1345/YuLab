# with open("TCGAMutationData.txt") as tcga_file:
#     next(tcga_file)
#     tcga_data = set()
#     for line in tcga_file:
#         line_data = line.split("\t")
#         tcga_data.add((line_data[0], line_data[2]))
#     print("TCGA Unique entries: " + str(len(tcga_data)))

with open("CosmicMutationData.txt") as cosmic_file:
    next(cosmic_file)
    cosmic_data = set()
    for line in cosmic_file:
        line_data = line.split("\t")
        line_data = (line_data[0], line_data[3])
        if line_data in cosmic_data:
            print(line, end="")
        cosmic_data.add(line_data)
    print("COSMIC Unique entries: " + str(len(cosmic_data)))
