"""Disprot is a database containing experimentally-determined disordered regions in proteins. This script submits
queries to obtain data for all human proteins from the website. disprot_human_ids.txt contains the IDs of all
human protein entries. Data are saved in ../Disordered Region Data/DisProt."""

from urllib.request import urlopen

# obtain IDs of proteins to query
ids_file = open("disprot_human_ids.txt")
ids = ids_file.read().split()
ids_file.close()

# perform queries and write to files
for x in ids:
    json_data = urlopen("http://www.disprot.org/ws/get/" + x).read().decode("ascii", "ignore")
    print(json_data)
    data_file = open("../Disordered Region Data/DisProt/" + x + ".json", "w")
    data_file.write(json_data)
    data_file.close()
