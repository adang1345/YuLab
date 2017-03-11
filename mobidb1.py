"""Generate raw file containing MobiDB entries to query"""

from urllib.request import urlopen
from urllib.error import HTTPError

# Do one query per protein length. We can obtain at most 100 proteins per query, so our search parameters are
# specific.
f = open("mobidb_human_entries.txt", "w")
for prot_length in range(1, 5460):
    try:
        address = "http://mobidb.bio.unipd.it/ws/search?q=organism:%22Homo%20sapiens%22%20AND%20length:[" + \
              str(prot_length) + "%20TO%20" + str(prot_length+9) + "]"
        json_data = urlopen(address).read().decode("ascii", "ignore")
        f.write(json_data + "\n")
    except HTTPError:
        print(prot_length)

# Do one query per range of 10 protein lengths to save time, as not many proteins are this large.
for prot_length in range(5460, 35000, 10):
    try:
        address = "http://mobidb.bio.unipd.it/ws/search?q=organism:%22Homo%20sapiens%22%20AND%20length:[" + \
              str(prot_length) + "%20TO%20" + str(prot_length+9) + "]"
        json_data = urlopen(address).read().decode("ascii", "ignore")
        f.write(json_data + "\n")
    except HTTPError:
        print(prot_length)

f.close()
