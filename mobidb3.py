"""Download data from MobiDB"""

from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from http.client import IncompleteRead

import threading
import datetime
import sys
import io

print(datetime.datetime.now())

num_threads = 100

# obtain IDs of proteins to query
ids_file = open("mobidb_human_ids.txt")
ids = ids_file.read().split()
ids_file.close()


# perform queries and write to files
def query(a, b):
    i = a
    while i < b:
        try:
            # write MobiDB data
            mobidb_data = urlopen("http://mobidb.bio.unipd.it/ws/entries/" + ids[i] + "/mobidb").read().decode("utf-8")
            with io.open("../Disordered Region Data/MobiDB/" + ids[i] + "_mobidb.json", "w", encoding="utf8") as mobidb_file:
                mobidb_file.write(mobidb_data)
            # write Pfam data
            pfam_data = urlopen("http://mobidb.bio.unipd.it/ws/entries/" + ids[i] + "/pfam").read().decode("utf-8")
            with io.open("../Disordered Region Data/MobiDB/" + ids[i] + "_pfam.json", "w", encoding="utf8") as pfam_file:
                pfam_file.write(pfam_data)
            # write Disorder data
            disorder_data = urlopen("http://mobidb.bio.unipd.it/ws/entries/" + ids[i] + "/disorder").read().decode("utf-8")
            with io.open("../Disordered Region Data/MobiDB/" + ids[i] + "_disorder.json", "w", encoding="utf8") as disorder_file:
                disorder_file.write(disorder_data)
            i += 1
        except (HTTPError, URLError, ConnectionResetError, IncompleteRead):
            print(ids[i])
            sys.stdout.flush()


# Divide up task between threads evenly based on number of threads specified at top of file.
threads = []
for a in range(num_threads):
    start = a * len(ids) // num_threads
    if a == num_threads-1:
        end = len(ids)
    else:
        end = start + len(ids) // num_threads
    t = threading.Thread(target=query, args=(start, end))
    threads.append(t)
    t.start()

for t in threads:
    t.join()

print(datetime.datetime.now())
