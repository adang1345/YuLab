"""Download FASTA sequences from UniProt."""

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
ids_file = open("uniprot_human_ids.txt")
ids = ids_file.read().split()
ids_file.close()


def query(a, b):
    """Perform queries and write files."""
    i = a
    while i < b:
        try:
            uniprot_data = urlopen("http://www.uniprot.org/uniprot/" + ids[i] + ".fasta").read().decode("utf-8")
            with io.open("../Mutation Control/UniProt/" + ids[i] + ".fasta", "w", encoding="utf8") as uniprot_file:
                uniprot_file.write(uniprot_data)
            i += 1
        except (HTTPError, URLError, ConnectionResetError, IncompleteRead, TimeoutError):
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
