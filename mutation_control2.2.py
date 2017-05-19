"""Download XML data from dbSNP. Remember to specify at the top of the file the number of threads and the size of each
query. Each downloaded file is named with the refSNP of the first mutation of the file.

Since I have had """

from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from http.client import IncompleteRead

import threading
import datetime
import io
import sys
import os

print(datetime.datetime.now())

num_threads = 20  # number of threads that perform web queries simultaneously
mutations_per_query = 100  # number of mutations to download per query

# obtain IDs of Reference SNP cluster reports already obtained
done_ids = set()
for f in os.listdir("../Mutation Control/dbSNP/XML"):
    done_ids.add(f.split(".")[0])
    with open("../Mutation Control/dbSNP/XML/" + f) as f2:
        for line in f2:
            try:
                rsid_start = line.index('<Rs rsId="') + len('<Rs rsId="')
                rsid_end = line.index('"', rsid_start)
                done_ids.add("rs" + line[rsid_start: rsid_end])
            except ValueError:
                if line != "\n":
                    print(f)
                    print(line)

# obtain IDs of Reference SNP cluster reports to query
ids_file = open("../Mutation Control/dbSNP/refseq_list.txt")
ids = [x for x in ids_file.read().split() if x not in done_ids]
ids_file.close()


def query(a, b):
    """Perform queries and write files."""
    i = a
    while i < b:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id="
        for x in range(i, min(i+mutations_per_query, b)):
            url += ids[x][2:] + ","
        url = url[:-1] + "&report=XML"  # remove final comma
        try:
            dbsnp_data = urlopen(url).read().decode("utf-8")
            data_start = dbsnp_data.index('<Rs rsId="')
            data_end = dbsnp_data.index("</ExchangeSet>")
            with io.open("../Mutation Control/dbSNP/XML/" + ids[i] + ".xml", "w", encoding="utf8") as dbsnp_file:
                dbsnp_file.write(dbsnp_data[data_start: data_end])
            i = i + mutations_per_query
        except (HTTPError, URLError, ConnectionResetError, IncompleteRead, TimeoutError):
            print(ids[i] + " failed, retrying")
            sys.stdout.flush()
        except ValueError:
            print("Skipping: " + url)
            sys.stdout.flush()
            i = i + mutations_per_query


# Divide up task between threads evenly based on number of threads specified at top of file.
threads = []
for a in range(num_threads):
    start = a * (len(ids) // num_threads)
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
