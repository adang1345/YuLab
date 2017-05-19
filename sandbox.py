ids = range(46382)
num_threads = 101

threads = []
for a in range(num_threads):
    start = a * (len(ids) // num_threads)
    if a == num_threads-1:
        end = len(ids)
    else:
        end = start + len(ids) // num_threads
    print(start, end)
