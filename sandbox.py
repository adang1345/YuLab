def range_to_set2(r):
    """Given a string r similar to the form "a;b-c", return the set {a, b, b+1, ... , c-1, c}"""
    r = [x.split("-") for x in r.split(";")]
    s = set()
    for x in r:
        if len(x) == 1:
            s.add(int(x[0]))
        else:
            s.update(list(range(int(x[0]), int(x[1])+1)))
    return s

print(range_to_set2("1-12;19-21;65"))