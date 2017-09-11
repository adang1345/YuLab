"""Common functions and constants to be reused in multiple scripts."""


amino_acid_abbrev = {'Cys':'C', 'Asp':'D', 'Ser':'S', 'Gln':'Q', 'Lys':'K', 'Ile':'I', 'Pro':'P', 'Thr':'T', 'Phe':'F',
                     'Asn':'N', 'Gly':'G', 'His':'H', 'Leu':'L', 'Arg':'R', 'Trp':'W', 'Ala':'A', 'Val':'V', 'Glu':'E',
                     'Tyr':'Y', 'Met':'M'}


def mutation_abbrev(m):
    """Given mutation m as a string that has 3-letter amino acid abbreviations, return the mutation using 1-letter
    amino acid abbreviations. Example: "Ser180Asn" becomes "S180N"."""
    start_position = -1
    end_position = -1
    for x in range(len(m)):
        if m[x].isdigit() and start_position == -1:
            start_position = x
        elif start_position != -1 and end_position == -1 and not m[x].isdigit():
            end_position = x
            break
    return amino_acid_abbrev[m[:start_position]] + m[start_position: end_position] + amino_acid_abbrev[m[end_position:]]


def uniprot2seq(uniprots):
    """Return a dictionary that maps each UniProtID of the list uniprots to its sequence, prepended with ">".
    This is done by searching the folder ../UniProt Sequences/ for the FASTA sequence. If a certain sequence is not
    found or if the file does not contain a valid FASTA sequence, then ignore it."""
    uniprot_seq = {}
    for uniprot in uniprots:
        try:
            with open("../UniProt Sequences/" + uniprot + ".fasta") as seq_file:
                seq = seq_file.read()
                seq = ">" + seq[seq.index("\n")+1:].replace("\n", "")
                uniprot_seq[uniprot] = seq
        except (FileNotFoundError, ValueError):
            pass
    return uniprot_seq


def head(filepath, n, start=0):
    """Return the first n lines of the filepath file, starting at 0-indexed line start. If file has fewer than n lines,
    then return all the file contents."""
    file = open(filepath)
    file_iter = iter(file)
    lines = []
    for _ in range(n + start):
        try:
            if start <= 0:
                lines.append(next(file_iter))
            else:
                next(file_iter)
                start -= 1
        except StopIteration:
            break
    file.close()
    return "".join(lines)
