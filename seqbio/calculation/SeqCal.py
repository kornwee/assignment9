def gcContent(seq):
    # G+C/(A+T+G+C)
    sum = countBase(seq, 'G') + countBase(seq, 'C')
    n = len(seq)
    n = float(n)

    return sum/n

def atContent(seq):
    # A+T/(A+T+G+C)
    sum = countBase(seq, 'A') + countBase(seq, 'T')
    n = len(seq)
    n = float(n)

    return sum/n

def countBase(seq, base):
    return seq.count(base.upper())

def countBasesDict(seq):
    basesM = {}
    for base in seq:
        basesM[base] = basesM.get(base, 0)+1
    return basesM