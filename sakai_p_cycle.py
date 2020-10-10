import itertools, random

def create_from_dcycle(d_cycles):
    tbl = {}
    for i, p in enumerate(d_cycles):
        tbl.setdefault(p, []).append(i)

    lst = [-1] * len(d_cycles)
    for p, idx_list in tbl.items():
        idx_list2 = list(idx_list) + [idx_list[0]]
        for i, idx in enumerate(idx_list):
            lst[idx] = idx_list2[i + 1]
    return lst

def gen_primes():
    D = {}
    q = 2
    while True:
        if q not in D:
            yield q
            D[q*q] = [q]
        else:
            for p in D[q]:
                D.setdefault(p+q, []).append(p)
            del D[q]
        q += 1


def random_cycle_permutation(n):
    primes = list(itertools.islice(gen_primes(), n))
    print("primes:      %s" % primes)

    d_cycles = []
    for p in primes:
        d_cycles.extend([p] * p)

    print("d_cycles:    %s" % d_cycles)

    random.shuffle(d_cycles)
    print("shuffled:    %s" % d_cycles)

    perm = create_from_dcycle(d_cycles)
    print("permutation: %s" % perm)
    return perm

random_cycle_permutation(33)

