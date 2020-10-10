import itertools, random, functools
import math
from datetime import datetime
import inspect
from sympy.combinatorics import Permutation

def random_cycle_permutation(n):
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

    primes = list(itertools.islice(gen_primes(), n))
    #print("primes:      %s" % primes)

    d_cycles = []
    for p in primes:
        d_cycles.extend([p] * p)

    #print("d_cycles:    %s" % d_cycles)

    random.shuffle(d_cycles)
    #print("shuffled:    %s" % d_cycles)

    perm = create_from_dcycle(d_cycles)
    #print("permutation: %s" % perm)
    return perm


# x  = Permutation(random_cycle_permutation(3))(1,2)(5,6,7)
# _A = Permutation(random_cycle_permutation(3))

# print("x =      " + str(x.cyclic_form))
# print("A =      " + str(_A.cyclic_form))
# print("xAx^-1 = " + str((x * _A * ~x).cyclic_form))

# invert of a in mod p
def mod_inv(a, p):
  a = a % p
  d = p
  x = 0
  s = 1
  while (a != 0):
    q = d // a

    r = d % a
    d = a
    a = r

    t = x - q * s
    x = s
    s = t
  return (x + p) % (p // d)


def cycle_pair(c_p1, c_p2, verbose = False):
    ret = []
    for c1 in c_p1:
        if verbose:
            print("---")
        c1_len = len(c1)
        c1 = Permutation([c1])
        c = c1
        for i in range(1, c1_len + 1):
            c_c = c.cyclic_form
            if verbose:
                print(c_c)
            if all([c in c_p2 for c in c_c]):
                ret.append((i % c1_len, c1_len))
                break
            c = c * c1
    return ret

def sim_con(m, i_p, verbose = False):
    ret = 0
    for i, p in i_p:
        _gcd = gcd(i,p)
        if i > 0 and _gcd > 1:
            i //= _gcd
            p //= _gcd
        _P = m//p
        _A = mod_inv(_P, p)
        if verbose:
            print("x = %d mod %d" % (i, p))
            # print("%d * %d = %d mod %d" % (_P, _A, _P*_A % p, p))
        ret = (ret + i * _P * _A) % m
    return ret

def brute(i_p, verbose = False):
    def lcm(iter):
        return functools.reduce(lambda a,b: (a*b)//gcd(a,b), iter, 1)
    m = lcm(map(lambda a: a[1], i_p))
    goal = list(map(lambda a: a[0], i_p))
    count = [1] * len(i_p)
    for x in range(1, m):
        if verbose:
            print(count, goal)
        if count == goal:
            return x
        for j,ip in enumerate(i_p):
            i, p = ip
            count[j] = (count[j] + 1) % p
    return m


def perm_log(x, base):
    i_p = cycle_pair(base.cyclic_form, x.full_cyclic_form)
    timestamp()

    print("cycle_pair: " + str(i_p))
    timestamp()

    m = functools.reduce(lambda a,b: (a*b), map(lambda a: a[1], i_p), 1)
    ret = sim_con(m, i_p, True  )
    # ret = brute(i_p)
    timestamp()

    return m, ret

start = datetime.now()
def timestamp():
    frame = inspect.stack()[1]
    #print("%s %s" % (frame.lineno, datetime.now()-start))

timestamp()

x  = Permutation(10,9,8,7)(6,5,4,3)(1,2)
timestamp()

n = 6
y = x**n
timestamp()

print("%s ^ %d = %s" % (x.cyclic_form, n, y.cyclic_form))
m, log = perm_log(y, x)
print("kotae awase %d == %d mod %d" % (n % m, log , m))
if n % m == log:
    print("yatta!")
else:
    print("baka")

print(mod_inv(3,5))

for n in range(1, 20):
     y = x**n
     print(y.array_form, y.full_cyclic_form)
