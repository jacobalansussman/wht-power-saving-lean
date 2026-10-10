"""pool1.py -- the sums of saved hierarchies, how invariant they are under the order-64 group, pool sizes.
usage: pool1.py <out pool.json> <hierarchy.json> ..."""
import sys, time
sys.dont_write_bytecode = True
from xlib import *
from collections import Counter
t0 = time.time()
A = Act("F8uu")
print("group elements:", A.n)
pool = set()
for path in sys.argv[2:]:
    h = load_h(path)
    nd = needed_set(h)
    syms = set(sym_of(h.content[u]) for u in nd if u >= V)
    inv = sum(1 for s in syms if A.orbit(s) <= syms)
    # is the use by targets invariant?
    fl = set((S, sym_of(h.content[u])) for S in range(V) for u in h.need[S])
    finv = 0
    for S, s in fl:
        im = A.images(s)
        if all((A.els[gi][0][S], im[gi]) in fl for gi in range(A.n)): finv += 1
    print(path, proxy(h), "sums %d, of which with the whole orbit in the design %d; uses (target, sum) %d, with whole orbit %d" % (len(syms), inv, len(fl), finv))
    pool |= syms
print("union of sums:", len(pool))
clo = set()
for s in pool: clo |= A.orbit(s)
print("closed under the group:", len(clo), " %.0fs" % (time.time() - t0))
orbs = {}
for s in clo:
    r = min(A.images(s)); orbs.setdefault(r, set()).add(s)
print("orbits:", len(orbs), "sizes", sorted(Counter(len(o) for o in orbs.values()).items()))
print("sources per sum:", sorted(Counter(len(bits(s[0])) for s in clo).items()))
print("targets per sum:", sorted(Counter(len(targets_of(s)) for s in clo).items()))
hio.dump_pool(sorted(clo), sys.argv[1])
print("%.0fs" % (time.time() - t0))
