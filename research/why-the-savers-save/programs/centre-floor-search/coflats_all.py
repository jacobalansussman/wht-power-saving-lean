#!/usr/bin/env python3
"""ctr-search / coflats_all.py: every hyperplane flat of the 112 lines with >= 32 lines, found from normals v in the
E8 lattice (Hamming coordinates: x in Z^8, x mod 2 in the extended Hamming code) of norm <= 24, deduplicated by
the flat itself.  Every hyperplane spanned by lines has a normal in the (self-dual) lattice; rank-7 root
subsystems with >= 16 roots in E8 have normals of small norm, the check below also confirms by random spanning
7-sets that no flat outside the list has >= 32 lines.  Output: co-flats (complements) of size <= 80, as bitmasks."""
import sys, os, itertools, json, random, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *
code = set()
for t in TET: code.add(sum(1 << (7 - i) for i in t))
code |= {0, 0xFF}
def in_lattice(x):
    w = sum((1 << (7 - i)) for i in range(8) if x[i] % 2)
    return w in code
R = np.array(LINES, dtype=np.int64)
flats = {}
rng = range(-4, 5)
for x in itertools.product(range(-4, 5), repeat=8):
    nx = sum(v * v for v in x)
    if nx == 0 or nx > 24: continue
    # canonical sign: first nonzero positive; skip non-primitive multiples later via flat dedup
    fnz = next(v for v in x if v != 0)
    if fnz < 0: continue
    if not in_lattice(x): continue
    v = np.array(x, dtype=np.int64)
    m = 0
    z = (R @ v) == 0
    for b in np.nonzero(z)[0].tolist(): m |= 1 << b
    c = bin(m).count('1')
    if c >= 32 and np.linalg.matrix_rank(R[z].astype(float)) == 7:
        if m not in flats: flats[m] = (c, x, nx)
sizes = sorted(set(v[0] for v in flats.values()), reverse=True)
import collections
print('hyperplane flats with >= 32 lines:', len(flats), 'size histogram', sorted(collections.Counter(v[0] for v in flats.values()).items(), reverse=True))
# random check: spanning 7-sets -> flat; flag any flat >= 32 not in the list
r = random.Random(3); miss = 0
Rf = R.astype(float)
for _ in range(30000):
    S = r.sample(range(N), 7)
    if np.linalg.matrix_rank(Rf[S]) < 7: continue
    u, s, vt = np.linalg.svd(Rf[S]); v = vt[-1]
    z = np.abs(Rf @ v) < 1e-9
    m = 0
    for b in np.nonzero(z)[0].tolist(): m |= 1 << b
    if bin(m).count('1') >= 32 and m not in flats: miss += 1
print('random spanning 7-sets giving a >=32-line flat not in the list:', miss)
full = (1 << N) - 1
cof = sorted(((full & ~m) for m in flats), key=lambda m: bin(m).count('1'))
pickle.dump(cof, open(os.path.join(HERE, 'out/coflats80.pkl'), 'wb'))
print('co-flats of size <= 80:', len(cof), 'sizes', sorted(collections.Counter(bin(m).count('1') for m in cof).items()))
