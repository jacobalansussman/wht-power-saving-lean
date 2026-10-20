"""km7.py (t-km): Kramer-Mesner search for C_p-invariant star decoders on P112 (copy gap 53-63), Hall stage.
Prescribed group G = <g> (g from gens.npz, order p, acting on ports and on F_2^9).  A G-invariant frame family =
one G-orbit of p frames of dim d (non-fixed flat, orbit representative) + a multiset of G-fixed flats, with total
star cost <= COST.  Necessary condition tested exactly: octad Hall (every one of the 1260 octads of P112 matches its
8 ports into 8 distinct frames containing them; G-invariance lets us test one octad per G-orbit).
usage: python -I km7.py p d COST   (writes km_p{p}_d{d}.json)"""
import sys, json, time, itertools
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
import numpy as np
from collections import Counter
p = int(sys.argv[1]); d = int(sys.argv[2]); COST = int(sys.argv[3]); KMIN = int(sys.argv[4]) if len(sys.argv) > 4 else 9
z = np.load('gens.npz'); perm = z['perm%d' % p].astype(np.int64)
inv = np.argsort(perm)
fixed = [(np.array(sorted(f)), dd) for f, dd in json.load(open('fixed_flats_g%d.json' % p))]
fixed = [(f, dd) for f, dd in fixed if dd < 9]
print('fixed flats (dim, ports):', [(dd, len(f)) for f, dd in fixed], flush=True)
# octad orbits under g
OC = [tuple(sorted(o)) for o in OCTP]; OCS = set(OC)
seen = set(); oreps = []
for o in OC:
    if o in seen: continue
    orb = set(); q = o
    for k in range(p):
        orb.add(q); q = tuple(sorted(perm[list(q)]))
    seen |= orb; oreps.append(o)
OR = np.array(oreps)
print('octads', len(OC), 'octad orbits', len(oreps), flush=True)
# all flats of dim d (distinct port sets whose span has dim d)
t0 = time.time()
ds, M = pool(d, minports=1)
M = np.unique(M, axis=0)
dims = np.array([f2dim([LAB[i] for i in np.nonzero(r)[0]]) for r in M])
M = M[dims == d]
print('flats of dim', d, ':', len(M), '%.0fs' % (time.time() - t0), flush=True)
# orbit representatives (non-fixed): canonical = lexicographically smallest packed mask among the p images
pk = lambda A: np.packbits(A, axis=1)
imgs = [M]
for k in range(1, p): imgs.append(imgs[-1][:, inv])      # port q in g^k F iff inv^k(q) in F
P8 = np.stack([pk(A) for A in imgs])                     # p x n x 14
keyv = [np.array([int.from_bytes(bytes(r), 'big') for r in P8[k]], dtype=object) for k in range(p)]
canon = np.array([min(keyv[k][i] for k in range(p)) for i in range(len(M))], dtype=object)
isfix = np.array([all(keyv[k][i] == keyv[0][i] for k in range(p)) for i in range(len(M))])
rep = np.array([keyv[0][i] == canon[i] and not isfix[i] for i in range(len(M))])
R = M[rep]
print('non-fixed orbit representatives:', len(R), ' fixed among them', int(isfix.sum()), flush=True)
# orbit port masks: W[r, q] = bitmask over the p orbit frames containing port q
W = np.zeros((len(R), NP), dtype=np.int64)
A = R.copy()
for k in range(p):
    W |= A.astype(np.int64) << k
    A = A[:, inv]
# fixed-flat combos (multisets) within budget
budget = COST - p * d
combos = []
def rec(start, cur, cost):
    if cur and len(cur) + p >= KMIN: combos.append(tuple(cur))
    for j in range(start, len(fixed)):
        if cost + fixed[j][1] <= budget: rec(j, cur + [j], cost + fixed[j][1])
rec(0, [], 0)
if p >= KMIN: combos.append(())
print('fixed combos within budget %d (K >= %d): %d' % (budget, KMIN, len(combos)), flush=True)
SUB = [[i for i in range(8) if q >> i & 1] for q in range(256)]
POP = np.array([bin(i).count('1') for i in range(1 << 16)], dtype=np.int64)
def deficits(Wc, octs):
    """Wc: n x 112 port bitmasks over frames; octs: m x 8.  returns n x m max Hall deficit"""
    pm = Wc[:, octs]                                    # n x m x 8
    N = np.zeros((256,) + pm.shape[:2], dtype=np.int64)
    best = np.full(pm.shape[:2], -8, dtype=np.int64)
    for q in range(1, 256):
        low = (q & -q).bit_length() - 1
        N[q] = N[q & (q - 1)] | pm[:, :, low]
        best = np.maximum(best, bin(q).count('1') - POP[N[q] & 0xFFFF] - POP[(N[q] >> 16) & 0xFFFF])
    return best
results = []
t0 = time.time()
for cb in combos:
    Wc = W.copy()
    for j, fi in enumerate(cb):
        bit = np.zeros(NP, dtype=np.int64); bit[fixed[fi][0]] = 1
        Wc |= bit[None, :] << (p + j)
    alive = np.arange(len(R)); mins = None
    # staged: octad reps in chunks, drop failing reps early
    total_def = np.zeros(len(R), dtype=np.int64)
    for a in range(0, len(OR), 30):
        if len(alive) == 0: break
        dd = np.zeros(len(alive), dtype=np.int64)
        for b in range(0, len(alive), 4000):
            idx = alive[b:b + 4000]
            df = deficits(Wc[idx], OR[a:a + 30]).clip(0).sum(1)
            dd[b:b + 4000] = df
        total_def[alive] += dd
        alive = alive[dd == 0]
    cost = p * d + sum(fixed[fi][1] for fi in cb)
    best_def = int(total_def.min()) if len(R) else None
    results.append(dict(fixed=[fixed[fi][1] for fi in cb], fixed_idx=list(cb), cost=cost, K=p + len(cb),
                        passing=len(alive), passing_reps=[np.nonzero(R[i])[0].tolist() for i in alive[:20]],
                        least_partial_deficit=best_def))
    print('combo fixed dims', [fixed[fi][1] for fi in cb], 'cost', cost, 'K', p + len(cb), 'Hall-passing reps', len(alive),
          'least (partial) deficit', best_def, '%.0fs' % (time.time() - t0), flush=True)
json.dump(dict(p=p, d=d, COST=COST, flats=int(len(M)), reps=int(len(R)), octad_orbits=len(oreps), results=results),
          open('km_p%d_d%d_c%d.json' % (p, d, COST), 'w'), indent=1)
