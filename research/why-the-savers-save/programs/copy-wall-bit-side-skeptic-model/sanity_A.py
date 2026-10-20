"""sanity_A.py (v-g-9hyp-model) -- adversarial sampling for Theorem A: W = W0 + <x> (dim 10) and W0 + <x, y> (dim 11)
with x, y = r(q) cut by a root-orthogonal 6-flat or other small flats (light vectors), exact min-weight basis (greedy
over all vectors of W, weight = Q-rank of support).  Sanity only; it cannot prove anything."""
import sys, random, collections; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import numpy as np, geo9
name = sys.argv[1]; N = int(sys.argv[2])
P, live, W0b = geo9.unit(name); v = len(live); LA = geo9.LA.astype(float); G = geo9.GR
W0all = [x for x in geo9.f2_span_all(W0b) if x]
flats = []
for l in range(120):
    for m in range(l + 1, 120):
        flats.append(sum(1 << t for t in range(v) if G[live[t], l] == 0 and G[live[t], m] == 0))
flats = list(set(flats))
cache = {}
def w(x):
    if x not in cache:
        idx = [live[t] for t in range(v) if (x >> t) & 1]
        cache[x] = int(np.linalg.matrix_rank(LA[idx])) if idx else 0
    return cache[x]
def minbasis(gens):
    allv = [x for x in geo9.f2_span_all(gens) if x]
    allv.sort(key=w); basis = {}; cost = 0; ws = []
    for x in allv:
        y = x
        while y:
            h = y.bit_length() - 1
            if h in basis: y ^= basis[h]
            else: basis[h] = y; break
        if y: cost += w(x); ws.append(w(x))
        if len(basis) == len(gens): break
    return cost, sorted(ws)
rng = random.Random(3)
def light():
    k = rng.random()
    q = rng.choice(W0all); F = rng.choice(flats)
    if k < 0.5: return q & F                     # piece of a label vector inside an E6/D6 flat
    if k < 0.8: return F & rng.getrandbits(v)    # random subset of a 6-flat
    return 1 << rng.randrange(v)                 # a single line
h10 = collections.Counter(); h11 = collections.Counter(); best = [99, 99]
for it in range(N):
    x = light(); y = light()
    if x == 0 or y == 0: continue
    c10, ws10 = minbasis(W0b + [x]); h10[c10] += 1
    if c10 < best[0]: best[0] = c10; print('dim10 new best', c10, ws10, flush=True)
    c11, ws11 = minbasis(W0b + [x, y])
    if len(ws11) == 11:
        h11[c11] += 1
        if c11 < best[1]: best[1] = c11; print('dim11 new best', c11, ws11, flush=True)
print(name, 'dim10 cost histogram', sorted(h10.items())); print(name, 'dim11 cost histogram', sorted(h11.items()))
