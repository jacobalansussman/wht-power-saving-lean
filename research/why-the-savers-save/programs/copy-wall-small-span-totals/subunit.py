"""subunit.py (d-R7): small case.  Sub-unit U = ports of P112 whose labels lie in a 7-dim subspace V of F_2^9
(36 ports: the largest).  (1) natural floor of U (linear totals; frames = subspaces of V, exhaustive: 29,212
subspaces); (2) for frame multisets below that floor: Hall on the octads inside U (annealing over subspaces of V),
then masked ALS (restarts).  usage: subunit.py h1 h2 restarts maxsecs"""
import sys, time, itertools, collections; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
from als import als
h1, h2, NR, MAXS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
t0 = time.time(); rng = np.random.default_rng(h1 * 1000 + h2)
Vann = [h1, h2]
U = np.array(lines_of_dual(Vann))
_O = ORTH[np.ix_(U, U)]
DROP = __import__('os').environ.get('DROPISO', '0') == '1'
if DROP: U = U[_O.sum(1) < len(U) - 1]          # drop ports orthogonal to every other port (they cost 1 each, apart)
nU = len(U)
Vb = f2basis([x for x in range(512) if par(x & h1) == 0 and par(x & h2) == 0])
assert len(Vb) == 7
OU = ORTH[np.ix_(U, U)]
uidx = {int(t): i for i, t in enumerate(U)}
OCTU = [tuple(uidx[t] for t in oc) for oc in OCTP if all(t in uidx for t in oc)]
if DROP:   # maximal orthogonal cliques inside U (size 7 after dropping the isolated port), padded to 8 with a dummy
    OCTU = sorted({tuple(sorted(uidx[t] for t in oc if t in uidx)) for oc in OCTP if sum(t in uidx for t in oc) == 7})
LAMU = LAM[U]; rU = qrank(LAMU)
print('sub-unit: duals (%d,%d), %d ports, root rank %d, octads inside %d, orth degree %s' % (h1, h2, nU, rU, len(OCTU), sorted(collections.Counter(OU.sum(1)).items())))
# subspaces of V: all subsets spanned; enumerate by span of combinations of basis (coordinates in F_2^7)
def subspaces_of_V():
    out = {}
    vecs = [sum(Vb[i] for i in range(7) if c >> i & 1) for c in range(128)]
    # enumerate subspaces of F_2^7 by reduced echelon form
    for d in range(1, 8):
        for piv in itertools.combinations(range(7), d):
            frees = [[j for j in range(p + 1, 7) if j not in piv] for p in piv]
            for bits in itertools.product(*[range(1 << len(f)) for f in frees]):
                rows = [(1 << p) | sum(1 << j for i, j in enumerate(f) if b >> i & 1) for p, f, b in zip(piv, frees, bits)]
                span = {0}
                for r in rows: span |= {x ^ r for x in span}
                labs = {vecs[c] for c in span}
                mask = np.array([LAB[t] in labs for t in U])
                out.setdefault(d, []).append(mask)
    return {d: np.array(v) for d, v in out.items()}
SUB = subspaces_of_V()
print('subspaces of V by dim:', {d: len(v) for d, v in SUB.items()}, '%.0fs' % (time.time() - t0))
OUTER = np.einsum('ti,tj->tij', LAMU, LAMU).reshape(nU, 64)
u_, s_, vt_ = np.linalg.svd(LAMU); SP = vt_[:rU]
acc = np.zeros((0, 8)); cost = 0; per = []
for d in range(1, 8):
    C = ~SUB[d]
    G = (C.astype(float) @ OUTER).reshape(-1, 8, 8)
    Gs = np.einsum('ai,nij,bj->nab', SP, G, SP)
    ev = np.linalg.eigvalsh(Gs); low = ev[:, 0] < 1e-8 * max(1.0, ev[:, -1].max())
    new = []
    for n in np.nonzero(low)[0]:
        w, v = np.linalg.eigh(Gs[n]); new.append(v[:, w < 1e-8 * max(1.0, w[-1])].T @ SP)
    before = qrank(acc) if len(acc) else 0
    if new: acc = np.vstack([acc] + new)
    after = qrank(acc) if len(acc) else 0
    cost += d * (after - before); per.append((d, after - before))
    if after >= rU: break
NATF = cost
print('natural floor of U: %d  (totals per frame dim: %s)' % (NATF, per), flush=True)
OCTA_U = np.array(OCTU)
def hall_def(S):
    K = S.shape[0]
    w = (S.astype(np.int64) << np.arange(K)[:, None]).sum(0); pm = w[OCTA_U]
    L = OCTA_U.shape[1]
    N = np.zeros((1 << L, len(OCTA_U)), dtype=np.int64)
    for q in range(1, 1 << L):
        low = (q & -q).bit_length() - 1; N[q] = N[q & (q - 1)] | pm[:, low]
    return (SUBSZ[:1 << L, None] - POPC[N]).max(0).clip(0).sum() + int((S.sum(0) == 0).sum())
def shapes_below(nat):
    out = []
    lo = int(__import__('os').environ.get('KLO', '8')); hi = int(__import__('os').environ.get('KHI', '11'))
    dl = int(__import__('os').environ.get('DLO', '3')); gap = int(__import__('os').environ.get('GAP', '6'))
    for K in range(lo, hi + 1):
        for combo in itertools.combinations_with_replacement(range(dl, 8), K):
            if sum(combo) < nat + int(__import__('os').environ.get('EQ', '0')) and sum(combo) >= nat - gap: out.append(combo)
    return out
shp = shapes_below(NATF)
print('shapes tried (cost in [floor-6, floor-1], K = 8..11, dims 3..7):', len(shp), flush=True)
best_overall = (1e9, None)
for dims in sorted(shp, key=lambda c: -sum(c)):
    if time.time() - t0 > MAXS: print('time box reached'); break
    cur = [rng.integers(len(SUB[d])) for d in dims]
    S = np.array([SUB[d][i] for d, i in zip(dims, cur)]); f = hall_def(S); T = 2.0
    for step in range(4000):
        j = rng.integers(len(dims)); d = dims[j]; new = rng.integers(len(SUB[d]))
        S2 = S.copy(); S2[j] = SUB[d][new]; f2 = hall_def(S2)
        if f2 <= f or rng.random() < np.exp((f - f2) / T): S = S2; f = f2
        T = max(0.05, T * 0.998)
        if f == 0: break
    if f: continue
    res = min(als(S, OU, iters=int(__import__('os').environ.get('IT', '800')), rng=rng)[0] for _ in range(NR))
    best_overall = min(best_overall, (res, dims))
    print('  shape %s cost %d: Hall OK; masked ALS best residual %.3e' % (dims, sum(dims), res), flush=True)
print('best below the natural floor:', best_overall, '%.0fs' % (time.time() - t0))
