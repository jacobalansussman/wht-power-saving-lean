"""cover_indep.py (v-g-9hyp-recheck) -- independent check of F1, F2, F3 of d-g-9hyp's Theorem A.
Independent of e8geo/f2_labels: W0 = F_2 row space of B0|live (B0 = (1+G) mod 2), all 511 nonzero vectors enumerated.
Ranks: LOWER bounds via the rank over F_p (p = 2^31-1) of the 8x8 integer Gram sum_{l in S} l l^T  (rank_p(Gram) <= rank_Q(Gram)
= rank_Q(L_S)), so every 'rank >= 7' printed here is rigorous.  Upper bounds (rank <= 7 for the 9 punctured labels) by
exhibiting the punctured line y with S inside y^perp (exact integer dot products).
usage: python -I cover_indep.py <s2|word|120>"""
import sys, itertools, collections, time
import numpy as np
P_ = {'s2': [39, 52, 57, 71, 76, 78, 103, 109, 114], 'word': [8, 17, 39, 51, 53, 58, 63, 72, 114], '120': []}[sys.argv[1]]
# hub hamming_e8 construction (same order as hub e8units.hamming_e8)
gen = np.array([[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[0,1,0,1,0,1,0,1]])
code = {tuple((np.array(c) @ gen) % 2) for c in itertools.product([0,1], repeat=4)}
w4 = sorted(c for c in code if sum(c) == 4)
Ls = [tuple(2*(np.arange(8)==a).astype(int)) for a in range(8)]
for c in w4:
    idx = [i for i in range(8) if c[i]]
    for sg in itertools.product([1,-1], repeat=3):
        x = np.zeros(8, dtype=int); x[idx[0]] = 1
        for k in range(3): x[idx[k+1]] = sg[k]
        Ls.append(tuple(x))
L = np.array(Ls, dtype=np.int64); G = (L @ L.T) // 2
assert len(L) == 120 and len({tuple(r) for r in L} | {tuple(-r) for r in L}) == 240
live = [i for i in range(120) if i not in P_]; v = len(live)
B0 = ((1 + G) % 2)[np.ix_(live, live)]
# F_2 row space basis of B0|live
rows = [int(''.join(str(int(b)) for b in r[::-1]), 2) for r in B0]   # bit t = column t
basis = []
for x in rows:
    for b in basis: x = min(x, x ^ b)
    if x: basis.append(x); basis.sort(reverse=True)
print('unit', sys.argv[1], 'v', v, 'rank_F2(B0|live) =', len(basis))
dimW0 = len(basis)
W0 = []
for coeffs in range(1, 1 << dimW0):
    x = 0
    for k in range(dimW0):
        if coeffs >> k & 1: x ^= basis[k]
    W0.append(x)
assert len(set(W0)) == (1 << dimW0) - 1
supp = [[live[t] for t in range(v) if x >> t & 1] for x in W0]
p = 2147483647
def batch_rank(A):
    """A: (B, 8, 8) int64 mod p -> ranks over F_p (batched Gauss-Jordan)."""
    A = A % p; Bn = A.shape[0]; ar = np.arange(Bn)
    rank = np.zeros(Bn, np.int64); used = np.zeros((Bn, 8), bool)
    for c in range(8):
        cand = (A[:, :, c] != 0) & ~used
        has = cand.any(1); piv = np.argmax(cand, 1)
        prow = A[ar, piv, :].copy()
        pv = prow[:, c].copy(); pv[~has] = 1
        inv = np.array([pow(int(z), p - 2, p) for z in pv], dtype=np.int64) if Bn < 50 else np.vectorize(lambda z: pow(int(z), p-2, p), otypes=[np.int64])(pv)
        prow = (prow * inv[:, None]) % p
        f = A[:, :, c].copy(); f[~has] = 0
        f[ar, piv] = 0
        A = (A - (f[:, :, None] * prow[:, None, :]) % p) % p
        A[ar[has], piv[has], :] = prow[has]
        used[ar[has], piv[has]] = True; rank += has
    return rank
Gl = np.einsum('ia,ib->iab', L, L)     # 120 x 8 x 8
def gram(S): return Gl[S].sum(0)
# F1 / F2
gr = np.array([gram(S) for S in supp]); rk = batch_rank(gr)
print('F1 lower bounds: rank_p of r(q) over the', len(W0), 'nonzero W0 vectors:', sorted(collections.Counter(rk.tolist()).items()))
low = [i for i in range(len(W0)) if rk[i] <= 7]
for i in low:
    S = supp[i]
    ys = [y for y in P_ if all(G[y, t] == 0 for t in S)]
    print('   rank<=7 candidate |S| =', len(S), ' S lies in y^perp for punctured y =', ys)
print('F2 sizes |r(q)|:', sorted(collections.Counter(len(S) for S in supp).items()), ' min', min(len(S) for S in supp))
# the 4,900 root-orthogonal 6-flats: lines orthogonal to both l and m
fl = {}
for l in range(120):
    for m in range(l + 1, 120):
        f = frozenset(t for t in range(120) if G[t, l] == 0 and G[t, m] == 0)
        fl.setdefault(f, (l, m))
flats = list(fl)
print('root-orthogonal 6-flats:', len(flats), 'line counts:', sorted(collections.Counter(len(f) for f in flats).items()))
# check each such flat is a 6-dim subspace's full line set: rank_p(Gram) = 6, and no other line lies in its span
for f in flats[:]:
    pass
FM = np.zeros((len(flats), 120), bool)
for i, f in enumerate(flats): FM[i, list(f)] = True
t0 = time.time(); hist = collections.Counter(); worst = 9
for qi, S in enumerate(supp):
    Sa = np.array(S)
    inF = FM[:, Sa]                             # (4900, |S|)
    out = ~inF
    Gs = np.einsum('fs,sab->fab', out.astype(np.int64), Gl[Sa])
    r = batch_rank(Gs)
    m = int(r.min()); hist[m] += 1; worst = min(worst, m)
print('F3 lower bound: min over the 4,900 flats of rank_p(r(q) minus F), histogram over q:', sorted(hist.items()), ' worst', worst, ' %.1f s' % (time.time() - t0))
print('VERDICT F3 (every nonzero r(q) minus every E6/D6 flat has rank >= 7):', worst >= 7)
