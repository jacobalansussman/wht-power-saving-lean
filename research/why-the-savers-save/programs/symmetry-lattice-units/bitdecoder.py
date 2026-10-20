"""Bit-side decoder for the BW16 unit (labels = the 2160 lines in Q^16, payload F_2), built from an ovoid of
O+(8,2): 9 pairwise anticommuting real Paulis P_1..P_9; every line is an eigenvector of exactly one of them, so the
18 eigenspaces V_{P,+-} (dimension 8, 120 lines each, isometric to E8) partition the lines.  In each block the E8 bit
decoder B0 = (1 + G_E8) mod 2 = U U^T (9 totals) is used; cross-block entries of B are 0, which is always allowed.
Checks: B fits the BW16 pattern exactly (diag 1, zero on every non-orthogonal pair), and the star cost = sum of
Q-ranks of the read sets.  Also: the F_2 rank of the assembled B."""
import sys, os, json, itertools, random, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import linesys

L = linesys.bw16(); v = len(L); N = 16
G = L @ L.T                            # norm 8 (support-8 states have entries +-1, support-2 entries +-2)
assert np.all(np.diag(G) == 8)

def pauli(a, b):
    """matrix of X^a Z^b on 4 qubits, basis index = bit pattern p (as in linesys), real."""
    M = np.zeros((N, N), dtype=np.int64)
    for p in range(N):
        sign = -1 if bin(p & b).count('1') % 2 else 1        # Z^b |p> = (-1)^{b.p} |p>
        M[p ^ a, p] = sign                                    # X^a |p> = |p ^ a>
    return M

paulis = []
for a in range(N):
    for b in range(N):
        if (a, b) == (0, 0): continue
        if bin(a & b).count('1') % 2: continue                 # X^a Z^b symmetric iff a.b even
        P = pauli(a, b)
        PL = L @ P.T                                           # rows = P psi
        plus = np.all(PL == L, axis=1); minus = np.all(PL == -L, axis=1)
        paulis.append(((a, b), plus, minus))
assert len(paulis) == 135, len(paulis)
cnt = [int(pl.sum() + mi.sum()) for _, pl, mi in paulis]
print('real Paulis:', len(paulis), 'lines that are eigenvectors of one Pauli:', set(cnt))

def anticomm(p, q):
    (a1, b1), (a2, b2) = p, q
    return (bin(a1 & b2).count('1') + bin(a2 & b1).count('1')) % 2 == 1

# find an ovoid: 9 pairwise anticommuting Paulis covering every line exactly once
idx = list(range(135)); best = None
def search(chosen, cov):
    global best
    if len(chosen) == 9:
        if np.all(cov == 1): best = list(chosen); return True
        return False
    start = chosen[-1] + 1 if chosen else 0
    for i in range(start, 135):
        if all(anticomm(paulis[i][0], paulis[j][0]) for j in chosen):
            c2 = cov + paulis[i][1] + paulis[i][2]
            if c2.max() > 1: continue
            if search(chosen + [i], c2): return True
    return False
search([], np.zeros(v, dtype=np.int64))  # takes ~2 min
assert best is not None
print('ovoid:', [paulis[i][0] for i in best])

sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-floor-checker')))
from copyfloor import f2_labels

def qrank(rows):
    if len(rows) == 0: return 0
    p = 2147483647
    M = np.array(rows, dtype=object) % p
    M = [list(map(int, r)) for r in M]
    rk = 0; nc = len(M[0])
    for c in range(nc):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % p), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; inv = pow(M[rk][c], p - 2, p)
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c] * inv % p; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[rk])]
        rk += 1
    r2 = np.linalg.matrix_rank(np.array(rows, dtype=float))
    assert r2 == rk
    return rk

terms = []; blocks = []
for i in best:
    for sel in (paulis[i][1], paulis[i][2]):
        B_idx = np.nonzero(sel)[0]; blocks.append(B_idx)
        Gb = G[np.ix_(B_idx, B_idx)] // 4                       # E8 scale: diag 2, off-diag 0, +-1 (block ips are 0, +-4)
        assert set(np.unique(Gb).tolist()) <= {-1, 0, 1, 2}
        B0 = ((1 + Gb) % 2).tolist()
        U = f2_labels(B0)                                       # 120 x 9, U U^T = B0 over F_2
        for c in range(len(U[0])):
            r = np.zeros(v, dtype=np.int64)
            for k, T in enumerate(B_idx): r[T] = U[k][c]
            terms.append((r, r.copy()))
print('blocks', len(blocks), 'sizes', sorted(set(len(b) for b in blocks)), 'totals', len(terms))
R = np.array([t[0] for t in terms]); Gm = np.array([t[1] for t in terms])
B = (Gm.T @ R) % 2
nonorth = (G != 0); np.fill_diagonal(nonorth, False)
fits = bool(np.all(np.diag(B) == 1) and not np.any(B[nonorth]))
costs = [qrank(L[np.nonzero(r)[0]].tolist()) for r, _ in terms]
cst = sum(costs)
# F_2 rank of B
rows = [int(''.join(map(str, row)), 2) for row in B.astype(int).tolist()]
basis = {}
for x in rows:
    while x:
        h = x.bit_length() - 1
        if h in basis: x ^= basis[h]
        else: basis[h] = x; break
out = dict(v=v, h=16, ovoid=[paulis[i][0] for i in best], totals=len(terms), fits=fits, cst=cst, x=cst / v,
           cost_per_total=sorted(set(costs)), F2_rank_B=len(basis))
print(out)
json.dump(out, open('bitdecoder_out.json', 'w'), indent=1)
np.savez_compressed('bitdecoder_terms.npz', R=R, Gm=Gm, L=L)
