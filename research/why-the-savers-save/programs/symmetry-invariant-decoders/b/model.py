"""model.py (d-R7) -- P112 in the u-ours/toy.py numbering (wht35 e8 harness: 9-point model, labels in F_2^9).
Ports P = E8 lines minus octad OCT[0] (the 8 lines (0,13,22,27,92,105,114,119)).  Orthogonal = Gram 0 = odd label
overlap.  W = 112 x 10 Lorentz root vectors (rank 8 on P); Lambda = row space of G[P,P].
Exports: P, NP (=112), LAB (labels of P as ints), W (np int), GP (Gram on P), ORTH (bool 112x112, off-diag),
OCTP (1260 octads inside P, as tuples of P-indices), lines_in(F_basis or ann) helpers, ann/subspace utilities."""
import sys, itertools
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), "../../e8-labels"))
import e8
import numpy as np

def par(x): return bin(x).count("1") & 1
V = e8.V
GRAM = np.array(e8.GRAM, dtype=np.int64)
def _cliques(adj, size, universe):
    out = []
    def rec(cur, cand):
        if len(cur) == size:
            out.append(tuple(cur)); return
        c = cand
        while c:
            low = c & -c; i = low.bit_length() - 1; c ^= low
            if bin(cand >> i).count("1") + len(cur) < size: break
            rec(cur + [i], cand & adj[i] & ~((low << 1) - 1))
    rec([], universe)
    return out
_OM = [sum(1 << T for T in range(V) if T != S and GRAM[S, T] == 0) for S in range(V)]
OCT_ALL = _cliques(_OM, 8, (1 << V) - 1)
assert len(OCT_ALL) == 2025
O0 = OCT_ALL[0]
P = [T for T in range(V) if T not in O0]
NP = len(P)
PIDX = {T: i for i, T in enumerate(P)}
LAB = [e8.PORTS[T] for T in P]
W = np.array([e8.VEC[T] for T in P], dtype=np.int64)
GP = GRAM[np.ix_(P, P)]
ORTH = (GP == 0)
np.fill_diagonal(ORTH, False)
OCTP = [tuple(PIDX[T] for T in oc) for oc in OCT_ALL if all(T in PIDX for T in oc)]
assert len(OCTP) == 1260
NORTH = [np.nonzero(ORTH[S])[0] for S in range(NP)]

def f2basis(xs):
    b = []
    for x in xs:
        for y in b: x = min(x, x ^ y)
        if x: b.append(x); b.sort(reverse=True)
    return b
def f2dim(xs): return len(f2basis(xs))
def ann(basis):
    """annihilator {h : h.b = 0 for all b} as a list of ints (including 0)"""
    return [h for h in range(512) if all(par(h & b) == 0 for b in basis)]
def lines_of_dual(hs):
    """ports of P whose label lies in the subspace {x : h.x = 0 for h in hs}"""
    return [i for i in range(NP) if all(par(LAB[i] & h) == 0 for h in hs)]
def lines_of_span(basis):
    sp = set([0])
    for b in basis:
        sp |= {x ^ b for x in sp}
    return [i for i in range(NP) if LAB[i] in sp]

def qrank(M, tol=1e-9):
    return int(np.linalg.matrix_rank(np.asarray(M, dtype=float), tol=tol)) if np.size(M) else 0

_u, _s, _vt = np.linalg.svd(GP.astype(float))
LAM = _u[:, :8] * _s[:8]                              # 112 x 8, columns span Lambda = rowspace(G[P,P])

def natural_test(B, tol=1e-6):
    """is there a diagonal D with Lambda D inside rowspace(B)?  Lambda D = {T -> d_T <z, w_T>}.
    Returns (numerical rank r of B, sigma_min/sigma_max of the linear system for d (small = solution exists),
    next sigma ratio, min|d_T|/max|d_T| of that solution)."""
    U, s, Vt = np.linalg.svd(B)
    r = int((s > tol * s[0]).sum())
    Q = Vt[r:]
    A = np.einsum('it,tc->ict', Q, LAM).reshape(-1, NP)
    u, s2, vt = np.linalg.svd(A, full_matrices=False)
    d = vt[-1]
    return r, float(s2[-1] / s2[0]), float(s2[-2] / s2[0]), float(np.abs(d).min() / np.abs(d).max())
