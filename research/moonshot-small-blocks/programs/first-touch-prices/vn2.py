# vn2.py -- w7-mixed3: f -> infinity (regular-representation) over-demand of TWO-TERM first entries D(b1 + c b2) into a
# two-row slot, hook certificate Y = N in V (x) Lambda^4 V (f = 5).  At level k the entry's constraint on E_R (x) V^(x)k is
# sum_b beta_b A_b (x) b^(x)k (A_b = pi_{Y'} rho(D b)|_{E_R}); normalised by 8^k its rank tends to the regular-rep rank
#   delta(c) = (1/|C|) sum_{omega in C} rank(A_1 + c omega A_2),
# C = the characters of <g> (g = b1^-1 b2) that occur in V^(x)k for the residue of k (all n-th roots if no power of g is a
# non-trivial scalar).  We take, for each pair, the residue class of k that is WORST for the program (largest delta), and
# minimise over c: candidates c = t*/omega for the drop points t* of the pencil A_1 + t A_2.
import sys, os, itertools, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from blk3 import D, RM, KM, rand_engine, rand_widened, mono, I8

S4 = list(itertools.combinations(range(8), 4))
U4 = np.array(S4)
BASIS = [(x, T) for x in range(8) for T in S4]
NIDX = np.array([i for i, (x, T) in enumerate(BASIS) if x in T]); OIDX = np.array([i for i, (x, T) in enumerate(BASIS) if x not in T])

def lam4(g):
    sub = g[U4[:, None, :, None], U4[None, :, None, :]]
    return np.linalg.det(sub)

def rho(g): return np.kron(g, lam4(g))

def rank(M, tol=1e-8):
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > tol * max(1.0, s[0])).sum())

def nullsp(M, tol=1e-9):
    u, s, vh = np.linalg.svd(M)
    r = int((s > tol * max(1.0, s[0])).sum())
    return vh[r:].conj().T

TGT = {'R': RM, 'K': KM}
_E = {}
def E_of(t):
    if t not in _E:
        R = rho(TGT[t])
        Z = nullsp(R[np.ix_(OIDX, NIDX)])            # coordinates in N
        _E[t] = Z
    return _E[t]

def A_of(b, t='R', z=1):
    R = rho(D[z] @ b)
    return R[np.ix_(OIDX, NIDX)] @ E_of(t)

def order_info(g, nmax=200):
    """n = order of g; j = least power with g^j scalar, zeta = that scalar"""
    P = np.eye(8, dtype=complex); j = None; zeta = None
    for e in range(1, nmax + 1):
        P = P @ g
        if j is None and np.allclose(P, P[0, 0] * np.eye(8), atol=1e-9):
            j, zeta = e, P[0, 0]
        if np.allclose(P, np.eye(8), atol=1e-9): return e, j, zeta
    return None, j, zeta

def char_classes(n, j, zeta):
    """list of character sets C_r (one per residue r of k mod n): omega^n = 1 with omega^j = zeta^k"""
    roots = [np.exp(2j * np.pi * a / n) for a in range(n)]
    out = []
    for k in range(n):
        target = zeta ** k
        out.append([w for w in roots if abs(w ** j - target) < 1e-8])
    return out

def delta_pair(b1, b2, t='R', z=1):
    g = np.linalg.inv(b1) @ b2
    n, j, zeta = order_info(g)
    if n is None: return None
    A1 = A_of(b1, t, z); A2 = A_of(b2, t, z)
    r1, r2 = rank(A1), rank(A2)
    rng = np.random.default_rng(0)
    tg = rng.normal() + 1j * rng.normal()
    rgen = rank(A1 + tg * A2)
    # drop points of the pencil: generalized eigenvalues of a random square compression
    P = rng.normal(size=(A1.shape[0], rgen)) + 1j * rng.normal(size=(A1.shape[0], rgen))
    Qm = rng.normal(size=(A1.shape[1], rgen)) + 1j * rng.normal(size=(A1.shape[1], rgen))
    from scipy.linalg import eigvals
    ev = eigvals(P.conj().T @ A1 @ Qm, -(P.conj().T @ A2 @ Qm))
    drops = []
    for e in ev:
        if not np.isfinite(e) or abs(e) > 1e8 or abs(e) < 1e-10: continue
        rr = rank(A1 + e * A2)
        if rr < rgen: drops.append((e, rr))
    # merge numerically equal drop points
    md = []
    for e, rr in drops:
        if not any(abs(e - f) < 1e-6 * max(1, abs(e)) for f, _ in md): md.append((e, rr))
    classes = char_classes(n, j, zeta)
    best = None
    for C in classes:
        if not C: continue
        cands = [None] + [e / w for (e, _) in md for w in C]
        bc = len(C) * rgen; bcand = None
        for c in cands:
            if c is None: continue
            tot = sum(rank(A1 + c * w * A2) for w in C)
            if tot < bc: bc, bcand = tot, c
        val = bc / len(C)
        if best is None or val > best[0]: best = (val, len(C), bcand)
    return dict(n=n, j=j, r1=r1, r2=r2, rgen=rgen, ndrops=len(md), drops=[(complex(e), r) for e, r in md][:6],
                delta=best[0], nchars=best[1], c=None if best[2] is None else complex(best[2]))

if __name__ == '__main__':
    kind = sys.argv[1]; npairs = int(sys.argv[2]); seed = int(sys.argv[3]); t = sys.argv[4] if len(sys.argv) > 4 else 'R'
    rng = np.random.default_rng(seed)
    from fam import POOL
    for it in range(npairs):
        if kind == 'engine': b1, b2 = rand_engine(rng), rand_engine(rng)
        elif kind == 'widened': b1, b2 = rand_widened(rng), rand_widened(rng)
        elif kind == 'pool1': b1 = I8; b2 = POOL[int(rng.integers(1, len(POOL)))]
        elif kind == 'pool2':
            b1 = POOL[int(rng.integers(len(POOL)))] @ POOL[int(rng.integers(len(POOL)))]
            b2 = b1 @ POOL[int(rng.integers(1, len(POOL)))] @ POOL[int(rng.integers(len(POOL)))]
        t0 = time.time()
        d = delta_pair(b1, b2, t)
        if d is None: continue
        d['kind'] = kind; d['t'] = t; d['secs'] = round(time.time() - t0, 2)
        d['drops'] = [(round(e.real, 6), round(e.imag, 6), r) for (e, r) in d['drops']]
        if d['c'] is not None: d['c'] = (round(d['c'].real, 6), round(d['c'].imag, 6))
        print(json.dumps(d), flush=True)
