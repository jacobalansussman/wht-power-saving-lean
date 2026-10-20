"""vrec.py (v-rm2-exact): skeptic's own re-implementation of the star recursion + pooled zero blocks + clique +
matching + |Y|-W tests for the per-array-target system at block 5, one move per direction, canonical targets.
Written from scratch on top of vf.py (own field / forms / parities). Not reading any rm2-exact table or certificate.

Facts (re-derived by hand in ../NOTE.md):
  T10  rank E[D_Q,D_Q] = |D_Q| - r_Q.
  STAR E[X,X] of nullity nu, X = P + Q_1..Q_g  =>  (g-2) rank E[P,P] <= sum_j rank E[P u Q_j] - 2|P| + 2 nu.
  CLQ  E[Y,Y] = 0  =>  |Y| <= W.          RKY  rank E[Y,Y] >= |Y| - W for every Y.
  MAT  rank E[D,D] <= max matching of D x D avoiding known zero cells.
usage: vrec.py [--maxd 16] [--rounds 3] [--allzero] [--truecheck] < lines "k1,k2,..."   (prints a verdict per line)
"""
import sys, itertools, time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
import vf

MASK, WT, PAR = vf.MASK, vf.WT, vf.PAR
NF = vf.NF
PAT = (PAR.astype(np.int64) << np.arange(31)[None, :]).sum(axis=1)
pats = sorted(set(PAT.tolist())); CLS = {p: i for i, p in enumerate(pats)}
VID = np.array([CLS[p] for p in PAT.tolist()], dtype=np.int64)
NV = len(pats)
CLSPAR = np.array([[(p >> k) & 1 for k in range(31)] for p in pats], dtype=np.int64)   # NV x 31

# ---- affine flats of F_2^5 as 32-bit point masks (bit p = point p) ----
def span(gens):
    S = {0}
    for g in gens: S |= {x ^ g for x in S}
    return S
FLATS = []
for d in range(5):
    subs = set()
    for gens in itertools.combinations(range(1, 32), d):
        S = span(gens)
        if len(S) == 1 << d: subs.add(sum(1 << x for x in S))
    fl = set()
    for V in subs:
        pts = [x for x in range(32) if V >> x & 1]
        for a in range(32): fl.add(sum(1 << (x ^ a) for x in pts))
    FLATS.append(sorted(fl))
assert [len(f) for f in FLATS] == [32, 496, 1240, 620, 62], [len(f) for f in FLATS]
FIDX = [{A: i for i, A in enumerate(fl)} for fl in FLATS]
UP = []
for d in range(4):
    g = (1 << (5 - d)) - 1; rows = []
    for A in FLATS[d]:
        pts = [x for x in range(32) if A >> x & 1]; a0 = pts[0]; Bs = set()
        for p in range(32):
            if A >> p & 1: continue
            w = p ^ a0
            Bs.add(A | sum(1 << (x ^ w) for x in pts))
        assert len(Bs) == g
        rows.append([FIDX[d + 1][B] for B in Bs])
    UP.append(np.array(rows, dtype=np.int64))
DIRMASK = [np.array([(A >> 1) & 0x7fffffff for A in fl], dtype=np.int64) for fl in FLATS]  # bit z-1 <-> direction z

def popc(x):
    x = np.asarray(x, dtype=np.int64); c = np.zeros_like(x)
    for b in range(31): c += (x >> b) & 1
    return c

_pre = {}
def pre(F):
    """target-independent: for every flat A and class v, min weight of a form G != F of class v with D_G >= D_F n A."""
    if F in _pre: return _pre[F]
    A = int(MASK[F]); dirs = [z for z in range(1, 32) if A >> (z - 1) & 1]; d = len(dirs)
    G = np.arange(1, NF); G = G[G != F]
    u = np.zeros(len(G), dtype=np.int64)
    for t, z in enumerate(dirs): u |= ((MASK[G] >> (z - 1)) & 1) << t
    T = np.full((1 << d, NV), 255, dtype=np.int16)
    np.minimum.at(T, (u, VID[G]), WT[G].astype(np.int16))
    for b in range(d):
        Tv = T.reshape(-1, 2, 1 << b, NV)
        np.minimum(Tv[:, 0], Tv[:, 1], out=Tv[:, 0])
    out = []
    for dd in range(5):
        S = DIRMASK[dd] & A
        idx = np.zeros(len(S), dtype=np.int64)
        for t, z in enumerate(dirs): idx |= ((S >> (z - 1)) & 1) << t
        out.append((T[idx].astype(np.int16), popc(S), S))
    _pre[F] = out
    return out

def matching(D, Z):
    """max matching of rows D x cols D (lists of directions) avoiding zero cells Z (set of (x,y))."""
    D = list(D); ix = {z: i for i, z in enumerate(D)}
    r, c = [], []
    for x in D:
        for y in D:
            if (x, y) not in Z: r.append(ix[x]); c.append(ix[y])
    if not r: return 0
    M = csr_matrix((np.ones(len(r)), (r, c)), shape=(len(D), len(D)))
    return int((maximum_bipartite_matching(M, perm_type='column') >= 0).sum())

def max_clique(V, adj):
    best = [[]]
    def bk(R, P, X):
        if not P and not X:
            if len(R) > len(best[0]): best[0] = R[:]
            return
        if len(R) + len(P) <= len(best[0]): return
        u = max(P | X, key=lambda q: len(adj[q] & P))
        for v in list(P - adj[u]):
            bk(R + [v], P & adj[v], X & adj[v]); P = P - {v}; X = X | {v}
    bk([], set(V), set())
    return best[0]

def run(ks, maxd=16, rounds=3, allzero=False, numax=0, trueE=None, log=None):
    W = len(ks)
    w = CLSPAR[:, ks].sum(axis=1)                       # odd arrays per class
    r = w[VID]; s = WT - r
    if (s[1:] < 0).any(): return ("D", int(np.nonzero(s[1:] < 0)[0][0]) + 1)
    if allzero: tests = [F for F in range(1, NF) if r[F] <= numax and WT[F] <= maxd]
    else: tests = [F for F in range(1, NF) if VID[F] == VID[0] and WT[F] <= maxd and r[F] <= numax]
    Z = set()          # known zero cells (x, y) of E
    src = {}
    for rnd in range(1, rounds + 1):
        newz = False
        for F in tests:
            nu = int(r[F]); P = pre(F)
            u = []
            for dd in range(5):
                mw, sz, S = P[dd]
                b = np.minimum(sz, (mw - w[None, :]).min(axis=1))
                u.append(b.astype(np.int64))
            if rnd > 1:
                for dd in range(5):
                    S = P[dd][2]
                    for i in np.nonzero(u[dd] > 0)[0]:
                        Sd = [z for z in range(1, 32) if int(S[i]) >> (z - 1) & 1]
                        if any((x, y) in Z for x in Sd for y in Sd):
                            mm = matching(Sd, Z)
                            if mm < u[dd][i]: u[dd][i] = mm
            for dd in range(3, -1, -1):
                g = (1 << (5 - dd)) - 1; sz = P[dd][1]
                Bu = u[dd + 1][UP[dd]]            # nfl x g
                tot = Bu.sum(axis=1); val = tot - 2 * sz + 2 * nu
                bad = np.nonzero((val < 0) & (sz > 0))[0]
                if len(bad): return ("Z", rnd, F, dd, int(bad[0]))
                nb = np.minimum(u[dd], np.minimum(val // (g - 2), Bu.min(axis=1)))
                u[dd] = np.where(sz > 0, nb, u[dd])
            for dd in range(5):
                S = P[dd][2]; sz = P[dd][1]
                for i in np.nonzero((u[dd] == 0) & (sz > 0))[0]:
                    Sd = [z for z in range(1, 32) if int(S[i]) >> (z - 1) & 1]
                    for x in Sd:
                        for y in Sd:
                            if (x, y) not in Z:
                                Z.add((x, y)); src[(x, y)] = (F, dd, int(S[i]), rnd); newz = True
            if trueE is not None:
                for dd in range(5):
                    S = P[dd][2]
                    for i in range(len(S)):
                        Sd = [z - 1 for z in range(1, 32) if int(S[i]) >> (z - 1) & 1]
                        if Sd:
                            tr_ = np.linalg.matrix_rank(trueE[np.ix_(Sd, Sd)])
                            if tr_ > u[dd][i]: return ("UNSOUND", F, dd, i, int(u[dd][i]), int(tr_))
        if trueE is not None:
            for (x, y) in Z:
                if abs(trueE[x - 1, y - 1]) > 1e-12: return ("UNSOUND-ZERO", x, y, src[(x, y)])
        V = [x for x in range(1, 32) if (x, x) in Z]
        adj = {x: {y for y in V if y != x and (x, y) in Z and (y, x) in Z} for x in V}
        Y = max_clique(V, adj)
        if len(Y) > W: return ("C", rnd, len(Y), sorted(Y))
        allm = list(range(1, 32))
        if matching(allm, Z) < 31 - W: return ("A", rnd)
        for Q in range(1, NF):
            Dq = [z for z in range(1, 32) if int(MASK[Q]) >> (z - 1) & 1]
            if not any((x, y) in Z for x in Dq for y in Dq): continue
            mm = matching(Dq, Z)
            if mm < s[Q]: return ("M", rnd, Q, mm, int(s[Q]))
        if not newz: break
    return ("P", len(Z), len(Y))

if __name__ == "__main__":
    args = sys.argv[1:]
    maxd = int(args[args.index("--maxd") + 1]) if "--maxd" in args else 16
    rounds = int(args[args.index("--rounds") + 1]) if "--rounds" in args else 3
    allzero = "--allzero" in args
    for line in sys.stdin:
        line = line.strip()
        if not line or line.startswith("#"): continue
        ks = [int(x) for x in line.split(",")]
        t0 = time.time()
        res = run(ks, maxd=maxd, rounds=rounds, allzero=allzero)
        print(line, res, "%.1fs" % (time.time() - t0), flush=True)
