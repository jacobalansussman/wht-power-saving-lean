"""vtests.py: independent exact tests for a multi-target commutative program shape (multiset of directions, list of
targets): S (span), D (deletion), C (tight-pair clique <= W), Z (zero-block), TR (term rank of E[D_F,D_F] with the
forced zeros). Forms: all 2^P binary quadratic forms (diagonal = linear terms)."""
import itertools
import numpy as np
from vfield import Model
class Tester:
    def __init__(s, m):
        s.m = m; s.Mo = Model(m); s.N = (1 << m) - 1
        P = len(s.Mo.mons); s.P = P
        coeffs = np.array(list(itertools.product([0, 1], repeat=P)), dtype=np.int64)
        zb = np.array([[(z >> j) & 1 for j in range(m)] for z in range(1, s.N + 1)])
        mon = np.stack([zb[:, i] * zb[:, j] for (i, j) in s.Mo.mons], axis=1)
        s.E = ((coeffs @ mon.T) % 2).astype(np.int8)          # forms x directions (column z-1)
        s.coeffs = coeffs
        s.parcache = {}
    def parity(s, k):
        if k not in s.parcache:
            T = s.Mo.target(k); tv = np.array([T[i][j] for (i, j) in s.Mo.mons])
            s.parcache[k] = ((s.coeffs @ tv) % 2).astype(np.int8)
        return s.parcache[k]
    def parT(s, T):
        tv = np.array([T[i][j] for (i, j) in s.Mo.mons]); return ((s.coeffs @ tv) % 2).astype(np.int8)
    def run(s, dirs, ks, want_detail=False, Ts=None):
        """dirs: list of direction ints (moves); ks: target labels (or Ts: explicit target matrices in the SAME
        coordinates as dirs). Returns dict of verdicts."""
        W = len(ks) if Ts is None else len(Ts); M = len(dirs)
        Emv = s.E[:, [z - 1 for z in dirs]].astype(np.int64)    # forms x moves
        D = Emv.sum(1); r = (sum(s.parity(k).astype(np.int64) for k in ks) if Ts is None
             else sum(s.parT(T).astype(np.int64) for T in Ts))
        res = {}
        res["S"] = not ((D == 0) & (r > 0)).any()
        res["D"] = not (D < r).any()
        if not (res["S"] and res["D"]): return res
        tight = (D == r) & (r > 0)
        Tm = Emv[tight]                                          # tight supports
        live = Tm.any(0)
        adj = (Tm.T @ Tm) > 0                                    # moves sharing a tight support (diag = live)
        # max clique (brute force over subsets of live moves, M small)
        L = [i for i in range(M) if live[i]]
        def maxclique(cand):
            best = []
            def rec(cl, rest):
                nonlocal best
                if len(cl) + len(rest) <= len(best): return
                if not rest:
                    best = list(cl); return
                v = rest[0]
                rec(cl + [v], [u for u in rest[1:] if adj[v, u]])
                rec(cl, rest[1:])
            rec([], list(cand)); return best
        cq = maxclique(L)
        res["clique"] = len(cq); res["C"] = len(cq) <= W
        res["ntight"] = int(tight.sum())
        # zero-block and term rank over every form
        zfail = 0; trfail = 0
        Zforced = adj.copy()                                     # forced-zero entries of E
        seenD = {}
        for f in range(Emv.shape[0]):
            if D[f] == 0: continue
            Dset = tuple(np.nonzero(Emv[f])[0])
            key = (Dset, int(r[f]))
            if key in seenD: continue
            seenD[key] = 1
            Y = maxclique([i for i in Dset if live[i]])
            if D[f] + r[f] < 2 * len(Y): zfail += 1
            # term rank: max matching in allowed entries
            allowed = {i: [j for j in Dset if not Zforced[i, j]] for i in Dset}
            match = {}
            def aug(i, vis):
                for j in allowed[i]:
                    if j in vis: continue
                    vis.add(j)
                    if j not in match or aug(match[j], vis):
                        match[j] = i; return True
                return False
            tr = sum(1 for i in Dset if aug(i, set()))
            if tr < D[f] - r[f]: trfail += 1
        res["Zfail_forms"] = zfail; res["Z"] = zfail == 0
        res["TRfail_forms"] = trfail; res["TR"] = trfail == 0
        return res
