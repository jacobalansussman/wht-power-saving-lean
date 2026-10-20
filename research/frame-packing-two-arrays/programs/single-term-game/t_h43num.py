#!/usr/bin/env python3
"""t_h43num.py: two-array units (no scratch, C[P] mixes; toylab Solver) in the HAMMING sub-game H(k,3) of the block-k
engine frame game: G = S_3^k (= Sp(2,2)^k inside Sp(2k,2)), P = Z_2^k, paid letters = one-coordinate moves, target at
level k.  Irreps exactly, as tensor products of S_3 irreps.  Letter sequences: every word in the k coordinates using
each coordinate, up to relabelling of coordinates and reversal.  Usage: t_h43num.py k M part nparts restarts"""
import sys, time, itertools, warnings
warnings.filterwarnings("ignore")
import numpy as np
from toylab import perm_group, cyc, FrameGame, Solver, search_lm2

def s3_irreps():
    B = np.array([[1, -1, 0], [1, 1, -2]], dtype=float).T
    B = B / np.linalg.norm(B, axis=0)
    def triv(p): return np.eye(1)
    def sign(p):
        inv = sum(1 for i in range(3) for j in range(i + 1, 3) if p[i] > p[j]); return np.eye(1) * (-1) ** inv
    def std(p):
        Pm = np.zeros((3, 3))
        for i in range(3): Pm[p[i], i] = 1
        return B.T @ Pm @ B
    return [triv, sign, std]

def build(k):
    n = 3 * k
    gens = []
    for i in range(k): gens += [cyc(n, [3 * i, 3 * i + 1]), cyc(n, [3 * i, 3 * i + 1, 3 * i + 2])]
    G = perm_group(gens, "H(%d,3)" % k)
    base = {3 * i for i in range(k)}
    P = [e for e, p in enumerate(G.elems) if {p[b] for b in base} == base]
    hreps = [G.idx[cyc(n, [3 * i, 3 * i + 1])] for i in range(k)]
    fg = FrameGame(G, P, hreps)
    ir = s3_irreps()
    reps = []
    for ks in itertools.product(range(3), repeat=k):
        mats = []
        for e in G.elems:
            M = np.eye(1)
            for i in range(k):
                comp = tuple(e[3 * i + j] - 3 * i for j in range(3))
                M = np.kron(M, ir[ks[i]](comp))
            mats.append(M)
        reps.append(np.array(mats, dtype=complex))
    # homomorphism check on random pairs
    rng = np.random.default_rng(0)
    for _ in range(20):
        a, b = rng.integers(G.n, size=2)
        for r in reps[:5] + reps[-3:]:
            assert np.allclose(r[G.m(a, b)], r[a] @ r[b])
    assert sum(r.shape[1] ** 2 for r in reps) == G.n
    return G, P, hreps, fg, reps

def words(k, M):
    """words in range(k) of length M using every coordinate, canonical under relabelling (first occurrences in order)
    and reversal"""
    out = []
    for w in itertools.product(range(k), repeat=M):
        if len(set(w)) < k: continue
        def canon(v):
            mp = {}; r = []
            for x in v:
                if x not in mp: mp[x] = len(mp)
                r.append(mp[x])
            return tuple(r)
        c = canon(w)
        if c != w: continue
        if canon(w[::-1]) < w: continue
        out.append(w)
    return out

if __name__ == "__main__":
    k, M, part, nparts, rs = (int(x) for x in sys.argv[1:6])
    G, P, hreps, fg, reps = build(k)
    T = [g for g in range(G.n) if fg.ell[g] == k][0]
    S = Solver(fg, reps, W=2)
    W = words(k, M)
    print("H(%d,3): |G|=%d |P|=%d frames=%d m=%d irreps=%d; M=%d: %d words, part %d/%d" % (
        k, G.n, len(P), fg.ncos, fg.diam, len(reps), M, len(W), part, nparts), flush=True)
    best = (np.inf, None)
    for i, w in enumerate(W):
        if i % nparts != part: continue
        t = time.time()
        res, Fs, allr = search_lm2(S, [hreps[c] for c in w], T, restarts=rs, seed=i, lm_iters=150)
        print("  %s  best %.3e  (%s)  %.0fs" % ("".join(map(str, w)), res, " ".join("%s:%.2f" % (a, b) for a, b in allr), time.time() - t), flush=True)
        if res < best[0]: best = (res, w)
        if res < 1e-10:
            np.save("unit_%s.npy" % "".join(map(str, w)), np.array(Fs)); print("  FOUND (residual %.2e)" % res, flush=True)
    print("RESULT k=%d M=%d part %d: best %.3e at %s" % (k, M, part, best[0], best[1]), flush=True)
