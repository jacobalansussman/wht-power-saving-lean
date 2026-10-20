#!/usr/bin/env python3
"""ctr-search / sas.py: subspace-assignment search for Problem D (O1, one-sided charge).
An element is (direction l in R^8, label subspace W_l of F_2^9); port b may use l iff u_b in W_l; port b is served iff
r_b in span{l : u_b in W_l}.  Cost = sum_l dim W_l (one direction carries one subspace; adding (l, W) to a present l
replaces W_l by W_l + W, which is Lemma D's merge).  Every served assignment gives an exact dictionary decoder
(dsearch.verify rebuilds it with rational coefficients and checks it on raw matrices).
Search: randomised greedy construction (best fixes per unit cost) + large-neighbourhood destroy/repair."""
import sys, os, json, time, random, itertools, argparse, pickle, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *
from dsearch import pool, verify, FLOOR_ONE_SIDED

R8 = np.array(LINES, dtype=float)

def span_list(basis):
    els = [0]
    for x in basis: els += [e ^ x for e in els]
    return els

def reduce_basis(vs):
    basis = []
    for v in vs:
        for x in basis: v = min(v, v ^ x)
        if v: basis.append(v)
    return basis

LABPOS = {u: b for b, u in enumerate(LABS)}

def port_mask(basis):
    m = 0
    for e in span_list(basis):
        b = LABPOS.get(e)
        if b is not None: m |= 1 << b
    return m

def candidate_subspaces(minports7=30, minports6=20, include_small=True):
    """hyperplanes, rich 7- and 6-dim subspaces, tetrad spans, singletons: list of (basis, dim, mask)"""
    cands = []
    # hyperplanes c^perp: basis of c^perp
    def perp_basis(cs):
        # basis of {x : x.c = 0 for c in cs}
        sol = [x for x in range(1, 512) if all(dot(x, c) == 0 for c in cs)]
        return reduce_basis(sol)
    for c in range(1, 512):
        bs = perp_basis([c]); m = port_mask(bs)
        if m: cands.append((bs, len(bs), m))
    seen = set()
    for c1 in range(1, 512):
        for c2 in range(c1 + 1, 512):
            key = tuple(sorted({c1, c2, c1 ^ c2}))
            if key in seen: continue
            seen.add(key)
            # quick count
            cnt = sum(1 for u in LABS if dot(u, c1) == 0 and dot(u, c2) == 0)
            if cnt >= minports7:
                bs = perp_basis([c1, c2]); cands.append((bs, 7, port_mask(bs)))
    if minports6 <= 20:
        # 6-dim with exactly 20 ports: label spans of 20-port sets found via 3 constraints (sampled exhaustively)
        seen3 = set()
        for c1 in range(1, 512):
            for c2 in range(c1 + 1, 512):
                for c3 in range(c2 + 1, 512):
                    pass
    full = list(range(9)); cands.append(([1 << i for i in range(9)], 9, (1 << N) - 1))
    if include_small:
        for t in range(14):
            bs = reduce_basis([LABS[b] for b in range(N) if TETOF[b] == t]); cands.append((bs, len(bs), port_mask(bs)))
        for b in range(N):
            cands.append(([LABS[b]], 1, 1 << b))
    return cands

class SAS:
    def __init__(s, P, cands, seed=0):
        s.P = [np.array(p, dtype=float) for p in P]; s.Pi = P; s.p = len(P)
        s.Pm = np.array(P, dtype=float)
        s.cands = cands
        s.by_port = [[i for i, (bs, d, m) in enumerate(cands) if m >> b & 1] for b in range(N)]
        s.rng = random.Random(seed)
        s.W = {}                    # direction index -> (basis, dim, mask)
        s.avail = [set() for _ in range(N)]
        s.proj = np.repeat(np.eye(8)[None], N, axis=0)       # projector onto span(avail)^perp
        s.res = R8.copy()                                    # residual of r_b
    def _update_port(s, b):
        A = sorted(s.avail[b])
        if A:
            M = s.Pm[A].T
            q, r = np.linalg.qr(M)
            rk = int(np.sum(np.abs(np.diag(r)) > 1e-9)) if r.size else 0
            # robust orthonormal basis of span
            u, sv, vt = np.linalg.svd(M, full_matrices=False)
            Q = u[:, sv > 1e-9]
            s.proj[b] = np.eye(8) - Q @ Q.T
        else:
            s.proj[b] = np.eye(8)
        s.res[b] = s.proj[b] @ R8[b]
    def served(s, b): return np.linalg.norm(s.res[b]) < 1e-9
    def unserved(s): return [b for b in range(N) if not s.served(b)]
    def cost(s): return sum(w[1] for w in s.W.values())
    def set_W(s, l, basis):
        old = s.W.get(l)
        oldm = old[2] if old else 0
        if basis:
            basis = reduce_basis(basis); m = port_mask(basis)
            s.W[l] = (basis, len(basis), m)
        else:
            m = 0; s.W.pop(l, None)
        ch = oldm ^ m
        while ch:
            b = (ch & -ch).bit_length() - 1; ch &= ch - 1
            if m >> b & 1: s.avail[b].add(l)
            else: s.avail[b].discard(l)
            s._update_port(b)
    def helped(s, l):
        """mask of unserved ports that become served if l were available to them"""
        rl = np.einsum('bij,j->bi', s.proj, s.P[l])
        num = np.abs(np.einsum('bi,bi->b', rl, s.res))
        den = np.linalg.norm(rl, axis=1) * np.linalg.norm(s.res, axis=1)
        ok = (den > 1e-9) & (np.abs(num - den) < 1e-7 * np.maximum(den, 1))
        m = 0
        for b in np.nonzero(ok)[0].tolist(): m |= 1 << b
        return m
    def progress(s, l):
        """mask of unserved ports whose span grows if l were available"""
        rl = np.einsum('bij,j->bi', s.proj, s.P[l])
        ok = (np.linalg.norm(rl, axis=1) > 1e-9) & (np.linalg.norm(s.res, axis=1) > 1e-9)
        m = 0
        for b in np.nonzero(ok)[0].tolist(): m |= 1 << b
        return m
    def greedy_step(s, b, topk=12):
        """add one element that serves (or advances) unserved port b, best ratio over candidates"""
        unm = 0
        for c in s.unserved(): unm |= 1 << c
        # candidate directions: those that help b, else those that advance b
        H = [(l, s.helped(l)) for l in range(s.p)]
        hl = [(l, h) for l, h in H if h >> b & 1]
        if not hl:
            hl = [(l, s.progress(l)) for l in range(s.p)]
            hl = [(l, h) for l, h in hl if h >> b & 1]
            partial = True
        else:
            partial = False
        best = None
        s.rng.shuffle(hl)
        for l, h in hl[:60]:
            cur = s.W.get(l)
            curb = cur[0] if cur else []; curd = cur[1] if cur else 0; curm = cur[2] if cur else 0
            for ci in s.by_port[b]:
                bs, d, m = s.cands[ci]
                if cur:
                    nb = reduce_basis(curb + bs); nd = len(nb)
                    if nd == curd: continue
                    nm = port_mask(nb) if nd != 9 else (1 << N) - 1
                else:
                    nb, nd, nm = bs, d, m
                gain = bin(h & nm & ~curm).count('1')
                inc = nd - curd
                score = (gain + (0.0 if not partial else 0.3)) / inc + s.rng.random() * 0.05
                if best is None or score > best[0]: best = (score, l, nb)
        if best is None: return False
        s.set_W(best[1], best[2])
        return True
    def greedy_global(s, noise=0.05, merge=True):
        """add the (direction, subspace) with the best drop of the residual potential sum_b |res_b|^2 per unit cost"""
        if not hasattr(s, 'CB'):
            s.CB = np.array([[(m >> b) & 1 for b in range(N)] for (bs, d, m) in s.cands], dtype=float)
            s.CD = np.array([d for (bs, d, m) in s.cands], dtype=float)
        best = None
        resn = np.einsum('bi,bi->b', s.res, s.res)
        for l in range(s.p):
            rl = np.einsum('bij,j->bi', s.proj, s.P[l])
            n2 = np.einsum('bi,bi->b', rl, rl)
            dotp = np.einsum('bi,bi->b', rl, s.res)
            red = np.where(n2 > 1e-12, dotp * dotp / np.maximum(n2, 1e-12), 0.0)
            if red.max() < 1e-9: continue
            cur = s.W.get(l)
            if cur is None:
                g = s.CB @ red
                sc = g / s.CD * (1 + noise * np.array([s.rng.random() for _ in range(1)]))
                i = int(np.argmax(sc))
                if best is None or sc[i] > best[0]: best = (sc[i], l, list(s.cands[i][0]))
            elif merge:
                curm = np.array([(cur[2] >> b) & 1 for b in range(N)], dtype=float)
                g = s.CB @ (red * (1 - curm))
                top = np.argsort(-g)[:8]
                for i in top.tolist():
                    if g[i] <= 1e-9: continue
                    nb = reduce_basis(cur[0] + s.cands[i][0]); inc = len(nb) - cur[1]
                    if inc <= 0: continue
                    sc = g[i] / inc
                    if best is None or sc > best[0]: best = (sc, l, nb)
        if best is None: return False
        s.set_W(best[1], best[2])
        return True
    def repair_global(s, maxsteps=300):
        for _ in range(maxsteps):
            if not s.unserved(): return True
            if not s.greedy_global(): return False
        return not s.unserved()
    def repair(s, maxsteps=400):
        for _ in range(maxsteps):
            un = s.unserved()
            if not un: return True
            b = s.rng.choice(un)
            if not s.greedy_step(b): return False
        return not s.unserved()
    def prune(s):
        """drop directions or shrink to a candidate subspace while every port stays served"""
        improved = True
        while improved:
            improved = False
            ls = list(s.W.keys()); s.rng.shuffle(ls)
            for l in ls:
                old = s.W[l][0]
                s.set_W(l, [])
                if not s.unserved(): improved = True; continue
                s.set_W(l, old)
        return s.cost()
    def snapshot(s): return {l: list(w[0]) for l, w in s.W.items()}
    def restore(s, snap):
        for l in list(s.W.keys()): s.set_W(l, [])
        for l, bs in snap.items(): s.set_W(l, bs)

def to_choice(S):
    """turn a served state into dsearch-style (P, reps, choice): each port's minimal subset of available directions"""
    reps, choice = [], []
    for b in range(N):
        A = sorted(S.avail[b])
        # greedy minimal subset spanning r_b
        F = list(A)
        for l in list(F):
            G = [x for x in F if x != l]
            if G:
                M = S.Pm[G].T
                u, sv, vt = np.linalg.svd(M, full_matrices=False); Q = u[:, sv > 1e-9]
                if np.linalg.norm(R8[b] - Q @ (Q.T @ R8[b])) < 1e-9: F = G
        reps.append([tuple(F)]); choice.append(0)
    return reps, choice

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--fams', default='FDR')
    ap.add_argument('--seeds', type=int, default=3)
    ap.add_argument('--lns', type=int, default=150)
    ap.add_argument('--destroy', type=int, default=3)
    ap.add_argument('--min7', type=int, default=30)
    ap.add_argument('--out', default=None)
    ap.add_argument('--calib', action='store_true', help='start from the coordinate decoder and check it is served')
    a = ap.parse_args()
    P = pool(a.fams)
    cache = os.path.join(HERE, f'out/cands9_{a.min7}.pkl')
    if os.path.exists(cache): cands = pickle.load(open(cache, 'rb'))
    else:
        cands = candidate_subspaces(minports7=a.min7); pickle.dump(cands, open(cache, 'wb'))
    print(f'pool {a.fams} {len(P)} directions; {len(cands)} candidate subspaces', flush=True)
    results = []
    t0 = time.time()
    for sd in range(a.seeds):
        S = SAS(P, cands, seed=sd)
        if a.calib:
            for i in range(8):
                l = P.index(tuple(1 if j == i else 0 for j in range(8)))
                S.set_W(l, reduce_basis([LABS[b] for b in range(N) if LINES[b][i] != 0]))
            print('calibration: coordinate decoder served:', not S.unserved(), 'cost', S.cost(), flush=True)
        else:
            ok = S.repair_global(); S.prune()
            print(f'seed {sd}: greedy served {ok} cost {S.cost()} ({time.time()-t0:.0f} s)', flush=True)
        best, bsnap = S.cost(), S.snapshot()
        for it in range(a.lns):
            ls = list(S.W.keys())
            for l in S.rng.sample(ls, min(a.destroy, len(ls))): S.set_W(l, [])
            ok = S.repair();
            if ok: S.prune()
            c = S.cost()
            if ok and c <= best:
                if c < best: print(f'   seed {sd} lns {it}: {c} ({time.time()-t0:.0f} s)', flush=True)
                best, bsnap = c, S.snapshot()
            else:
                S.restore(bsnap)
        S.restore(bsnap)
        reps, choice = to_choice(S)
        v = verify(P, reps, choice, f'sas {a.fams} seed {sd}')
        print(f'seed {sd}: best {best}; verified exact {v["exact"]} n {v["n"]} pair {v["pair"]} star {v["star_gather"]} '
              f'{v.get("REFUSED","")} terms {sorted(v["star_terms"], reverse=True)} ({time.time()-t0:.0f} s)', flush=True)
        results.append(v)
    if a.out: json.dump(results, open(a.out, 'w'), indent=1, default=str)
