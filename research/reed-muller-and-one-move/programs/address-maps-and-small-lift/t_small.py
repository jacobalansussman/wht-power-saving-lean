"""t_small.py (rm-theory): exhaustive exact pruning of SMALL multi-target slice programs (the small lift of Theorem B).

For W arrays at block m with M moves: every M-multiset of directions up to GL(m,2); every W-multiset of targets
congruent to I (invertible, non-alternating symmetric). Exact necessary tests, in order:
  (S)  every target lies in the F_2-span of the moves' rank-ones z z^T;
  (D)  deletion: |D_Q| >= r_Q for every form Q   (|D_Q| = moves with Q(z) = 1, with multiplicity; r_Q = odd arrays);
  (C)  tight-pair cliques: tight Q (|D_Q| = r_Q >= 1) have K[D,D] = I + Gamma_D, so a set of live moves that are
       pairwise in a common tight support has size <= W  (ex2-families Theorem 1; the multi-target proof is the same,
       see NOTE section 8).
A survivor of all three is printed (it is NOT a program: it goes to the numerical search and then to exact checks).
By Theorem B(a), an exact program for a survivor is a unit W in M of the no-phase game with any address group H
containing maps h_l with h_l h_l^T = T_l.

Usage: python3 t_small.py m W M [maxprint]
       python3 t_small.py m W list:z1,z2,...      (one multiset only; for calibration)
"""
import sys, itertools, time
import numpy as np, os
SKIP_INDEP = os.environ.get('SKIP_INDEP', '0') == '1'
TERMRANK = os.environ.get('TERMRANK', '0') == '1'
DUMP = os.environ.get('DUMP', ''); DUMPCAP = int(os.environ.get('DUMPCAP', '40'))
from t_free import forms, rank2, gl_elements


def targets_of(m, pairs):
    P = len(pairs); out = []
    for bits in range(1 << P):
        S = [[0] * m for _ in range(m)]
        for idx, (i, j) in enumerate(pairs):
            if (bits >> idx) & 1: S[i][j] = S[j][i] = 1
        if not any(S[i][i] for i in range(m)): continue
        if rank2([sum(S[i][j] << j for j in range(m)) for i in range(m)]) != m: continue
        out.append(np.array([S[i][j] for (i, j) in pairs], dtype=np.int64))
    return np.stack(out)


def max_clique(adj, verts):
    best = 0
    verts = list(verts)
    def rec(R, Pset):
        nonlocal best
        if not Pset:
            best = max(best, len(R)); return
        if len(R) + len(Pset) <= best: return
        for v in list(Pset):
            rec(R + [v], [u for u in Pset if u in adj[v] and u > v])
    rec([], verts)
    return best


def clique_ok(ms, vals, cnt, r, W):
    tight = np.where((cnt == r) & (r >= 1))[0]
    if len(tight) == 0: return True, 0, 0
    M = len(ms)
    adj = {j: set() for j in range(M)}
    live = set()
    Zv = vals[:, list(ms)]                       # forms x M : Q(z_j)
    for q in tight:
        Dm = [j for j in range(M) if Zv[q, j]]
        live |= set(Dm)
        for a in Dm:
            for b in Dm:
                if a != b: adj[a].add(b)
    c = max_clique(adj, sorted(live))
    return c <= W, len(tight), c


def clique_bits(adjm, verts):
    """max clique size; adjm: list of int bitmasks; verts: bitmask of allowed vertices"""
    best = 0
    def rec(R, P):
        nonlocal best
        if P == 0:
            best = max(best, R); return
        if R + bin(P).count("1") <= best: return
        while P:
            v = (P & -P).bit_length() - 1
            P &= P - 1
            rec(R + 1, P & adjm[v])
    rec(0, verts)
    return best


def leaf_tests(Zv, cnt, r, W):
    """(C) tight-pair cliques <= W and (Z) zero-block test: for every form F and every pairwise-joined set Y of
    live moves inside D_F (so E[Y,Y] = 0), |D_F| + r_F >= 2|Y|  (rank E[D_F,D_F] = |D_F| - r_F <= 2(|D_F| - |Y|)).
    Returns (pass, ntight, clique, zero_block_violation_or_None)."""
    tight = (cnt == r) & (r >= 1)
    nt = int(tight.sum())
    if nt == 0:
        return (True, 0, 0, None)
    T = Zv[tight].astype(np.int32)
    A = (T.T @ T) > 0
    M = Zv.shape[1]
    live = [bool(A[j, j]) for j in range(M)]
    adjm = [sum(1 << u for u in range(M) if u != v and A[v, u]) for v in range(M)]
    livemask = sum(1 << v for v in range(M) if live[v])
    c = clique_bits(adjm, livemask)
    if c > W:
        return (False, nt, c, None)
    # zero-block test: only forms with cnt + r < 2c can fail
    cand = np.where((cnt + r < 2 * c) & (cnt > 0))[0]
    for f in cand:
        dm = 0
        for j in np.nonzero(Zv[f])[0]:
            if live[j]: dm |= 1 << int(j)
        w = clique_bits(adjm, dm)
        if 2 * w > int(cnt[f]) + int(r[f]):
            return (False, nt, c, (int(f), int(cnt[f]), int(r[f]), w))
    if TERMRANK:
        # term-rank test: E[a,b] = 0 whenever a, b lie in a common tight support (a = b allowed: live moves), so
        # rank E[D_F, D_F] <= maximum matching avoiding those cells; it must be >= |D_F| - r_F.
        for f in np.where(cnt > 0)[0]:
            Dm = [int(j) for j in np.nonzero(Zv[f])[0]]
            if not any(live[j] for j in Dm): continue
            need = int(cnt[f]) - int(r[f])
            allowed = {a: [b for b in Dm if not (A[a, b])] for a in Dm}
            match = {}
            def aug(a, seen):
                for b in allowed[a]:
                    if b in seen: continue
                    seen.add(b)
                    if b not in match or aug(match[b], seen):
                        match[b] = a; return True
                return False
            tr = 0
            for a in Dm:
                if aug(a, set()): tr += 1
            if tr < need:
                return (False, nt, c, (int(f), int(cnt[f]), int(r[f]), -tr))
    return (True, nt, c, None)


def run(m, W, Mspec, maxprint=30):
    t0 = time.time()
    pairs, mon, coeffs, vals = forms(m)
    nn = 1 << m
    Tm = targets_of(m, pairs)
    par = (coeffs @ Tm.T % 2).T.astype(np.int16)          # t x forms
    nt = len(Tm)
    if Mspec.startswith("list:"):
        classes = [(tuple(sorted(int(x) for x in Mspec[5:].split(","))), 1)]
        M = len(classes[0][0])
    else:
        M = int(Mspec)
        G = gl_elements(m)
        A = np.zeros((len(G), nn), dtype=np.int64)
        for gi, h in enumerate(G):
            for z in range(nn):
                r = 0; i = 0; zz = z
                while zz:
                    if zz & 1: r ^= h[i]
                    zz >>= 1; i += 1
                A[gi, z] = r
        pw = (nn ** np.arange(M)[::-1]).astype(np.int64)
        seen = set(); classes = []
        for ms in itertools.combinations_with_replacement(range(1, nn), M):
            key = 0
            for z in ms: key = key * nn + z
            if key in seen: continue
            imgs = np.sort(A[:, list(ms)], axis=1)
            keys = set((imgs @ pw).tolist())
            seen |= keys; classes.append((ms, len(keys)))
        del seen
    print(f"m={m} W={W} M={M}: {nt} targets; {len(classes)} classes of direction multisets ({time.time()-t0:.1f} s)",
          flush=True)
    nS = nD = nC = 0; printed = 0; surv_classes = 0; nI = [0]; nZ = [0]; exC = [0]; exZ = [0]
    for ms, osz in classes:
        if SKIP_INDEP:
            dist = sorted(set(ms))
            ro = [int("".join(str(int(x)) for x in mon[z][::-1]), 2) for z in dist]
            if rank2(ro) == len(dist) and M < m * W:
                nI[0] += 1; continue          # free-direction lemma (proved, NOTE 8.5): deletion excludes these
        cnt = np.zeros(vals.shape[0], dtype=np.int16)
        for z in ms: cnt += vals[:, z].astype(np.int16)
        Z = cnt == 0
        inV = np.where(~(par[:, Z].astype(bool)).any(axis=1))[0]
        if len(inV) == 0: continue
        nS += 1
        sub = par[inV].astype(np.int16)                    # k x forms
        k = len(inV)
        Zv = vals[:, list(ms)].astype(bool)                # forms x M
        cls_surv = []; anyD = [False]; nleaf = [0]
        # depth-first over W-multisets of targets (indices nondecreasing), deletion-pruned at every level
        def dfs(tup, s, depth):
            a0 = tup[-1] if tup else 0
            cand = (s[None, :] if s is not None else 0) + sub[a0:]
            ok = np.where((cnt[None, :] >= cand).all(axis=1))[0]
            for o in ok:
                t2 = tup + (a0 + int(o),)
                if depth + 1 == W:
                    anyD[0] = True; nleaf[0] += 1
                    res = leaf_tests(Zv, cnt, cand[o], W)
                    if res[0]: cls_surv.append((t2,) + res[1:])
                    elif res[3] is None: exC[0] += 1
                    else:
                        exZ[0] += 1
                        if exZ[0] <= 3: print(f"    zero-block exclusion: dirs {ms} targets {[int(inV[a]) for a in t2]}: form {res[3][0]} |D| = {res[3][1]}, r = {res[3][2]}, zero clique {res[3][3]}", flush=True)
                else:
                    dfs(t2, cand[o], depth + 1)
        dfs((), None, 0)
        if anyD[0]: nD += 1
        nZ[0] += sum(1 for c in cls_surv if True)
        if cls_surv:
            nC += 1; surv_classes += 1
            if DUMP:
                import random as _r
                pick = cls_surv if len(cls_surv) <= DUMPCAP else _r.Random(len(cls_surv)).sample(cls_surv, DUMPCAP)
                with open(DUMP, "a") as fh:
                    for c in pick:
                        fh.write(",".join(map(str, ms)) + "|" + ",".join(str(int(inV[a])) for a in c[0]) + "\n")
            dist = sorted(set(ms))
            ro = [int("".join(str(int(x)) for x in mon[z][::-1]), 2) for z in dist]
            indep = rank2(ro) == len(dist)
            if printed < maxprint:
                printed += 1
                ex = cls_surv[0]
                tl = ",".join(str(int(inV[a])) for a in ex[0])
                ntf = sum(1 for c in cls_surv if c[1] == 0)
                print(f"  SURVIVOR dirs {ms} (orbit {osz}; distinct dirs independent: {indep}): {len(cls_surv)} target "
                      f"{W}-multisets pass S+D+C+Z ({ntf} with no tight form); example targets idx {tl} "
                      f"(tight {ex[1]}, clique {ex[2]})", flush=True)
                if Mspec.startswith("list:"):
                    for tup, ntt, c, _ in cls_surv[:10]:
                        print("     targets", ",".join(str(int(inV[a])) for a in tup), "tight", ntt, "clique", c)
    print(f"  classes with some target in span: {nS}; with a target {W}-multiset passing deletion: {nD}; "
          f"passing cliques and the zero-block test too: {nC}")
    print(f"  target tuples passing deletion but excluded: by cliques {exC[0]}, by the zero-block test {exZ[0]}")
    print(f"  classes skipped by the free-direction lemma (independent distinct directions, M < mW): {nI[0]}")
    print(f"  done in {time.time()-t0:.1f} s")


if __name__ == "__main__":
    m, W = int(sys.argv[1]), int(sys.argv[2])
    run(m, W, sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 30)
