"""rc_ilp.py (rm2-exact): piece + term-rank bound on rank E[D_F, D_F], optimised exactly by an ILP (HiGHS).

Bound: E[D_F,D_F] = sum over chosen disjoint pieces R x C (each inside D_G x D_G for a form G, so rank <= s_G)
       + the leftover cells, whose rank is at most the minimum number of rows + columns covering them (term rank).
If the ILP optimum is < s_F = |D_F| - r_F (the exact rank, Theorem 10 per-array form), no program exists.
Candidate pieces from each other form G with A = D_F n D_G: the square A x A and every split R x (A \ R), kept
only when s_G is below the piece's trivial cost.  Violations are re-verified independently (rc_run.verify-style).

usage: rc_ilp.py <case> <which> [maxforms] [smax]
   which: even (forms with r_F = 0), all, or a comma list of form ids
"""
import sys, time, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from xbase import Field, Forms, System
from rc_run import build, second_parity

def present_mask(fo, sy):
    return sum(1 << (z - 1) for z in fo.dirs if sy.counts.get(z, 0) > 0)

def pieces_for(fo, sy, F, smax, maxsplit=10):
    nd = len(fo.dirs); pm = present_mask(fo, sy); DF = int(fo.mask[F]) & pm; c = sy.counts
    cw = np.array([c.get(z, 0) for z in range(1, nd + 1)], dtype=np.int64)
    A = (fo.mask & DF).astype(np.int64)
    bits = ((A[:, None] >> np.arange(nd)[None, :]) & 1)
    wA = bits @ cw
    ok = (sy.s < wA) & (sy.s <= smax) & (A != 0)
    ok[F] = False
    idx = np.nonzero(ok)[0]
    # dedupe by A keeping the smallest s
    order = idx[np.lexsort((sy.s[idx], A[idx]))]
    best = {}
    seenA = {}
    for G in order:
        a = int(A[G])
        if a in seenA: continue
        seenA[a] = (int(sy.s[G]), int(G))
    w = lambda msk: int(sum(cw[z - 1] for z in range(1, nd + 1) if msk >> (z - 1) & 1))
    for a, (sG, G) in seenA.items():
        if best.get((a, a), (10 ** 9,))[0] > sG: best[(a, a)] = (sG, G)
        el = [z for z in range(1, nd + 1) if a >> (z - 1) & 1]
        if len(el) > maxsplit: continue
        for k in range(1, len(el)):
            for Rl in itertools.combinations(el, k):
                R = sum(1 << (z - 1) for z in Rl); C = a & ~R
                if sG < min(w(R), w(C)):
                    if best.get((R, C), (10 ** 9,))[0] > sG: best[(R, C)] = (sG, G)
    return best

def solve(fo, sy, F, smax, tlimit=60):
    nd = len(fo.dirs); DF = int(fo.mask[F]) & present_mask(fo, sy); c = sy.counts
    D = [z for z in range(1, nd + 1) if DF >> (z - 1) & 1]; d = len(D); pos = {z: i for i, z in enumerate(D)}
    P = pieces_for(fo, sy, F, smax)
    plist = list(P.items())
    npc = len(plist); nv = npc + 2 * d
    cost = np.array([v[0] for (_, v) in plist] + [c[z] for z in D] * 2, dtype=float)
    ncell = d * d
    Acov = lil_matrix((ncell, nv)); Adis = lil_matrix((ncell, nv))
    for j, ((R, C), _) in enumerate(plist):
        for x in D:
            if R >> (x - 1) & 1:
                for y in D:
                    if C >> (y - 1) & 1:
                        cell = pos[x] * d + pos[y]; Acov[cell, j] = 1; Adis[cell, j] = 1
    for i in range(d):
        for k in range(d):
            Acov[i * d + k, npc + i] = 1; Acov[i * d + k, npc + d + k] = 1
    cons = [LinearConstraint(Acov.tocsr(), lb=np.ones(ncell), ub=np.inf)]
    if npc: cons.append(LinearConstraint(Adis.tocsr(), lb=-np.inf, ub=np.ones(ncell)))
    res = milp(cost, constraints=cons, integrality=np.ones(nv), bounds=Bounds(0, 1),
               options={"time_limit": tlimit, "disp": False})
    if res.x is None: return None, npc, None
    x = np.round(res.x).astype(int)
    chosen = [(plist[j][0][0], plist[j][0][1], plist[j][1][0], plist[j][1][1]) for j in range(npc) if x[j]]
    rows = [D[i] for i in range(d) if x[npc + i]]; cols = [D[k] for k in range(d) if x[npc + d + k]]
    val = int(round(res.fun))
    return val, npc, (chosen, rows, cols, res.status)

def verify(Fd, fo, sy, targets, F, chosen, rows, cols):
    """independent check: pieces disjoint, inside D_F, inside D_G, bound = s_G recomputed by field evaluation;
    leftover cells covered by the chosen rows/columns; total < s_F."""
    nd = len(fo.dirs); cnt = sy.counts
    def Qv(q, z):
        v = 0
        for p, (j, k) in enumerate(fo.PAIRS):
            if (q >> p) & 1: v ^= ((z >> j) & 1) & ((z >> k) & 1)
        return v
    supp = lambda q: {z for z in range(1, nd + 1) if cnt.get(z, 0) and Qv(q, z)}
    slack = lambda q: sum(cnt[z] for z in supp(q)) - sum(second_parity(Fd, fo, q, T) for T in targets)
    DF = supp(F); used = set(); tot = 0
    for (R, C, sG, G) in chosen:
        Rs = {z for z in range(1, nd + 1) if R >> (z - 1) & 1}; Cs = {z for z in range(1, nd + 1) if C >> (z - 1) & 1}
        assert Rs <= DF and Cs <= DF and (Rs | Cs) <= supp(G)
        assert slack(G) == sG
        for cell in itertools.product(Rs, Cs):
            assert cell not in used; used.add(cell)
        tot += sG
    for cell in itertools.product(DF, DF):
        if cell not in used: assert cell[0] in rows or cell[1] in cols, "uncovered leftover cell"
    tot += sum(cnt[z] for z in rows) + sum(cnt[z] for z in cols)
    sF = slack(F)
    assert tot < sF
    return tot, sF

def main():
    case = sys.argv[1]; which = sys.argv[2]
    maxforms = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
    smax = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9
    m, targets, counts, expected = build(case)
    Fd = Field(m); fo = Forms(m); sy = System(fo, targets, counts)
    if len(sy.violations()):
        print("CASE %s: excluded by deletion" % case); return
    present = sum(1 << i for i, z in enumerate(fo.dirs) if sy.counts.get(z, 0) > 0)
    if which == "even": cand = [F for F in range(1, fo.nf) if sy.r[F] == 0 and fo.mask[F] & present]
    elif which == "all": cand = [F for F in range(1, fo.nf) if sy.s[F] > 0 and fo.mask[F] & present]
    else: cand = [int(x) for x in which.split(",")]
    cand.sort(key=lambda F: (-(sy.s[F] / max(1, sy.d[F])), sy.d[F]))
    cand = cand[:maxforms]
    t0 = time.time(); nviol = 0; hist = {}; first = None
    for F in cand:
        val, npc, info = solve(fo, sy, F, smax)
        if val is None: hist["fail"] = hist.get("fail", 0) + 1; continue
        mg = val - int(sy.s[F]); hist[mg] = hist.get(mg, 0) + 1
        if val < sy.s[F]:
            chosen, rows, cols, st = info
            tot, sF = verify(Fd, fo, sy, targets, F, chosen, rows, cols); nviol += 1
            if first is None: first = (F, tot, sF, len(chosen), len(rows) + len(cols))
    print("CASE %s expected=%s W=%d M=%d: ILP on %d forms (%s) in %.1fs; VERIFIED violations %d%s; margin hist %s"
          % (case, expected, sy.W, sum(sy.counts.values()), len(cand), which, time.time() - t0, nviol,
             (" first F=%d bound %d < s_F %d (%d pieces + %d lines)" % first) if first else "", sorted(hist.items(), key=str)[:8]))

if __name__ == "__main__":
    main()
