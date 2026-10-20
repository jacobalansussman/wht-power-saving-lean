#!/usr/bin/env python3
"""vunits.py (v-p6-prove): rebuild every unit named in p6-prove's PROOF.md with my own engine (veng.py, true model),
words re-typed from the PROOF.md prose (not imported).  For each: formal identity in C[G] (true model), M,
raw complex arrays at several f, letter legitimacy (free letters fix L0; paid letters are dir z), and the hypotheses:
  single-term at every moment; one Sp label per entry at every moment; Theorem C hypothesis with h = level and with
  h = m - d(., L_T) (every MOVE OUTPUT at one height); relative Sp labels of two-input registers (holonomy);
  my own Bellman-Ford height filter (existence of ANY height system) and the best quantitative Theorem-C bound.
Usage: python3 vunits.py [fmax_m2 fmax_m3 fmax_m4]
"""
import sys
import numpy as np
from fractions import Fraction as Fr
from veng import Eng, run, raw_check, Q

def words(E):
    m = E.m
    out = []
    plain = lambda W: [("move", r, 1 << j) for j in range(m) for r in range(W)]
    out.append(("plain", 2, 2, plain(2)))
    # gauge exchange: move x0 by e1; c0 += c1; e2..em on both; c0 -= c1; e1 on c1
    w = [("move", 0, 1), ("mix", 0, [(0, 1, None), (1, 1, None)])]
    for j in range(1, m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
    w += [("mix", 0, [(0, 1, None), (1, -1, None)]), ("move", 1, 1)]
    out.append(("gaugeX", 2, 2, w))
    if m >= 3:
        s = E.addr([2, 4, 1] + [1 << j for j in range(3, m)])
        s2 = E.addr([4, 1, 2] + [1 << j for j in range(3, m)])
        h = Fr(1, 2)
        w = [("mix", 0, [(0, 1, None), (1, 1, None)]),                       # c0 = x0 + x1
             ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s)])]          # c1 = x0 - s x1 (repeated source c1)
        for j in range(m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        w += [("mix", 1, [(0, h, None), (1, -h, None), (0, -h, s), (1, h, s), (0, h, s2), (1, -h, s2)]),
              ("mix", 0, [(0, 1, None), (1, -1, None)])]
        out.append(("C%dO3t" % (2 * m), 2, 2, w))
    # SPLIT1 (W = 1; registers c=0, t=1, s=2)
    w = [("mix", 1, [(0, 1, None)]), ("mix", 2, [(0, 1, None)]), ("move", 2, 1),
         ("mix", 0, [(0, 1, None), (2, 1, None)]), ("move", 0, 1), ("mix", 0, [(0, 1, None), (1, -1, "X1")])]
    w += [("move", 0, 1 << j) for j in range(1, m)]
    out.append(("SPLIT1", 1, 3, w))
    # HOLO3 (W = 3, no scratch)
    if m >= 2:
        sw = E.addr([2, 1] + [1 << j for j in range(2, m)])
        w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("free", 1, sw), ("mix", 1, [(1, 1, None), (2, 1, None)]),
             ("mix", 2, [(2, 1, None), (0, 1, None)])] + plain(3)
        w += [("mix", 2, [(2, 1, None), (0, -1, None)]), ("mix", 1, [(1, 1, None), (2, -1, None)]), ("free", 1, sw),
              ("mix", 0, [(0, 1, None), (1, -1, None)])]
        out.append(("HOLO3", 3, 3, w))
    if m == 1:
        # t_scalar word: s := c; move s; c := c + s; t := copy of x (made first), then X Z X Z on t (= -I); c := c + t
        Z2 = ["P1_0", "P1_0"]       # (P1_0)^2 = diag(1, -1) = Z
        w = [("mix", 2, [(0, 1, None)]), ("mix", 1, [(0, 1, None)]), ("move", 1, 1), ("mix", 0, [(0, 1, None), (1, 1, None)])]
        w += [("free", 2, n) for n in ["X1"] + Z2 + ["X1"] + Z2]
        w += [("mix", 0, [(0, 1, None), (2, 1, None)])]
        out.append(("tscalar", 1, 3, w))
    return out

def analyse(E, name, W, R, word, fmax, rng):
    ok, M, hist = run(E, W, R, word)
    L0 = E.L0(); LT = E.pos(E.K)
    single = True; onelabel = True
    for C in hist:
        for c in C:
            per = {}
            for (l, g) in c: per.setdefault(l, []).append(g)
            for l, gs in per.items():
                if len(gs) > 1: single = False
                if len({tuple(E.symp(g)) for g in gs}) > 1: onelabel = False
    # Theorem C hypothesis
    hC_lv = hC_co = True; outs = []
    rel = {}
    for i, st in enumerate(word):
        C = hist[i + 1]
        if st[0] == "move":
            c = C[st[1]]
            lv = {E.level(g) for (l, g) in c}; co = {E.m - E.d(E.pos(g), LT) for (l, g) in c}
            if len(lv) > 1: hC_lv = False
            if len(co) > 1: hC_co = False
            outs.append({(l, E.pos(g)) for (l, g) in c})
        for r, c in enumerate(C):
            per = {}
            for (l, g) in c: per.setdefault(l, []).append(g)
            ls = sorted(per)
            for a in ls:
                for b in ls:
                    if a < b and len(per[a]) == 1 and len(per[b]) == 1:
                        # relative Sp label gbar_a^-1 gbar_b
                        ga, gb = per[a][0], per[b][0]
                        rk = tuple(E.symp(E.mulk(E.dag(ga), gb)))
                        rel.setdefault((a, b), set()).add(rk)
    holo = {k: len(v) for k, v in rel.items()}
    feas, bound = heightfilter(E, outs, W)
    raw = {}
    for f in range(1, fmax + 1):
        raw[f] = raw_check(E, W, R, word, f, rng)
    print("== %-7s m=%d W=%d R=%d  unit(true C[G])=%s  M=%d (Wm=%d)" % (name, E.m, W, R, ok, M, W * E.m))
    print("   raw arrays max err by f: %s" % {f: float("%.1e" % v) for f, v in raw.items()})
    print("   single-term always=%s  one Sp label per entry always=%s  distinct relative labels per input pair=%s"
          % (single, onelabel, holo))
    print("   ThmC hyp: level-homog=%s  colevel-homog=%s;  my height filter: system exists=%s, best bound M>=%d"
          % (hC_lv, hC_co, feas, bound))
    return ok

def heightfilter(E, outs, W):
    """independent Bellman-Ford: variables h_l(x) for x in positions used + L0, LT; constraints
    h(L0)=0, |h(x)-h(y)| <= d(x,y) per input, equal heights inside each move output.  Feasibility with h(LT)=m, and
    the max of sum h_l(LT) without that equality (pointwise max = shortest paths from the zero node)."""
    L0 = E.L0(); LT = E.pos(E.K)
    nodes = {("Z",)}
    for l in range(W): nodes |= {(l, L0), (l, LT)}
    for o in outs: nodes |= o
    nodes = list(nodes); ix = {v: i for i, v in enumerate(nodes)}
    base = []
    def eq(u, v, c): base.append((ix[u], ix[v], c)); base.append((ix[v], ix[u], -c))
    for l in range(W): eq(("Z",), (l, L0), 0)
    for a in nodes:
        for b in nodes:
            if a != b and a[0] != "Z" and b[0] != "Z" and a[0] == b[0]:
                base.append((ix[a], ix[b], E.d(a[1], b[1])))
    for o in outs:
        o = list(o)
        for a, b in zip(o, o[1:]): eq(a, b, 0)
    def sp(edges):
        INF = 10 ** 9; n = len(nodes); dist = [INF] * n; dist[ix[("Z",)]] = 0
        for _ in range(n + 1):
            ch = False
            for (u, v, w) in edges:
                if dist[u] < INF and dist[u] + w < dist[v]: dist[v] = dist[u] + w; ch = True
            if not ch: return dist
        return None
    dist = sp(base)
    bound = sum(dist[ix[(l, LT)]] for l in range(W))
    ext = list(base)
    for l in range(W): ext += [(ix[("Z",)], ix[(l, LT)], E.m), (ix[(l, LT)], ix[("Z",)], -E.m)]
    feas = sp(ext) is not None
    return feas, bound

if __name__ == "__main__":
    fm = [int(x) for x in sys.argv[1:]] or [8, 8, 6, 4]
    rng = np.random.default_rng(7)
    for m, fmax in zip((1, 2, 3, 4), fm):
        E = Eng(m)
        for (name, W, R, w) in words(E):
            analyse(E, name, W, R, w, fmax, rng)
        sys.stdout.flush()
