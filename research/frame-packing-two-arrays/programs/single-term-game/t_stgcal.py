#!/usr/bin/env python3
"""t_stgcal.py: calibrate stg.py's position-level closure against the actual module V in C[G]^2.

Random single-term plays in a toy frame game.  Every generator is an explicit vector of C[G]^2 (random complex
coefficients); V = C-span of the left P-translates of the generators and of C[P]^2.  At every step compare
  (a) stg's predicted available positions (permissive closure on the position graph), with
  (b) the true availability: e_(l,g) in V for g in the coset;
and every composite: the vector is FOUND by linear algebra (a vector of V supported on (1,g_u), (2,g') for some g' in
the target coset), and its twist is compared with the twist stg predicts from the path."""
import sys, random
import numpy as np
from toylab import *
from sweep import GAMES
from stg import *
from t_stg_toys import paid_classes

def closure_nocycle(g, A1, A2, lins):
    N = g.N
    while True:
        par = {}
        def f(a):
            while par.get(a, a) != a: a = par[a]
            return a
        edges = []
        for (t, xs) in lins:
            for x in bits(xs):
                u, w = x, N + g.tw[t][x]; edges.append((u, w))
                ru, rw = f(u), f(w)
                if ru != rw: par[ru] = rw
        V = {}
        for (u, w) in edges:
            V.setdefault(f(u), set()).update((u, w))
        nA1, nA2 = A1, A2
        for r, vs in V.items():
            if any(((A1 >> v) & 1) if v < N else ((A2 >> (v - N)) & 1) for v in vs):
                for v in vs:
                    if v < N: nA1 |= 1 << v
                    else: nA2 |= 1 << (v - N)
        if (nA1, nA2) == (A1, A2): return A1, A2, lins
        A1, A2 = nA1, nA2

KW = 3
def run(name, ci, plays=20, steps=8, seed=1, verbose=False):
    G, P, _ = GAMES[name]()
    hr = paid_classes(G, P)[ci]
    g = from_toylab(G, P, hr, name=name)
    n = G.n; fg = g.fg
    H = fg.H
    rng = random.Random(seed); nrng = np.random.default_rng(seed)
    mism = {"avail": 0, "permissive_only": 0, "missed": 0, "comp_twist": 0, "comp_absent": 0}
    checks = 0
    stats = {"moves": {}, "unlock_by_cycle": 0, "release_jump": 0, "nocycle_rule_would_miss": 0}
    coset_elems = {}
    for e in range(n): coset_elems.setdefault(g.elem_pos[e], []).append(e)
    for play in range(plays):
        gens = []                        # explicit vectors (2n)
        def tr(v, h):                    # left multiply by group element h
            w = np.zeros(2 * n, complex)
            for l in range(2):
                for e in np.nonzero(v[l * n:(l + 1) * n])[0]:
                    w[l * n + G.m(h, int(e))] += v[l * n + e]
            return w
        def span_basis():
            rows = []
            for l in range(2):
                for p in P:
                    v = np.zeros(2 * n, complex); v[l * n + p] = 1; rows.append(v)
            for v in gens:
                for p in P: rows.append(tr(v, p))
            A = np.array(rows)
            u, s, vh = np.linalg.svd(A, full_matrices=False)
            r = int((s > 1e-9 * s[0]).sum())
            return vh[:r]
        def inV(Bs, v):
            res = v - Bs.conj().T @ (Bs @ v) if False else v - (v @ Bs.conj().T) @ Bs
            return np.linalg.norm(res) < 1e-8 * max(1, np.linalg.norm(v))
        A1 = A2 = 1                      # abstract state
        lins = []                        # [twist, xs]
        pairs = []                       # explicit pair generators: (vec, lineage id, g1, g2)
        for st in range(steps):
            Bs = span_basis()
            # true availability
            tA = [0, 0]
            for l in range(2):
                for x, els in coset_elems.items():
                    v = np.zeros(2 * n, complex); v[l * n + els[0]] = 1
                    if inV(Bs, v): tA[l] |= 1 << x
            pA1, pA2, _ = closure(g, A1, A2, tuple((t, xs) for (t, xs) in lins))
            # control: closure WITHOUT the cycle rule (only loops unlock) must miss whenever a cycle unlocked something
            g._nocyc = True
            qA1, qA2, _ = closure_nocycle(g, A1, A2, tuple((t, xs) for (t, xs) in lins))
            if (qA1, qA2) != (pA1, pA2):
                stats["unlock_by_cycle"] += 1
                if (qA1, qA2) != (tA[0], tA[1]): stats["nocycle_rule_would_miss"] += 1
            if bin(pA1).count("1") + bin(pA2).count("1") > bin(A1).count("1") + bin(A2).count("1"): stats["release_jump"] += 1
            checks += 1
            if (pA1, pA2) != (tA[0], tA[1]):
                over = (tA[0] & ~pA1) | (tA[1] & ~pA2)
                if over: mism["missed"] += 1
                else: mism["permissive_only"] += 1
                if verbose: print("mismatch", name, play, st, bin(pA1), bin(tA[0]), bin(pA2), bin(tA[1]))
            A1, A2 = pA1, pA2
            # choose a move
            opts = ["s1", "s2"]
            locked = [i for i, (v, lid, g1, g2) in enumerate(pairs) if not ((A1 >> g.elem_pos[g1]) & 1)]
            if locked: opts += ["c", "c"]
            if A1 != 1 or A2 != 1: opts += ["p", "p"]
            else: opts += ["p"]
            comps = components_tree_paths(g, A1, A2, tuple((t, xs & ~A1) for (t, xs) in lins if xs & ~A1))
            if comps: opts += ["k"] * KW
            kind = rng.choice(opts)
            stats["moves"][kind] = stats["moves"].get(kind, 0) + 1
            h = rng.choice(H)
            if kind in ("s1", "s2"):
                l = int(kind[1]) - 1
                xs_av = [x for x in range(g.N) if ((A1 if l == 0 else A2) >> x) & 1]
                x = rng.choice(xs_av); e = rng.choice(coset_elems[x])
                v = np.zeros(2 * n, complex); v[l * n + e] = 1
                w = tr(v, h); gens.append(w)
                x2 = g.elem_pos[G.m(h, e)]
                if l == 0: A1 |= 1 << x2
                else: A2 |= 1 << x2
            elif kind == "c":
                i = rng.choice(locked); v, lid, g1, g2 = pairs[i]
                w = tr(v, h); gens.append(w)
                n1, n2 = G.m(h, g1), G.m(h, g2)
                pairs.append((w, lid, n1, n2)); lins[lid][1] |= 1 << g.elem_pos[n1]
            elif kind == "p":
                x0 = rng.choice([x for x in range(g.N) if (A1 >> x) & 1])
                y0 = rng.choice([y for y in range(g.N) if (A2 >> y) & 1])
                g1 = rng.choice(coset_elems[x0]); g2 = rng.choice(coset_elems[y0])
                v = np.zeros(2 * n, complex); v[g1] = nrng.normal() + 1j * nrng.normal(); v[n + g2] = nrng.normal() + 1j * nrng.normal()
                w = tr(v, h); gens.append(w)
                n1, n2 = G.m(h, g1), G.m(h, g2)
                t = g.elem_tw[G.m(int(G.inv[n1]), n2)]
                lins.append([t, 1 << g.elem_pos[n1]]); pairs.append((w, len(lins) - 1, n1, n2))
            else:
                (t, u) = rng.choice(sorted(comps))
                # find the composite by linear algebra: a vector of V supported on (1, gu) and (2, g') with g' in coset w
                wpos = g.tw[t][u]
                gu = coset_elems[u][0]
                found = None
                for gp in coset_elems[wpos]:
                    a = np.zeros(2 * n, complex); a[gu] = 1
                    b = np.zeros(2 * n, complex); b[n + gp] = 1
                    # intersection of span{a, b} with V: solve [a b] c in V  <=> residual matrix rank
                    R = np.array([a - (a @ Bs.conj().T) @ Bs, b - (b @ Bs.conj().T) @ Bs]).T
                    u_, s_, vh_ = np.linalg.svd(R)
                    if s_[-1] < 1e-8:
                        c = vh_[-1].conj()
                        if abs(c[0]) > 1e-6 and abs(c[1]) > 1e-6: found = (gp, c[0] * a + c[1] * b); break
                if found is None:
                    mism["comp_absent"] += 1
                    continue
                gp, vec = found
                tw_true = g.elem_tw[G.m(int(G.inv[gu]), gp)]
                if tw_true != t: mism["comp_twist"] += 1
                w = tr(vec, h); gens.append(w)
                n1, n2 = G.m(h, gu), G.m(h, gp)
                lins.append([tw_true, 1 << g.elem_pos[n1]]); pairs.append((w, len(lins) - 1, n1, n2))
    print("%s class %d: %d plays x %d steps, %d availability checks; mismatches %s; %s" % (name, ci, plays, steps, checks, mism, stats), flush=True)

if __name__ == "__main__":
    for spec in sys.argv[1:]:
        nm, ci, pl, stp = spec.split(":"); run(nm, int(ci), plays=int(pl), steps=int(stp))
