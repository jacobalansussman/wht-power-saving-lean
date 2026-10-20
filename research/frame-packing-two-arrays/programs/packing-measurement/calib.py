#!/usr/bin/env python3
"""fp-measure calibration.
(a) raw-gauge ledger (P, B, D, Z, I0) of the known small units, with sl-combinatorics' exact engine ledger.py
    (imported read-only): plain W = 1..3, gauge exchange, the order-3 controls C6O3 / C8O3, the +-1 halving controls,
    the packing gadget, at blocks 2-4.
(b) the frame-model level formula used by fm_ledger.py / oa_ledger.py, checked exactly in the engine: for every pair
    of subspaces U, S of F_2^m that have an orthonormal basis (m = 2, 3, 4), Phi_U = product of dir z over an
    orthonormal basis of U, and the raw label Phi_U Phi_S^-1 has level dim U + dim S - 2 dim(U cap S); plus the
    Lagrangian form d(F_U, F_S) = m - dim(F_U cap F_S) for EVERY pair of subspaces (degenerate ones included).
(c) the gauge-exchange unit replayed in the frame model gives the engine's raw numbers.
usage: python3 calib.py   (numpy needed)"""
import sys, os, itertools, random, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ledger"))
from ledger import *          # Engine, run, plain_word, Fr, gauge_trivial

out = dict(a=[], b={}, c=[])
rng = random.Random(7)


def row(name, E, W, R, word):
    r = run(E, W, R, word)
    d = dict(case=name, m=E.m, W=W, M=r["M"], unit=r["unit"], P=r["P"], B=r.get("B"), D=r["D"], Z=r["Z"],
             I0=r["I0"], maxlab=r["maxlab"], Delta=W * E.m - r["M"], identity=r.get("unit_ok"))
    print("%-22s m=%d W=%d M=%2d unit=%s | P=%d B=%s D=%d Z=%d I0=%d | Delta=%d = P-B-2D-Z-I0: %s" %
          (name, E.m, W, r["M"], r["unit"], r["P"], r.get("B"), r["D"], r["Z"], r["I0"], d["Delta"], r.get("unit_ok")))
    out["a"].append(d)


print("== (a) raw-gauge ledgers of the known small units (ledger.py engine, exact)")
for m in (2, 3, 4):
    E = Engine(m)
    for W in (1, 2, 3):
        row("plain", E, W, W, plain_word(m, W, rng))
    w = [("move", 0, 1), ("mix", 0, [(0, 1, None), (1, 1, None)])]
    for j in range(1, m):
        w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
    w += [("mix", 0, [(0, 1, None), (1, -1, None)]), ("move", 1, 1)]
    row("gauge exchange", E, 2, 2, w)
    if m >= 3:
        s = E.addr([2, 4, 1] + [1 << j for j in range(3, m)]); s2 = E.addr([4, 1, 2] + [1 << j for j in range(3, m)])
        w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s)])]
        for j in range(m):
            w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        h = Fr(1, 2)
        w += [("mix", 1, [(0, h, None), (1, -h, None), (0, -h, s), (1, h, s), (0, h, s2), (1, -h, s2)]),
              ("mix", 0, [(0, 1, None), (1, -1, None)])]
        row("order-3 control C%dO3" % (2 * m), E, 2, 2, w)
        w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -2, None)])]
        for j in range(m):
            w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        w += [("mix", 1, [(0, h, None), (1, -h, None)]), ("mix", 0, [(0, 1, None), (1, -1, None)])]
        row("+-1 halving C%dH" % (2 * m), E, 2, 2, w)
    w = [("mix", 2, [(1, 1, None)]), ("mix", 0, [(0, 1, None), (1, 1, None)]), ("move", 0, 1),
         ("mix", 0, [(0, 1, None), (2, -1, None)]), ("move", 0, 2), ("move", 2, 2),
         ("mix", 0, [(0, 1, None), (2, 1, None)]), ("move", 1, 1), ("move", 1, 2),
         ("mix", 0, [(0, 1, None), (1, -1, None)])]
    for j in range(2, m):
        w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
    row("packing gadget", E, 2, 3, w)

print("\n== (b) the level formula of the frame model")
def pc(x): return bin(x).count("1")
def dt(a, b): return pc(a & b) & 1
def span(vs):
    S = {0}
    for v_ in vs:
        S |= {x ^ v_ for x in S}
    return frozenset(S)
def all_subspaces(m):
    seen = set()
    for k in range(m + 1):
        for vs in itertools.combinations(range(1, 1 << m), k):
            S = span(vs)
            if len(S) == 1 << k:
                seen.add(S)
    return sorted(seen, key=lambda S: (len(S), sorted(S)))
def dim(S): return len(S).bit_length() - 1
def onb(S):
    """an orthonormal basis (odd weights, pairwise orthogonal) of S, or None"""
    k = dim(S)
    odd = [x for x in S if x and pc(x) & 1]
    for vs in itertools.combinations(odd, k):
        if all(dt(a, b) == 0 for a, b in itertools.combinations(vs, 2)) and span(vs) == S:
            return vs
    return None
def lagr(S, m):
    """F_S = {(c | a + c) : c in S, a in S^perp} as a set of 2m-bit ints (x part high, z part low)"""
    perp = [a for a in range(1 << m) if all(dt(a, c) == 0 for c in S)]
    return frozenset((c << m) | (a ^ c) for c in S for a in perp)
for m in (2, 3, 4):
    E = Engine(m)
    subs = all_subspaces(m)
    good = [(S, onb(S)) for S in subs]
    good = [(S, b) for S, b in good if b is not None]
    phi = {}
    for S, b in good:
        gid = E.ID
        for z in b:
            gid = E.lmul("D%d" % z, gid)[0]
        phi[S] = gid
    bad = 0
    for (U, _), (S, _) in itertools.product(good, repeat=2):
        g = E.mul(phi[U], E.inv(phi[S])[0])[0]
        lev = E.level(E.syms[g])
        if lev != dim(U) + dim(S) - 2 * dim(U & S):
            bad += 1
    badL = 0
    for U, S in itertools.product(subs, repeat=2):
        FU, FS = lagr(U, m), lagr(S, m)
        dL = m - (len(FU & FS).bit_length() - 1)
        if dL != dim(U) + dim(S) - 2 * dim(U & S):
            badL += 1
    print("m=%d: %d subspaces, %d with an orthonormal basis; engine levels of Phi_U Phi_S^-1: %d pairs, %d mismatches;"
          " Lagrangian distances, all %d pairs: %d mismatches" % (m, len(subs), len(good), len(good) ** 2, bad,
                                                                   len(subs) ** 2, badL))
    out["b"][m] = dict(subspaces=len(subs), onb=len(good), engine_pairs=len(good) ** 2, engine_bad=bad,
                       lagr_pairs=len(subs) ** 2, lagr_bad=badL)

print("\n== (c) gauge exchange in the frame model (coherent frames; input 1 twisted by Phi_<e1>)")
for m in (2, 3, 4):
    # frames as subspaces (bitmask sets); classes: input 0 source 0, input 1 source <e1>
    def lev(U, S): return dim(U) + dim(S) - 2 * dim(U & S)
    S0, S1 = span([]), span([1])
    P = Bm = D = M = 0
    r0 = span([]); r1 = span([1])        # coherent frames: r0 at 0, r1 sits at its source <e1>
    c0 = [S0]; c1 = [S1]
    def climb(U, z, cls):
        global P, D, M
        U2 = span(list(U) + [z]) if z not in U else frozenset(x for x in U if True)
        return U2
    # move r0 along e1: 0 -> <e1>
    U = span([1]); M += 1; D += sum(1 for S in c0 if lev(U, S) < lev(r0, S)); P += len(c0) - 1; r0 = U
    # mix r0 += r1 at frame <e1>
    before = sum(lev(r0, S) for S in c0); c0 = [S0, S1]; Bm -= sum(lev(r0, S) for S in c0) - before
    for j in range(1, m):
        for which in (0, 1):
            Ucur = r0 if which == 0 else r1
            cls = c0 if which == 0 else c1
            U = span(list(Ucur) + [1 << j])
            M += 1; P += len(cls) - 1; D += sum(1 for S in cls if lev(U, S) < lev(Ucur, S))
            if which == 0: r0 = U
            else: r1 = U
    before = sum(lev(r0, S) for S in c0); c0 = [S0]; Bm -= sum(lev(r0, S) for S in c0) - before
    # move r1 along e1: full -> <e1>^perp (a descent of the frame); its label S1 ascends
    U = span([1 << j for j in range(1, m)])
    M += 1; P += len(c1) - 1; D += sum(1 for S in c1 if lev(U, S) < lev(r1, S)); r1 = U
    psi_end = sum(lev(r0, S) for S in c0) + sum(lev(r1, S) for S in c1)
    print("m=%d: frame model M=%d P=%d B=%d D=%d; ends at levels %d (r0) and %d (r1); engine (a): M=2m, P=B=m-1, D=0"
          % (m, M, P, Bm, D, sum(lev(r0, S) for S in c0), sum(lev(r1, S) for S in c1)))
    out["c"].append(dict(m=m, M=M, P=P, B=Bm, D=D, psi_end=psi_end))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "calib.json"), "w"), indent=1,
          default=str)
