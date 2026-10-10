#!/usr/bin/env python3
"""xsub.py -- device "x array as carrier before the scatter", applied to an accepted E8 certificate.
1. every helper-to-helper gate of phase B that does not depend on a delivery is moved to phase A (before the totals);
2. a helper that only ever holds one source x_T (loaded on line(T)) and is only read by such gates is deleted; the
   gates read x_T itself, which climbs through the same frames on its way to the full frame;
3. (swap) if such a helper is the host of an in-place sum whose other part dies there, the other part becomes the host.
usage: xsub.py <in cert> <out cert> [noswap] [nototals]"""
import sys, json; sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__))); import e8
from fractions import Fraction as Fr
from collections import defaultdict

def load_gates(c):
    fr = [tuple(f) for f in c["frames"]]
    def conv(g):
        if g[0] == "out":
            return ["out", fr[g[1]], g[2], [[t, Fr(n, d)] for t, n, d in g[3]], list(g[4]) if len(g) > 4 else []]
        return ["in", fr[g[1]], g[2], [[s, Fr(n, d)] for s, n, d in g[3]]]
    return [conv(g) for g in c["A"]], [conv(g) for g in c["B"]]

def regs(g):
    if g[0] == "out": return [g[2]] + [t for t, _ in g[3]] + list(g[4])
    return [g[2]] + [s for s, _ in g[3]]

def reads(g):          # registers read
    return [g[2]] if g[0] == "out" else [s for s, _ in g[3]]
def writes(g):
    return [t for t, _ in g[3]] if g[0] == "out" else [g[2]]

def main(src, dst, swap=True, verbose=True):
    c = e8.load(src)
    v, R = c["v"], c["R"]
    A, B = load_gates(c)
    isx = lambda r: r < v
    isy = lambda r: v <= r < 2 * v
    ish = lambda r: r >= 2 * v
    ret = [s for k, s, f in c["ret"]]
    # ---- 1. move rules
    nload = 0
    while nload < len(A) and A[nload][0] == "out" and isx(A[nload][2]) and len(A[nload][1]) == 1: nload += 1
    loads, totg = A[:nload], A[nload:]
    assert all(g[0] == "in" and g[2] in ret for g in totg), "phase A is not loads + totals"
    taint = set(); moved = []; stay = []
    for g in B:
        rg = regs(g)
        if any(isy(r) for r in rg) or any(r in taint for r in rg) or any(isx(r) for r in rg):
            stay.append(g); taint.update(rg)
        else:
            moved.append(g)
    if verbose: print("rules moved to phase A: %d of %d B gates; staying %d" % (len(moved), len(B), len(stay)))
    mid = moved
    # ---- index
    def occurrences(gl):
        occ = defaultdict(list)
        for i, g in enumerate(gl):
            for r in set(regs(g)): occ[r].append(i)
        return occ
    occ_mid = occurrences(mid); occ_stay = occurrences(stay)
    loadof = {}
    for g in loads:
        for t, co in g[3]:
            assert co == 1
            loadof[t] = g[2]
    xframes = {T: [(-1, e8.line(T))] for T in range(v)}      # (gate index in mid, frame) where x_T is named
    dead = set(); nsub = nswap = 0
    def x_ok(T, extra):
        ev = sorted(xframes[T] + extra)
        for (i, f), (j, g) in zip(ev, ev[1:]):
            if not (f == g or (len(f) < len(g) and e8.inside(f, g))): return False
        return True
    # candidates in order of helper number
    for rnd, h in [(0, h) for h in sorted(loadof)] + [(1, h) for h in sorted(loadof)]:
        if h in dead: continue
        T = loadof[h]
        if h in ret or h in occ_stay: 
            # h is used after the scatter: its pure period must end inside phase A by a swap; skip unless swap below
            pass
        idx = occ_mid.get(h, [])
        # pure period: gates of mid where h is only read, up to its first write
        firstw = None
        for i in idx:
            if h in writes(mid[i]): firstw = i; break
        pure = [i for i in idx if firstw is None or i < firstw]
        if firstw is None:
            if rnd == 1 or NOPLAIN: continue
            if h in occ_stay or h in ret: continue          # read after the scatter while still a single source: x is at FULL then
            ev = [(i, mid[i][1]) for i in pure]
            if not x_ok(T, ev): continue
            for i in pure:
                g = mid[i]
                if g[0] == "out": g[2] = T
                else:
                    for e in g[3]:
                        if e[0] == h: e[0] = T
            xframes[T] += ev; dead.add(h); nsub += 1
            continue
        if not swap or rnd == 0: continue
        G = mid[firstw]
        if G[0] != "in" or G[2] != h: continue
        # a partner that dies at G
        part = None
        for s, co in G[3]:
            if ish(s) and s not in ret and s not in occ_stay and occ_mid[s][-1] == firstw and abs(co) == 1 and [x for x, _ in G[3]].count(s) == 1:
                part = (s, co); break
        if part is None: continue
        q1, c1 = part
        ev = [(i, mid[i][1]) for i in pure] + [(firstw, G[1])]
        if not x_ok(T, ev): continue
        for i in pure:
            g = mid[i]
            if g[0] == "out": g[2] = T
            else:
                for e in g[3]:
                    if e[0] == h: e[0] = T
        # new G: q1 += (1/c1) x_T + sum (cj/c1) qj    ->  q1 = (1/c1) * (old content of h after G)
        mid[firstw] = ["in", G[1], q1, [[T, 1 / c1]] + [[s, co / c1] for s, co in G[3] if s != q1]]
        # later occurrences of h: h held Z, now q1 holds Z / c1
        def fix(g):
            if g[0] == "out":
                if g[2] == h:
                    g[2] = q1; g[3] = [[t, co * c1] for t, co in g[3]]
                else:
                    g[3] = [[q1, co / c1] if t == h else [t, co] for t, co in g[3]]
                    g[4] = [q1 if r == h else r for r in g[4]]
            else:
                if g[2] == h:
                    g[2] = q1; g[3] = [[s, co / c1] for s, co in g[3]]
                else:
                    g[3] = [[q1, co * c1] if s == h else [s, co] for s, co in g[3]]
        for i in idx:
            if i > firstw: fix(mid[i])
        for i in occ_stay.get(h, []): fix(stay[i])
        # update indices
        occ_mid[q1] = sorted(set(occ_mid[q1] + [i for i in idx if i > firstw]))
        if h in occ_stay: occ_stay[q1] = sorted(set(occ_stay.get(q1, []) + occ_stay[h])); del occ_stay[h]
        occ_mid[h] = []
        xframes[T] += ev; dead.add(h); nswap += 1
    if verbose: print("helpers replaced by their x array: %d plain, %d by host swap" % (nsub, nswap))
    # ---- rebuild
    keep = [r for r in range(2 * v, 2 * v + R) if r not in dead]
    ren = {r: r for r in range(2 * v)}
    for i, r in enumerate(keep): ren[r] = 2 * v + i
    def out(g):
        if g[0] == "out":
            tg = [[ren[t], co.numerator, co.denominator] for t, co in g[3] if t not in dead]
            return ["out", list(g[1]), ren[g[2]], tg, [ren[r] for r in g[4]]]
        return ["in", list(g[1]), ren[g[2]], [[ren[s], co.numerator, co.denominator] for s, co in g[3]]]
    newA = [out(g) for g in loads]
    newA = [g for g in newA if g[3]]                     # a load gate with no target left is dropped
    newA += [out(g) for g in mid] + [out(g) for g in totg]
    newB = [out(g) for g in stay]
    D = dict(h=c["h"], ports=c["ports"], R=len(keep), A=newA, B=newB, ret=[ren[s] for s in ret], scat=c["scat"],
             name="xsub.py on %s: rules before the scatter, %d single-source helpers replaced by the x array itself (%d by host swap)" % (src.split("/")[-1], nsub + nswap, nswap))
    cert = e8.write_cert(D, dst)
    return cert, nsub, nswap

NOPLAIN = False
if __name__ == "__main__":
    a = sys.argv[1:]
    NOPLAIN = "noplain" in a
    cert, nsub, nswap = main(a[0], a[1], swap="noswap" not in a)
    r = e8.quick(cert)
    print({k: r.get(k) for k in ("accepted_ref", "accepted_mirror", "ref_msg", "mirror_msg", "R", "N", "cst", "whole_block", "per_rank", "adds", "most_regs_one_gate")})
    print("ref stats den/maxabs:", r.get("ref_stats", {}).get("den"), r.get("ref_stats", {}).get("maxabs"))
    P = e8.profile(cert)
    print("x patterns", sorted(P["x"].items(), key=lambda kv: -kv[1])[:8])
