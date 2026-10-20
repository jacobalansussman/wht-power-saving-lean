#!/usr/bin/env python3
"""vunits.py (v-fp-capacity): rebuild fp-capacity's calibration units with my own engine (vengine.py) and test
Theorem A (Sp labels, raw gauge, level potential), Theorem A' (image frames: every level step crossed >= W times by
live frames; Delta <= P_img), Corollary A'2 (co-image frames, gauge K^-1), and P* (exact min over relative gauges, W=2).
usage: python3 vunits.py seed n_random"""
import sys, json, random, itertools
import numpy as np
from vengine import Eng, run, check_unit, raw_arrays, analyse, moves_report

def words(E, rng, nrand):
    m = E.m; D = E.dirm; e = [1 << j for j in range(m)]
    out = []
    out.append(("plain W=1", 1, 1, [("move", 0, D(z)) for z in e], {0: 0}))
    out.append(("plain W=2", 2, 2, [("move", r, D(z)) for r in (0, 1) for z in e], {0: 0, 1: 1}))
    H2 = [[1, 1], [1, -1]]; H2i = [[.5, .5], [.5, -.5]]
    out.append(("pair trick", 2, 2, [("gate", [0, 1], H2)] + [("move", r, D(z)) for z in e for r in (0, 1)] + [("gate", [0, 1], H2i)], {0: 0, 1: 1}))
    w = [("move", 0, D(e[0])), ("gate", [0, 1], H2)] + [("move", r, D(z)) for z in e[1:] for r in (0, 1)] + [("gate", [0, 1], H2i), ("move", 1, D(e[0]))]
    out.append(("crossing pair", 2, 2, w, {0: 0, 1: 1}))
    if m >= 3:
        S = E.addr([1, 2, 0] + list(range(3, m))); S2 = S @ S
        assert np.allclose(S @ S2, E.I) and np.allclose(S @ E.Kmat, E.Kmat @ S)
        w = [("mix", 1, [(0, 1, None), (1, -1, S)]), ("mix", 0, [(0, 1, None), (0, 1, S2), (1, -1, S2)])]
        w += [("move", r, D(z)) for z in e for r in (0, 1)]
        w += [("mix", 1, [(0, .5, None), (1, -.5, None), (0, -.5, S), (1, .5, S), (0, .5, S2), (1, -.5, S2)]),
              ("mix", 0, [(0, 1, None), (1, -1, None)])]
        out.append(("order-3 control C%dO3" % (2 * m), 2, 2, w, {0: 0, 1: 1}))
        # random order-3 packers: random move interleaving, random extra free letters inside
        for i in range(nrand):
            perm = list(range(m)); rng.shuffle(perm)
            ph = E.phase(rng.randrange(1, 1 << m))
            mv = [("move", r, D(z)) for z in e for r in (0, 1)]; rng.shuffle(mv)
            w2 = w[:2] + mv + w[-2:]
            out.append(("random order-3 interleave %d" % i, 2, 2, w2, {0: 0, 1: 1}))
    # random crossing units (no saving, live image-frame packing)
    for i in range(nrand):
        W = rng.choice([2, 3]) if m >= 3 else 2
        R = W
        csz = rng.randrange(1, m)
        Cset = rng.sample(e, csz)
        rest = [z for z in e if z not in Cset]
        pre = [[z for z in rest if rng.random() < .5] for a in range(W)]
        suf = [[z for z in rest if z not in pre[a]] for a in range(W)]
        # gate: product of elementary Gaussian-integer matrices (exact inverse)
        G = np.eye(W, dtype=complex); els = []
        for _ in range(3):
            a, b = rng.sample(range(W), 2); c = complex(rng.choice([-1, 1, 2, -2]), rng.choice([0, 1, -1]))
            Eab = np.eye(W, dtype=complex); Eab[a, b] = c; els.append(Eab); G = Eab @ G
        Gi = np.linalg.inv(G)
        assert np.allclose(np.round(Gi * 4) / 4, Gi); Gi = np.round(Gi * 4) / 4
        w = [("move", a, D(z)) for a in range(W) for z in pre[a]]
        w += [("gate", list(range(W)), G.tolist())]
        w += [("move", a, D(z)) for z in Cset for a in range(W)]
        w += [("gate", list(range(W)), Gi.tolist())]
        w += [("move", a, D(z)) for a in range(W) for z in suf[a]]
        out.append(("random crossing W=%d |C|=%d #%d" % (W, csz, i), W, R, w, {a: a for a in range(W)}))
    return out

def analyse_unit(E, name, W, R, w, deliver, rng, fmax):
    m = E.m
    snaps = run(E, W, R, w)
    rec = dict(m=m, name=name, W=W, M=sum(1 for s in w if s[0] == "move"))
    rec["Delta"] = W * m - rec["M"]
    rec["formal_identity"] = check_unit(E, W, snaps, deliver)
    rec["raw_err_f"] = {f: raw_arrays(E, W, R, w, deliver, f, np.random.default_rng(f)) for f in range(1, fmax + 1)}
    L0 = E.L0(); L1base = [E.apply(E.sp(E.Kk), b) for b in L0]
    Lfull = E.frame(tuple(1 << j for j in range(2 * m)), L0)
    lev = lambda F: m - (len(F & Lfull).bit_length() - 1)
    viol = []
    # Sp labels, raw gauge, level potential
    cls = lambda g, l: E.sp(g)
    nodes, live = analyse(E, W, R, w, snaps, deliver, cls)
    levsp = lambda S: lev(E.frame(S, L0))
    rep = moves_report(E, W, R, w, snaps, nodes, live, cls, levsp)
    rec["P_live_sp_raw"] = sum(max(k - 1, 0) for k, _ in rep)
    rec["P_asc_sp_raw"] = sum(sum(1 for a, b in st if b == a + 1) - 1 for k, st in rep)
    if rec["Delta"] > rec["P_asc_sp_raw"]: viol.append("Delta > P_asc(sp)")
    # image frames, raw gauge
    clsI = lambda g, l: E.frame(E.sp(g), L0)
    nodes, live = analyse(E, W, R, w, snaps, deliver, clsI)
    rep = moves_report(E, W, R, w, snaps, nodes, live, clsI, lev)
    cross = [0] * m
    for k, st in rep:
        for a, b in st:
            if b == a + 1: cross[a] += 1
    rec["img_cross_per_level"] = cross
    rec["P_img"] = sum(sum(1 for a, b in st if b == a + 1) - 1 for k, st in rep)
    rec["P_live_img"] = sum(max(k - 1, 0) for k, _ in rep)
    rec["moves_raising_2plus_img"] = sum(1 for k, st in rep if sum(1 for a, b in st if b == a + 1) >= 2)
    rec["moves_lifting_2_off_L0"] = sum(1 for k, st in rep if sum(1 for a, b in st if a == 0 and b == 1) >= 2)
    if min(cross) < W: viol.append("img level step crossed < W")
    if rec["Delta"] > rec["P_img"]: viol.append("Delta > P_img")
    # co-image frames: class gbar L1 (gauge K^-1, H = B); co-level = d(gbar L1, L0) runs m -> 0
    clsC = lambda g, l: E.frame(E.sp(g), L1base)
    nodes, live = analyse(E, W, R, w, snaps, deliver, clsC)
    rep = moves_report(E, W, R, w, snaps, nodes, live, clsC, lev)
    crossC = [0] * m
    for k, st in rep:
        for a, b in st:
            if b == a - 1: crossC[b] += 1
    rec["coimg_cross_per_colevel"] = crossC
    if min(crossC) < W: viol.append("co-img step crossed < W")
    # pair-Levi classes: class of g = (gbar L0, gbar L1); f = d(F, L0) - d(F', L0); weight max(0, df)/2
    clsO = lambda g, l: (E.frame(E.sp(g), L0), E.frame(E.sp(g), L1base))
    fO = lambda k: lev(k[0]) - lev(k[1])
    nodes, live = analyse(E, W, R, w, snaps, deliver, clsO)
    rep = moves_report(E, W, R, w, snaps, nodes, live, clsO, fO)
    wsum = sum(max(0, b - a) / 2 for k, st in rep for a, b in st)
    rec["opp_weight_sum"] = wsum; rec["P_opp_w"] = sum(sum(max(0, b - a) / 2 for a, b in st) - 1 for k, st in rep)
    if wsum < W * m: viol.append("opp weights < Wm")
    if rec["Delta"] > rec["P_opp_w"] + 1e-12: viol.append("Delta > P_opp_w")
    # P* (W = 2): exact min over relative gauges h = (1, rho)
    if W == 2:
        cand = set()
        for t in range(len(snaps)):
            for r in range(R):
                g1 = [g for (g, l) in snaps[t][r] if l == 0]; g2 = [g for (g, l) in snaps[t][r] if l == 1]
                for a in g1:
                    for b in g2: cand.add(E.comp(E.spinv(b), E.sp(a)))
        # a generic relative gauge: random Sp element outside cand
        gen = None
        rr = random.Random(5)
        letters = [E.dirm(z) for z in range(1, E.n)] + [E.phase(z) for z in range(1, E.n)]
        while gen is None:
            M = E.I
            for _ in range(12): M = letters[rr.randrange(len(letters))] @ M
            Sg = E.sp_of_mat(M)
            if Sg not in cand: gen = Sg
        ident = tuple(1 << j for j in range(2 * m))
        vals = {}
        for rho in list(cand) + [gen]:
            gauge = [ident, rho]
            c2 = lambda g, l, gauge=gauge: E.comp(E.sp(g), gauge[l])
            nodes, live = analyse(E, W, R, w, snaps, deliver, c2)
            rep = moves_report(E, W, R, w, snaps, nodes, live, c2, None)
            Pl = sum(max(k - 1, 0) for k, _ in rep)
            tot = sum(k for k, _ in rep)
            if tot < 2 * m: viol.append("sum k_j < 2m in a gauge")
            vals[rho] = Pl
        rec["P_star"] = min(vals.values()); rec["P_live_generic"] = vals[gen]; rec["n_rel_gauges"] = len(vals)
        if rec["Delta"] > rec["P_star"]: viol.append("Delta > P*")
    rec["violations"] = viol
    return rec

if __name__ == "__main__":
    seed = int(sys.argv[1]); nr = int(sys.argv[2])
    rng = random.Random(seed); allrec = []
    for m, fmax in ((2, 8), (3, 6), (4, 4)):
        E = Eng(m)
        for (name, W, R, w, dl) in words(E, rng, nr):
            rec = analyse_unit(E, name, W, R, w, dl, rng, fmax)
            allrec.append(rec)
            print(json.dumps(rec), flush=True)
    nv = sum(len(r["violations"]) for r in allrec)
    nf = sum(1 for r in allrec if not r["formal_identity"])
    mx = max(max(r["raw_err_f"].values()) for r in allrec)
    print("SUMMARY units=%d formal_identity_failures=%d max_raw_err=%.2e violations=%d" % (len(allrec), nf, mx, nv))
    json.dump(allrec, open("vunits_s%d.json" % seed, "w"), indent=1)
