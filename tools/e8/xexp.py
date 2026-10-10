#!/usr/bin/env python3
"""xexp.py -- device "a helper that lives only before the scatter is replaced by direct reads of the x arrays".
Input: a certificate written by xsub.py (phase A = loads on the lines, helper gates, the totals at the full frame).
A helper d that is never named after the scatter and is no total is ELIMINATED: every gate that read d reads instead
the x arrays d held at that moment (they climb to the gate's frame); the gates that wrote d are dropped.
Greedy: candidates by size of what they hold; a candidate is kept if the exact figure (e8.figure) grows.
Then (totals): each total helper takes over the life of one remaining such helper; its junk is cancelled at the full
frame, where all x arrays stand.
usage: xexp.py <in> <out> [maxn] [nototals]
NOTE: the greedy stops after TLIM = 420 seconds of wall-clock time (the recorded run needs about 16)."""
import sys, json, time; sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__))); import e8
from fractions import Fraction as Fr
from collections import defaultdict

SEED = None; MODE = "tier"; PASSES = 2; TLIM = 420; VERB = True; LIFT = False
def main(src, dst, maxn=6, do_totals=True, minden=2):
    c = e8.load(src)
    v, R, h = c["v"], c["R"], c["h"]
    fr = [tuple(f) for f in c["frames"]]
    def conv(g):
        if g[0] == "out":
            assert not (len(g) > 4 and g[4]), "extras not handled"
            return ["out", fr[g[1]], g[2], [[t, Fr(n, d)] for t, n, d in g[3]]]
        return ["in", fr[g[1]], g[2], [[s, Fr(n, d)] for s, n, d in g[3]]]
    A = [conv(g) for g in c["A"]]; B = [conv(g) for g in c["B"]]
    isx = lambda r: r < v
    ret = [s for k, s, f in c["ret"]]; retset = set(ret)
    nload = 0
    while nload < len(A) and A[nload][0] == "out" and isx(A[nload][2]) and len(A[nload][1]) == 1: nload += 1
    ntot = 0
    while A[len(A) - 1 - ntot][0] == "in" and A[len(A) - 1 - ntot][2] in retset: ntot += 1
    loads, mid, totg = A[:nload], A[nload:len(A) - ntot], A[len(A) - ntot:]
    inB = set()
    for g in B:
        inB.add(g[2]); inB.update(r for r, _ in g[3])
    # ---- replay contents, snapshots of every read
    cont = defaultdict(dict)
    for T in range(v): cont[T] = {T: Fr(1)}
    def addc(t, snap, co):
        ct = cont[t]
        for k, x in snap.items():
            y = ct.get(k, 0) + co * x
            if y: ct[k] = y
            else: ct.pop(k, None)
    for g in loads:
        for t, co in g[3]: addc(t, cont[g[2]], co)
    snaps = []                              # per mid gate: {source register: content at the read}
    for g in mid:
        sn = {}
        if g[0] == "out":
            sn[g[2]] = dict(cont[g[2]])
            for t, co in g[3]: addc(t, sn[g[2]], co)
        else:
            for s, co in g[3]:
                sn[s] = dict(cont[s]); addc(g[2], sn[s], co)
        snaps.append(sn)
    endcont = {r: dict(cont[r]) for r in range(2 * v, 2 * v + R)}
    gone = set()
    def build(gone, host=None):
        host = host or {}
        """gate lists with the registers in `gone` eliminated (cascade included); returns (design dict, gone) or None"""
        gone = set(gone)
        while True:
            readers = defaultdict(int)
            for g in mid:
                if g[0] == "out":
                    if any(t not in gone for t, _ in g[3]): readers[g[2]] += 1
                elif g[2] not in gone:
                    for s, _ in g[3]: readers[s] += 1
            more = [r for r in range(2 * v, 2 * v + R) if r not in gone and r not in retset and r not in inB and not readers[r]]
            if not more: break
            gone.update(more)
        newmid = []
        for g, sn in zip(mid, snaps):
            F = g[1]
            if g[0] == "in":
                if g[2] in gone: continue
                acc = {}; order = []
                for s, co in g[3]:
                    if s in gone:
                        for T, x in sn[s].items():
                            if T not in acc: order.append(T)
                            acc[T] = acc.get(T, 0) + co * x
                    else:
                        if s not in acc: order.append(s)
                        acc[s] = acc.get(s, 0) + co
                srcs = [[s, acc[s]] for s in order if acc[s]]
                if srcs: newmid.append(["in", F, g[2], srcs])
            else:
                tg = [[t, co] for t, co in g[3] if t not in gone]
                if not tg: continue
                if g[2] in gone:
                    for t, co in tg:
                        newmid.append(["in", F, t, [[T, co * x] for T, x in sn[g[2]].items()]])
                else:
                    newmid.append(["out", F, g[2], tg])
        keep = [r for r in range(2 * v, 2 * v + R) if r not in gone and r not in host]
        ren = {r: r for r in range(2 * v)}
        for i, r in enumerate(keep): ren[r] = 2 * v + i
        for dd, tk in host.items(): ren[dd] = ren[tk]
        tg2 = []
        for g in totg:
            junk = [endcont[dd] for dd, tk in host.items() if tk == g[2]]
            if junk:
                acc = {}; order = []
                for sx, co in g[3]:
                    acc[sx] = acc.get(sx, 0) + co; order.append(sx)
                for T, x in junk[0].items():
                    if T not in acc: order.append(T)
                    acc[T] = acc.get(T, 0) - x
                g = ["in", g[1], g[2], [[sx, acc[sx]] for sx in order if acc[sx]]]
            tg2.append(g)
        def out(g):
            if g[0] == "out":
                return ["out", list(g[1]), ren[g[2]], [[ren[t], co.numerator, co.denominator] for t, co in g[3] if t not in gone], []]
            return ["in", list(g[1]), ren[g[2]], [[ren[s], co.numerator, co.denominator] for s, co in g[3]]]
        newA = [g for g in (out(g) for g in loads) if g[3]] + [out(g) for g in newmid] + [out(g) for g in tg2]
        D = dict(h=h, ports=c["ports"], R=len(keep), A=newA, B=[out(g) for g in B], ret=[ren[s] for s in ret], scat=c["scat"], name="")
        return D, gone
    def lift(D):
        """a gate whose registers do not all stand inside its frame is moved up to the join (frames only matter to the
        label side; the contents do not change)"""
        cur = {}
        for T in range(v): cur[T] = e8.line(T)
        nl = 0
        for ph in "AB":
            for g in D[ph]:
                F = tuple(g[1]); J = F
                rs = [g[2]] + [e[0] for e in g[3]]
                for r in rs:
                    f = cur.get(r, ())
                    if f and f != J and not e8.inside(f, J): J = e8.join(J, f)
                if J != F:
                    g[1] = list(J); nl += 1
                for r in rs: cur[r] = J
        return nl
    def score(D):
        if LIFT: lift(D)
        try:
            cert = e8.build_cert(D, strict=True)
        except e8.DesignError:
            return None, None
        f = e8.figure(e8.inv_hist(cert), h, v, cert["R"], cert["cst"])
        return f["whole_block"], cert
    D0, gone = build(set())
    best, cert = score(D0)
    if VERB: print("start: R=%d figure=%s (cascade removed %d)" % (D0["R"], best, len(gone)), flush=True)
    # candidates: helpers living only in phase A
    cand = []
    for d in range(2 * v, 2 * v + R):
        if d in retset or d in inB or d in gone: continue
        sz = 0; ok = True
        for g, sn in zip(mid, snaps):
            if d in sn:
                sz = max(sz, len(sn[d]))
                if any(x.denominator > minden for x in sn[d].values()): ok = False
        if ok and 0 < sz <= maxn: cand.append((sz, d))
    cand.sort()
    if SEED is not None:
        import random
        rng = random.Random(SEED)
        if MODE == "tier":
            cand = [(sz, d) for sz, _, d in sorted((sz, rng.random(), d) for sz, d in cand)]
        elif MODE == "soft":
            cand = [(sz, d) for _, sz, d in sorted((sz + 2.5 * rng.random(), sz, d) for sz, d in cand)]
        else:
            rng.shuffle(cand)
    if VERB: print("candidates (only before the scatter, at most %d sources when read): %d" % (maxn, len(cand)), flush=True)
    t0 = time.time(); nacc = 0
    rej = {"frames": 0, "nogain": 0}
    for rnd in range(PASSES):
        acc0 = nacc
        for sz, d in cand:
            if d in gone: continue
            D1, g1 = build(gone | {d})
            f1, c1 = score(D1)
            if f1 is not None and f1 > best:
                best, gone, cert, nacc = f1, g1, c1, nacc + 1
            elif rnd == 0:
                rej["frames" if f1 is None else "nogain"] += 1
            if time.time() - t0 > TLIM: print("time limit in greedy"); break
        if nacc == acc0: break
    if VERB: print("rejected in the first pass:", rej)
    print("eliminated %d helpers (accepted %d candidates): R=%d figure=%d  [%.0f s]" % (len(gone), nacc, cert["R"], best, time.time() - t0), flush=True)
    host = {}
    if do_totals:
        rest = sorted((len(endcont[d]), d) for d in range(2 * v, 2 * v + R)
                      if d not in retset and d not in inB and d not in gone and all(x.denominator == 1 for x in endcont[d].values()))
        for tk in ret:
            for _, d in rest:
                if d in host: continue
                h2 = dict(host); h2[d] = tk
                D1, g1 = build(gone, h2)
                if g1 != gone: continue
                f1, c1 = score(D1)
                if f1 is not None and f1 > best:
                    best, cert, host = f1, c1, h2
                    break
        print("totals that first serve as a helper before the scatter: %d: R=%d figure=%d" % (len(host), cert["R"], best), flush=True)
    cert["derived_from"] = "xexp.py on %s: %d helpers that lived only before the scatter replaced by direct reads of x arrays" % (src.split("/")[-1], len(gone))
    e8.dump(cert, dst)
    return cert

if __name__ == "__main__":
    a = sys.argv[1:]
    maxn = int(a[2]) if len(a) > 2 else 6
    LIFT = "lift" in a
    for w in a:
        if w.startswith("seed="): SEED = int(w[5:])
        if w.startswith("mode="): MODE = w[5:]
    cert = main(a[0], a[1], maxn, do_totals="nototals" not in a)
    r = e8.quick(cert)
    print({k: r.get(k) for k in ("accepted_ref", "accepted_mirror", "ref_msg", "mirror_msg", "R", "N", "cst", "whole_block", "per_rank", "adds", "most_regs_one_gate")})
    print("ref stats den/maxabs:", r.get("ref_stats", {}).get("den"), r.get("ref_stats", {}).get("maxabs"))
    P = e8.profile(cert)
    print("x patterns", sorted(P["x"].items(), key=lambda kv: -kv[1])[:10])
