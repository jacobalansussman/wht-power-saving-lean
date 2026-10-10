#!/usr/bin/env python3
"""reuse2.py: an E8 design (a saved hierarchy and the list of dropped target stops) with SLOT REUSE: a copy of a sum is
written on top of a helper that is already dead ("donor"), instead of into a new helper.  The donor's stale content
then travels with the copy into every sum built from it; every target that ends up with some of it reads the donor
ONCE at the moment the donor dies, with the opposite sign.  Price: one stop of each such target at the donor's frame.

CREDIT.  Both ideas come from outside circuits (hub = CrocSwap/integer-mult-bounds):
  * slot reuse with a compensating read (the "birth-cut"): jamesyc, hub #124, and eumemic, hub #143;
  * the late pairing, that is, which dead helper a new value takes, chosen late: Chafik Boukhalfa (account chafreaky),
    hub #200 and #233 and pull request #2 of this repository, where the choice is a maximum-weight one.  Here it is the
    greedy rule of make_plan() below (the donor that died latest and still fits scores best); his program, his pairs
    and his certificate are not used.

usage: reuse2.py <out name> [maxk=<n>] [minbenefit=<x>] [gen] [join] [passes=<n>] [seeds=<n>] [tlim=<s>] [prune]
Reads hierarchy.json and drop.json in the working directory, writes <out name>.json.gz and <out name>.plan.json there.
In-process checks only (e8.quick).  NOTE: the loop over trial plans stops after tlim seconds of wall-clock time
(default 400; the recorded run needs about 20), so a much slower machine could end with another plan."""
import sys, time, math, json, os
from collections import Counter
sys.dont_write_bytecode = True
from lib import *
import comp2, hio
OUTD = os.getcwd() + "/"
outname = sys.argv[1]
opt = dict(a.split("=") for a in sys.argv[2:] if "=" in a)
flags = set(a for a in sys.argv[2:] if "=" not in a)
MAXK = int(opt.get("maxk", 4)); MINB = float(opt.get("minbenefit", 0.5)); GEN = "gen" in flags; JOIN = "join" in flags; EARLY = "early" in flags; PASSES = int(opt.get("passes", 1))
t0 = time.time()
gd = hio.load_h(OUTD + "hierarchy.json")
class H: pass
h = H(); h.content, h.frame, h.kids, h.need, h.twin = gd["content"], gd["frame"], gd["kids"], gd["need"], gd["twin"]
drop = set(map(tuple, json.load(open(OUTD + "drop.json"))))
f = lambda r: r * math.log(45.0 / r) if r > 0 else 0.0

def weight(c):
    inv = e8.inv_hist(c)
    return 5 * sum(n * r * math.log(45.0 / r) for r, n in inv.items())

def build(plan, name, trace=False):
    st = {}
    d = comp2.compile_design(h, name=name, ylow="all", stats=st, drop=drop, plan=plan, trace=trace)
    return e8.build_cert(d), st

c0, st0 = build({}, "trace", trace=True)
w0 = weight(c0)
deaths, recips, ylf = st0["deaths"], st0["recips"], st0["ylow_frames"]
aff = {r["key"]: {S: k for S, k in st0["kappa"].get(("rec", r["key"]), {}).items() if k} for r in recips}
den = {r["key"]: st0["maxden"].get(("rec", r["key"]), 1) for r in recips}
print("baseline: R %d weight %.1f; dead helpers %d %s; copies %d %s" % (st0["R"], w0, len(deaths), sorted(Counter(len(x["F"]) for x in deaths).items()),
      len(recips), sorted(Counter(len(x["F"]) for x in recips).items())))
print("targets touched by a copy (how many copies): %s;  largest denominator of the stale part: %s" % (
    sorted(Counter(len(aff[r["key"]]) for r in recips).items()), sorted(Counter(den.values()).items())))
sys.stdout.flush()

def nested(a, b):
    return a == b or (len(a) < len(b) and e8.inside(a, b))

def chain_ok(frames, F):
    for G in frames:
        if G == F: continue
        if len(G) == len(F): return False
        if not (e8.inside(G, F) if len(G) < len(F) else e8.inside(F, G)): return False
    return True

def stop_cost(frames, F):
    """weight of one more stop of a target at F, given its present low stops"""
    if F in frames: return 0.0
    dims = sorted([0] + [len(G) for G in frames] + [8])
    dq = len(F)
    a = max(x for x in dims if x < dq); b = min(x for x in dims if x > dq)
    return f(dq - a) + f(b - dq) - f(b - a)

import random
by_target = {}
def make_plan(maxk, minb, order="k", seed=None, fut=0.0):
    rnd = random.Random(seed)
    if not by_target:
        for r in recips:
            for S in aff[r["key"]]: by_target.setdefault(S, []).append(r)
    ych = {S: [tuple(G) for G in ylf[S]] for S in range(V)}
    used = set(); plan = {}; ben = {}
    key = (lambda r: (len(aff[r["key"]]), -len(r["F"]))) if order == "k" else (lambda r: (-len(r["F"]), len(aff[r["key"]])))
    rl = sorted(recips, key=key)
    if seed is not None:
        rl = sorted(recips, key=lambda r: (len(aff[r["key"]]), rnd.random()))
    for r in rl * PASSES:
        if r["key"] in plan: continue
        A = aff[r["key"]]
        if len(A) > maxk or den[r["key"]] > 2 or any(Fr(k).denominator > 2 for k in A.values()): continue
        best = None
        cands_G = set([r["F"]]) if GEN else set()
        if GEN:
            for S in A:
                for G in ych[S]:
                    if nested(G, r["F"]): cands_G.add(G)
        for x in deaths:
            if x["key"] in used or abs(x["scale"]) != 1 or x["pos"] > r["pos"] or not nested(x["F"], r["F"]): continue
            if x["key"] == r["key"]: continue
            dq, dF = len(x["F"]), len(r["F"])
            gl = [x["F"]] + [G for G in cands_G if G != x["F"] and nested(x["F"], G)]
            if JOIN:
                js = set()
                for S in A:
                    for Gs in ych[S]:
                        J = tuple(e8.join(x["F"], Gs))
                        if nested(J, r["F"]): js.add(J)
                tops = [max((Gs for Gs in ych[S] if len(Gs) < dF), key=len, default=None) for S in A]
                J = tuple(e8.join(x["F"], *[t for t in tops if t is not None]))
                if nested(J, r["F"]): js.add(J)
                gl += [J for J in js if J not in gl]
            for G in gl:
                if not all(chain_ok(ych[S], G) for S in A): continue
                g = len(G)
                b = f(9 - dq) + f(dF) - f(dF - dq) - sum(stop_cost(ych[S], G) for S in A) - (f(g - dq) + f(dF - g) - f(dF - dq))
                if b <= minb: continue
                b2 = b
                if fut:
                    b2 += fut * sum(1 for S in A for r2 in by_target[S] if r2["key"] != r["key"] and r2["key"] not in plan and len(aff[r2["key"]]) <= maxk and nested(G, r2["F"]))
                if seed is not None: b2 += rnd.random() * 2.0
                if EARLY: b2 = -len(x["F"]) * 100 + b * 0.01   # control: the donor that died LOWEST (the opposite of the late choice)
                if best is None or b2 > best[0]: best = (b2, x, G)
        if best is not None:
            b, x, G = best
            plan[r["key"]] = (x["key"], G); used.add(x["key"]); ben[r["key"]] = b
            for S in A:
                if G not in ych[S]: ych[S].append(G)
    return plan, ben

def accepted(plan):
    c, st = build(plan, "reuse2.py: hierarchy + %d helpers reused with birth-cut compensation" % len(plan))
    q = e8.quick(c, rate=False)
    return (q["accepted_ref"] and q["accepted_mirror"]), c, st, q

def sift(plan):
    """largest sub-plan the two checkers accept, by splitting (the reuses are nearly independent)"""
    ok, c, st, q = accepted(plan)
    if ok or len(plan) <= 1:
        if not ok: print("   refused single reuse:", list(plan.items()), q["ref_msg"][:80], "|", q["mirror_msg"][:80])
        return plan if ok else {}
    ks = sorted(plan); a = {k: plan[k] for k in ks[:len(ks) // 2]}; b = {k: plan[k] for k in ks[len(ks) // 2:]}
    pa, pb = sift(a), sift(b)
    both = dict(pa); both.update(pb)
    if accepted(both)[0]: return both
    out = dict(pa)                       # add the reuses of the second half one by one
    for k in sorted(pb):
        t = dict(out); t[k] = pb[k]
        if accepted(t)[0]: out = t
    return out

best = None
trials = [("k", mk, None, ft) for mk in sorted(set([2, MAXK])) for ft in (0.0, 1.0, 3.0, 6.0)]
trials += [("k", MAXK, sd, ft) for sd in range(int(opt.get("seeds", 0))) for ft in (1.0, 3.0)]
for order, maxk, seed, fut in trials:
    if time.time() - t0 > float(opt.get("tlim", 400)): break
    if True:
        plan, ben = make_plan(maxk, MINB, order, seed, fut)
        ok, c, st, q = accepted(plan)
        if not ok:
            print("plan order=%s maxk=%d: %d reuses REFUSED as a whole (%s | %s): sifting" % (order, maxk, len(plan), q["ref_msg"][:60], q["mirror_msg"][:60]))
            plan = sift(plan)
            ok, c, st, q = accepted(plan)
        w = weight(c)
        if seed is None or best is None or w < best[0]:
            print("plan maxk=%d seed=%s fut=%s: reuses %3d accepted %s R %d weight %.1f (%+.1f) y blocks %s  %.0fs" % (maxk, seed, fut, len(plan), ok, st["R"], w, w - w0, c["blocks"]["y"], time.time() - t0))
        sys.stdout.flush()
        if ok and (best is None or w < best[0]):
            best = (w, plan, c, st)
w, plan, c, st = best
if "prune" in flags:
    for rk in sorted(plan):
        p2 = dict(plan); p2.pop(rk)
        ok2, c2, st2, q2 = accepted(p2)
        if ok2 and weight(c2) < w - 1e-9:
            plan, w, c, st = p2, weight(c2), c2, st2
    print("after pruning: reuses %d R %d weight %.1f  %.0fs" % (len(plan), st["R"], w, time.time() - t0))
q = e8.quick(c)
print("best:", {k: q[k] for k in ("accepted_ref", "accepted_mirror", "ref_msg", "mirror_msg", "R", "cst", "whole_block")}, "reuses", len(plan))
print("blocks:", c["blocks"])
if q["accepted_ref"] and q["accepted_mirror"]:
    e8.dump(c, OUTD + outname + ".json.gz")
    json.dump([[list(k), [list(x), list(G)]] for k, (x, G) in plan.items()], open(OUTD + outname + ".plan.json", "w"))
    print("wrote", OUTD + outname + ".json.gz", "most registers in one gate", q.get("most_regs_one_gate"), "adds", q.get("adds"))
print("%.0fs" % (time.time() - t0))
