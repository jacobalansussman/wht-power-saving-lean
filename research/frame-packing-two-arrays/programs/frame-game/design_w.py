#!/usr/bin/env python3
"""design_w.py (fp-build): the W-input version of design.py: ONE trunk move packs W lineages.

Skeleton (W in W(m-1) + 1; at W = 3, m = 4 this is 3 in 10, the brief's target): input l has a prefix register
R<l> (a_l moves; input 0 plain e_0.., the others plain directions in random orders), the trunk Rt := sum_l p_l R<l>
(p_0 = 1, p_l random free charts with which lineage l ascends at the trunk move and stays geodesic to K), ONE trunk
move (P = W - 1 on it), then W branches A<l>, each completing lineage l (A0 plain, the others random geodesic
completions).  Repairs exactly as in design.py (landing injections, Lemma 3.2), end solves for all W outputs, and the
linearised psi screen at the end.  Usage: design_w.py W m trials depth seed
"""
import sys, random, time
from fpa import *
import design as D

class ProgW(D.Prog):
    def __init__(s, S, m, W, rng):
        s.S, s.m, s.W, s.rng = S, m, W, rng
        s.a = [rng.randrange(m) for _ in range(W)]
        s.dup = 0; s.tlen = 1
        s.a1, s.a2 = s.a[0], s.a[1]
        s.inj = []
        s.ok = s._skeleton()
    def _skeleton(s):
        S, m, W, rng = s.S, s.m, s.W, s.rng
        s.regs = {}
        orders = []
        for l in range(W):
            o = list(range(m))
            if l: rng.shuffle(o)
            orders.append(o)
            s.regs["R%d" % l] = [(1 << o[k]) for k in range(s.a[l])]
        zt = 1 << orders[0][s.a[0]] if s.a[0] < m else None
        s.regs["Rt"] = [1 << s.a[0]]
        s.p = [S.I]
        hs = []
        g = S.I
        for u in s.regs["R0"]: g = S.mul(S.transv(u), g)
        hs.append(S.mul(S.T(1 << s.a[0]), g))
        for l in range(1, W):
            g = S.I
            for u in s.regs["R%d" % l]: g = S.mul(S.transv(u), g)
            for tries in range(5000):
                p = S.rand_P(rng)
                h = S.mul(S.T(1 << s.a[0]), S.mul(p, g))
                if S.level(h) == s.a[l] + 1 and S.level(S.mul(S.K, S.inv(h))) == m - s.a[l] - 1: break
            else:
                return False
            s.p.append(p); hs.append(h)
        s.regs["A0"] = [(1 << k) for k in range(s.a[0] + 1, m)]
        allu = [u for u in range(1, 1 << (2 * m)) if (u & ((1 << m) - 1))]
        for l in range(1, W):
            h = hs[l]; w = []
            need = S.level(S.mul(S.K, S.inv(h)))
            while need > 0:
                rng.shuffle(allu)
                for u in allu:
                    h2 = S.mul(S.transv(u), h)
                    if S.level(S.mul(S.K, S.inv(h2))) == need - 1:
                        w.append(u); h = h2; need -= 1; break
                else:
                    return False
            if len(w) != m - 1 - s.a[l]: return False
            s.regs["A%d" % l] = w
        br = ["A%d" % l for l in range(W)]; rng.shuffle(br)
        s.order = ["R%d" % l for l in range(W)] + ["Rt"] + br
        return True
    def run(s):
        S = s.S; W = s.W
        states = [inp(S, l) for l in range(W)]
        out = {}; avail = {}
        for reg in s.order:
            if reg.startswith("R") and reg != "Rt":
                c = inp(S, int(reg[1:]))
            elif reg == "Rt":
                c = {}
                for l in range(W):
                    cl = states[out[("R%d" % l, s.a[l] - 1)]] if s.a[l] else inp(S, l)
                    c = cadd((1, c), (1, lmul(S, s.p[l], cl)))
            else:
                c = dict(states[out[("Rt", 0)]])
            for k, u in enumerate(s.regs[reg]):
                avail[(reg, k)] = len(states)
                for (rg, sl, src, mu, co) in s.inj:
                    if rg == reg and sl == k:
                        c = cadd((1, c), (co, lmul(S, mu, states[src])))
                c = lmul(S, S.transv(u), c)
                states.append(c); out[(reg, k)] = len(states) - 1
        s.states, s.out, s.avail = states, out, avail
        return states
    def evaluate(s):
        S = s.S; st = s.run()
        return [end_solve(S, st, {(l, S.K): Fr(1)}, cap_cols=6000, cap_rows=20000) for l in range(s.W)]

if __name__ == "__main__":
    W, m, trials, depth, seed = map(int, sys.argv[1:6])
    S = Sp(m); rng = random.Random(seed); tally = {}; t0 = time.time()
    print(f"design_w.py W={W} m={m} trials={trials} depth={depth} seed={seed} (target {W} in {W * (m - 1) + 1})", flush=True)
    for tr in range(trials):
        P = ProgW(S, m, W, rng)
        if not P.ok: tally["skel-fail"] = tally.get("skel-fail", 0) + 1; continue
        for step in range(depth + 1):
            res = P.evaluate()
            if all(r[0] == "yes" for r in res):
                print("FRAME-GAME UNIT CANDIDATE", W, m, P.a, P.inj, flush=True); break
            orph = P.orphan_list(res)
            if step == depth or not orph: break
            if P.try_repair(rng.choice(orph)) is None: break
        reasons = tuple((r[0] if r[0] != "no" else ("no:" + r[1].get("reason", "?")[:12])) for r in res)
        num = list(P.inj)
        P.inj = [(rg, sl, src, mu, Poly.var(i)) for i, (rg, sl, src, mu, co) in enumerate(num)]
        P.run()
        pl = "".join("1" if psi_lin_screen(S, P.states, l) else "0" for l in range(W))
        P.inj = num
        key = (tuple(P.a), reasons + ("psilin:" + pl,))
        tally[key] = tally.get(key, 0) + 1
    print(f"W={W} m={m}: {trials} designs in {time.time() - t0:.0f}s", flush=True)
    for k in sorted(tally, key=str): print("  ", k, tally[k], flush=True)
