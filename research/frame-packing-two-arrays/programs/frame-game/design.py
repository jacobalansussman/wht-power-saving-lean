#!/usr/bin/env python3
"""design.py (fp-build): packing-first repair search for two-array units in the frame game.

Skeleton (2 in 2m - 1, one shared move): y1 prefix register R1 (a1 plain moves), y2 prefix register R2 (a2 moves,
plain directions in random order), trunk Rt := R1 + p2 R2 then ONE move (shared by both lineages), branch A (y1's
remaining m - 1 - a1 plain moves), branch B (y2's remaining m - 1 - a2 moves, a random geodesic completion to K).
Optional extra skeleton: 'dup' = number of B's first moves that copy A's first moves (two branches with a common
segment), which the analysis in NOTE.md section 3 singles out.

Packing pattern = injections: before a move of a register, add coef * mu * (a stored state), mu a free chart chosen
so that one walker of the source LANDS on the position of an orphan row (a junk term with a single producer) at a
chosen later output of that register (the metric condition d(pos source walker, X) = paid length of the segment).
Every stored state is an input or a move output; the end solve uses every stored state with C[P] coefficients.

Usage: design.py m trials depth seed [mode]
"""
import sys, random, time, json
from fpa import *

PSILIN = "psilin" in sys.argv

class Prog:
    def __init__(s, S, m, rng, a1=None, a2=None, dup=0, tlen=1):
        s.S, s.m, s.rng = S, m, rng
        s.tlen = tlen
        s.a1 = rng.randrange(m - tlen + 1) if a1 is None else a1
        s.a2 = rng.randrange(m - tlen + 1) if a2 is None else a2
        s.dup = dup
        s.inj = []          # (reg, slot, src_state, mu, coef)
        s.ok = s._skeleton()
    def _skeleton(s):
        S, m, rng = s.S, s.m, s.rng
        s.order2 = list(range(m)); rng.shuffle(s.order2)
        # registers: name -> list of move vectors u (paid T_u)
        s.regs = {"R1": [(1 << k) for k in range(s.a1)],
                  "R2": [(1 << s.order2[k]) for k in range(s.a2)],
                  "Rt": [(1 << (s.a1 + j)) for j in range(s.tlen)],
                  "A": [(1 << k) for k in range(s.a1 + s.tlen, m)]}
        # y2's element after the trunk, for a random relative chart p2
        g2 = S.I
        for u in s.regs["R2"]: g2 = S.mul(S.transv(u), g2)
        for tries in range(5000):
            p2 = S.rand_P(rng)
            h = S.mul(p2, g2); good = True
            for j in range(s.tlen):
                h = S.mul(S.T(1 << (s.a1 + j)), h)
                if S.level(h) != s.a2 + 1 + j or S.level(S.mul(S.K, S.inv(h))) != m - s.a2 - 1 - j:
                    good = False; break
            if good:
                break
        else:
            return False
        s.p2 = p2
        # branch B: optionally copy A's first `dup` moves, then complete greedily
        wB = []
        for k in range(min(s.dup, len(s.regs["A"]))):
            u = s.regs["A"][k]
            h2 = S.mul(S.transv(u), h)
            if S.level(S.mul(S.K, S.inv(h2))) != S.level(S.mul(S.K, S.inv(h))) - 1: break
            wB.append(u); h = h2
        need = S.level(S.mul(S.K, S.inv(h)))
        allu = [u for u in range(1, 1 << (2 * m)) if (u & ((1 << m) - 1))]
        while need > 0:
            rng.shuffle(allu)
            for u in allu:
                h2 = S.mul(S.transv(u), h)
                if S.level(S.mul(S.K, S.inv(h2))) == need - 1:
                    wB.append(u); h = h2; need -= 1; break
            else:
                return False
        if len(wB) != m - s.tlen - s.a2: return False
        s.regs["B"] = wB
        s.order = ["R1", "R2", "Rt", "B", "A"] if rng.random() < 0.5 else ["R1", "R2", "Rt", "A", "B"]
        return True
    def run(s):
        """returns states, index map (reg, k) -> state index of the output of the k-th move, and the global time
        (number of states existing) at each slot (reg, k)"""
        S = s.S
        states = [inp(S, 0), inp(S, 1)]
        out = {}; avail = {}
        cont = {}
        for reg in s.order:
            if reg == "R1": c = inp(S, 0)
            elif reg == "R2": c = inp(S, 1)
            elif reg == "Rt":
                c1 = states[out[("R1", s.a1 - 1)]] if s.a1 else inp(S, 0)
                c2 = states[out[("R2", s.a2 - 1)]] if s.a2 else inp(S, 1)
                c = cadd((1, c1), (1, lmul(S, s.p2, c2)))
            else:
                c = dict(states[out[("Rt", s.tlen - 1)]])
            for k, u in enumerate(s.regs[reg]):
                avail[(reg, k)] = len(states)
                for (rg, sl, src, mu, co) in s.inj:
                    if rg == reg and sl == k:
                        c = cadd((1, c), (co, lmul(S, mu, states[src])))
                c = lmul(S, S.transv(u), c)
                states.append(c); out[(reg, k)] = len(states) - 1
        s.states, s.out, s.avail = states, out, avail
        return states
    def deliver(s, l):
        return s.out[("A", len(s.regs["A"]) - 1)] if l == 0 else s.out[("B", len(s.regs["B"]) - 1)]
    def evaluate(s):
        S = s.S
        st = s.run()
        res = []
        for l in (0, 1):
            r = end_solve(S, st, {(l, S.K): Fr(1)}, cap_cols=6000, cap_rows=20000)
            res.append(r)
        return res
    def orphan_list(s, res=None):
        S = s.S; st = s.states; outl = []
        if res is not None:
            for r in res:
                if r[0] == "no" and "singles" in r[1]:
                    sing = r[1]["singles"]
                    if sing:
                        dmin = sing[0][0]
                        for (dd, row) in sing:
                            if dd <= dmin + 1:
                                outl.append(dict(row=row, pos=S.pos(row[1]), depth=dd))
            return outl
        for l in (0, 1):
            for o in orphans(S, st, {(l, S.K): Fr(1)}):
                if o["nprod"] == 1: outl.append(o)
        return outl
    def try_repair(s, orph):
        """pick a register slot and a source walker that can land on orph['pos'] at a later output; add it"""
        S, rng = s.S, s.rng
        l, X = orph["row"][0], orph["pos"]
        cands = []
        for reg in s.order:
            n = len(s.regs[reg])
            for k in range(n):
                for t in range(k, n):
                    w = S.I
                    for u in s.regs[reg][k:t + 1]: w = S.mul(S.transv(u), w)
                    lw = S.level(w)
                    for src in range(s.avail[(reg, k)]):
                        for (l2, g), v in s.states[src].items():
                            if l2 != l: continue
                            cands.append((reg, k, t, src, g, w, lw))
        rng.shuffle(cands)
        for (reg, k, t, src, g, w, lw) in cands[:400]:
            # metric precheck via a representative of X
            if s._dist_pos(S.pos(g), X) != lw: continue
            mu = find_chart_c(S, w, g, X)
            if mu is None: continue
            co = Fr(rng.choice([1, -1, 2, -2, 3])) / rng.choice([1, 2, 3])
            s.inj.append((reg, k, src, mu, co))
            return (reg, k, t, src)
        return None
    def _dist_pos(s, X, Y):
        # dim(X cap Y) via rank of the union
        S = s.S
        return S.rank_rows(list(X) + list(Y)) - S.m

def search(m, trials, depth, seed, dup_mode=False, log=None, tlen=1):
    S = Sp(m); rng = random.Random(seed)
    best = None; tally = {}
    t0 = time.time()
    for tr in range(trials):
        dup = rng.randrange(0, m) if dup_mode else 0
        P = Prog(S, m, rng, dup=dup, tlen=tlen)
        if not P.ok: tally["skel-fail"] = tally.get("skel-fail", 0) + 1; continue
        hist = []
        for step in range(depth + 1):
            res = P.evaluate()
            st = [r[0] for r in res]
            orph = P.orphan_list(res)
            hist.append((st, len(orph), sum(len(x) for x in P.states)))
            if st == ["yes", "yes"]:
                print("FRAME-GAME UNIT CANDIDATE", m, P.a1, P.a2, P.inj, flush=True)
                return P, tally
            if not orph or step == depth: break
            o = rng.choice(orph)
            if P.try_repair(o) is None: break
        reasons = tuple((r[0] if r[0] != "no" else ("no:" + r[1].get("reason", "?")[:12])) for r in res)
        if PSILIN:
            num = list(P.inj)
            P.inj = [(rg, sl, src, mu, Poly.var(i)) for i, (rg, sl, src, mu, co) in enumerate(num)]
            P.run()
            pl = tuple(psi_lin_screen(S, P.states, l) for l in (0, 1))
            P.inj = num; P.run()
            reasons = reasons + ("psilin:%d%d" % pl,)
        key = (P.a1, P.a2, P.dup, reasons)
        tally[key] = tally.get(key, 0) + 1
        score = (sum(1 for x in hist[-1][0] if x == "yes"), -hist[-1][1])
        if best is None or score > best[0]:
            best = (score, P.a1, P.a2, P.dup, len(P.inj), hist)
    print(f"m={m}: {trials} designs, depth {depth}, {time.time()-t0:.0f}s; best {best[:5] if best else None}", flush=True)
    if best: print("   best history (statuses, orphans, total terms):", best[5], flush=True)
    for k in sorted(tally, key=str): print("  ", k, tally[k], flush=True)
    return None, tally

if __name__ == "__main__":
    m, trials, depth, seed = map(int, sys.argv[1:5])
    dup_mode = "dup" in sys.argv[5:]
    tlen = 2 if "t2" in sys.argv[5:] else 1
    print(f"design.py m={m} trials={trials} depth={depth} seed={seed} dup={dup_mode} trunk length={tlen} "
          f"(target 2 in {2 * m - tlen})", flush=True)
    search(m, trials, depth, seed, dup_mode, tlen=tlen)
