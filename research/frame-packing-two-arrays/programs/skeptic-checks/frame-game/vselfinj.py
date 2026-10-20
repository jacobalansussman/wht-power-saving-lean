"""vselfinj.py (v-fp-build): how often are fp-build's repairs null? Generate designs with fp-build's own design.py
(their search, run unchanged), then classify every injection with my own vfg.py arithmetic:
  rescale  = the injected content is supported inside the register's current content and proportional to it
  partial  = supported inside the current content but not proportional (a reweighting, e.g. of the trunk mix)
  genuine  = brings at least one walker the register does not already hold."""
import sys, random, time, collections
import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)); sys.path.insert(0, _os.path.join(_HERE, "..", "..", "frame-game"))
import fpa, design as D
from vfg import FG, cadd
m, trials, depth, seed = map(int, sys.argv[1:5]); dup = "dup" in sys.argv[5:]
S = fpa.Sp(m); V = FG(m); rng = random.Random(seed)
tot = collections.Counter(); per_design = []
t0 = time.time()
for tr in range(trials):
    P = D.Prog(S, m, rng, dup=(rng.randrange(m) if dup else 0))
    if not P.ok: continue
    for step in range(depth + 1):
        res = P.evaluate(); orph = P.orphan_list(res)
        if step == depth or not orph or [r[0] for r in res] == ["yes", "yes"]: break
        if P.try_repair(rng.choice(orph)) is None: break
    # replay with my arithmetic and classify
    p2 = V.from_rows(P.p2)
    inj = [(rg, sl, src, V.from_rows(mu), float(co)) for (rg, sl, src, mu, co) in P.inj]
    states = [{(0, V.I): 1.0}, {(1, V.I): 1.0}]; out = {}; cls = []
    for reg in P.order:
        if reg == "R1": cur = dict(states[0])
        elif reg == "R2": cur = dict(states[1])
        elif reg == "Rt":
            c1 = states[out[("R1", P.a1 - 1)]] if P.a1 else states[0]
            c2 = states[out[("R2", P.a2 - 1)]] if P.a2 else states[1]
            cur = cadd((1, c1), (1, V.lmul(p2, c2)))
        else:
            cur = dict(states[out[("Rt", P.tlen - 1)]])
        for k, u in enumerate(P.regs[reg]):
            for (rg, sl, src, mu, co) in inj:
                if rg == reg and sl == k:
                    add = V.lmul(mu, states[src])
                    ka = {kk for kk, v in add.items() if v != 0}; kc = {kk for kk, v in cur.items() if v != 0}
                    if ka <= kc:
                        ratios = {round(add[kk] / cur[kk], 9) for kk in ka}
                        cls.append("rescale" if (len(ratios) == 1 and ka == kc) else "partial")
                    else:
                        cls.append("genuine")
                    cur = cadd((1, cur), (co, add))
            cur = V.lmul(V.tv(u), cur)
            states.append(cur); out[(reg, k)] = len(states) - 1
    c = collections.Counter(cls); tot.update(c); per_design.append(c)
n = len(per_design)
print(f"m={m} seed={seed} dup={dup} depth={depth}: {n} designs in {time.time()-t0:.0f}s")
print("  injections by class:", dict(tot), " total", sum(tot.values()))
print("  designs with >= 1 injection:", sum(1 for c in per_design if sum(c.values())))
print("  genuine landings per design: histogram", sorted(collections.Counter(c["genuine"] for c in per_design).items()))
print("  null (rescale) injections per design: histogram", sorted(collections.Counter(c["rescale"] for c in per_design).items()))
