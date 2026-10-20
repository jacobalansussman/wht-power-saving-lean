#!/usr/bin/env python3
"""growth.py (pk2-orbits): does a CLOSED (orbit) repair stop the junk growth that one-at-a-time landings show?
Same skeleton (fp-build design.Prog, 2 in 2m - 1), two repair rules from the same start:
  (a) fp-build's single landing (design.Prog.try_repair), as in fp-build's runs;
  (b) an orbit repair (orb.Orbit.add_orbit) aimed at the same kind of orphan row, coefficients drawn rationals.
After every repair: exact end solves, fp-build's orphan count (singleton rows near the target in the closure; the
measure behind fp-build's "2, 4, 7, 9, 13, 20, 48, 64"), number of terms over all stored states.
Usage: growth.py m trials depth seed [kinds]
"""
import sys, os, random, copy, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "heights")); sys.path.insert(0, os.path.join(HERE, "..", "frame-game")); sys.path.insert(0, HERE)
from fpa import *
import design as D
import orb

def status(S, P):
    res = P.evaluate()
    st = "".join(r[0][0] if r[0] != "yes" else "Y" for r in res)
    orph = P.orphan_list(res)
    defined = any(r[0] == "no" and r[1].get("singles") for r in res)
    return res, st, (len(orph) if defined else None), sum(len(x) for x in P.states), orph

def numeric_orbits(P):
    # replace symbolic orbit coefficients by drawn rationals (deterministic per index)
    out = []
    for (rg, sl, src, mu, co) in P.inj:
        if isinstance(co, Poly):
            (mono, _), = co.t.items()
            i = mono[0][0]
            r = random.Random(1000 + i)
            co = Fr(r.choice([1, -1, 2, -2, 3, -3, 5])) / r.choice([1, 2, 3, 7])
        out.append((rg, sl, src, mu, co))
    P.inj = out

def main():
    m, trials, depth, seed = map(int, sys.argv[1:5])
    kinds = sys.argv[5].split(",") if len(sys.argv) > 5 else ["stab", "cycle", "twof", "perm", "rand"]
    S = Sp(m); rng = random.Random(seed)
    print(f"growth.py m={m} trials={trials} depth={depth} seed={seed} kinds={kinds}", flush=True)
    agg = {"single": [[] for _ in range(depth + 1)], "orbit": [[] for _ in range(depth + 1)]}
    for tr in range(trials):
        P0 = D.Prog(S, m, rng, dup=0)
        if not P0.ok: continue
        for rule in ("single", "orbit"):
            P = copy.deepcopy(P0); P.rng = random.Random(seed * 1000 + tr)
            O = orb.Orbit(S, P, 2, P.rng, kinds)
            seq = []
            for step in range(depth + 1):
                res, st, no, nt, orph = status(S, P)
                seq.append((st, no, nt))
                if no is not None: agg[rule][step].append(no)
                if "Y" * 2 == st or step == depth or not orph: break
                o = P.rng.choice(orph)
                if rule == "single":
                    ok = P.try_repair(o) is not None
                else:
                    ok = O.add_orbit(o["row"], [x["row"] for x in orph])
                    numeric_orbits(P)
                if not ok: break
            print(f"trial {tr} a=({P0.a1},{P0.a2}) {rule:6s}: " + " ".join(f"{s}/{'-' if o is None else o}/{t}" for (s, o, t) in seq)
                  + (f"  orbits={O.log}" if rule == "orbit" else ""), flush=True)
    print("mean orphan count per repair step (designs still running at that step):")
    for rule in ("single", "orbit"):
        print("  ", rule, [("%.1f (n=%d)" % (sum(v) / len(v), len(v))) if v else "-" for v in agg[rule]], flush=True)

if __name__ == "__main__":
    main()
