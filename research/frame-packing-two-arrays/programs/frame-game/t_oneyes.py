"""t_oneyes.py: find designs where one output closes, and re-verify that solution term by term (verify_combo)."""
import random
from fpa import *
import design as D
S = Sp(3); rng = random.Random(52); found = 0
for tr in range(300):
    P = D.Prog(S, 3, rng)
    if not P.ok: continue
    for step in range(9):
        res = P.evaluate()
        ys = [l for l, r in enumerate(res) if r[0] == "yes"]
        if ys:
            for l in ys:
                ok, diff = verify_combo(S, P.states, res[l][1]["solution"], {(l, S.K): Fr(1)})
                print(f"design {tr}: output {l} exact; re-verified term by term: {ok}; other output: "
                      f"{res[1 - l][0]} {res[1 - l][1].get('reason', '')}; landings {len(P.inj)}", flush=True)
            found += 1; break
        orph = P.orphan_list(res)
        if not orph or P.try_repair(rng.choice(orph)) is None: break
    if found >= 4: break
