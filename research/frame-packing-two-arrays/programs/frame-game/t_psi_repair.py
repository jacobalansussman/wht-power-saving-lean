"""t_psi_repair.py: run the repair loop (design.py) and record, after every step, the exact end-solve status and the
psi screen; question: when repairs remove every singleton row, is the remaining failure visible at position level?"""
import sys, random
from fpa import *
import design as D
m, trials, depth, seed = map(int, sys.argv[1:5])
dup = len(sys.argv) > 5
S = Sp(m); rng = random.Random(seed)
tally = {}
for tr in range(trials):
    P = D.Prog(S, m, rng, dup=(rng.randrange(m) if dup else 0))
    if not P.ok: continue
    for step in range(depth + 1):
        res = P.evaluate()
        st = tuple(r[0] if r[0] != "no" else "no:" + r[1]["reason"][:6] for r in res)
        ps = tuple(psi_screen(S, P.states, l) for l in (0, 1))
        orph = P.orphan_list(res)
        if step == depth or not orph or ("yes", "yes") == tuple(r[0] for r in res):
            break
        if P.try_repair(rng.choice(orph)) is None: break
    k = (st, ps, len(P.inj))
    tally[k] = tally.get(k, 0) + 1
for k in sorted(tally, key=str): print(k, tally[k])
