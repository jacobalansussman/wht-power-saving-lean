"""t_split.py: for repaired designs whose LINEARISED psi screen passes for both outputs, run the exact chart-level end
solve on the SPLIT states (each state x monomial in the injection symbols its own state with its own free coefficient):
a relaxation of 'some choice of injection scalars works'.  If the split solve fails, the design is dead for every
choice of scalars."""
import sys, random
from fpa import *
import design as D
from t_psilin import split_states
m, trials, depth, seed = map(int, sys.argv[1:5])
dup = len(sys.argv) > 5 and sys.argv[5] == "dup"
S = Sp(m); rng = random.Random(seed)
tally = {}; found = []
for tr in range(trials):
    P = D.Prog(S, m, rng, dup=(rng.randrange(m) if dup else 0))
    if not P.ok: continue
    for step in range(depth + 1):
        res = P.evaluate()
        orph = P.orphan_list(res)
        if step == depth or not orph or ("yes", "yes") == tuple(r[0] for r in res): break
        if P.try_repair(rng.choice(orph)) is None: break
    num = list(P.inj)
    P.inj = [(rg, sl, src, mu, Poly.var(i)) for i, (rg, sl, src, mu, co) in enumerate(num)]
    P.run()
    ps = tuple(psi_lin_screen(S, P.states, l) for l in (0, 1))
    if ps != (True, True):
        tally[("psi-lin fails",)] = tally.get(("psi-lin fails",), 0) + 1; continue
    sp = split_states(P.states)
    r = [end_solve(S, sp, {(l, S.K): Fr(1)}, cap_cols=8000, cap_rows=30000) for l in (0, 1)]
    k = ("psi-lin passes", tuple(x[0] if x[0] != "no" else "no:" + x[1]["reason"][:6] for x in r), len(num))
    tally[k] = tally.get(k, 0) + 1
    if all(x[0] == "yes" for x in r): found.append((P.a1, P.a2, P.dup, num))
print(f"m={m} trials={trials} depth={depth} seed={seed} dup={dup}")
for k in sorted(tally, key=str): print("  ", k, tally[k])
print("split relaxation passes for both outputs:", len(found))
