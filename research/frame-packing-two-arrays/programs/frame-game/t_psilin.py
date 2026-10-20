"""t_psilin.py: repaired single-trunk designs, re-run with every injection scalar a SYMBOL; the linearised psi screen
(each state x monomial its own column) decides whether ANY choice of injection scalars could pass psi."""
import sys, random
from fpa import *
import design as D
if __name__ == '__main__':
    m, trials, depth, seed = map(int, sys.argv[1:5])
    dup = len(sys.argv) > 5
    S = Sp(m); rng = random.Random(seed)
    tally = {}
    for tr in range(trials):
        P = D.Prog(S, m, rng, dup=(rng.randrange(m) if dup else 0))
        if not P.ok: continue
        for step in range(depth + 1):
            res = P.evaluate()
            orph = P.orphan_list(res)
            if step == depth or not orph or ("yes", "yes") == tuple(r[0] for r in res): break
            if P.try_repair(rng.choice(orph)) is None: break
        # symbolic re-run
        num = list(P.inj)
        P.inj = [(rg, sl, src, mu, Poly.var(i)) for i, (rg, sl, src, mu, co) in enumerate(num)]
        P.run()
        ps = tuple(psi_lin_screen(S, P.states, l) for l in (0, 1))
        P.inj = num; P.run()
        ps0 = tuple(psi_screen(S, P.states, l) for l in (0, 1))
        k = (ps0, ps, len(num))
        tally[k] = tally.get(k, 0) + 1
    print(f"m={m} trials={trials} depth={depth} seed={seed} dup={dup}: (psi with the drawn scalars, linearised psi, #injections)")
    for k in sorted(tally, key=str): print("  ", k, tally[k])
    

def split_states(states):
    out = []
    for st in states:
        bymono = {}
        for k, v in st.items():
            pv = v if isinstance(v, Poly) else Poly.const(v)
            for mono, co in pv.t.items():
                bymono.setdefault(mono, {})[k] = co
        for mono, d in bymono.items():
            d = {k: v for k, v in d.items() if v != 0}
            if d: out.append(d)
    return out
if __name__ == "__main__" and "--split" in sys.argv:
    pass
