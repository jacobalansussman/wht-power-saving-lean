"""t_hit2.py: continue repairing the one design that passed the bilinear psi screen (hit_m4_s102_4.pkl), now at the
CHART level: find the rows that die in the structural peel, land a walker on one, re-optimise all landing scalars
for psi (bilinear), and re-check the chart level.  Several random continuations.  Numerical psi; structural peel exact."""
import sys, random, pickle, time
import numpy as np
from scipy.optimize import minimize
from fpa import *
import design as D
from t_psibil import psi_parts, make_res
from t_bilin import closure_system

def rebuild(S, d, rng):
    P = object.__new__(D.Prog)
    P.S, P.m, P.rng, P.tlen = S, d["m"], rng, 1
    P.a1, P.a2, P.dup, P.regs, P.order, P.p2 = d["a1"], d["a2"], d["dup"], d["regs"], d["order"], d["p2"]
    P.inj = [(rg, sl, src, mu, Fr(1)) for (rg, sl, src, mu) in d["inj"]]
    P.ok = True
    return P

def psi_opt(S, P, rng, starts=8, maxit=300, c0=None):
    num = list(P.inj); k = len(num)
    P.inj = [(rg, sl, src, mu, Poly.var(i)) for i, (rg, sl, src, mu, co) in enumerate(num)]
    P.run(); sym = P.states; parts = psi_parts(S, sym); P.inj = num
    f = make_res(S, parts, k)
    inits = [np.array([rng.gauss(0, 1.5) for _ in range(2 * k)]) for _ in range(starts)]
    if c0 is not None:
        inits[0] = np.concatenate([np.real(c0), np.imag(c0)])
    best = min((minimize(f, v, method="Nelder-Mead", options=dict(maxiter=maxit * k, xatol=1e-12, fatol=1e-15))
                for v in inits), key=lambda r: r.fun)
    return best.fun, best.x[:k] + 1j * best.x[k:], sym

def chart_status(S, sym):
    out = []
    for l in (0, 1):
        s = closure_system(S, sym, {(l, S.K): Fr(1)}, cap=20000)
        out.append("cap" if s is None else ("dead" if s == "dead" else f"live {s[0]}x{s[1]}"))
    return out

d = pickle.load(open(sys.argv[1], "rb"))
S = Sp(d["m"]); seed = int(sys.argv[2]); conts = int(sys.argv[3]); steps = int(sys.argv[4])
t0 = time.time()
for ci in range(conts):
    rng = random.Random(seed * 1000 + ci)
    P = rebuild(S, d, rng); c = np.array(d["c"])
    hist = []
    for st in range(steps):
        # chart-level singleton rows (exact end solve with generic Fraction scalars)
        P.inj = [(rg, sl, src, mu, Fr(rng.randint(1, 9), rng.randint(1, 9))) for (rg, sl, src, mu, co) in P.inj]
        res = P.evaluate()
        orph = P.orphan_list(res)
        if not orph: hist.append(("no singles", [r[0] for r in res])); break
        if P.try_repair(rng.choice(orph)) is None: hist.append(("no landing",)); break
        cc = np.concatenate([c, [1.0 + 0j]])
        r, c, sym = psi_opt(S, P, rng, c0=cc)
        cs = chart_status(S, sym)
        hist.append((len(P.inj), round(r, 4), cs))
        print(f"  continuation {ci} step {st}: landings {len(P.inj)}, psi residual {r:.3e}, chart {cs}", flush=True)
        if r > 1e-4: break
        if all(x.startswith("live") for x in cs):
            print("  CHART-LEVEL SURVIVOR (psi passes, peel leaves the target alive)", flush=True)
            pickle.dump(dict(d, inj=[(rg, sl, src, mu) for (rg, sl, src, mu, co) in P.inj], c=c),
                        open(f"hit2_{seed}_{ci}_{st}.pkl", "wb"))
print(f"done in {time.time() - t0:.0f}s", flush=True)
