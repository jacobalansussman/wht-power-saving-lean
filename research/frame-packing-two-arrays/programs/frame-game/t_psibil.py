"""t_psibil.py (fp-build): the psi screen with the injection scalars as UNKNOWNS (the true bilinear position-level
problem).  For repaired single-trunk designs whose linearised psi screen passes, minimise over complex scalars c the
relative least-squares residual of  sum_r x_r psi(S_r(c)) = (l, L1)  for both outputs at once (variable projection:
x is solved exactly for each c; Nelder-Mead over c from many starts).  A unit needs residual 0.  If the minimum is 0,
the scalars found are re-checked with the exact chart-level end solve (rationalised).  Numerical evidence only."""
import sys, random, time
import numpy as np
from scipy.optimize import minimize
from fpa import *
import design as D

def psi_parts(S, states):
    """per state: dict monomial -> dict (l, pos) -> coef (Fraction)"""
    out = []
    for st in states:
        bym = {}
        for (l, g), v in st.items():
            pv = v if isinstance(v, Poly) else Poly.const(v)
            key = (l, S.pos(g))
            for mono, co in pv.t.items():
                d = bym.setdefault(mono, {}); d[key] = d.get(key, 0) + co
        out.append(bym)
    return out

def make_res(S, parts, k):
    keys = sorted({kk for bym in parts for d in bym.values() for kk in d} | {(0, S.pos(S.K)), (1, S.pos(S.K))}, key=str)
    ki = {kk: i for i, kk in enumerate(keys)}
    nr, ns = len(keys), len(parts)
    mats = {}
    for r, bym in enumerate(parts):
        for mono, d in bym.items():
            M = mats.setdefault(mono, np.zeros((nr, ns), complex))
            for kk, v in d.items(): M[ki[kk], r] += float(v)
    b = [np.zeros(nr, complex) for _ in (0, 1)]
    for l in (0, 1): b[l][ki[(l, S.pos(S.K))]] = 1
    def A_of(c):
        A = np.zeros((nr, ns), complex)
        for mono, M in mats.items():
            v = 1 + 0j
            for i, e in mono: v *= c[i] ** e
            A += v * M
        return A
    def res(v):
        c = v[:k] + 1j * v[k:]
        A = A_of(c); tot = 0
        for l in (0, 1):
            x, *_ = np.linalg.lstsq(A, b[l], rcond=None)
            tot += float(np.linalg.norm(A @ x - b[l]))
        return tot
    return res

if __name__ == "__main__":
    m, trials, depth, seed = map(int, sys.argv[1:5])
    dup = "dup" in sys.argv[5:]; tlen = 2 if "t2" in sys.argv[5:] else 1
    STARTS = next((int(a[6:]) for a in sys.argv[5:] if a.startswith("starts")), 20)
    MAXIT = next((int(a[5:]) for a in sys.argv[5:] if a.startswith("maxit")), 500)
    S = Sp(m); rng = random.Random(seed); t0 = time.time()
    res_all = []
    for tr in range(trials):
        P = D.Prog(S, m, rng, dup=(rng.randrange(m) if dup else 0), tlen=tlen)
        if not P.ok: continue
        for step in range(depth + 1):
            res = P.evaluate()
            orph = P.orphan_list(res)
            if step == depth or not orph or ("yes", "yes") == tuple(r[0] for r in res): break
            if P.try_repair(rng.choice(orph)) is None: break
        k = len(P.inj)
        if k == 0: continue
        num = list(P.inj)
        P.inj = [(rg, sl, src, mu, Poly.var(i)) for i, (rg, sl, src, mu, co) in enumerate(num)]
        P.run()
        if not all(psi_lin_screen(S, P.states, l) for l in (0, 1)):
            P.inj = num; continue
        parts = psi_parts(S, P.states)
        P.inj = num
        f = make_res(S, parts, k)
        best = min((minimize(f, np.array([rng.gauss(0, 1.5) for _ in range(2 * k)]), method="Nelder-Mead",
                             options=dict(maxiter=MAXIT * k, xatol=1e-12, fatol=1e-15)).fun for _ in range(STARTS)))
        res_all.append(best)
        print(f"  design a1={P.a1} a2={P.a2} dup={P.dup} injections={k}: min psi residual over scalars {best:.3e}", flush=True)
    print(f"m={m} tlen={tlen} dup={dup}: {len(res_all)} designs passed linearised psi; min residuals "
          f"{sorted(round(x, 4) for x in res_all)} ({time.time()-t0:.0f}s)", flush=True)
