"""t_hit.py: reproduce t_psibil's run with the same arguments; for every design whose bilinear psi residual falls
below 1e-4, save the design and the optimal scalars, then test it at the CHART level: closure + peel on the symbolic
states (Poly coefficients), then least squares with the scalars plugged in (both outputs).  Also report which
scalars are near 0 or large (degenerate optima)."""
import sys, random, time, pickle
import numpy as np
from scipy.optimize import minimize
from fpa import *
import design as D
from t_psibil import psi_parts, make_res
from t_bilin import closure_system, mono_val
import scipy.sparse as sps
from scipy.sparse.linalg import lsqr

def sparse_residual(nr, nc, ent, trow, k):
    b = np.zeros(nr, complex)
    for ri, v in trow.items(): b[ri] = float(v)
    by = {}
    for (ri, ci, mono), v in ent.items(): by.setdefault(mono, []).append((ri, ci, float(v)))
    mats = [(mono, sps.csr_matrix((np.array([v for _, _, v in L], complex), ([r for r, _, _ in L], [c for _, c, _ in L])),
                                  shape=(nr, nc))) for mono, L in by.items()]
    def res(c):
        A = None
        for mono, M in mats:
            t = mono_val(mono, c) * M
            A = t if A is None else A + t
        x = lsqr(A, b, atol=1e-14, btol=1e-14, iter_lim=20000)[0]
        return float(np.linalg.norm(A @ x - b) / np.linalg.norm(b))
    return res

m, trials, depth, seed = map(int, sys.argv[1:5])
dup = "dup" in sys.argv[5:]; tlen = 2 if "t2" in sys.argv[5:] else 1
STARTS = next((int(a[6:]) for a in sys.argv[5:] if a.startswith("starts")), 20)
MAXIT = next((int(a[5:]) for a in sys.argv[5:] if a.startswith("maxit")), 500)
STOP = next((int(a[4:]) for a in sys.argv[5:] if a.startswith("stop")), 10 ** 9)
S = Sp(m); rng = random.Random(seed); t0 = time.time()
nprint = 0
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
    symstates = P.states
    P.inj = num
    f = make_res(S, parts, k)
    runs = [minimize(f, np.array([rng.gauss(0, 1.5) for _ in range(2 * k)]), method="Nelder-Mead",
                     options=dict(maxiter=MAXIT * k, xatol=1e-12, fatol=1e-15)) for _ in range(STARTS)]
    best = min(runs, key=lambda r: r.fun)
    nprint += 1
    print(f"  design #{nprint} a1={P.a1} a2={P.a2} dup={P.dup} injections={k}: psi residual {best.fun:.3e}", flush=True)
    if best.fun < 1e-4:
        c = best.x[:k] + 1j * best.x[k:]
        print("    scalars |c_i|:", np.round(np.abs(c), 4).tolist(), flush=True)
        print("    injections (register, slot, source state):", [(rg, sl, src) for (rg, sl, src, mu, co) in num], flush=True)
        pickle.dump(dict(m=m, a1=P.a1, a2=P.a2, dup=P.dup, regs=P.regs, order=P.order, p2=P.p2,
                         inj=[(rg, sl, src, mu) for (rg, sl, src, mu, co) in num], c=c, psi=best.fun),
                    open(f"hit_m{m}_s{seed}_{nprint}.pkl", "wb"))
        # chart level, with these scalars
        out = []
        for l in (0, 1):
            sysm = closure_system(S, symstates, {(l, S.K): Fr(1)}, cap=20000)
            if sysm is None: out.append("cap"); continue
            if sysm == "dead": out.append("dead after peel"); continue
            rf = sparse_residual(*sysm, k)
            out.append(f"live {sysm[0]}x{sysm[1]}, residual {rf(c):.3e}")
        print("    chart level with these scalars:", out, flush=True)
    if nprint >= STOP: break
print(f"done in {time.time() - t0:.0f}s", flush=True)
