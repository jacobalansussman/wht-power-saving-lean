"""t_bilin.py (fp-build): the TRUE scalar problem for repaired designs whose split relaxation passes.
Unknowns: the injection scalars c (k of them) and the readout coefficients beta_(r,p) on the closure columns of the
split system.  A(c) = sum over monomials mu of mu(c) A_mu (A_mu the columns of the split states of monomial mu,
re-indexed by (state r, p)).  For each c the best beta is least squares; we minimise the relative residual over c
(variable projection, Nelder-Mead from many starts, complex c).  A unit needs residual 0 for both outputs with the
SAME c.  Numerical: a residual bounded away from 0 is evidence, not proof."""
import sys, random, time
import numpy as np
from scipy.optimize import minimize
from fpa import *
import design as D
from t_psilin import split_states

def split_with_tags(states):
    out = []
    for r, st in enumerate(states):
        bymono = {}
        for k, v in st.items():
            pv = v if isinstance(v, Poly) else Poly.const(v)
            for mono, co in pv.t.items():
                bymono.setdefault(mono, {})[k] = co
        for mono, d in bymono.items():
            d = {k: v for k, v in d.items() if v != 0}
            if d: out.append((r, mono, d))
    return out

def closure_system(S, states, target, cap=6000):
    """closure + structural peel on the SYMBOLIC states (Poly coefficients); returns the live system with Poly
    entries: rows, cols, ent[(ri, ci)] = Poly"""
    walk = {}
    for r, st in enumerate(states):
        for (l, g), v in st.items():
            walk.setdefault((l, S.pos(g)), []).append((r, g))
    rows = {}; cols = {}; colent = []; rowprod = []
    queue = []
    def add_row(k):
        if k not in rows:
            rows[k] = len(rows); rowprod.append(set()); queue.append(k)
        return rows[k]
    for k in target: add_row(k)
    while queue:
        if len(cols) > cap: return None
        key = queue.pop(); l, g = key
        for (r, h) in walk.get((l, S.pos(g)), []):
            p = S.mul(g, S.inv(h)); ck = (r, p)
            if ck in cols: continue
            cols[ck] = len(cols); e = {}
            for (l2, g2), v in states[r].items():
                ri = add_row((l2, S.mul(p, g2)))
                pv = v if isinstance(v, Poly) else Poly.const(v)
                e[ri] = e.get(ri, Poly()) + pv
            e = {a: b for a, b in e.items() if b.t}
            colent.append(e)
            for ri in e: rowprod[ri].add(cols[ck])
    tgt = {rows[k] for k in target}
    nr, nc = len(rows), len(cols)
    alive = [True] * nc; cnt = [len(rowprod[i]) for i in range(nr)]
    stack = [i for i in range(nr) if cnt[i] == 1 and i not in tgt]
    while stack:
        i = stack.pop()
        if cnt[i] != 1: continue
        c = next(cc for cc in rowprod[i] if alive[cc]); alive[c] = False
        for ri in colent[c]:
            cnt[ri] -= 1
            if cnt[ri] == 1 and ri not in tgt: stack.append(ri)
    if any(cnt[i] == 0 for i in tgt): return "dead"
    live = [c for c in range(nc) if alive[c]]
    lrows = sorted({ri for c in live for ri in colent[c]} | tgt)
    rix = {r: i for i, r in enumerate(lrows)}
    ent = {}
    for j, c in enumerate(live):
        for ri, pv in colent[c].items():
            for mono, co in pv.t.items():
                ent[(rix[ri], j, mono)] = ent.get((rix[ri], j, mono), 0) + co
    trow = {rix[rows[k]]: v for k, v in target.items()}
    return len(lrows), len(live), ent, trow

def mono_val(mono, c):
    v = 1.0 + 0j
    for i, e in mono: v *= c[i] ** e
    return v

def residual_fn(nr, nc, ent, trow, k):
    b = np.zeros(nr, complex)
    for ri, v in trow.items(): b[ri] = float(v)
    by_mono = {}
    for (ri, ci, mono), v in ent.items():
        by_mono.setdefault(mono, []).append((ri, ci, float(v)))
    mats = []
    for mono, lst in by_mono.items():
        M = np.zeros((nr, nc), complex)
        for ri, ci, v in lst: M[ri, ci] += v
        mats.append((mono, M))
    def res(c):
        A = sum(mono_val(mono, c) * M for mono, M in mats)
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        return float(np.linalg.norm(A @ x - b) / np.linalg.norm(b))
    return res

def solve_design(S, P, rng, starts=10):
    num = list(P.inj); k = len(num)
    P.inj = [(rg, sl, src, mu, Poly.var(i)) for i, (rg, sl, src, mu, co) in enumerate(num)]
    P.run(); states = P.states
    P.inj = num
    fns = []; sizes = []
    for l in (0, 1):
        tgt = {(l, S.K): Fr(1)}
        sysm = closure_system(S, states, tgt)
        if sysm is None: return "cap"
        if sysm == "dead": return "dead"
        sizes.append(sysm[:2])
        fns.append(residual_fn(*sysm, k))
    def total(v):
        c = v[:k] + 1j * v[k:]
        return fns[0](c) + fns[1](c)
    best = None
    for s_ in range(starts):
        v0 = np.array([rng.gauss(0, 1.5) for _ in range(2 * k)])
        r = minimize(total, v0, method="Nelder-Mead", options=dict(maxiter=300 * k, xatol=1e-10, fatol=1e-14))
        if best is None or r.fun < best[0]: best = (r.fun, r.x)
    c = best[1][:k] + 1j * best[1][k:]
    return best[0], fns[0](c), fns[1](c), k, sizes

if __name__ == "__main__":
    m, trials, depth, seed = map(int, sys.argv[1:5])
    dup = len(sys.argv) > 5 and sys.argv[5] == "dup"
    S = Sp(m); rng = random.Random(seed)
    out = []
    t0 = time.time()
    for tr in range(trials):
        P = D.Prog(S, m, rng, dup=(rng.randrange(m) if dup else 0))
        if not P.ok: continue
        for step in range(depth + 1):
            res = P.evaluate()
            orph = P.orphan_list(res)
            if step == depth or not orph or ("yes", "yes") == tuple(r[0] for r in res): break
            if P.try_repair(rng.choice(orph)) is None: break
        if not P.inj: continue
        num = list(P.inj)
        P.inj = [(rg, sl, src, mu, Poly.var(i)) for i, (rg, sl, src, mu, co) in enumerate(num)]
        P.run()
        ok = all(psi_lin_screen(S, P.states, l) for l in (0, 1))
        P.inj = num
        if not ok: continue
        r = solve_design(S, P, rng)
        if isinstance(r, str): print("  ", r, flush=True); continue
        print(f"  design a1={P.a1} a2={P.a2} dup={P.dup} k={r[3]} live systems {r[4]}: best residual sum {r[0]:.3e} (outputs {r[1]:.3e}, {r[2]:.3e})", flush=True)
        out.append(r[0])
    print(f"m={m}: {len(out)} designs solved for scalars in {time.time()-t0:.0f}s; min residual sum "
          f"{min(out) if out else None}", flush=True)
