#!/usr/bin/env python3
"""bilin.py (pk2-orbits): EVIDENCE ONLY (numerical).  For a replayed orbit design that passes the linearised psi
screen (psilin) but fails psi at drawn scalars (psinum), search for special orbit scalars c with which the position
screen passes: minimise || sum_k beta_k psi(S_k(c)) - [L1] e_l || over (beta, c) by alternating least squares
(designs whose orbit coefficients are all of degree 1 only, so the problem is bilinear), random restarts.
Reports the smallest residual and the size of beta at it (a residual that only falls as beta blows up is the
degenerate limit fp-build's skeptic found, not a pass).
Usage: bilin.py mode m depth seed trial [kinds] [restarts]"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "heights")); sys.path.insert(0, os.path.join(HERE, "..", "frame-game")); sys.path.insert(0, HERE)
import numpy as np
from fpa import *
import orb, replay

def psi_poly(S, st):
    out = {}
    for (l, g), v in st.items():
        k = (l, S.pos(g)); pv = v if isinstance(v, Poly) else Poly.const(v)
        out[k] = out.get(k, Poly()) + pv
    return out

def main():
    mode, m, depth, seed, trial = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    kinds = sys.argv[6].split(",") if len(sys.argv) > 6 and sys.argv[6] != "all" else ["stab", "cycle", "twof", "perm", "rand"]
    restarts = int(sys.argv[7]) if len(sys.argv) > 7 else 30
    S, P, O, res, st = replay.build(mode, m, depth, seed, trial, kinds)
    W = 3 if mode == "three" else 2
    print(f"design {mode} m={m} seed={seed} trial={trial}: psilin={res['psilin']} psinum={res.get('psinum')} "
          f"exact={res['exact']} orbits={O.log}", flush=True)
    vecs = [psi_poly(S, x) for x in st]
    nv = max(O.nvar, 1)
    rows = sorted({k for v in vecs for k in v} | {(l, S.pos(S.K)) for l in range(W)}, key=str)
    ri = {k: i for i, k in enumerate(rows)}
    R, C = len(rows), len(vecs)
    terms = []                                     # (row, col, coef, ((i, e), ...))
    for j, v in enumerate(vecs):
        for k, p in v.items():
            for mono, co in p.t.items(): terms.append((ri[k], j, float(co), mono))
    def A_and_J(c):
        A = np.zeros((R, C)); dA = np.zeros((nv, R, C))
        for (r_, j, co, mono) in terms:
            val = co
            for (i, e) in mono: val *= c[i] ** e
            A[r_, j] += val
            for (i, e) in mono:
                d = co * e * c[i] ** (e - 1)
                for (i2, e2) in mono:
                    if i2 != i: d *= c[i2] ** e2
                dA[i, r_, j] += d
        return A, dA
    deg = max((sum(e for (_, e) in mono) for (_, _, _, mono) in terms), default=0)
    print(f"   psi system: {R} rows x {C} states, {nv} orbit scalars, max degree {deg}", flush=True)
    rng = np.random.default_rng(seed * 100 + trial)
    T = np.zeros((R, W))
    for l in range(W): T[ri[(l, S.pos(S.K))], l] = 1.0
    def resid_all(c):
        A, _ = A_and_J(c); Bt = np.linalg.lstsq(A, T, rcond=None)[0]
        return np.linalg.norm(A @ Bt - T, axis=0), Bt
    np.seterr(all="ignore")
    for which in [[l] for l in range(W)] + [list(range(W))]:
        best = (9e9, None, None, None)
        for r in range(restarts):
            c = rng.normal(size=nv); A, _ = A_and_J(c); Bt = np.linalg.lstsq(A, T, rcond=None)[0]; lam = 1e-2
            for it in range(150):            # Levenberg-Marquardt on (beta_l for l in which, c), one c for all of them
                A, dA = A_and_J(c)
                res_ = np.concatenate([A @ Bt[:, l] - T[:, l] for l in which])
                blocks = []
                for q, l in enumerate(which):
                    row = [np.zeros((R, C))] * len(which); row = list(row); row[q] = A
                    Jc = np.stack([dA[i] @ Bt[:, l] for i in range(nv)], axis=1)
                    blocks.append(np.concatenate(row + [Jc], axis=1))
                J = np.concatenate(blocks, axis=0)
                if not np.all(np.isfinite(J)): break
                g = J.T @ res_; H = J.T @ J
                try: step = np.linalg.solve(H + lam * np.eye(H.shape[0]), -g)
                except np.linalg.LinAlgError: break
                Bt2 = Bt.copy()
                for q, l in enumerate(which): Bt2[:, l] = Bt[:, l] + step[q * C:(q + 1) * C]
                c2 = c + step[len(which) * C:]
                A2, _ = A_and_J(c2)
                new = np.linalg.norm(np.concatenate([A2 @ Bt2[:, l] - T[:, l] for l in which]))
                if new < np.linalg.norm(res_): Bt, c, lam = Bt2, c2, lam / 3
                else: lam *= 4
            rr, Bt_ = resid_all(c)
            tot = float(np.linalg.norm(rr[which]))
            if np.isfinite(tot) and tot < best[0]:
                best = (tot, float(np.linalg.norm(Bt_[:, which])), float(np.linalg.norm(c)), [float(x) for x in rr[which]])
        print(f"   outputs {which}: smallest psi residual with ONE choice of orbit scalars: {best[0]:.3g} "
              f"(per output {['%.2g' % x for x in (best[3] or [])]}, |beta| {best[1] if best[1] is None else round(best[1], 3)}, "
              f"|c| {best[2] if best[2] is None else '%.3g' % best[2]}); a pass needs 0", flush=True)

if __name__ == "__main__":
    main()
