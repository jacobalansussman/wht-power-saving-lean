"""Coxeter carousel (accumulator) on real arrays (track wildcard, wht44/ideas, 16 Oct 2026).

Unit: r carrier arrays s_1..s_r (scratch, start at 0).  Round t:
  1. inject:  s_i += x(t)                                  (free)
  2. step:    u_i = lam_i * c(s_i),  c = sigma o C_1       (one paid move per carrier)
  3. read:    y(t) = sum_i alpha_i u_i                     (free)
  4. drain:   s_i = u_i - lam_i^m * y(t)                   (free)
With sum_i alpha_i lam_i^j = delta_{j,m} (j = 1..m) the output is y(t) = C^(x)m x(t-m+1) exactly, and every
carrier holds each live input at m different frames of the Petrie polygon: one paid move advances m inputs.
Rate = r paid moves per input.  Calibration: r = m works (plain rate).  The claim tested: r < m is impossible
(the scalar Vandermonde system that the trivial representation forces has no solution).
"""
import numpy as np
from coxeter import apply_dir, sigma_arr

C = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])

def step(arr, m, f):
    return sigma_arr(apply_dir(arr, C, 0, m, f), m)

def target(arr, m, f):
    for j in range(m):
        arr = apply_dir(arr, C, j, m, f)
    return arr

def run(m, f, lam, alpha, T, rng):
    shape = (2 ** f,) * m
    xs = [rng.standard_normal(shape) + 1j * rng.standard_normal(shape) for _ in range(T)]
    r = len(lam)
    s = [np.zeros(shape, complex) for _ in range(r)]
    err, moves, outs = 0.0, 0, 0
    for t in range(T + m - 1):
        x = xs[t] if t < T else np.zeros(shape, complex)
        for i in range(r):
            s[i] = s[i] + x
        u = [lam[i] * step(s[i], m, f) for i in range(r)]
        moves += r
        y = sum(alpha[i] * u[i] for i in range(r))
        k = t - m + 1
        if 0 <= k < T:
            err = max(err, np.max(np.abs(y - target(xs[k], m, f))))
            outs += 1
        s = [u[i] - lam[i] ** m * y for i in range(r)]
    return err, moves, outs

def vandermonde_solve(lam, m):
    V = np.array([[l ** j for l in lam] for j in range(1, m + 1)], dtype=complex)
    e = np.zeros(m, complex)
    e[-1] = 1
    sol, res, rk, _ = np.linalg.lstsq(V, e, rcond=None)
    return sol, np.linalg.norm(V @ sol - e), rk

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    out = []
    for m in (2, 3, 4):
        for f in (1, 2):
            lam = list(np.arange(1, m + 1) * 1.0)
            alpha, resid, rk = vandermonde_solve(lam, m)
            err, moves, outs = run(m, f, lam, alpha, T=10, rng=rng)
            out.append(f"m={m} f={f} r={m}: Vandermonde residual {resid:.1e}; stream of 10 inputs: "
                       f"{outs} outputs, max error {err:.1e}, paid moves {moves} = {moves/10:.1f} per input "
                       f"(r per input in the limit)")
    out.append("# r = m - 1: least-squares residual of the age conditions over random distinct lam (min of 2000)")
    for m in (2, 3, 4, 5):
        best = np.inf
        for _ in range(2000):
            lam = list(rng.standard_normal(m - 1) * rng.choice([0.3, 1, 3]) + 1j * rng.standard_normal(m - 1))
            _, resid, rk = vandermonde_solve(lam, m)
            best = min(best, resid)
        out.append(f"m={m} r={m-1}: smallest residual {best:.3f} (exact zero is impossible: an order-(m-1) "
                   f"recurrence with m-1 leading zeros vanishes)")
    txt = "\n".join(out)
    print(txt)
    open("out_carousel.txt", "w").write(txt + "\n")
