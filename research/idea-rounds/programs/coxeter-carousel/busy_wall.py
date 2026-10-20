"""Busy-register wall, checked on scalar shadows (track wildcard, wht44/ideas, 16 Oct 2026).

Claim [proved in NOTE.md section 5]: in a time-invariant cell with r registers and p inputs per round, where every
register gets exactly one paid move in every layer (d layers per round), inputs are injected on arrival and outputs
drain the registers, correctness forces r >= p(delta+1) with d(delta+1) >= m, so the rate d*r/p >= m.

The proof applies the trivial representation of the frame group (every group element -> 1), which turns the cell's
feedback-free Markov parameters a_k = alpha^T L (E L)^k beta into scalars; their block Hankel matrix factors through
r states.  This script checks the linear-algebra step on random scalar shadows:
  (1) for random (alpha, L, E, beta) with r states, rank of the Hankel matrix [a_{i+j}] (i, j <= delta) is <= r;
  (2) the required pattern (zeros below lag delta, identity at lag delta) has Hankel rank p(delta+1), so it is
      reachable by least squares only when r >= p(delta+1): we minimise the pattern residual over random starts
      for r = p(delta+1) - 1 and r = p(delta+1).
"""
import numpy as np
from scipy.optimize import least_squares

def markov(alpha, L, E, beta, K):
    M = E @ L
    out, V = [], beta.copy()
    for k in range(K):
        out.append(alpha.T @ L @ V)
        V = M @ V
    return out

def hankel_rank(seq, delta, p):
    H = np.block([[seq[i + j] for j in range(delta + 1)] for i in range(delta + 1)])
    s = np.linalg.svd(H, compute_uv=False)
    return int(np.sum(s > 1e-9 * max(1.0, s[0])))

def pattern_residual(r, p, delta, rng, starts=30):
    best = np.inf
    nvar = r * p + r + r * r + r * p
    def unpack(v):
        i = 0
        alpha = v[i:i + r * p].reshape(r, p); i += r * p
        L = np.diag(v[i:i + r]); i += r
        E = v[i:i + r * r].reshape(r, r); i += r * r
        beta = v[i:i + r * p].reshape(r, p)
        return alpha, L, E, beta
    def resid(v):
        seq = markov(*unpack(v), delta + 1)
        res = [seq[k].ravel() for k in range(delta)] + [(seq[delta] - np.eye(p)).ravel()]
        return np.concatenate(res)
    for _ in range(starts):
        sol = least_squares(resid, rng.standard_normal(nvar), method="trf", max_nfev=3000)
        best = min(best, np.linalg.norm(resid(sol.x)))
    return best

if __name__ == "__main__":
    rng = np.random.default_rng(11)
    out = ["# (1) Hankel rank of random shadows never exceeds r"]
    for r, p, delta in ((2, 1, 4), (3, 1, 5), (3, 2, 3), (5, 2, 4)):
        worst = 0
        for _ in range(200):
            alpha = rng.standard_normal((r, p)); beta = rng.standard_normal((r, p))
            L = np.diag(rng.standard_normal(r)); E = rng.standard_normal((r, r))
            worst = max(worst, hankel_rank(markov(alpha, L, E, beta, 2 * delta + 1), delta, p))
        out.append(f"r={r} p={p} delta={delta}: max Hankel rank over 200 random shadows = {worst} (<= r)")
    out.append("# (2) residual of the target pattern (a_k = 0 for k < delta, a_delta = I_p), best of 30 starts")
    for p, delta in ((1, 2), (1, 3), (2, 1), (2, 2)):
        need = p * (delta + 1)
        for r in (need - 1, need):
            out.append(f"p={p} delta={delta} r={r} (need {need}): residual {pattern_residual(r, p, delta, rng):.2e}")
    txt = "\n".join(out)
    print(txt)
    open("out_busy_wall.txt", "w").write(txt + "\n")
