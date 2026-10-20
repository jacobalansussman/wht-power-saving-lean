# recur.py (bound-check/rl-prove, 14 Oct 2026).  Cost of the full schedule of PROOF.md Theorem 1G on small rules:
#   T(n) = min( n,  sum_src w(s) T(n - s) + sum_exit w(e) T(e) )      (free passes F = 0; T(1) = 1, T(0) = 0)
# where, for cutoff L = n, w = the dilated lockstep tree's sources normalised to total weight 1 and the exits carry the
# same normalisation (exit jobs are new top calls Do(e); size-1 exits are leaves).  Prints T(n)/n^tau_c and
# T(n)/(n^tau_c ln(n)^kappa) for the kappa of the proof, and (optionally) the fluid-LP optimum of meta/ex-laser
# grouped_lp.py (a lower bound for every fluid schedule, order ignored) for comparison.
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import log
from tree_flow import check, tau_c

def eta_min(tau, m, steps=20000):
    return min((1 - (1 - t) ** tau) / t ** tau for t in [1 / m + (1 - 1 / m) * i / steps for i in range(steps + 1)])

def run(name, Wt, rounds, nmax):
    P = {}
    for r in rounds: P[r] = P.get(r, 0) + 1
    tau = tau_c(Wt, P); m = max(Wt); kmin = min(Wt); rmax = max(P)
    theta = log(rmax) / log(kmin) if rmax > 1 else 0.0
    eta = eta_min(tau, m)
    kappa = 0 if theta == 0 else max(0, -(-log(2 / eta) // log(1 / theta)))
    T = [0.0, 1.0]
    for n in range(2, nmax + 1):
        best = float(n)
        if n >= m + 1:
            ch = check(Wt, P, n, tau)
            S = ch['sources']
            # source distribution by size (exact counts -> floats)
            from tree_flow import Tree
            val = 0.0
            # re-derive per-size source weights (check() only keeps the set): recompute cheaply
            srcw = source_weights(Wt, P, n)
            tot = sum(srcw.values())
            val += sum(w / tot * T[n - s] for s, w in srcw.items())
            val += sum(x / tot * T[e] for e, x in ch['ex_sizes'].items())
            best = min(best, val)
        T.append(best)
    out = [f"{name}: tau_c={tau:.6f}, theta={theta:.4f}, eta={eta:.4f}, kappa={kappa}"]
    for n in (50, 100, 200, 400, 800, 1600, 3200, 6400):
        if n <= nmax:
            out.append(f"   n={n}: T={T[n]:.3f}  T/n^tau={T[n]/n**tau:.4f}  T/(n^tau ln^kappa n)={T[n]/n**tau/max(1,log(n))**kappa:.4f}")
    mx = max(T[n] / n ** tau for n in range(2, nmax + 1))
    out.append(f"   max_n<={nmax} T/n^tau = {mx:.4f}")
    return T, tau, "\n".join(out)

def source_weights(Wt, P, L):
    from tree_flow import Tree
    T = Tree(Wt, P, L); G = T.G
    w = {}
    for b, x in T.demand(0, 1).items():
        s = T.c(b); w[s] = w.get(s, 0) + x
    N = {(0, 1): 1}
    for d in range(0, G):
        for (dd, rho), n in list(N.items()):
            if dd != d: continue
            for r, pr in P.items():
                rc = rho * r
                N[(d + 1, rc)] = N.get((d + 1, rc), 0) + n * pr
                for b, x in T.demand(d + 1, rc).items():
                    if not T.inR(b, rho):
                        s = rc * T.c(b); w[s] = w.get(s, 0) + n * pr * x
    return w

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 800
    rules = [("A: targets {4:1,3:2}, rounds [1,1,1,1,2,2]", {4: 1, 3: 2}, [1, 1, 1, 1, 2, 2]),
             ("B: targets {7:1,4:1}, rounds [1,1,1,1,2,2,2]", {7: 1, 4: 1}, [1, 1, 1, 1, 2, 2, 2]),
             ("C: targets {5:1,4:3}, rounds [1]*6+[2]*2+[3]*2", {5: 1, 4: 3}, [1] * 6 + [2] * 2 + [3] * 2),
             ("D (rank-1 control = Theorem 1): targets {3:1,2:4}, rounds [1]*10", {3: 1, 2: 4}, [1] * 10)]
    for name, Wt, rl in rules:
        T, tau, txt = run(name, Wt, rl, nmax)
        print(txt, flush=True)
