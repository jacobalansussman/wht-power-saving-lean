# check_E.py (rl-vprove-compute). Checks, for every n in [n0, Nmax], the hypothesis (E) of rl-prove Lemma A in the
# form used for Theorem 1G (Corollary A3): psi = max(1, ln)^kappa, beta = eta_m/2,
#   LHS = sum_e x(e) e^tau psi(e) + psi(n) sum_s w(s)(n-s)^tau   <=   RHS = psi(n)(n^tau - beta Phi_hat),
# on the ACTUAL flows (my merged computation, cost_fit.Rule). Also: Lemma E unit bound (units/S <= Phi_hat/(W kmin^tau
# (1-kmin^-tau))), the exit bound ln e <= theta ln n, and the psi = 1 version (A4 predicts it FAILS everywhere).
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import log, ceil, exp
from cost_fit import Rule, RULES

def main(which, Nmax, step=1):
    name, Wt, rl = RULES[which]; R = Rule(Wt, rl); tau = R.tau; m = R.m; kmin = min(R.K); rmax = max(R.Pr)
    theta = log(rmax) / log(kmin)
    eta = min((1 - (1 - t) ** tau) / t ** tau for t in [1 / m + (1 - 1 / m) * i / 50000 for i in range(50001)])
    kappa = ceil(log(2 / eta) / log(1 / theta)); beta = eta / 2
    psi = lambda x: max(1.0, log(x)) ** kappa if x > 0 else 1.0
    n0 = max(m, ceil(exp(1 / theta)))
    fails = 0; worst = -1e9; plain_ok = 0; unitbad = 0; exitbad = 0; checked = 0; worst_n = None
    for n in range(n0, Nmax + 1, step):
        f = R.flow(n); S = sum(f['src'].values())
        phat = sum(x * s ** tau for s, x in f['src'].items()) / S
        lhs = sum(x / S * e ** tau * psi(e) for e, x in f['ex'].items()) + psi(n) * sum(x / S * (n - s) ** tau for s, x in f['src'].items())
        rhs = psi(n) * (n ** tau - beta * phat)
        r = (lhs - rhs) / (psi(n) * n ** tau)
        if r > worst: worst, worst_n = r, n
        if lhs > rhs * (1 + 1e-12): fails += 1
        # psi = 1 version (A4: must fail since exits carry all the Phi)
        lhs1 = sum(x / S * e ** tau for e, x in f['ex'].items()) + sum(x / S * (n - s) ** tau for s, x in f['src'].items())
        if lhs1 <= n ** tau - 1e-9 * phat: plain_ok += 1
        if f['units'] / S > phat / (R.W * kmin ** tau * (1 - kmin ** -tau)) * (1 + 1e-12): unitbad += 1
        if max(log(e) for e in f['ex']) > theta * log(n) + 1e-9: exitbad += 1
        checked += 1
    print(f"{name}: tau={tau:.6f} theta={theta:.4f} eta_m={eta:.4f} kappa={kappa} beta={beta:.4f} n0={n0}; checked {checked} n in [{n0},{Nmax}] step {step}: "
          f"(E) with psi=ln^kappa FAILS at {fails} n (max (LHS-RHS)/(psi n^tau) = {worst:.3e} at n={worst_n}); "
          f"(E) with psi=1 holds at {plain_ok} n (A4 predicts 0); Lemma E unit bound violated at {unitbad} n; exit bound violated at {exitbad} n", flush=True)

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 1)
