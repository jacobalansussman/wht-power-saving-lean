# indep_tree.py (rl-v-logic, skeptic of rl-prove Theorem 1G). Written from PROOF.md section 3 alone, without
# importing or copying rl-prove's tree_flow.py. Exact integer balance of the dilated lockstep tree, merged by
# (depth, rho) with multiplicities, on a given mixed rule; Phi bookkeeping in Decimal; exit/source size bounds.
# Usage: python3 -I -B indep_tree.py
import sys
from decimal import Decimal as D, getcontext
from math import comb, log
getcontext().prec = 50

def taus(Wt, P):
    lo, hi = D(0), D(1)
    pw = lambda x, t: (D(x).ln() * t).exp()
    for _ in range(200):
        t = (lo + hi) / 2
        if sum(n * pw(r, t) for r, n in P.items()) > sum(w * pw(k, t) for k, w in Wt.items()): lo = t
        else: hi = t
    return (lo + hi) / 2

def run(Wt, P, L):
    K = sorted(Wt); m = K[-1]; assert max(P) < K[0]
    tau = taus(Wt, P); pw = lambda x: (D(x).ln() * tau).exp()
    def c(a):
        v = 1
        for k, x in zip(K, a): v *= k ** x
        return v
    def u(a):  # multinomial(|a|; a) prod W_k^a_k
        v, s = 1, 0
        for k, x in zip(K, a): s += x; v *= comb(s, x) * Wt[k] ** x
        return v
    def pts(rho, g):  # a in R(rho) with |a| = g  (2 or 3 target ranks)
        out = []
        def rec(i, a, left):
            if i == len(K) - 1:
                b = a + [left]
                if m * rho * c(b) <= L: out.append(tuple(b))
                return
            for x in range(left + 1): rec(i + 1, a + [x], left - x)
        rec(0, [], g); return out
    G = 0
    while pts(1, G + 1): G += 1
    # root demand = sources
    N = {(0, 1): 1}; bad = 0; src_n = 0; src_phi = D(0); ex_n = 0; ex_phi = D(0); smin = None; smax = 0; emax = 0
    def demand(rho, g):
        need = {}
        for a in pts(rho, g):
            ua = u(a)
            for i, k in enumerate(K):
                b = list(a); b[i] += 1; b = tuple(b); need[b] = need.get(b, 0) + Wt[k] * ua
        return need
    for b, x in demand(1, G).items():
        s = c(b)
        if m * s <= L: bad += 1
        if not (L / m < s <= L): bad += 1
        src_n += x; src_phi += x * pw(s); smin = s if smin is None else min(smin, s); smax = max(smax, s)
    for d in range(G + 1):
        for (dd, rho), mult in list(N.items()):
            if dd != d: continue
            units = {a: u(a) for a in pts(rho, G - d)}
            if d == G:
                for a, ua in units.items():
                    assert sum(a) == 0
                    for r, pr in P.items():
                        e = rho * r; ex_n += mult * pr * ua; ex_phi += mult * pr * ua * pw(e); emax = max(emax, e)
                continue
            for r, pr in P.items():
                rc = rho * r; N[(d + 1, rc)] = N.get((d + 1, rc), 0) + mult * pr
                need = demand(rc, G - d - 1)
                for b in set(need) | set(units):
                    sup, dem = units.get(b, 0), need.get(b, 0)
                    if m * rho * c(b) <= L:            # parent point: must balance exactly
                        if sup != dem: bad += 1
                    else:
                        if sup: bad += 1
                        if dem:
                            s = rc * c(b)
                            if not (r * L / m < s <= L): bad += 1
                            src_n += mult * pr * dem; src_phi += mult * pr * dem * pw(s)
                            smin = min(smin, s); smax = max(smax, s)
    theta = log(max(P)) / log(K[0])
    return dict(L=L, G=G, tau=float(tau), bad=bad, sources=src_n, src_range=(smin, smax), L_over_m=L / m,
                phi_ratio=float(ex_phi / src_phi), exit_max=emax, exit_bound_L_theta=L ** theta,
                classes=sum(N.values()))

if __name__ == '__main__':
    T74 = ({25: 22748, 73: 43556}, {1: 1560432, 2: 336336, 3: 278724, 4: 84912, 5: 54612, 8: 5328})
    toyA = ({3: 2, 4: 1}, {1: 4, 2: 2})
    for name, (Wt, P) in [('toyA', toyA), ('T74-mixed', T74)]:
        for L in ([40, 150, 1000, 10**5] if name == 'toyA' else [10**5, 10**6, 10**8]):
            r = run(Wt, P, L)
            print(name, r, flush=True)
