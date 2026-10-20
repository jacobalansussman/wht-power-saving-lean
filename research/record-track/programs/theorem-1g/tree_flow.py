# tree_flow.py (bound-check/rl-prove, 14 Oct 2026).  The DILATED LOCKSTEP TREE of PROOF.md section 3 (Theorem 1G),
# for a mixed rule with GROUPED paid blocks: targets {k: W_k} (k in K), rounds r_1..r_P (P_r rounds of rank r), and
# the separation hypothesis r_max < k_min.  Exact integers / fractions.
#   u(a) = multinomial(|a|; a) prod_k W_k^a_k  (a in N^K),  c(a) = prod_k k^a_k,  m = max K,
#   R(rho) = {a : m rho c(a) <= L},  G = max{|a| : a in R(1)}.
#   Class s (a round sequence of length d <= G, rho_s = product of its ranks): u(a) units at column rho_s c(a) for
#   every a in R(rho_s) with |a| = G - d.  Round j of class s feeds class s+(j) only.
# Checks (merged by (d, rho) with path multiplicities N(d, rho)):
#   (1) exact balance: at every point b in R(rho') the child's demand equals the parent's round output;
#   (2) every other demand is a SOURCE of size in (r_j L/m, L] (root: (L/m, L]);
#   (3) exits only at generation 0: P jobs of sizes rho*r_j per class, all <= r_max * (L/m)^theta, theta = ln r_max/ln k_min;
#   (4) Phi conservation at tau_c: Phi(sources) = Phi(exits);  (5) unit count <= Phi(sources)/(W k_min^tau (1-k_min^-tau)).
# Literal lockstep simulation (lockstep_literal): explicit recursion over round sequences, a pool of issued jobs per
# class, every class must find all its inputs issued (parent round j + sources) and consume every job it is given.
import sys
from fractions import Fraction as Fr
from math import factorial, log, prod
from itertools import product as iprod

def tau_c(Wt, P):
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if sum(w * k ** mid for k, w in Wt.items()) > sum(p * r ** mid for r, p in P.items()): hi = mid
        else: lo = mid
    return (lo + hi) / 2

def lattice(Wt, L, rho):
    """points a in R(rho), as tuples over sorted K, grouped by |a|"""
    K = sorted(Wt); m = max(K)
    pts = {}
    def rec(i, a, c):
        if i == len(K):
            pts.setdefault(sum(a), []).append(tuple(a)); return
        k = K[i]; e = 0; cc = c
        while m * rho * cc <= L:
            rec(i + 1, a + [e], cc); e += 1; cc *= k
    rec(0, [], 1)
    return pts

def mult(a):
    return factorial(sum(a)) // prod(factorial(x) for x in a)

class Tree:
    def __init__(self, Wt, P, L):
        self.Wt, self.P, self.L = dict(Wt), dict(P), L
        self.K = sorted(Wt); self.m = max(self.K)
        assert max(P) < min(self.K), "separation r_max < k_min fails"
        self.Wv = [Wt[k] for k in self.K]
        self.lat = {}
        self.G = max(lattice(Wt, L, 1))
    def R(self, rho):
        if rho not in self.lat: self.lat[rho] = lattice(self.Wt, self.L, rho)
        return self.lat[rho]
    def u(self, a): return mult(a) * prod(w ** x for w, x in zip(self.Wv, a))
    def c(self, a): return prod(k ** x for k, x in zip(self.K, a))
    def inR(self, a, rho): return self.m * rho * self.c(a) <= self.L
    def units(self, d, rho):
        return {a: self.u(a) for a in self.R(rho).get(self.G - d, [])}
    def demand(self, d, rho):
        need = {}
        for a, ua in self.units(d, rho).items():
            for i, k in enumerate(self.K):
                b = list(a); b[i] += 1; b = tuple(b)
                need[b] = need.get(b, 0) + self.Wt[k] * ua
        return need

def check(Wt, P, L, tau=None):
    T = Tree(Wt, P, L); tau = tau_c(Wt, P) if tau is None else tau
    G, m, kmin, rmax = T.G, T.m, min(T.K), max(P)
    theta = log(rmax) / log(kmin) if rmax > 1 else 0.0
    N = {(0, 1): 1}
    src_w, src_phi, ex_w, ex_phi, units, bad = 0, 0.0, 0, 0.0, 0, []
    src_sizes, ex_sizes = set(), {}
    # root
    for b, x in T.demand(0, 1).items():
        s = T.c(b)
        if T.inR(b, 1): bad.append(('root demand inside R(1)', b))
        if not (L / m < s <= L): bad.append(('root source size', s))
        src_w += x; src_phi += x * s ** tau; src_sizes.add(s)
    for d in range(0, G + 1):
        layer = {k: v for k, v in N.items() if k[0] == d}
        for (_, rho), n in layer.items():
            us = T.units(d, rho)
            units += n * sum(us.values())
            if d == G:                                    # generation 0: exits
                for a, ua in us.items():
                    assert a == tuple(0 for _ in T.K)
                    for r, pr in P.items():
                        e = rho * r
                        ex_w += n * pr * ua; ex_phi += n * pr * ua * e ** tau
                        ex_sizes[e] = ex_sizes.get(e, 0) + n * pr * ua
                continue
            for r, pr in P.items():                       # children via a round of rank r
                rc = rho * r
                N[(d + 1, rc)] = N.get((d + 1, rc), 0) + n * pr
                need = T.demand(d + 1, rc)
                supply = {a: ua for a, ua in us.items()}  # parent round output at point a, size rc*c(a)
                for b in set(need) | set(supply):
                    sup = supply.get(b, 0); dem = need.get(b, 0)
                    if T.inR(b, rho):
                        if sup != dem: bad.append(('unbalanced', d, rho, r, b, sup, dem))
                    else:
                        if sup: bad.append(('supply outside R', b))
                        s = rc * T.c(b)
                        if dem and not (r * L / m < s <= L): bad.append(('source size', s, r))
                        src_w += n * pr * dem; src_phi += n * pr * dem * s ** tau
                        if dem: src_sizes.add(s)
    emax = max(ex_sizes)
    bound = rmax * (L / m) ** theta if rmax > 1 else 1
    ubound = src_phi / (sum(Wt.values()) * kmin ** tau * (1 - kmin ** -tau))
    return dict(L=L, G=G, tau=tau, unbalanced_or_bad=len(bad), bad=bad[:3], sources=src_w, Phi_src=src_phi,
                exits=ex_w, Phi_exit=ex_phi, Phi_ratio=ex_phi / src_phi, src_size_range=(min(src_sizes), max(src_sizes)),
                exit_max=emax, exit_bound=bound, units=units, unit_bound=ubound, classes=sum(N.values()),
                ex_sizes=ex_sizes, src_sizes=sorted(src_sizes))

def lockstep_literal(Wt, P_list, L):
    """Explicit recursion over round sequences (P_list = ranks of rounds 1..P in word order)."""
    Pd = {}
    for r in P_list: Pd[r] = Pd.get(r, 0) + 1
    T = Tree(Wt, Pd, L); G = T.G
    st = dict(classes=0, src=0, exits={}, issued_before_start=True)
    def process(path, rho, fluid):
        d = len(path); st['classes'] += 1
        us = T.units(d, rho)
        need = T.demand(d, rho)
        # fluid = parent's round output (points in R(parent rho)); top it up with sources, then it must equal need
        for b, x in need.items():
            have = fluid.get(b, 0)
            if have > x: raise AssertionError(('surplus', path, b, have, x))
            if have < x: st['src'] += x - have          # source from the top split (always available)
        for b in fluid:
            if b not in need and fluid[b]: raise AssertionError(('unconsumed', path, b))
        for j, r in enumerate(P_list):                  # rounds in word order; round j+1 only after round j completes
            out = dict(us)                              # one job of rank r per unit, at point a, size rho*r*c(a)
            if d == G:
                for a, ua in out.items():
                    e = rho * r; st['exits'][e] = st['exits'].get(e, 0) + ua   # handled by a new top call Do(e)
            else:
                process(path + (j,), rho * r, out)
    process((), 1, {})
    return st, T

if __name__ == "__main__":
    rules = [("A: targets {4:1,3:2}, rounds [1,1,1,1,2,2]", {4: 1, 3: 2}, [1, 1, 1, 1, 2, 2]),
             ("B: pure {7:2}+idle 3 -> targets {7:1,4:1}, rounds [1,1,1,1,2,2,2]", {7: 1, 4: 1}, [1, 1, 1, 1, 2, 2, 2]),
             ("C: targets {5:1,4:3}, rounds [1]*6+[2]*2+[3]*2", {5: 1, 4: 3}, [1] * 6 + [2] * 2 + [3] * 2)]
    print("(A) literal lockstep (explicit round sequences):", flush=True)
    for name, Wt, rl in rules:
        for L in (40, 150):
            Pd = {}
            for r in rl: Pd[r] = Pd.get(r, 0) + 1
            ref = check(Wt, Pd, L)
            st, T = lockstep_literal(Wt, rl, L)
            ok = (st['src'] == ref['sources'] and st['exits'] == ref['ex_sizes'])
            print(f"  {name}, L={L}: G={T.G}, classes {st['classes']}, every class exactly fed, sources {st['src']} "
                  f"(= merged count: {ok}), exits by size {dict(sorted(st['exits'].items()))}", flush=True)
    print("(B) merged checks at larger L:", flush=True)
    for name, Wt, rl in rules:
        Pd = {}
        for r in rl: Pd[r] = Pd.get(r, 0) + 1
        for L in (1000, 10 ** 5, 10 ** 7):
            ch = check(Wt, Pd, L)
            print(f"  {name}, L={L}: tau_c={ch['tau']:.6f} G={ch['G']} classes={ch['classes']} bad={ch['unbalanced_or_bad']} "
                  f"sources in {ch['src_size_range']} (L/m={L/max(Wt):.1f}); Phi_exit/Phi_src={ch['Phi_ratio']:.12f}; "
                  f"max exit {ch['exit_max']} <= bound {ch['exit_bound']:.1f}; units {float(ch['units']):.4g} <= "
                  f"bound {ch['unit_bound']:.4g}", flush=True)
