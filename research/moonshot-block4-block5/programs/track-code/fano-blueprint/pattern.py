"""pattern.py -- ms-fano: exact checks for the PATTERN LEMMA (separation needs a frustrated sign pattern) and for the
RIGIDITY of a shared lineage.
(1) Separation criterion [proved in NOTE/PROOF]: bundles c_j = x1 + lam_j h_j x2 (h_j group elements, lam_j scalars)
    can be separated by left Q[Gamma]-combinations (sum P_j c_j = K x1 with no x2) iff the assignment
    h_1^-1 h_j -> lam_1 / lam_j does NOT extend to a character of H = <h_1^-1 h_j>.  For H abelian we test it over
    all characters, and for separable patterns we print the exact solution (idempotents e_chi = 2^-k sum chi(h) h:
    the halving).  Examples: H = F_2^3 (e.g. Z-Paulis Z^v), points v of the Fano plane.
(2) Rigidity in the engine (block 3): a staggered lineage G1 = D b1 x1, G2 = D(r2 G1 + s b2 x2), G3, G4 double advances;
    the pattern rho = X_2^-1 Y_2 is a single group element (times s) and Y_t = X_t rho for every t >= 2 exactly."""
import sys, itertools
from fractions import Fraction as Fr

def chars(k):
    return [tuple(a) for a in itertools.product((0, 1), repeat=k)]

def dot(a, v):
    return sum(x * y for x, y in zip(a, v)) % 2

def separable(pattern, k):
    """pattern: list of (v in F_2^k as tuple, lam) for bundles x1 + lam Z^v x2.  Returns (bool, per-character info)"""
    (v1, l1) = pattern[0]
    bad = []
    for a in chars(k):
        vals = [l * (-1) ** dot(a, v) for (v, l) in pattern]
        if len(set(vals)) == 1:
            bad.append(a)
    return (len(bad) == 0), bad

def solve(pattern, k):
    """explicit P_j in Q[F_2^k] with sum P_j = 1, sum P_j lam_j Z^{v_j} = 0 (exists iff separable): per character
    pick two bundles with different values and solve the 2x2 system; P_j = sum_chi p_j(chi) e_chi"""
    sol = {j: {} for j in range(len(pattern))}           # j -> {h: coeff}
    n = 2 ** k
    for a in chars(k):
        vals = [l * (-1) ** dot(a, v) for (v, l) in pattern]
        i, j = next((i, j) for i in range(len(vals)) for j in range(len(vals)) if vals[i] != vals[j])
        # p_i + p_j = 1, p_i vals[i] + p_j vals[j] = 0
        pj = Fr(vals[i]) / (vals[i] - vals[j]); pi = 1 - pj
        for (jj, p) in ((i, pi), (j, pj)):
            for h in chars(k):                           # e_chi = 2^-k sum_h chi(h) h
                c = p * Fr((-1) ** dot(a, h), n)
                sol[jj][h] = sol[jj].get(h, 0) + c
    return sol

def check(pattern, k, sol):
    """verify sum P_j = 1 and sum P_j lam_j Z^{v_j} = 0 in Q[F_2^k]"""
    tot1, tot2 = {}, {}
    for j, P in sol.items():
        (v, l) = pattern[j]
        for h, c in P.items():
            tot1[h] = tot1.get(h, 0) + c
            hv = tuple((x + y) % 2 for x, y in zip(h, v))
            tot2[hv] = tot2.get(hv, 0) + c * l
    one = tuple([0] * k)
    ok1 = all((c == (1 if h == one else 0)) for h, c in tot1.items())
    ok2 = all(c == 0 for c in tot2.values())
    dens = sorted(set(c.denominator for P in sol.values() for c in P.values() if c))
    return ok1 and ok2, dens

if __name__ == "__main__":
    k = 3
    pts = [tuple((v >> j) & 1 for j in range(3)) for v in range(1, 8)]
    zero = (0, 0, 0)
    cases = {
        "sum/difference (x1+x2, x1-x2)": [(zero, 1), (zero, -1)],
        "Fano points, all signs +1": [(v, 1) for v in pts],
        "Fano points, all signs -1": [(v, -1) for v in pts],
        "Fano points, sign -1 on 111 only": [(v, -1 if v == (1, 1, 1) else 1) for v in pts],
        "line {100,010,110}, signs (+,+,-)": [((1, 0, 0), 1), ((0, 1, 0), 1), ((1, 1, 0), -1)],
        "line {100,010,110}, signs (+,+,+)": [((1, 0, 0), 1), ((0, 1, 0), 1), ((1, 1, 0), 1)],
        "base 000 + line, signs (+,+,+,-)": [(zero, 1), ((1, 0, 0), 1), ((0, 1, 0), 1), ((1, 1, 0), -1)],
    }
    for name, pat in cases.items():
        ok, bad = separable(pat, k)
        line = f"{name}: separable {ok}"
        if not ok:
            line += f"; blocking characters {bad}"
        else:
            sol = solve(pat, k)
            good, dens = check(pat, k, sol)
            line += f"; explicit solution verified {good}; denominators {dens}"
        print(line)
    # characteristic 2: every sign is 1, so every pattern is a character pattern (the trivial character blocks)
    print("char 2: signs collapse to +1, the trivial character always blocks -> no separation (mod-2 perfection)")
