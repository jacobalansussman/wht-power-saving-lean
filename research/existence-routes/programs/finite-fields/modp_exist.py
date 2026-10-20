"""modp_exist.py -- ex2-modp: the existence pipeline "find mod p, lift p-adically, reconstruct, verify exactly"
(PROOF.md, Theorem 1 + Proposition 4) for systems linear in the beta variables (names b*, c*), alphas (a*) shared.

  python3 modp_exist.py PREFIX_ARRAY0 [PREFIX_ARRAY1 ...] [--tries N] [--box B] [--seed S]

Each PREFIX has PREFIX.p1.ms and PREFIX.p2.ms (nag-rerun exports of one array's full-system core mod two primes).
Steps:
 0. rebuild the system over Q (rational reconstruction mod p1, checked mod p2); join the arrays (alphas shared).
 1. search an F_p1 point: alpha = small integers in [-B, B] (so the lift of alpha is exact), beta from the linear
    system mod p1.
 2. Hensel-lift beta with alpha fixed (Prop. 4: unique when J restricted to beta is injective; free beta coords
    fixed to 0 otherwise), to precision p1^8.
 3. rational reconstruction of beta; 4. EXACT verification of every equation of the joint system over Q.
Only step 4 is the proof; steps 1-3 are the search.  Prints one verdict line.
"""
import sys, random, argparse
from fractions import Fraction
sys.path.insert(0, __import__('os').path.dirname(__file__))
from zsys import *

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('prefixes', nargs='+'); ap.add_argument('--tries', type=int, default=200)
    ap.add_argument('--box', type=int, default=3); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--K', type=int, default=8)
    a = ap.parse_args()
    systems = []
    for pf in a.prefixes:
        if open(pf + '.p1.ms').read().startswith('#'):
            print(f"{pf}: EMPTY CORE (every beta eliminated by constant pivots: solvable for every alpha)"); continue
        systems.append(load_Q(pf))
    if not systems:
        print("EXISTS [every array's core is empty: solvable for every alpha; constant-pivot elimination is exact]"); return 0
    names, polys = joint(systems)
    n = len(names)
    p = 1073741789
    ai, bi = split_linear(names, polys)
    rng = random.Random(a.seed)
    print(f"system: {len(polys)} equations, {n} unknowns ({len(ai)} alphas shared, {len(bi)} betas), arrays {len(systems)}")
    for t in range(a.tries):
        alpha = {i: rng.randint(-a.box, a.box) for i in ai}
        A, r = beta_system(names, polys, ai, bi, alpha, p)
        sol = solve_mod(A, r, p)
        if sol is None:
            continue
        xbar = [0] * n
        for i in ai: xbar[i] = alpha[i]          # the exact integer (fixed during the lift)
        for j, b in enumerate(bi): xbar[b] = sol[j]
        intpolys = [red(P, p ** a.K) for P in polys]        # integer polys (coefficients in Z_(p), as residues)
        # Prop. 4 with the alphas fixed; fix also beta coordinates outside a column basis of J|beta
        J = jac(polys, n)
        Jb = [[ev(red(J[q][b], p), xbar, p) for b in bi] for q in range(len(polys))]
        _, piv = rref_mod(Jb, p) if Jb else ([], [])
        # piv = beta columns forming a basis of the column space of J|beta at xbar
        basis = [bi[j] for j in piv]
        fixed = [i for i in range(n) if i not in basis]
        for i in fixed:
            if i in bi: xbar[i] = 0
        # recompute a point with the non-basis betas at 0
        sol2 = solve_mod([[row[bi.index(b)] for b in basis] for row in A], r, p)
        if sol2 is None:
            continue
        for j, b in enumerate(basis): xbar[b] = sol2[j]
        st = hensel(intpolys, n, xbar, p, a.K, fixed=fixed)
        if st[0] != 'lifted':
            print(f"try {t}: alpha {[alpha[i] for i in ai]}: {st[0]} at step {st[1]}")
            continue
        x, pk = st[1], st[2]
        xq = []
        ok = True
        for i in range(n):
            if i in ai:
                xq.append(Fraction(alpha[i]))
            else:
                rr = ratrec(x[i], pk)
                if rr is None: ok = False; break
                xq.append(rr)
        if not ok:
            print(f"try {t}: lifted but no rational reconstruction (algebraic irrational beta?)"); continue
        bad = [q for q, P in enumerate(polys) if ev(P, xq) != 0]
        if bad:
            print(f"try {t}: reconstruction fails exact check on {len(bad)} equations"); continue
        print(f"EXISTS [exact rational point, all {len(polys)} equations verified over Q] try {t}: "
              + ", ".join(f"{names[i]}={xq[i]}" for i in range(n)))
        return 0
    print(f"NO CERTIFICATE after {a.tries} tries (box {a.box})")
    return 1

if __name__ == '__main__':
    sys.exit(main())
