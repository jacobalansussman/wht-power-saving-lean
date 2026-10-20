"""spurious_scan.py -- ex2-modp, REFUSAL calibration: for a system known to have no complex solution, list the small
primes p at which it nevertheless has F_p points ("spurious points"), and test each such point with Proposition 4
(Hensel with the Jacobian frozen at the point): if the Jacobian is injective there, an obstruction PROVES that the point
has no p-adic lift; otherwise the point is reported as singular (no conclusion from one path).

  python3 spurious_scan.py PREFIX [PREFIX2 ...] [--primes 2,3,5,...] [--cap N] [--K K]

PREFIXes are joined (alphas shared, betas per array), as for a unit.  Systems must be linear in the betas.
For each prime: enumerate every alpha in F_p^k when p^k <= cap (else a random sample of cap alphas), solve the beta
system mod p; the F_p points over that alpha form an affine space (particular solution + kernel).
"""
import sys, os, random, argparse, itertools
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zsys import *

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('prefixes', nargs='+'); ap.add_argument('--primes', default='2,3,5,7,11,13,17,19,23,29,31')
    ap.add_argument('--cap', type=int, default=60000); ap.add_argument('--K', type=int, default=6)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--maxtest', type=int, default=20)
    a = ap.parse_args()
    systems = []
    for pf in a.prefixes:
        if open(pf + '.p1.ms').read().startswith('#'):
            print(f"{pf}: EMPTY CORE"); continue
        systems.append(load_Q(pf))
    names, polys = joint(systems); n = len(names)
    ai, bi = split_linear(names, polys)
    dens = set()
    for P in polys:
        for c in P.values():
            d = Fraction(c).denominator
            for q in range(2, 64):
                if d % q == 0: dens.add(q)
    print(f"{' + '.join(os.path.basename(x) for x in a.prefixes)}: {len(polys)} eqs, {len(ai)} alphas, {len(bi)} betas; "
          f"coefficient denominators divisible by {sorted(dens) or 'none'}")
    rng = random.Random(a.seed)
    J = jac(polys, n)
    for p in [int(x) for x in a.primes.split(',')]:
        if any(p % q == 0 for q in dens if q > 1 and all(q % r for r in range(2, q))):
            print(f"  p={p}: skipped (divides a denominator)"); continue
        total = p ** len(ai)
        if total <= a.cap:
            alphas = itertools.product(range(p), repeat=len(ai)); mode = 'all'
        else:
            alphas = (tuple(rng.randrange(p) for _ in ai) for _ in range(a.cap)); mode = f'sample {a.cap}'
        npts = 0; tested = {'obstructed': 0, 'lifted': 0}
        first = None
        for al in alphas:
            alpha = dict(zip(ai, al))
            A, r = beta_system(names, polys, ai, bi, alpha, p)
            sol = solve_mod(A, r, p)
            if sol is None: continue
            npts += 1
            if sum(tested.values()) >= a.maxtest: continue
            x = [0] * n
            for i in ai: x[i] = alpha[i]
            for j, b in enumerate(bi): x[b] = sol[j]
            inj, rk, nf = injective_on([red(P, p) for P in polys], n, x, p)
            if not inj:
                # tangent space > 0: fix the coordinates outside a column basis of J (a slice) and lift along it
                Jb = [[ev(red(J[q][i], p), x, p) for i in range(n)] for q in range(len(polys))]
                _, piv = rref_mod(Jb, p)
                fixed = [i for i in range(n) if i not in piv]
                st = hensel([red(P, p ** a.K) for P in polys], n, x, p, a.K, fixed=fixed)
                key = 'slice-obstructed' if st[0] == 'obstructed' else 'slice-lifted'
                tested[key] = tested.get(key, 0) + 1
                if first is None: first = (al, f'tangent dim {nf - rk}; {key}' + (f' at step {st[1]}' if st[0] == 'obstructed' else ''))
                continue
            st = hensel([red(P, p ** a.K) for P in polys], n, x, p, a.K)
            tested['obstructed' if st[0] == 'obstructed' else 'lifted'] += 1
            if first is None: first = (al, f"{st[0]} at step {st[1]}" if st[0] != 'lifted' else 'lifted to p^K')
        print(f"  p={p}: alphas {mode}: {npts} alpha values carry F_p points; Prop.4 on up to {a.maxtest}: {tested}"
              + (f"; first: alpha={first[0]} {first[1]}" if first else ''), flush=True)

if __name__ == '__main__':
    main()
