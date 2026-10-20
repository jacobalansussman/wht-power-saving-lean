"""special_exist.py -- ex2-modp: existence certificates in SPECIAL POSITION for systems linear in beta.

  python3 special_exist.py PREFIX [--gb1 F --gb2 F] [--tries N] [--box B] [--seed S]

PREFIX.p1.ms / PREFIX.p2.ms: one system mod two primes (nag-rerun / msolve format).  GB files (msolve output, same
primes; default PREFIX.p1.gb / PREFIX.p2.gb) are used ONLY to locate the special position: their alpha-only elements
are reconstructed over Q (rational reconstruction mod p1, checked mod p2) and solved for a rational alpha by random
small-integer assignment + exact rational roots of the polynomials that become univariate.  Then beta is solved
EXACTLY over Q.  The verdict EXISTS is printed only after every equation of the original system (rebuilt over Q) has
been checked exactly at the rational point; nothing from the GB is trusted for the verdict (PROOF.md Theorem 1).
"""
import sys, os, random, argparse, math
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zsys import *

def parse_gb(path):
    txt = open(path).read()
    names = None
    for line in txt.split('\n'):
        if line.startswith('#variable order:'):
            names = [v.strip() for v in line.split(':', 1)[1].split(',')]
        if line.startswith('#field characteristic:'):
            p = int(line.split(':')[1])
    body = txt[txt.index('['):].strip().rstrip(':').strip()[1:-1]
    polys = [parse_poly(c.strip(), names, p) for c in body.split(',') if c.strip()]
    return names, p, polys

def rat_roots(P, i):
    """rational roots of a univariate polynomial (in variable i, exact Fraction coefficients)"""
    co = {}
    for k, c in P.items():
        co[k[i]] = co.get(k[i], 0) + c
    d = max(co); low = min(co)
    roots = set()
    if low > 0:
        roots.add(Fraction(0))
    co = {e - low: c for e, c in co.items()}
    d -= low
    if d == 0:
        return roots
    L = 1
    for c in co.values():
        L = L * c.denominator // math.gcd(L, c.denominator)
    ic = {e: int(c * L) for e, c in co.items()}
    a0, ad = abs(ic.get(0, 0)), abs(ic[d])
    if a0 == 0:
        return roots
    def divs(N):
        out = set()
        for q in range(1, int(N ** 0.5) + 1):
            if N % q == 0: out |= {q, N // q}
        return out
    if a0 > 10 ** 8 or ad > 10 ** 8:
        return roots
    for u in divs(a0):
        for v in divs(ad):
            for s in (1, -1):
                r = Fraction(s * u, v)
                if sum(c * r ** e for e, c in co.items()) == 0:
                    roots.add(r)
    return roots

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('prefix'); ap.add_argument('--gb1'); ap.add_argument('--gb2')
    ap.add_argument('--tries', type=int, default=3000); ap.add_argument('--box', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--nogb', action='store_true')
    a = ap.parse_args()
    names, polys = load_Q(a.prefix)
    n = len(names)
    ai, bi = split_linear(names, polys)
    G = []
    if not a.nogb:
        g1 = parse_gb(a.gb1 or a.prefix + '.p1.gb'); g2 = parse_gb(a.gb2 or a.prefix + '.p2.gb')
        gn, GQ = to_Q(g1, g2)
        assert gn == names, 'GB variable order differs'
        G = [P for P in GQ if all(sum(k[b] for b in bi) == 0 for k in P)]
    print(f"{os.path.basename(a.prefix)}: {len(polys)} eqs, {len(ai)} alphas, {len(bi)} betas; "
          f"alpha-only GB elements used as hints: {len(G)}", flush=True)
    rng = random.Random(a.seed)
    stats = {'dead': 0, 'beta-inconsistent': 0}
    for t in range(a.tries):
        val = {}
        Gs = [dict(P) for P in G]
        order = ai[:]; rng.shuffle(order)
        dead = False
        while len(val) < len(ai) and not dead:
            # propagate: any hint that is univariate -> pick one of its rational roots
            prog = True
            while prog and not dead:
                prog = False
                for P in Gs:
                    if not P: continue
                    vs = {i for k in P for i in range(n) if k[i]}
                    if not vs:
                        dead = True; break
                    if len(vs) == 1:
                        i = vs.pop(); rs = sorted(rat_roots(P, i))
                        if not rs:
                            dead = True; break
                        val[i] = rng.choice(rs)
                        Gs = [subst(Q, i, val[i]) for Q in Gs]
                        prog = True; break
            if dead: break
            free = [i for i in order if i not in val]
            if not free: break
            i = free[0]
            val[i] = Fraction(rng.randint(-a.box, a.box))
            Gs = [subst(Q, i, val[i]) for Q in Gs]
        if dead or any(Q for Q in Gs):
            stats['dead'] += 1; continue
        A, r = beta_system(names, polys, ai, bi, val)
        sol = solve_Q(A, r)
        if sol is None:
            stats['beta-inconsistent'] += 1; continue
        x = [Fraction(0)] * n
        for i in ai: x[i] = val[i]
        for j, b in enumerate(bi): x[b] = sol[j]
        bad = [q for q, P in enumerate(polys) if ev(P, x) != 0]
        if bad:
            print('internal: exact check failed', bad); continue
        print(f"EXISTS [exact rational point; all {len(polys)} equations verified over Q] try {t}: "
              + ", ".join(f"{names[i]}={x[i]}" for i in range(n) if x[i] != 0) + "  (all other unknowns 0)")
        return 0
    print(f"NO CERTIFICATE after {a.tries} tries: {stats}")
    return 1

if __name__ == '__main__':
    sys.exit(main())
