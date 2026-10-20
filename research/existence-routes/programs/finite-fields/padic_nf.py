"""padic_nf.py -- ex2-modp: padic_exist.py generalised to number fields of degree up to --dmax (default 6).

  python3 padic_nf.py PREFIX [--primes 10009,10007,...] [--K 30] [--tries N] [--dmax D]

Same search as padic_exist.py (alpha hints from the modular Groebner bases, F_p point, Hensel lift of the forced
coordinates), then: the first irrational forced coordinate theta gets its minimal polynomial (degree d <= dmax, the
first d for which LLL on (1, v, ..., v^d) gives a small exact relation mod p^K); every other forced coordinate is
written in the basis 1, theta, ..., theta^(d-1) by LLL; then beta is solved EXACTLY in K = Q(theta) and every equation
is checked exactly in K.  Only that check can print EXISTS (PROOF.md Theorem 1).  Tries the primes in turn (a prime
where the needed field does not split carries no point: that is the census of Proposition 5)."""
import sys, os, random, argparse
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zsys import *
from special_exist import parse_gb
from padic_exist import attempt, lift_forced, relation

class NF:
    m = [Fraction(0), Fraction(1)]          # monic modulus, low degree first (default: Q itself, t = 0)
    def __init__(s, c):
        D = len(NF.m) - 1
        if isinstance(c, list): s.c = (c + [Fraction(0)] * D)[:D]
        else: s.c = [Fraction(c)] + [Fraction(0)] * (D - 1)
    @staticmethod
    def _red(p):
        m = NF.m; D = len(m) - 1; p = list(p)
        for k in range(len(p) - 1, D - 1, -1):
            c = p[k]
            if c:
                for j in range(D + 1): p[k - D + j] -= c * m[j]
        return (p + [Fraction(0)] * D)[:D]
    def __add__(s, o): o = o if isinstance(o, NF) else NF(o); return NF([a + b for a, b in zip(s.c, o.c)])
    __radd__ = __add__
    def __neg__(s): return NF([-a for a in s.c])
    def __sub__(s, o): return s + (-(o if isinstance(o, NF) else NF(o)))
    def __mul__(s, o):
        o = o if isinstance(o, NF) else NF(o)
        p = [Fraction(0)] * (2 * len(s.c))
        for i, a in enumerate(s.c):
            if a:
                for j, b in enumerate(o.c):
                    if b: p[i + j] += a * b
        return NF(NF._red(p))
    __rmul__ = __mul__
    def __pow__(s, k):
        r = NF(1)
        for _ in range(k): r = r * s
        return r
    def iszero(s): return all(a == 0 for a in s.c)
    def inv(s):
        # extended Euclid in Q[t]: find u with u * s = 1 mod m
        def trim(a):
            a = list(a)
            while a and a[-1] == 0: a.pop()
            return a
        def divmod_(a, b):
            a = trim(a); b = trim(b); q = [Fraction(0)] * max(len(a) - len(b) + 1, 1)
            while len(a) >= len(b) and a:
                c = a[-1] / b[-1]; k = len(a) - len(b); q[k] = c
                for j in range(len(b)): a[k + j] -= c * b[j]
                a = trim(a)
            return trim(q), a
        def sub(a, b):
            n = max(len(a), len(b)); a = a + [Fraction(0)] * (n - len(a)); b = b + [Fraction(0)] * (n - len(b))
            return trim([x - y for x, y in zip(a, b)])
        def mulp(a, b):
            if not a or not b: return []
            r = [Fraction(0)] * (len(a) + len(b) - 1)
            for i, x in enumerate(a):
                for j, y in enumerate(b): r[i + j] += x * y
            return trim(r)
        r0, r1 = trim(NF.m), trim(s.c); u0, u1 = [], [Fraction(1)]
        while r1:
            q, r = divmod_(r0, r1)
            r0, r1 = r1, r
            u0, u1 = u1, sub(u0, mulp(q, u1))
        if len(r0) != 1:
            raise ZeroDivisionError('modulus not irreducible (zero divisor)')
        return NF([c / r0[0] for c in u0] if u0 else [Fraction(0)])
    def __truediv__(s, o): o = o if isinstance(o, NF) else NF(o); return s * o.inv()
    def __repr__(s):
        terms = [f"{c}" + (f"*t^{i}" if i > 1 else ("*t" if i == 1 else "")) for i, c in enumerate(s.c) if c]
        return "(" + " + ".join(terms) + ")" if len(terms) > 1 else (terms[0] if terms else "0")

def evNF(P, x):
    tot = NF(0)
    for k, c in P.items():
        t = NF(c)
        for xi, e in zip(x, k):
            if e: t = t * (xi ** e)
        tot = tot + t
    return tot

def solveNF(M, rhs):
    n = len(M[0]) if M else 0
    A = [list(row) + [b] for row, b in zip(M, rhs)]; piv = []; r = 0
    for c in range(n + 1):
        pr = next((i for i in range(r, len(A)) if not A[i][c].iszero()), None)
        if pr is None: continue
        if c == n: return None
        A[r], A[pr] = A[pr], A[r]; inv = A[r][c].inv(); A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and not A[i][c].iszero():
                f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
    x = [NF(0)] * n
    for i, c in enumerate(piv): x[c] = A[i][n]
    return x

def reconstruct_nf(X, forced_idx, M, dmax, bound=10 ** 12):
    theta = None; d = 1; out = {}; minpoly = None
    for i in forced_idx:
        v = X[i] % M
        if theta is None:
            for dd in range(1, dmax + 1):
                c = relation([pow(v, k, M) for k in range(1, dd + 1)], M)
                if c[dd] != 0 and max(abs(t) for t in c) < bound and sum(ck * pow(v, k, M) for k, ck in enumerate(c)) % M == 0:
                    break
            else:
                return None
            if dd == 1:
                out[i] = [Fraction(-c[0], c[1])]; continue
            theta = v; d = dd; minpoly = [Fraction(ck, c[dd]) for ck in c]
            out[i] = [Fraction(0), Fraction(1)]; continue
        c = relation([pow(theta, k, M) for k in range(1, d)] + [v], M)
        if c[d] == 0 or max(abs(t) for t in c) > bound or (sum(c[k] * pow(theta, k, M) for k in range(d)) + c[d] * v) % M:
            return None
        out[i] = [Fraction(-c[k], c[d]) for k in range(d)]
    return minpoly, out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('prefix'); ap.add_argument('--primes', default='10009,10007,10037,10039,10061,10067,10069,10079')
    ap.add_argument('--K', type=int, default=30); ap.add_argument('--tries', type=int, default=200)
    ap.add_argument('--box', type=int, default=2); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--dmax', type=int, default=6)
    a = ap.parse_args()
    names, polys = load_Q(a.prefix); n = len(names)
    ai, bi = split_linear(names, polys)
    gn, GQ = to_Q(parse_gb(a.prefix + '.p1.gb'), parse_gb(a.prefix + '.p2.gb'))
    G = [P for P in GQ if all(sum(k[b] for b in bi) == 0 for k in P)]
    print(f"{os.path.basename(a.prefix)}: {len(polys)} eqs, {len(ai)} alphas, {len(bi)} betas, {len(G)} hints", flush=True)
    rng = random.Random(a.seed)
    for p in [int(x) for x in a.primes.split(',')]:
        stats = {'dead': 0, 'beta': 0, 'norecon': 0, 'exactfail': 0}
        for t in range(a.tries):
            r = attempt(names, polys, ai, bi, G, p, rng, a.box)
            if r is None: stats['dead'] += 1; continue
            if r == 'beta': stats['beta'] += 1; continue
            x, forced = r
            X, M = lift_forced(x, forced, p, a.K)
            rec = reconstruct_nf(X, [i for i, _ in forced], M, a.dmax)
            if rec is None: stats['norecon'] += 1; continue
            minpoly, vals = rec
            NF.m = minpoly if minpoly else [Fraction(0), Fraction(1)]
            try:
                xK = [NF(0)] * n
                for i in ai: xK[i] = NF(vals[i]) if i in vals else NF(x[i])
                A = [[NF(0)] * len(bi) for _ in polys]; rhs = [NF(0)] * len(polys); pos = {b: j for j, b in enumerate(bi)}
                for q, P in enumerate(polys):
                    for k, c in P.items():
                        tt = NF(c)
                        for i in ai:
                            if k[i]: tt = tt * (xK[i] ** k[i])
                        bj = next((b for b in bi if k[b]), None)
                        if bj is None: rhs[q] = rhs[q] - tt
                        else: A[q][pos[bj]] = A[q][pos[bj]] + tt
                sol = solveNF(A, rhs)
            except ZeroDivisionError:
                stats['exactfail'] += 1; continue
            if sol is None: stats['exactfail'] += 1; continue
            for j, b in enumerate(bi): xK[b] = sol[j]
            if any(not evNF(P, xK).iszero() for P in polys): stats['exactfail'] += 1; continue
            fld = "Q" if not minpoly else "Q(t), " + " + ".join(f"({c})*t^{k}" for k, c in enumerate(minpoly) if c) + " = 0"
            print(f"EXISTS [exact point over {fld}; all {len(polys)} equations verified exactly] p = {p}, try {t}:")
            print("   " + ", ".join(f"{names[i]}={xK[i]}" for i in range(n) if not xK[i].iszero()) + "  (others 0)")
            return 0
        print(f"   p = {p}: no certificate in {a.tries} tries {stats}", flush=True)
    print("NO CERTIFICATE")
    return 1

if __name__ == '__main__':
    sys.exit(main())
