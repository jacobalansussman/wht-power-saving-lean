"""padic_exist.py -- ex2-modp: existence in special position via F_p points, Hensel lifting, algebraic reconstruction
and EXACT verification over a number field (PROOF.md Theorem 1, Propositions 4 and 5).

  python3 padic_exist.py PREFIX [--p P] [--K K] [--tries N] [--box B] [--seed S] [--census N_PRIMES]

PREFIX.p{1,2}.ms : the system mod two primes (rebuilt over Q); PREFIX.p{1,2}.gb : msolve GBs whose alpha-only
elements serve as HINTS (search only; never trusted).  Systems must be linear in the betas (b*, c*).

Search at the prime P (default 10009 = 1 mod 12, chosen as the first prime above 10^4 that splits in Q(i) and
Q(sqrt(-3)); --census reports which primes carry F_p points, which is how a field is read off):
 1. assign alphas in random order: a coordinate forced by a hint that became univariate gets a root mod P of that
    hint; any other coordinate gets a random small integer (kept exact).  A hint that becomes a nonzero constant
    kills the try.
 2. the beta system must be consistent mod P.
 3. Hensel: every forced coordinate is lifted to Z/P^K as a simple root of its (triangular) hint, in order.
 4. reconstruction: the first irrational forced coordinate theta gets a minimal polynomial of degree <= 2 by LLL;
    every other forced coordinate is written as (c0 + c1 theta)/c2 by LLL.
 5. EXACT: in K = Q(theta), recompute the forced coordinates, solve the beta system exactly over K, and evaluate
    EVERY equation of the original system exactly.  Only this step can print EXISTS.
"""
import sys, os, random, argparse
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zsys import *
from special_exist import parse_gb

# ------------------------------------------------------------ univariate polynomials mod p (lists, low degree first)
def ptrim(a):
    while a and a[-1] == 0: a.pop()
    return a
def pmod(a, b, p):
    a = a[:]; inv = pow(b[-1], -1, p)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p; s = len(a) - len(b)
        for i, x in enumerate(b): a[s + i] = (a[s + i] - c * x) % p
        ptrim(a)
    return a
def pmul(a, b, p):
    if not a or not b: return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] = (r[i + j] + x * y) % p
    return ptrim(r)
def ppow(base, e, mod, p):
    r = [1]; b = pmod(base, mod, p)
    while e:
        if e & 1: r = pmod(pmul(r, b, p), mod, p)
        b = pmod(pmul(b, b, p), mod, p); e >>= 1
    return r
def pgcd(a, b, p):
    a, b = ptrim(a[:]), ptrim(b[:])
    while b: a, b = b, pmod(a, b, p)
    if a:
        inv = pow(a[-1], -1, p); a = [x * inv % p for x in a]
    return a
def roots_mod(f, p, rng):
    f = ptrim([x % p for x in f])
    if len(f) <= 1: return []
    if p < 60:
        return [x for x in range(p) if sum(c * pow(x, i, p) for i, c in enumerate(f)) % p == 0]
    xp = ppow([0, 1], p, f, p)
    xp = xp + [0] * (2 - len(xp)) if len(xp) < 2 else xp
    h = xp[:]; h[1] = (h[1] - 1) % p; ptrim(h)
    g = pgcd(f, h, p) if h else f[:]
    out = []
    def split(g):
        if len(g) <= 1: return
        if len(g) == 2:
            out.append((-g[0]) * pow(g[1], -1, p) % p); return
        while True:
            a = rng.randrange(p)
            t = ppow([a, 1], (p - 1) // 2, g, p)
            t = t + [0] * (1 - len(t)) if not t else t
            t[0] = (t[0] - 1) % p; ptrim(t)
            d = pgcd(g, t, p) if t else g
            if 1 < len(d) < len(g):
                split(d); q = g[:]
                # q = g / d
                quo = []; r = g[:]; inv = pow(d[-1], -1, p)
                while len(r) >= len(d) and r:
                    c = r[-1] * inv % p; s = len(r) - len(d); quo.append((s, c))
                    for i, x in enumerate(d): r[s + i] = (r[s + i] - c * x) % p
                    ptrim(r)
                qq = [0] * (max(s for s, _ in quo) + 1)
                for s, c in quo: qq[s] = c
                split(ptrim(qq)); return
    split(g)
    return sorted(set(out))

# ------------------------------------------------------------ LLL (small dimension, exact)
def lll(B, delta=Fraction(3, 4)):
    B = [list(map(int, b)) for b in B]; n = len(B)
    def dot(u, v): return sum(x * y for x, y in zip(u, v))
    def gso(B):
        Bs = []; mu = [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = dot(B[i], Bs[j]) / dot(Bs[j], Bs[j]) if dot(Bs[j], Bs[j]) else Fraction(0)
                v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu
    Bs, mu = gso(B); k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [a - q * b for a, b in zip(B[k], B[j])]
                Bs, mu = gso(B)
        if dot(Bs[k], Bs[k]) >= (delta - mu[k][k - 1] ** 2) * dot(Bs[k - 1], Bs[k - 1]):
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu = gso(B); k = max(k - 1, 1)
    return B

def relation(vals, M):
    """small integer (c0, c1, ..., ck) with c0 + sum c_i vals_i = 0 mod M (vals: k residues)"""
    k = len(vals)
    B = [[M] + [0] * k] + [[(-v) % M] + [1 if j == i else 0 for j in range(k)] for i, v in enumerate(vals)]
    R = lll(B)
    return min(R, key=lambda r: sum(x * x for x in r))

# ------------------------------------------------------------ exact quadratic field Q(theta), theta^2 = -e theta - f
class QF:
    e = Fraction(0); f = Fraction(0)
    def __init__(s, a, b=0): s.a = Fraction(a); s.b = Fraction(b)
    def __add__(s, o): o = o if isinstance(o, QF) else QF(o); return QF(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __neg__(s): return QF(-s.a, -s.b)
    def __sub__(s, o): return s + (-(o if isinstance(o, QF) else QF(o)))
    def __rsub__(s, o): return QF(o) - s
    def __mul__(s, o):
        o = o if isinstance(o, QF) else QF(o)
        bb = s.b * o.b   # theta^2 = -e theta - f
        return QF(s.a * o.a - bb * QF.f, s.a * o.b + s.b * o.a - bb * QF.e)
    __rmul__ = __mul__
    def norm(s): return s.a * s.a - s.a * s.b * QF.e + s.b * s.b * QF.f
    def inv(s):
        N = s.norm(); assert N != 0
        return QF((s.a - s.b * QF.e) / N, -s.b / N)
    def __truediv__(s, o): o = o if isinstance(o, QF) else QF(o); return s * o.inv()
    def __pow__(s, k):
        r = QF(1)
        for _ in range(k): r = r * s
        return r
    def iszero(s): return s.a == 0 and s.b == 0
    def __eq__(s, o): o = o if isinstance(o, QF) else QF(o); return s.a == o.a and s.b == o.b
    def __repr__(s):
        if s.b == 0: return str(s.a)
        return f"({s.a}{'+' if s.b >= 0 else '-'}{abs(s.b)}*t)"

def evK(P, x):
    tot = QF(0)
    for k, c in P.items():
        t = QF(c)
        for xi, e in zip(x, k):
            if e: t = t * (xi ** e)
        tot = tot + t
    return tot

def solveK(M, rhs):
    n = len(M[0]) if M else 0
    A = [list(row) + [b] for row, b in zip(M, rhs)]
    piv = []; r = 0
    for c in range(n + 1):
        pr = next((i for i in range(r, len(A)) if not A[i][c].iszero()), None)
        if pr is None: continue
        if c == n: return None
        A[r], A[pr] = A[pr], A[r]
        inv = A[r][c].inv(); A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and not A[i][c].iszero():
                f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
    x = [QF(0)] * n
    for i, c in enumerate(piv): x[c] = A[i][n]
    return x

# ------------------------------------------------------------ the search
def univ_coeffs(P, i, x, m):
    """P with all variables but i evaluated at x (mod m) -> coefficient list in x_i mod m"""
    co = {}
    for k, c in P.items():
        t = coef_mod(c, m)
        for j, e in enumerate(k):
            if j != i and e: t = t * pow(x[j], e, m) % m
        co[k[i]] = (co.get(k[i], 0) + t) % m
    d = max(co) if co else 0
    return [co.get(e, 0) for e in range(d + 1)]

def attempt(names, polys, ai, bi, G, p, rng, box):
    n = len(names)
    x = [None] * n; forced = []          # forced: (i, hint poly)
    order = ai[:]; rng.shuffle(order)
    while any(x[i] is None for i in ai):
        prog = True
        while prog:
            prog = False
            for P in G:
                vs = {j for k in P for j in range(n) if k[j] and x[j] is None}
                if not vs:
                    if ev_partial(P, x, p) != 0: return None
                    continue
                if len(vs) == 1:
                    i = vs.pop()
                    xs = [v if v is not None else 0 for v in x]
                    co = univ_coeffs(P, i, xs, p)
                    if not any(co[1:]):
                        if co[0] % p: return None
                        continue
                    rs = roots_mod(co, p, rng)
                    # keep only simple roots (needed for Hensel)
                    dco = [(e * c) % p for e, c in enumerate(co)][1:]
                    rs = [r for r in rs if sum(c * pow(r, e, p) for e, c in enumerate(dco)) % p]
                    if not rs: return None
                    x[i] = rng.choice(rs); forced.append((i, P)); prog = True; break
        free = [i for i in order if x[i] is None]
        if not free: break
        x[free[0]] = rng.randint(-box, box)
    for P in G:
        if ev_partial(P, x, p) != 0: return None
    alpha = {i: x[i] % p for i in ai}
    A, r = beta_system(names, polys, ai, bi, alpha, p)
    if solve_mod(A, r, p) is None: return 'beta'
    return x, forced

def ev_partial(P, x, p):
    xs = [v if v is not None else 0 for v in x]
    return ev(red(P, p), xs, p)

def lift_forced(x, forced, p, K):
    M = p ** K; X = list(x); forced_idx = [i for i, _ in forced]
    for i, P in forced:
        co = univ_coeffs(P, i, X, M)
        r = X[i] % p; mod = p
        dco = [(e * c) for e, c in enumerate(co)][1:]
        while mod < M:
            mod2 = min(mod * mod, M)
            fv = sum(c * pow(r, e, mod2) for e, c in enumerate(co)) % mod2
            dv = sum(c * pow(r, e, mod2) for e, c in enumerate(dco)) % mod2
            r = (r - fv * pow(dv, -1, mod2)) % mod2; mod = mod2
        X[i] = r
    return X, M

def reconstruct(X, forced_idx, M):
    """returns (e, f, {i: (a, b)}) meaning x_i = a + b theta, theta^2 + e theta + f = 0; or None"""
    theta = None; out = {}; e = f = Fraction(0)
    for i in forced_idx:
        v = X[i]
        c = relation([v], M)                       # c0 + c1 v = 0 ?
        if c[1] != 0 and abs(c[0]) < 10 ** 12 and abs(c[1]) < 10 ** 12 and (c[0] + c[1] * v) % M == 0:
            out[i] = (Fraction(-c[0], c[1]), Fraction(0)); continue
        if theta is None:
            c = relation([v, v * v % M], M)        # c0 + c1 v + c2 v^2 = 0
            if c[2] == 0 or max(abs(t) for t in c) > 10 ** 12 or (c[0] + c[1] * v + c[2] * v * v) % M:
                return None
            theta = v; e = Fraction(c[1], c[2]); f = Fraction(c[0], c[2])
            out[i] = (Fraction(0), Fraction(1)); continue
        c = relation([theta, v], M)                # c0 + c1 theta + c2 v = 0
        if c[2] == 0 or max(abs(t) for t in c) > 10 ** 12 or (c[0] + c[1] * theta + c[2] * v) % M:
            return None
        out[i] = (Fraction(-c[0], c[2]), Fraction(-c[1], c[2]))
    return e, f, out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('prefix'); ap.add_argument('--p', type=int, default=10009); ap.add_argument('--K', type=int, default=24)
    ap.add_argument('--tries', type=int, default=400); ap.add_argument('--box', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--census', type=int, default=0)
    ap.add_argument('--census_tries', type=int, default=60)
    a = ap.parse_args()
    names, polys = load_Q(a.prefix); n = len(names)
    ai, bi = split_linear(names, polys)
    g1 = parse_gb(a.prefix + '.p1.gb'); g2 = parse_gb(a.prefix + '.p2.gb'); gn, GQ = to_Q(g1, g2)
    G = [P for P in GQ if all(sum(k[b] for b in bi) == 0 for k in P)]
    tag = os.path.basename(a.prefix)
    print(f"{tag}: {len(polys)} eqs, {len(ai)} alphas, {len(bi)} betas, {len(G)} alpha-only hints", flush=True)
    rng = random.Random(a.seed)
    if a.census:
        ps = []; q = 101
        while len(ps) < a.census:
            if all(q % d for d in range(2, int(q ** 0.5) + 1)): ps.append(q)
            q += 2
        res = []
        for q in ps:
            hit = any(isinstance(attempt(names, polys, ai, bi, G, q, rng, a.box), tuple) for _ in range(a.census_tries))
            res.append((q, hit))
        yes = [q for q, h in res if h]; no = [q for q, h in res if not h]
        print(f"census over {len(ps)} primes 101..{ps[-1]} ({a.census_tries} tries each): points at {len(yes)}")
        for mod in (3, 4, 8, 12):
            print(f"   residues mod {mod}: with points {sorted(set(q % mod for q in yes))}  without {sorted(set(q % mod for q in no))}")
    stats = {'dead': 0, 'beta': 0, 'norecon': 0, 'exactfail': 0}
    for t in range(a.tries):
        r = attempt(names, polys, ai, bi, G, a.p, rng, a.box)
        if r is None: stats['dead'] += 1; continue
        if r == 'beta': stats['beta'] += 1; continue
        x, forced = r
        X, M = lift_forced(x, forced, a.p, a.K)
        rec = reconstruct(X, [i for i, _ in forced], M)
        if rec is None: stats['norecon'] += 1; continue
        e, f, vals = rec
        QF.e, QF.f = e, f
        xK = [QF(0)] * n
        for i in ai:
            xK[i] = QF(*vals[i]) if i in vals else QF(x[i])
        # exact beta over K
        A = [[QF(0)] * len(bi) for _ in polys]; rhs = [QF(0)] * len(polys); pos = {b: j for j, b in enumerate(bi)}
        for q, P in enumerate(polys):
            for k, c in P.items():
                tt = QF(c)
                for i in ai:
                    if k[i]: tt = tt * (xK[i] ** k[i])
                bj = next((b for b in bi if k[b]), None)
                if bj is None: rhs[q] = rhs[q] - tt
                else: A[q][pos[bj]] = A[q][pos[bj]] + tt
        sol = solveK(A, rhs)
        if sol is None: stats['exactfail'] += 1; continue
        for j, b in enumerate(bi): xK[b] = sol[j]
        bad = [q for q, P in enumerate(polys) if not evK(P, xK).iszero()]
        if bad: stats['exactfail'] += 1; continue
        fld = "Q" if all(v.b == 0 for v in xK) else f"Q(t), t^2 + ({e}) t + ({f}) = 0"
        print(f"EXISTS [exact point over {fld}; all {len(polys)} equations verified exactly] try {t}, p = {a.p}:")
        print("   " + ", ".join(f"{names[i]}={xK[i]}" for i in range(n) if not xK[i].iszero()) + "  (others 0)")
        return 0
    print(f"NO CERTIFICATE after {a.tries} tries at p = {a.p}: {stats}")
    return 1

if __name__ == '__main__':
    sys.exit(main())
