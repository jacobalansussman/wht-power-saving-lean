"""allparam.py (w8-b5sparse): close the "special free parameters" gap of a point exclusion, exactly.

Input: a pickled truncated d = 3 system (alpha3 truncated mode), a localised subsystem (localize.py: rows within BFS
distance k of the seed, every column touching them; any set of complete rows is a valid test), a V4 root choice.
On that root choice alpha_k = const_k * t^{e_k} (t = the free parameters of V4, here two of them), so every entry of the
subsystem is a Laurent polynomial in t.  Clearing denominators column by column (multiplying a column by a monomial,
harmless since t != 0) gives a matrix M(t) = [N(t) | e] with polynomial entries.

If the system is consistent at some t (all t_i != 0), then rank [N(t)|e] <= rank N(t) <= r, the generic rank of N.  So
every (r+1)-minor of [N(t)|e] vanishes there; in particular f_j(t) = det(L_j M(t) R_j) for any fixed matrices L_j
((r+1) x rows) and R_j (cols+1 x (r+1)).  With two free parameters (t1, t2) we take three random pairs (L_j, R_j), compute
the resultants R_12(t1) = Res_t2(f_1, f_2) and R_13(t1) = Res_t2(f_1, f_3) by evaluation (Sylvester determinants over
F_p) and interpolation, and their gcd g(t1).  Every consistent t has g(t1) = 0 (or a vanishing leading coefficient in
t2, checked separately).  If g = c * t1^m, no t with t1 t2 != 0 is consistent: the component is excluded for EVERY
value of the free parameters (exact over F_p; the usual caveat that p is not one of finitely many bad primes).
Otherwise the roots of g in F_p are listed and tested directly.
usage: allparam.py sys.pkl loc.pkl rootchoice"""
import sys, pickle, random, itertools, time
import flint
import alpha3, v4torus
import modsolve as MS

P = alpha3.P2
D = pickle.load(open(sys.argv[1], 'rb')); LOC = pickle.load(open(sys.argv[2], 'rb')); ci0 = int(sys.argv[3])
l = LOC['l']; Rset = LOC['Rset']
nal = D['nal']; T = D.get('T'); V = D.get('V')
cl, seedcols, cap = D['sys'][l]
rows, cols, ent, kl = cl
klass, kterms = kl
# alpha_k = const_k * t^{e_k}
if T is not None:
    ch = v4torus.root_choices(T)
    roots = list(itertools.product(*ch))[ci0]
    pa = v4torus.alpha_parts(T, roots)
    nfree = T['N'] - T['k']
else:
    # two-column case: alpha_k = base(comp) * const: free parameters = the non-constant components
    a, b = V['sols'][ci0]
    z = alpha3.zeta(V['N']); s1 = pow(z, a, P); s2 = pow(z, b, P)
    comps = sorted(set(rt for rt in V['comp'].values() if rt != 0))
    pa = {}
    for k in range(1, nal + 1):
        if k not in V['pot']:
            continue
        C, E = V['pot'][k]; rt = V['comp'][k]
        x = (C.numerator % P) * pow(C.denominator % P, P - 2, P) % P
        x = x * pow(s1, E[0] % (P - 1), P) % P * pow(s2, E[1] % (P - 1), P) % P
        pa[k] = (x, tuple(int(rt == c_) for c_ in comps))
    nfree = len(comps)
assert nfree == 2, "this script handles exactly two free parameters (got %d)" % nfree
free_alpha = [k for k in range(1, nal + 1) if k not in pa]
# entries of the subsystem: {(row, col): {exponent: coeff}}
sub = {}
rng = random.Random(7 + l)
for (r, c, cf) in ent:
    if r + 1 not in Rset:
        continue
    d = sub.setdefault((r + 1, c), {})
    for k, v in cf.items():
        if k == 0:
            ex, cc = (0, 0), 1
        elif k in pa:
            cc, ex = pa[k]
        else:
            raise SystemExit("alpha_%d is free (not on V4's component): not handled" % k)
        x = MS.fval(v, P) * cc % P
        d[ex] = (d.get(ex, 0) + x) % P
# normalisation row 0 (seed K-row), constant entries
seedrows = set(rk for c in seedcols for (rk, cf) in kterms.get(c, []))
lam = {}
for c in klass:
    for (rk, cf) in kterms.get(c, []):
        if rk not in seedrows:
            continue
        lam.setdefault(rk, rng.randrange(1, P))
        d = sub.setdefault((0, c), {})
        for k, v in cf.items():
            assert k == 0
            d[(0, 0)] = (d.get((0, 0), 0) + lam[rk] * MS.fval(v, P)) % P
R_ = sorted(set(r for (r, c) in sub)); C_ = sorted(set(c for (r, c) in sub))
ri = {r: i for i, r in enumerate(R_)}; cj = {c: j for j, c in enumerate(C_)}
nR, nC = len(R_), len(C_)
# clear denominators per column: shift by the minimum exponent in each variable
shift = {}
for c in C_:
    exs = [ex for (r, cc), d in sub.items() if cc == c for ex in d]
    shift[c] = (min(e[0] for e in exs), min(e[1] for e in exs))
poly = {}   # (i, j) -> {(e1, e2): coeff}
for (r, c), d in sub.items():
    s = shift[c]
    poly[(ri[r], cj[c])] = {(e[0] - s[0], e[1] - s[1]): v for e, v in d.items() if v}
deg = [0, 0]
for d in poly.values():
    for e in d:
        deg[0] = max(deg[0], e[0]); deg[1] = max(deg[1], e[1])
print("subsystem: %d rows x %d columns (+ rhs); entry degrees in (t1, t2) up to %s; monomials used %s" % (
    nR, nC, deg, sorted(set(e for d in poly.values() for e in d))), flush=True)


def Mat(t1, t2, with_e=True):
    M = flint.nmod_mat(nR, nC + (1 if with_e else 0), P)
    for (i, j), d in poly.items():
        x = 0
        for (e1, e2), v in d.items():
            x = (x + v * pow(t1, e1, P) * pow(t2, e2, P)) % P
        M[i, j] = x
    if with_e:
        M[ri[0], nC] = 1
    return M


rr = random.Random(11)
pts = [(rr.randrange(1, P), rr.randrange(1, P)) for _ in range(5)]
rN = max(Mat(a, b, False).rank() for a, b in pts); rM = max(Mat(a, b, True).rank() for a, b in pts)
print("generic ranks: N %d, [N|e] %d" % (rN, rM), flush=True)
if rM == rN:
    raise SystemExit("subsystem generically consistent: nothing to do")
r = rN + 1                      # minors of size r vanish wherever the system is consistent
LR = []
for j in range(3):
    L = flint.nmod_mat(r, nR, [rr.randrange(P) for _ in range(r * nR)], P)
    Rm = flint.nmod_mat(nC + 1, r, [rr.randrange(P) for _ in range(r * (nC + 1))], P)
    LR.append((L, Rm))


def fval(j, t1, t2):
    L, Rm = LR[j]
    return (L * Mat(t1, t2) * Rm).det()


dmax2 = r * deg[1]              # degree of f_j in t2 at most r * (max entry degree in t2)
dmax1 = r * deg[0]


def newton(xs, ys):
    """coefficients (low to high) of the interpolating polynomial over F_p (Newton divided differences)"""
    n = len(xs); c = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) * pow((xs[i] - xs[i - j]) % P, P - 2, P) % P
    poly = [0] * n
    poly[0] = c[n - 1]
    deg = 0
    for i in range(n - 2, -1, -1):
        # poly = poly * (x - xs[i]) + c[i]
        new = [0] * n
        for k in range(deg + 1):
            new[k + 1] = (new[k + 1] + poly[k]) % P
            new[k] = (new[k] - xs[i] * poly[k]) % P
        new[0] = (new[0] + c[i]) % P
        poly = new; deg += 1
    return poly


def upoly_t2(j, t1):
    """coefficient list (low to high) of f_j(t1, .) in t2, interpolated at dmax2 + 1 points"""
    xs = list(range(1, dmax2 + 2))
    ys = [int(fval(j, t1, x)) for x in xs]
    return newton(xs, ys)


def resultant(a, b, m, n):
    """FORMAL resultant of polynomials with coefficient lists a, b (low to high) of formal degrees m, n"""
    a = (list(a) + [0] * (m + 1))[:m + 1][::-1]; b = (list(b) + [0] * (n + 1))[:n + 1][::-1]
    S = flint.nmod_mat(m + n, m + n, P)
    for i in range(n):
        for k, v in enumerate(a):
            S[i, i + k] = v
    for i in range(m):
        for k, v in enumerate(b):
            S[n + i, i + k] = v
    return int(S.det())


def degree(cl):
    d = len(cl) - 1
    while d >= 0 and cl[d] == 0:
        d -= 1
    return d


t0 = time.time()
f_test = [upoly_t2(j, 12345) for j in range(3)]
dg = [degree(f) for f in f_test]
print("f_j degrees in t2 at a generic t1: %s (bound %d); %.0fs" % (dg, dmax2, time.time() - t0), flush=True)
# degrees in t1 at a generic t2 (by symmetry of the construction)
def upoly_t1(j, t2):
    xs = list(range(1, dmax1 + 2))
    return newton(xs, [int(fval(j, x, t2)) for x in xs])
dg1 = [degree(upoly_t1(j, 54321)) for j in range(3)]
print("f_j degrees in t1 at a generic t2: %s" % dg1, flush=True)
# bivariate f_j: grid interpolation
import sympy
t1s, t2s = sympy.symbols('t1 t2')
F = []
for j in range(3):
    xs1 = list(range(1, dg1[j] + 2))
    cols_ = [upoly_t2(j, x) for x in xs1]                  # for each t1 value: coefficients in t2 (low..high)
    coef = {}
    for e2 in range(dg[j] + 1):
        c1 = newton(xs1, [cols_[i][e2] if e2 < len(cols_[i]) else 0 for i in range(len(xs1))])
        for e1, v in enumerate(c1):
            if v:
                coef[(e1, e2)] = v
    expr = sum(sympy.Integer(v if v <= P // 2 else v - P) * t1s ** e1 * t2s ** e2 for (e1, e2), v in coef.items())
    Fj = sympy.Poly(expr, t1s, t2s, modulus=P)
    # check at a fresh point
    a_, b_ = 777, 9191
    assert int(Fj.eval({t1s: a_, t2s: b_})) % P == int(fval(j, a_, b_)) % P, "bivariate interpolation check failed"
    F.append(Fj)
print("bivariate f_j built and checked (%.0fs)" % (time.time() - t0), flush=True)
G = sympy.gcd(F[0], F[1]); G = sympy.gcd(G, F[2])
print("gcd of f_1, f_2, f_3: total degree %d" % G.total_degree(), flush=True)
print("common part G =", G.as_expr(), flush=True)
fl = [(G.as_expr(), 1)]
mono = len(G.terms()) == 1
print("G is a monomial: %s" % mono, flush=True)
Q = [sympy.div(Fj, G)[0] for Fj in F]
r12 = sympy.resultant(Q[0], Q[1], t2s); r13 = sympy.resultant(Q[0], Q[2], t2s)
gg = sympy.gcd(sympy.Poly(r12, t1s, modulus=P), sympy.Poly(r13, t1s, modulus=P))
print("after removing the common part: gcd of resultants in t1 has degree %d: %s" % (gg.degree(), str(gg)[:200]), flush=True)
ggm = len(gg.terms()) == 1
if mono and ggm:
    print("VERDICT: EXCLUDED for every value of the free parameters (t1 t2 != 0): the common factor of the minors is the "
          "monomial %s and the remaining resultants share only %s" % (G.as_expr(), gg.as_expr()), flush=True)
else:
    print("VERDICT: inconclusive (common factor %s monomial: %s; resultant gcd monomial: %s)" % (G.as_expr(), mono, ggm), flush=True)
pickle.dump(dict(G=str(G.as_expr()), factors=[(str(f_), m_) for f_, m_ in fl], gg=str(gg.as_expr())), open('allparam_out_%s_c%d.pkl' % (sys.argv[2].split('/')[-1][:-4], ci0), 'wb'))
print("done (%.0fs)" % (time.time() - t0), flush=True)
