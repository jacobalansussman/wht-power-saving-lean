# certnum.py -- certified numerics for overdetermined polynomial systems with exact (Gaussian-rational) coefficients.
# ex2-certnum, 14 Oct 2026.  Ball arithmetic: python-flint (arb/acb), rigorous outward rounding.
#
# Tools
#   PolySys                        polynomial system F: C^n -> C^N, exact coefficients (Fraction pairs)
#   lm_solve                       damped Gauss-Newton from random starts (numpy, NOT rigorous; only proposes points)
#   choose_square                  rows / free columns from a pivoted QR of the Jacobian at a numerical point
#   krawczyk                       Krawczyk test on a square subsystem: K(X) in int(X)  =>  exactly one root in X
#   certify_point                  square subsystem + Krawczyk + enclosures of the leftover rows, then EXACT CLOSURE:
#                                  simplest Gaussian rationals inside the certified box, every row of F checked exactly
#   exclude_box                    rigorous exclusion: 0 not in an enclosure of some row, or of the dual combination
#                                  y.F (y = conj F(x0)), on a box  =>  F has no zero in that box
#   dual_certificate               global refusal: an exact vector y with sum_i y_i f_i == nonzero constant (degree 0),
#                                  read off a numerical near miss and rounded; or NulLA multipliers of degree d (exact)
import json, math, random, sys, time, itertools
from fractions import Fraction as Fr
import numpy as np
from flint import arb, acb, acb_mat, arb_mat, fmpq, fmpz, ctx, fmpq_mat

# ---------------------------------------------------------------- exact Gaussian rationals as (Fr, Fr)
def gmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def gadd(a, b): return (a[0] + b[0], a[1] + b[1])
GZERO = (Fr(0), Fr(0)); GONE = (Fr(1), Fr(0))

def fr_of_float(x): return Fr(x)          # exact binary value of a float

class PolySys:
    """polys[k] = dict {exponent tuple (len n): (Fr re, Fr im)}"""
    def __init__(self, n, polys, names=None, info=None):
        self.n = n; self.polys = [dict((tuple(e), (Fr(c[0]), Fr(c[1]))) for e, c in p.items()) for p in polys]
        self.N = len(self.polys); self.names = names or ['x%d' % i for i in range(n)]; self.info = info or {}
        self._np = None

    @classmethod
    def from_a2nag_json(cls, path):
        d = json.load(open(path)); n = d['nvar']; polys = []
        for eq in d['eqs']:
            p = {}
            for vars_, (re, im) in eq:
                e = [0] * n
                for v in vars_: e[v] += 1
                e = tuple(e); c = (fr_of_float(re), fr_of_float(im))
                p[e] = gadd(p.get(e, GZERO), c)
            polys.append({e: c for e, c in p.items() if c != GZERO})
        return cls(n, polys, info=d.get('info'))

    def substitute(self, fixed):
        """fix variables: fixed = {var: (Fr, Fr)}; returns (new PolySys in the remaining vars, list of kept var ids)"""
        keep = [v for v in range(self.n) if v not in fixed]; out = []
        for p in self.polys:
            q = {}
            for e, c in p.items():
                cc = c
                for v, val in fixed.items():
                    for _ in range(e[v]): cc = gmul(cc, val)
                ee = tuple(e[v] for v in keep)
                q[ee] = gadd(q.get(ee, GZERO), cc)
            out.append({e: c for e, c in q.items() if c != GZERO})
        return PolySys(len(keep), out, [self.names[v] for v in keep], self.info), keep

    def subsystem(self, rows):
        return PolySys(self.n, [self.polys[r] for r in rows], self.names, self.info)

    # ---------------- numpy (float) evaluation
    def _prep(self):
        if self._np is not None: return
        self._np = []
        for p in self.polys:
            if p:
                E = np.array(list(p.keys()), dtype=int).reshape(len(p), self.n)
                C = np.array([complex(float(c[0]), float(c[1])) for c in p.values()])
            else:
                E = np.zeros((0, self.n), dtype=int); C = np.zeros(0, complex)
            self._np.append((E, C))
    def eval_np(self, x):
        self._prep(); x = np.asarray(x, complex)
        return np.array([C @ np.prod(x[None, :] ** E, axis=1) if len(C) else 0j for E, C in self._np])
    def jac_np(self, x):
        self._prep(); x = np.asarray(x, complex); J = np.zeros((self.N, self.n), complex)
        for k, (E, C) in enumerate(self._np):
            if not len(C): continue
            for v in range(self.n):
                m = E[:, v] > 0
                if not m.any(): continue
                E2 = E[m].copy(); c2 = C[m] * E2[:, v]; E2[:, v] -= 1
                J[k, v] = c2 @ np.prod(x[None, :] ** E2, axis=1)
        return J

    # ---------------- ball evaluation (rigorous)
    @staticmethod
    def _acb_c(c): return acb(arb(fmpq(c[0].numerator, c[0].denominator)), arb(fmpq(c[1].numerator, c[1].denominator)))
    def _pows(self, X):
        maxe = [0] * self.n
        for p in self.polys:
            for e in p:
                for v in range(self.n): maxe[v] = max(maxe[v], e[v])
        P = []
        for v in range(self.n):
            row = [acb(1)]
            for k in range(maxe[v]): row.append(row[-1] * X[v])
            P.append(row)
        return P
    def eval_acb(self, X, rows=None):
        P = self._pows(X); out = []
        for k in (range(self.N) if rows is None else rows):
            s = acb(0)
            for e, c in self.polys[k].items():
                t = self._acb_c(c)
                for v in range(self.n):
                    if e[v]: t = t * P[v][e[v]]
                s += t
            out.append(s)
        return out
    def jac_acb(self, X, rows=None):
        P = self._pows(X); rows = range(self.N) if rows is None else rows; out = []
        for k in rows:
            r = [acb(0)] * self.n
            for e, c in self.polys[k].items():
                cc = self._acb_c(c)
                for v in range(self.n):
                    if not e[v]: continue
                    t = cc * e[v]
                    for w in range(self.n):
                        ew = e[w] - (1 if w == v else 0)
                        if ew: t = t * P[w][ew]
                    r[v] = r[v] + t
            out.append(r)
        return acb_mat(out) if out else None
    # ---------------- exact evaluation
    def eval_exact(self, xq, rows=None):
        out = []
        for k in (range(self.N) if rows is None else rows):
            s = GZERO
            for e, c in self.polys[k].items():
                t = c
                for v in range(self.n):
                    for _ in range(e[v]): t = gmul(t, xq[v])
                s = gadd(s, t)
            out.append(s)
        return out

# ---------------------------------------------------------------- numerical solver (proposes points only)
def lm_solve(F, x0, iters=400, tol=1e-28, lam=1e-3, fixed_mask=None):
    x = np.array(x0, complex); r = F.eval_np(x); f = float(np.vdot(r, r).real)
    for it in range(iters):
        J = F.jac_np(x)
        if fixed_mask is not None: J[:, fixed_mask] = 0
        A = np.concatenate([J, math.sqrt(lam) * np.eye(F.n)]); b = np.concatenate([-r, np.zeros(F.n)])
        dx = np.linalg.lstsq(A, b, rcond=None)[0]
        x2 = x + dx; r2 = F.eval_np(x2); f2 = float(np.vdot(r2, r2).real)
        if f2 < f: x, r, f = x2, r2, f2; lam = max(lam / 4, 1e-15)
        else: lam *= 8
        if f < tol or lam > 1e12: break
    return x, f

def numrank(J, gap=1e-7):
    s = np.linalg.svd(J, compute_uv=False)
    if not len(s) or s[0] == 0: return 0, s
    r = int((s > gap * s[0]).sum()); return r, s

def pivoted_qr_cols(M, k):
    """greedy column selection (Businger-Golub) -> list of k column indices"""
    M = np.array(M, complex); chosen = []; R = M.copy()
    for _ in range(k):
        norms = np.linalg.norm(R, axis=0)
        for c in chosen: norms[c] = -1
        j = int(np.argmax(norms)); chosen.append(j)
        q = R[:, j] / np.linalg.norm(R[:, j]); R = R - np.outer(q, q.conj() @ R)
    return chosen

# ---------------------------------------------------------------- simplest rationals
def simplest_between(lo, hi):
    """simplest fraction in the closed interval [lo, hi] (Fractions), Stern-Brocot / continued fractions"""
    if lo > hi: lo, hi = hi, lo
    if lo <= 0 <= hi: return Fr(0)
    if hi < 0: return -simplest_between(-hi, -lo)
    fl = math.floor(lo)
    if Fr(fl) == lo: return Fr(fl)
    if fl + 1 <= hi: return Fr(fl + 1)
    # lo, hi in (fl, fl+1)
    sub = simplest_between(1 / (hi - fl), 1 / (lo - fl))
    return fl + 1 / sub

def arb_to_interval(a):
    """rigorous [lo, hi] as Fractions from an arb ball"""
    m, e = a.mid().man_exp(); r_m, r_e = a.rad().man_exp() if a.rad() != 0 else (0, 0)
    mid = Fr(int(m)) * (Fr(2) ** int(e)); rad = Fr(int(r_m)) * (Fr(2) ** int(r_e))
    return mid - rad, mid + rad

def arb_mid_fr(a):
    m, e = a.mid().man_exp(); return Fr(int(m)) * (Fr(2) ** int(e))

# ---------------------------------------------------------------- Krawczyk
def box(y, rho):
    """complex box around midpoint vector y (list of complex/acb), half-width rho (float) in re and im"""
    out = []
    for v in y:
        v = acb(v) if not isinstance(v, acb) else v
        out.append(acb(arb(v.real.mid(), rho), arb(v.imag.mid(), rho)))
    return out

def newton_flint(G, y, steps=8):
    """Newton on a square system with acb midpoints (precision ctx.prec); returns midpoint list"""
    y = [acb(v) for v in y]
    for _ in range(steps):
        g = G.eval_acb(y); J = G.jac_acb(y)
        try:
            dy = J.solve(acb_mat([[v] for v in g]), nonstop=True)
        except Exception:
            Jm = acb_mat([[J[i, j].mid() for j in range(G.n)] for i in range(G.n)])
            dy = Jm.solve(acb_mat([[v.mid()] for v in g]))
        y = [acb((y[i] - dy[i, 0]).mid()) for i in range(G.n)]
    return y

def krawczyk(G, y, rho):
    """G square PolySys (n x n). y midpoints (acb, exact mids). Returns (ok, K box list)"""
    n = G.n; X = box(y, rho)
    Jy = G.jac_acb([acb(v.mid()) for v in y])
    Jm = acb_mat([[Jy[i, j].mid() for j in range(n)] for i in range(n)])
    try:
        C = Jm.inv()
    except Exception:
        return False, None
    C = acb_mat([[C[i, j].mid() for j in range(n)] for i in range(n)])     # any matrix is rigorous
    gy = G.eval_acb([acb(v.mid()) for v in y]); JX = G.jac_acb(X)
    I = acb_mat([[1 if i == j else 0 for j in range(n)] for i in range(n)])
    D = acb_mat([[X[i] - acb(y[i].mid())] for i in range(n)])
    K = acb_mat([[acb(y[i].mid())] for i in range(n)]) - C * acb_mat([[v] for v in gy]) + (I - C * JX) * D
    ok = all(X[i].contains_interior(K[i, 0]) for i in range(n))
    return ok, [K[i, 0] for i in range(n)]

# ---------------------------------------------------------------- certify an (approximate) solution of F
def certify_point(F, x0, prec=256, rhos=(1e-30, 1e-20, 1e-12, 1e-8), fix_free=True, round_den=64, verbose=True,
                  rank_gap=1e-7, fixed_override=None, prefer_fix=None):
    """Returns dict with keys: rank, rows, free (fixed coords), krawczyk ok, leftover enclosures, exact closure."""
    ctx.prec = prec; rep = {}
    x0 = np.array(x0, complex); J = F.jac_np(x0); r, s = numrank(J, rank_gap)
    rep['rank'] = r; rep['n'] = F.n; rep['N'] = F.N; rep['sv_gap'] = (float(s[r - 1] / s[0]), float(s[r] / s[0]) if r < len(s) else 0.0)
    rows = pivoted_qr_cols(J.T, r)                     # r independent rows
    Jc = J[rows, :].copy()
    if prefer_fix: Jc[:, prefer_fix] *= 1e-3                 # bias: these unknowns are fixed (sliced) first
    cols = pivoted_qr_cols(Jc, r)                      # r well-conditioned columns (determined unknowns)
    free = [v for v in range(F.n) if v not in cols]
    # coordinate slice: fix the free unknowns at simple Gaussian rationals near x0
    fixed = {}
    for v in free:
        if fixed_override and v in fixed_override: fixed[v] = fixed_override[v]; continue
        fixed[v] = (Fr(round(x0[v].real * round_den), round_den), Fr(round(x0[v].imag * round_den), round_den))
    Fs, keep = F.substitute(fixed)
    # re-solve the sliced system numerically (moved the free coords)
    xs = np.array([x0[v] for v in keep]); xs, fval = lm_solve(Fs, xs, iters=200)
    rep['resid_after_slice'] = fval
    Js = Fs.jac_np(xs); r2, _ = numrank(Js, rank_gap); rep['rank_sliced'] = r2
    rows2 = pivoted_qr_cols(Js.T, r2)
    G = Fs.subsystem(rows2); others = [k for k in range(Fs.N) if k not in rows2]
    rep['square'] = (G.N, G.n); rep['fixed'] = {F.names[v]: (str(c[0]), str(c[1])) for v, c in fixed.items()}
    if G.N != G.n:
        rep['krawczyk'] = 'not square after slicing (rank %d, %d unknowns): singular point' % (r2, Fs.n); return rep
    y = newton_flint(G, [complex(v) for v in xs], steps=10)
    ok = False; Kb = None
    for rho in rhos:
        ok, Kb = krawczyk(G, y, rho)
        if ok: rep['rho'] = rho; break
    rep['krawczyk'] = bool(ok)
    if not ok: return rep
    # leftover rows on the certified box
    encl = Fs.eval_acb(Kb, others); mx = 0.0; excl = []
    for k, e in zip(others, encl):
        if not e.contains(acb(0)): excl.append(k)
        mx = max(mx, float(abs(e).upper()))
    rep['leftover_rows'] = len(others); rep['leftover_max_abs_upper'] = mx; rep['leftover_excluding_zero'] = excl
    if excl:
        rep['verdict'] = 'REFUSED in box: the square subsystem has exactly one root there and row(s) %s are nonzero at it' % excl[:5]
        return rep
    # exact closure: simplest Gaussian rationals in the certified box
    xq = []
    for v in Kb:
        lo, hi = arb_to_interval(v.real); a = simplest_between(lo, hi)
        lo, hi = arb_to_interval(v.imag); b = simplest_between(lo, hi)
        xq.append((a, b))
    vals = Fs.eval_exact(xq)
    nz = [k for k, v in enumerate(vals) if v != GZERO]
    rep['exact_candidate_maxden'] = max(max(a.denominator, b.denominator) for a, b in xq)
    rep['exact_nonzero_rows'] = len(nz)
    if not nz:
        full = [None] * F.n
        for v, c in fixed.items(): full[v] = c
        for i, v in enumerate(keep): full[v] = xq[i]
        assert all(c == GZERO for c in F.eval_exact(full))
        rep['verdict'] = 'EXISTS (exact): a Gaussian-rational point in the certified box solves every row of F exactly'
        rep['exact_point'] = {F.names[v]: (str(full[v][0]), str(full[v][1])) for v in range(F.n)}
    else:
        rep['verdict'] = ('UNDECIDED: square subsystem root certified; %d leftover rows contain 0 in their enclosures '
                          '(max %.1e) but no simple Gaussian rational closes them' % (len(others), mx))
    return rep

# ---------------------------------------------------------------- exclusion
def exclude_box(F, x0, rho, y=None, prec=128):
    """rigorous: True if F has no zero in the box x0 +- rho (re and im).  Tests every row and the combination y.F"""
    ctx.prec = prec; X = box([complex(v) for v in x0], rho)
    vals = F.eval_acb(X)
    for k, v in enumerate(vals):
        if not v.contains(acb(0)): return True, 'row %d' % k
    if y is not None:
        # mean-value form of g = sum y_k f_k:  g(X) in g(x0) + grad g(X) (X - x0)
        y = [acb(complex(c)) for c in y]
        X0 = [acb(complex(v)) for v in x0]
        g0 = sum((yk * fk for yk, fk in zip(y, F.eval_acb(X0))), acb(0))
        JX = F.jac_acb(X)
        grad = [sum((y[k] * JX[k, j] for k in range(F.N)), acb(0)) for j in range(F.n)]
        g = g0 + sum((grad[j] * (X[j] - X0[j]) for j in range(F.n)), acb(0))
        if not g.contains(acb(0)): return True, 'dual combination y.F (mean-value form)'
    return False, None

def exclusion_radius(F, x0, y=None, lo_exp=-14, hi_exp=3, prec=128):
    """largest rho = 2^k (k in [lo_exp*3.32, hi_exp*3.32]) with a certified exclusion of the box x0 +- rho"""
    best = None; how = None
    for k in range(int(lo_exp * 3.33), int(hi_exp * 3.33) + 1):
        rho = 2.0 ** k; ok, h = exclude_box(F, x0, rho, y, prec)
        if ok: best, how = rho, h
        else: break
    return best, how

# ---------------------------------------------------------------- global refusal
def coeff_matrix(F, rows=None, mult_monos=None):
    """columns = rows (times multiplier monomials), entries = coefficients on all output monomials (exact)"""
    rows = list(range(F.N)) if rows is None else rows
    mult_monos = mult_monos or [tuple([0] * F.n)]
    cols = []; monos = {}
    for k in rows:
        for mm in mult_monos:
            col = {}
            for e, c in F.polys[k].items():
                ee = tuple(a + b for a, b in zip(e, mm))
                if ee not in monos: monos[ee] = len(monos)
                col[monos[ee]] = c
            cols.append(((k, mm), col))
    return cols, monos

def dual_from_point(F, x0, round_den_max=10**6):
    """degree-0 dual certificate read off a numerical near miss: y = conj(F(x0)), rescaled, rounded, verified exactly.
       Returns (ok, y_exact, constant) with sum y_k f_k == constant != 0 identically (=> V(F) empty over C)."""
    r = F.eval_np(np.array(x0, complex)); y = np.conj(r)
    if np.abs(y).max() == 0: return False, None, None
    # rescale so that the combination's constant term is 1 (numerically)
    cols, monos = coeff_matrix(F); z = tuple([0] * F.n)
    c0 = sum(y[k] * complex(float(col.get(monos.get(z, -1), GZERO)[0]), float(col.get(monos.get(z, -1), GZERO)[1]))
             for (k, _), col in cols if z in monos)
    if abs(c0) < 1e-12: return False, None, None
    y = y / c0
    for den in [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 128, 256, 1024, 10**4, round_den_max]:
        yq = [(Fr(round(v.real * den), den), Fr(round(v.imag * den), den)) for v in y]
        ok, const = check_combination(F, yq)
        if ok: return True, yq, const
    return False, None, None

def check_combination(F, yq):
    tot = {}
    for k, c in enumerate(yq):
        if c == GZERO: continue
        for e, a in F.polys[k].items():
            tot[e] = gadd(tot.get(e, GZERO), gmul(c, a))
    tot = {e: c for e, c in tot.items() if c != GZERO}
    z = tuple([0] * F.n)
    if set(tot.keys()) == {z}: return True, tot[z]
    return False, None

def monos_upto(n, d):
    out = []
    for deg in range(d + 1):
        for comb in itertools.combinations_with_replacement(range(n), deg):
            e = [0] * n
            for v in comb: e[v] += 1
            out.append(tuple(e))
    return out

def nulla(F, d, real_coeffs=True, var_subset=None):
    """exact Nullstellensatz certificate sum_k c_k f_k = 1 with deg c_k <= d, by exact linear algebra over Q
       (coefficients assumed real rational). var_subset: multiplier monomials only in these variables.
       Returns (found, info)."""
    n = F.n
    if var_subset is None: mm = monos_upto(n, d)
    else:
        base = monos_upto(len(var_subset), d); mm = []
        for e in base:
            ee = [0] * n
            for i, v in enumerate(var_subset): ee[v] = e[i]
            mm.append(tuple(ee))
    cols, monos = coeff_matrix(F, None, mm)
    z = tuple([0] * n)
    if z not in monos: monos[z] = len(monos)
    nr, nc = len(monos), len(cols)
    if nr * nc > 4e7: return None, 'too large %d x %d' % (nr, nc)
    cplx = any(c[1] != 0 for _, col in cols for c in col.values())
    if not cplx:
        A = fmpq_mat(nr, nc); b = fmpq_mat(nr, 1)
        for j, (_, col) in enumerate(cols):
            for i, c in col.items(): A[i, j] = fmpq(c[0].numerator, c[0].denominator)
        b[monos[z], 0] = 1
    else:   # Gaussian coefficients: realify (rows: monomial x re/im; columns: multiplier x re/im part)
        if 4 * nr * nc > 4e7: return None, 'too large %d x %d (realified)' % (2 * nr, 2 * nc)
        A = fmpq_mat(2 * nr, 2 * nc); b = fmpq_mat(2 * nr, 1)
        for j, (_, col) in enumerate(cols):
            for i, c in col.items():
                re = fmpq(c[0].numerator, c[0].denominator); im = fmpq(c[1].numerator, c[1].denominator)
                if c[0] != 0: A[i, 2 * j] = re; A[nr + i, 2 * j + 1] = re
                if c[1] != 0: A[nr + i, 2 * j] = im; A[i, 2 * j + 1] = -im
        b[monos[z], 0] = 1
        nr, nc = 2 * nr, 2 * nc
    # consistency of A c = b  <=>  rank [A b] = rank A
    rA = A.rank(); Ab = fmpq_mat(nr, nc + 1)
    for i in range(nr):
        for j in range(nc): Ab[i, j] = A[i, j]
        Ab[i, nc] = b[i, 0]
    rAb = Ab.rank()
    return (rAb == rA), 'matrix %d x %d, rank %d / %d' % (nr, nc, rA, rAb)

# ---------------------------------------------------------------- Route S: certified root + local syzygy closure
# If z* is the unique root of the square subsystem G in the Krawczyk box X, and for a leftover row f there are polynomials
# u, q_i with   u * f = sum_i q_i * g_i   identically (exact, over Q)   and   0 not in u(X),
# then u(z*) f(z*) = 0 with u(z*) != 0, so f(z*) = 0.  Doing this for every leftover row proves that z* is an EXACT zero
# of the whole overdetermined F, with no need to recognise z* (it may be irrational).  At a nonsingular isolated root
# such (u, q) always exist (f lies in the local ideal (G)_z* = m_z*), so the only limit is the degree searched.
def _poly_mul_mono(p, mm):
    return {tuple(a + b for a, b in zip(e, mm)): c for e, c in p.items()}

def syzygy_row(Fs, rows, j, Kb, du=1, dq=1, prec=128):
    from flint import fmpz_mat
    n = Fs.n; mu = monos_upto(n, du); mq = monos_upto(n, dq)
    cols = []
    f = Fs.polys[j]
    for m in mu: cols.append(('u', m, _poly_mul_mono(f, m)))
    for i in rows:
        for m in mq: cols.append(('q', (i, m), {e: (-c[0], -c[1]) for e, c in _poly_mul_mono(Fs.polys[i], m).items()}))
    monos = {}
    for _, _, p in cols:
        for e in p:
            if e not in monos: monos[e] = len(monos)
    nr, nc = 2 * len(monos), 2 * len(cols)          # realified: rows (monomial, re/im); columns (unknown, re/im part)
    if nr * nc > 3e7: return None, 'too large %d x %d' % (nr, nc)
    A = fmpz_mat(nr, nc); colL = []
    M0 = len(monos)
    for k, (_, _, p) in enumerate(cols):
        L = 1
        for c in p.values():
            for t in c: L = L * t.denominator // math.gcd(L, t.denominator)
        for e, c in p.items():
            re, im = int(c[0] * L), int(c[1] * L); r0 = monos[e]
            # real part a of the unknown: a*(re + i im);  imaginary part b: i b (re + i im) = -b im + i b re
            if re: A[r0, 2 * k] = re; A[M0 + r0, 2 * k + 1] = re
            if im: A[M0 + r0, 2 * k] = im; A[r0, 2 * k + 1] = -im
        colL.append(L)
    Nm, nul = A.nullspace()
    if nul == 0: return False, 'no syzygy at degree (%d,%d) [%d x %d]' % (du, dq, nr, nc)
    ctx.prec = prec
    P = Fs._pows(Kb)
    nu = len(mu)
    for b in range(nul):
        u = acb(0); nz = False
        for a in range(nu):
            cr, ci = int(Nm[2 * a, b]), int(Nm[2 * a + 1, b])
            if cr or ci:
                nz = True; t = acb(arb(fmpq(cr, colL[a])), arb(fmpq(ci, colL[a])))
                for v in range(n):
                    if mu[a][v]: t = t * P[v][mu[a][v]]
                u += t
        if nz and not u.contains(acb(0)):
            # independent exact re-check of the identity u*f - sum q_i g_i == 0 with Fractions (guards the matrix assembly)
            def coef(k):
                return (Fr(int(Nm[2 * k, b]), colL[k]), Fr(int(Nm[2 * k + 1, b]), colL[k]))
            tot = {}
            for k, (kind, key, _) in enumerate(cols):
                c = coef(k)
                if c == GZERO: continue
                if kind == 'u': src, mm = f, key
                else: src, mm = Fs.polys[key[0]], key[1]
                sgn = (Fr(1), Fr(0)) if kind == 'u' else (Fr(-1), Fr(0))
                for e, a in src.items():
                    ee = tuple(x + y for x, y in zip(e, mm))
                    tot[ee] = gadd(tot.get(ee, GZERO), gmul(gmul(c, sgn), a))
            if any(v != GZERO for v in tot.values()):
                return None, 'BUG: identity fails the exact re-check'
            return True, 'u*f = sum q_i g_i with deg u <= %d, deg q <= %d, u nonzero on the box [%d x %d, nullity %d; identity re-checked exactly]' % (du, dq, nr, nc, nul)
    return False, 'syzygies exist (nullity %d) but every u vanishes somewhere on the box [%d x %d]' % (nul, nr, nc)

def add_linear_slices(F, x0, k, rng, den=8):
    """k random integer linear equations through a Gaussian-rational point near x0 (real parts only if F is real)"""
    polys = list(F.polys)
    for _ in range(k):
        c = [int(v) for v in rng.integers(-3, 4, size=F.n)]
        p = {}
        const = GZERO
        for v in range(F.n):
            if c[v]:
                e = [0] * F.n; e[v] = 1; p[tuple(e)] = (Fr(c[v]), Fr(0))
                xv = (Fr(round(x0[v].real * den), den), Fr(round(x0[v].imag * den), den))
                const = gadd(const, gmul((Fr(-c[v]), Fr(0)), xv))
        z = tuple([0] * F.n)
        if const != GZERO: p[z] = const
        polys.append(p)
    return PolySys(F.n, polys, F.names, F.info)

def certify_point_S(F, x0, prec=256, rhos=(1e-30, 1e-20, 1e-12), rank_gap=1e-7, fix_vars=None, round_den=64,
                    degrees=((0, 0), (0, 1), (1, 1), (0, 2), (1, 2)), slice_mode='coord', rng=None):
    """Krawczyk on a square subsystem, then Route S for every leftover row.  fix_vars: which unknowns to slice
    (coordinate slices at simple rationals near x0); default: the complement of a pivoted column set."""
    ctx.prec = prec; rep = {}
    x0 = np.array(x0, complex); J = F.jac_np(x0); r, s = numrank(J, rank_gap)
    rep['rank'] = r
    if slice_mode == 'linear':
        N0 = F.N; F = add_linear_slices(F, x0, F.n - r, rng or np.random.default_rng(3)); fix_vars = []
        rep['slices'] = 'linear (%d random integer forms)' % (F.n - r)
    if fix_vars is None:
        rows = pivoted_qr_cols(J.T, r); cols = pivoted_qr_cols(J[rows, :], r)
        fix_vars = [v for v in range(F.n) if v not in cols]
    fixed = {v: (Fr(round(x0[v].real * round_den), round_den), Fr(round(x0[v].imag * round_den), round_den)) for v in fix_vars}
    Fs, keep = F.substitute(fixed)
    xs = np.array([x0[v] for v in keep]); xs, fval = lm_solve(Fs, xs, iters=300)
    rep['resid_after_slice'] = fval
    Js = Fs.jac_np(xs); r2, _ = numrank(Js, rank_gap)
    if r2 != Fs.n:
        rep['verdict'] = 'sliced point not isolated/nonsingular (rank %d of %d)' % (r2, Fs.n); return rep
    if slice_mode == 'linear':
        sl = list(range(N0, Fs.N)); J0 = Js[:N0]; rr = pivoted_qr_cols(J0.T, r2 - len(sl)) if r2 > len(sl) else []
        # keep the slices; complete with rows of F (pivoted on the Jacobian projected off the slice directions)
        Q, _ = np.linalg.qr(Js[sl].T); P = np.eye(Fs.n) - Q @ Q.conj().T
        rr = pivoted_qr_cols((J0 @ P).T, Fs.n - len(sl)); rows2 = rr + sl
    else:
        rows2 = pivoted_qr_cols(Js.T, r2)
    G = Fs.subsystem(rows2); others = [k for k in range(Fs.N) if k not in rows2 and (slice_mode != 'linear' or k < N0)]
    y = newton_flint(G, [complex(v) for v in xs], steps=10)
    ok = False
    for rho in rhos:
        ok, Kb = krawczyk(G, y, rho)
        if ok: break
    rep['krawczyk'] = bool(ok); rep['square'] = (G.N, G.n); rep['leftover'] = len(others)
    if not ok: return rep
    # is the root rational?  (only reported, not used)
    xq = []
    for v in Kb:
        lo, hi = arb_to_interval(v.real); a = simplest_between(lo, hi); lo, hi = arb_to_interval(v.imag); b = simplest_between(lo, hi)
        xq.append((a, b))
    rep['simplest_rational_closes'] = all(c == GZERO for c in Fs.eval_exact(xq))
    rep['root_sample'] = [str(Kb[i].real.mid())[:30] for i in range(min(3, len(Kb)))]
    closed = 0; failed = []
    for j in others:
        e = Fs.eval_acb(Kb, [j])[0]
        if not e.contains(acb(0)):
            rep['verdict'] = 'REFUSED in box: leftover row %d is nonzero at the unique root of the square subsystem' % j; return rep
        done = False
        for du, dq in degrees:
            okS, msg = syzygy_row(Fs, rows2, j, Kb, du, dq)
            if okS: done = True; break
        if done: closed += 1
        else: failed.append((j, msg))
    rep['route_S_closed'] = '%d of %d leftover rows' % (closed, len(others))
    if not failed:
        rep['verdict'] = 'EXISTS (Route S): the certified root of the square subsystem is an exact zero of every row of F'
    else:
        rep['verdict'] = 'UNDECIDED by Route S at these degrees: %d rows not closed, e.g. %s' % (len(failed), failed[0])
    return rep
