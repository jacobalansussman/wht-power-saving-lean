#!/usr/bin/env python3
"""fpa.py (fp-build, packing round, 15 Oct 2026): the FRAME GAME A = C[Sp(2m,2)] with C[P] free, exact.

Why this game (NOTE.md section 1, PROOF.md): the engine group G has the Pauli group N (with the scalars of G) as a
normal subgroup inside the free group B, and G/N = Sp(2m,2).  The central idempotents of C[N] lie in C[B], so the
engine game splits into blocks with a COMMON SKELETON; the Pauli-trivial block is exactly this game (e3-rmbuild
Lemma N), and the faithful block (odd numbers of digit columns) is M_n(C) (x) C[Sp(2m,2)], into which every program
of this game lifts unchanged (g -> 1 (x) g).  So every engine unit gives a program here with the same W and M, and
every program here is right on the blocks of odd f.

Conventions.  Vectors (a|b) in F_2^(2m), a = x-part (bits 0..m-1), b = z-part (bits m..2m-1).  A matrix is a tuple
of 2m row bitmasks; (M v)_i = parity(M[i] & v).  Symplectic form w((a|b),(a'|b')) = a.b' + b.a'.
  T_z  (paid; image of dir z)   : (a|b) -> (a + (z.b) z | b)
  addr(Mt) (free; address map)  : (a|b) -> (Mt a | Mt^-T b)
  ph(S) (free; phase letters)   : (a|b) -> (a | b + S a),  S symmetric
  P = stabiliser of L0 = {(0|b)} = {[[A, 0], [C, A^-T]] : A^T C symmetric}; level(g) = rank of the x-block on L0
  (= rank of the top-right m x m block) = paid length (sl-combinatorics Lemma 1.1(d)).
  K = T_e1 ... T_em = [[I, I], [0, I]].
Contents: dict {(l, g): Fraction}.  Positions: g^-1 L0 (two terms can be cancelled by a free factor iff they have
the same input and the same position, i.e. g h^-1 in P).
"""
import itertools, random
from fractions import Fraction as Fr

def par(x): return bin(x).count("1") & 1

class Sp:
    def __init__(s, m):
        s.m = m; s.n2 = 2 * m
        s.I = tuple(1 << i for i in range(2 * m))
        s.Tcache = {}
        s.K = s.T_prod([1 << j for j in range(m)])
    # ---------------------------------------------------------------- basic matrix ops (row bitmasks)
    def mul(s, A, B):
        out = []
        for row in A:
            acc = 0; j = 0; r = row
            while r:
                if r & 1: acc ^= B[j]
                r >>= 1; j += 1
            out.append(acc)
        return tuple(out)
    def apply(s, A, v):
        return sum(par(A[i] & v) << i for i in range(s.n2))
    def transpose(s, A):
        n = s.n2
        return tuple(sum(((A[j] >> i) & 1) << j for j in range(n)) for i in range(n))
    def inv(s, A):
        # symplectic inverse: J A^T J with J = [[0, I], [I, 0]]
        m = s.m; At = s.transpose(A)
        def swap(v): return ((v & ((1 << m) - 1)) << m) | (v >> m)
        rows = [swap(At[(i + m) % (2 * m)]) for i in range(2 * m)]
        return tuple(rows)
    def is_symp(s, A):
        m = s.m
        J = tuple((1 << (i + m)) if i < m else (1 << (i - m)) for i in range(2 * m))
        return s.mul(s.mul(s.transpose(A), J), A) == J
    # ---------------------------------------------------------------- letters
    def T(s, z):
        if z not in s.Tcache:
            m = s.m
            rows = []
            for i in range(m):
                rows.append((1 << i) | ((z << m) if (z >> i) & 1 else 0))
            for i in range(m): rows.append(1 << (m + i))
            s.Tcache[z] = tuple(rows)
        return s.Tcache[z]
    def T_prod(s, zs):
        A = tuple(1 << i for i in range(2 * s.m))
        for z in zs: A = s.mul(s.T(z), A)
        return A
    def addr(s, cols):
        """address map x -> sum x_k cols[k] (cols = images of e_k as ints); returns [[M,0],[0,M^-T]]"""
        m = s.m
        Mrows = [sum(((cols[k] >> i) & 1) << k for k in range(m)) for i in range(m)]   # row i of M
        # M^-T: invert M over F2 then transpose
        Minv = s._inv_m(Mrows)
        MinvT = [sum(((Minv[k] >> i) & 1) << k for k in range(m)) for i in range(m)]
        rows = [Mrows[i] for i in range(m)] + [MinvT[i] << m for i in range(m)]
        return tuple(rows)
    def _inv_m(s, rows):
        m = s.m
        A = list(rows); B = [1 << i for i in range(m)]
        for c in range(m):
            p = next(i for i in range(c, m) if (A[i] >> c) & 1)
            A[c], A[p] = A[p], A[c]; B[c], B[p] = B[p], B[c]
            for i in range(m):
                if i != c and (A[i] >> c) & 1: A[i] ^= A[c]; B[i] ^= B[c]
        return B
    def ph(s, S):
        """S as list of m row bitmasks (symmetric): (a|b) -> (a | b + S a)"""
        m = s.m
        rows = [1 << i for i in range(m)] + [(1 << (m + i)) | S[i] for i in range(m)]
        return tuple(rows)
    def ph_v(s, v):
        """rank-one phase S = v v^T"""
        return s.ph([v if (v >> i) & 1 else 0 for i in range(s.m)])
    # ---------------------------------------------------------------- geometry
    def xblock(s, A):
        m = s.m
        return [(A[i] >> m) & ((1 << m) - 1) for i in range(m)]
    def rank_rows(s, rows):
        rows = [r for r in rows if r]; rk = 0
        rows = list(rows)
        while rows:
            p = rows.pop()
            if not p: continue
            rk += 1
            lb = p & -p
            rows = [r ^ p if r & lb else r for r in rows]
            rows = [r for r in rows if r]
        return rk
    def level(s, A): return s.rank_rows(s.xblock(A))
    def in_P(s, A): return all(r == 0 for r in s.xblock(A))
    def pos(s, A):
        """position A^-1 L0 as a canonical RREF tuple of basis vectors"""
        Ai = s.inv(A); m = s.m
        vecs = [s.apply(Ai, 1 << (m + j)) for j in range(m)]
        return s._rref(vecs)
    def _rref(s, vecs):
        rows = [v for v in vecs if v]
        out = []
        for bit in reversed(range(2 * s.m)):
            piv = next((r for r in rows if (r >> bit) & 1), None)
            if piv is None: continue
            rows = [r ^ piv if (r >> bit) & 1 else r for r in rows if r != piv]
            out = [o ^ piv if (o >> bit) & 1 else o for o in out]
            out.append(piv)
            rows = [r for r in rows if r]
        return tuple(sorted(out))
    def dist(s, A, B):
        """distance between positions of A and B = level(A B^-1)"""
        return s.level(s.mul(A, s.inv(B)))
    # ---------------------------------------------------------------- random free letters
    def rand_P(s, rng, levi=True, phase=True):
        m = s.m
        A = s.I
        if levi:
            while True:
                cols = [rng.randrange(1, 1 << m) for _ in range(m)]
                if s.rank_rows(cols) == m: break
            A = s.mul(s.addr(cols), A)
        if phase:
            S = [0] * m
            for i in range(m):
                for j in range(i, m):
                    if rng.random() < 0.5:
                        S[i] |= 1 << j; S[j] |= 1 << i
            A = s.mul(s.ph(S), A)
        return A

# -------------------------------------------------------------------- contents
def cadd(*pairs):
    out = {}
    for co, c in pairs:
        for k, v in c.items():
            out[k] = out.get(k, 0) + co * v
    return {k: v for k, v in out.items() if v != 0}
def lmul(S, A, c):
    out = {}
    for (l, g), v in c.items():
        k = (l, S.mul(A, g)); out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}
def inp(S, l): return {(l, S.I): Fr(1)}

# -------------------------------------------------------------------- the exact end solve (closure + peel + solve)
def end_solve(S, states, target, cap_cols=20000, cap_rows=60000, want_solution=True, extra_target=None):
    """Is target (a content dict) in the C[P]-span of the states?  Unknowns beta_r(p) for (r, p).  Exact over Q.
    Returns (status, info): status in {"yes", "no", "cap"}; info has the solution or the reason."""
    # index walkers by (input, position)
    walk = {}
    for r, st in enumerate(states):
        for (l, g), v in st.items():
            walk.setdefault((l, S.pos(g)), []).append((r, g))
    posc = {}
    def P_of(g):
        if g not in posc: posc[g] = S.pos(g)
        return posc[g]
    rows = {}            # (l, g) -> index
    cols = {}            # (r, p) -> index
    colent = []          # list of dict row -> coef
    rowprod = []         # row -> set of cols
    queue = []
    def add_row(key):
        if key not in rows:
            rows[key] = len(rows); rowprod.append(set()); queue.append(key)
        return rows[key]
    for k in target: add_row(k)
    while queue:
        if len(cols) > cap_cols or len(rows) > cap_rows:
            return "cap", dict(rows=len(rows), cols=len(cols))
        key = queue.pop()
        l, g = key
        for (r, h) in walk.get((l, P_of(g)), []):
            p = S.mul(g, S.inv(h))
            assert S.in_P(p)
            ck = (r, p)
            if ck in cols: continue
            ci = len(cols); cols[ck] = ci
            ent = {}
            for (l2, g2), v in states[r].items():
                ri = add_row((l2, S.mul(p, g2)))
                ent[ri] = ent.get(ri, 0) + v
            ent = {a: b for a, b in ent.items() if b != 0}
            colent.append(ent)
            for ri in ent: rowprod[ri].add(ci)
    # peel: a non-target row with exactly one live producer kills it; target rows with none -> "no"
    nr, nc = len(rows), len(cols)
    tgt = {rows[k]: v for k, v in target.items()}
    alive = [True] * nc
    live_count = [len(rowprod[i]) for i in range(nr)]
    stack = [i for i in range(nr) if live_count[i] == 1 and i not in tgt]
    while stack:
        i = stack.pop()
        if live_count[i] != 1: continue
        c = next(cc for cc in rowprod[i] if alive[cc])
        alive[c] = False
        for ri in colent[c]:
            live_count[ri] -= 1
            if live_count[ri] == 1 and ri not in tgt: stack.append(ri)
    for i in tgt:
        if live_count[i] == 0:
            inv_rows = {v: k for k, v in rows.items()}
            # singleton rows of the closure (before the peel), with their BFS depth from the target
            depth = {i0: 0 for i0 in tgt}; frontier = list(tgt); seen_c = set()
            while frontier:
                nf = []
                for ri in frontier:
                    for c in rowprod[ri]:
                        if c in seen_c: continue
                        seen_c.add(c)
                        for r2 in colent[c]:
                            if r2 not in depth: depth[r2] = depth[ri] + 1; nf.append(r2)
                frontier = nf
            singles = sorted([(depth.get(ri, 99), inv_rows[ri]) for ri in range(nr)
                              if len(rowprod[ri]) == 1 and ri not in tgt], key=lambda x: x[0])
            return "no", dict(reason="target row unproduced after peel", rows=nr, cols=nc,
                              alive=sum(alive), singles=singles)
    live_cols = [c for c in range(nc) if alive[c]]
    live_rows = sorted({ri for c in live_cols for ri in colent[c]} | set(tgt))
    # exact solve: sum_c x_c colent[c] = tgt on live rows (rows of dead columns only are 0 anyway)
    ok, sol = solve_exact(live_rows, live_cols, colent, tgt)
    info = dict(rows=nr, cols=nc, live_cols=len(live_cols), live_rows=len(live_rows))
    if not ok:
        return "no", dict(info, reason="linear system inconsistent")
    inv_cols = {v: k for k, v in cols.items()}
    info["solution"] = [(inv_cols[c][0], inv_cols[c][1], x) for c, x in sol.items() if x != 0]
    return "yes", info

P_MOD = 2147483647
def solve_exact(live_rows, live_cols, colent, tgt):
    """mod-p elimination to find a consistent pivot set, then exact Fraction solve on the pivots; verify."""
    ridx = {r: i for i, r in enumerate(live_rows)}
    n = len(live_cols)
    # build sparse rows: equation per live row: sum_c a_{rc} x_c = t_r
    eqs = [dict() for _ in live_rows]
    for j, c in enumerate(live_cols):
        for ri, v in colent[c].items():
            if ri in ridx: eqs[ridx[ri]][j] = v
    rhs = [tgt.get(r, Fr(0)) for r in live_rows]
    # Gaussian elimination over Q with Fractions (sparse dicts); fine for systems up to a few thousand
    piv = {}  # col -> (eq dict, rhs)
    order = []
    for e, b in zip(eqs, rhs):
        e = dict(e); b = Fr(b)
        for cpiv in order:
            if cpiv in e:
                f = e[cpiv]; pe, pb = piv[cpiv]
                for k, v in pe.items():
                    nv = e.get(k, 0) - f * v
                    if nv == 0: e.pop(k, None)
                    else: e[k] = nv
                b -= f * pb
        e = {k: v for k, v in e.items() if v != 0}
        if not e:
            if b != 0: return False, None
            continue
        cp = min(e)
        f = e[cp]
        e = {k: v / f for k, v in e.items()}; b = b / f
        # eliminate cp from existing pivots (keep reduced)
        for c2 in order:
            pe, pb = piv[c2]
            if cp in pe:
                g = pe[cp]
                for k, v in e.items():
                    nv = pe.get(k, 0) - g * v
                    if nv == 0: pe.pop(k, None)
                    else: pe[k] = nv
                piv[c2] = (pe, pb - g * b)
        piv[cp] = (e, b); order.append(cp)
    sol = {}
    for cp in order:
        pe, pb = piv[cp]
        sol[live_cols[cp]] = pb     # free variables = 0
    return True, sol

def verify_combo(S, states, sol, target):
    acc = {}
    for (r, p, x) in sol:
        acc = cadd((1, acc), (x, lmul(S, p, states[r])))
    diff = cadd((1, acc), (-1, target))
    return len(diff) == 0, diff

# -------------------------------------------------------------------- transvections, charts, orphans
def transv(S, u):
    """T_u (x) = x + w(x,u) u"""
    m = S.m; ua, ub = u & ((1 << m) - 1), u >> m
    Ju = ub | (ua << m)
    return tuple(((1 << i) ^ (Ju if (u >> i) & 1 else 0)) for i in range(2 * m))
Sp.transv = transv

def find_chart(S, w, g, X, rng, tries=200000):
    """mu in P with pos(w mu g) = X (w: group element, g: the source walker's element), or None"""
    if S.dist(S.I, S.inv(w)) != S.level(w): pass
    # necessary: d(pos g, X) == level(w)
    for _ in range(tries):
        mu = S.rand_P(rng)
        if S.pos(S.mul(w, S.mul(mu, g))) == X:
            return mu
    return None

def orphans(S, states, target, cap_cols=20000):
    """Rows that appear in columns producing the target and have only one producer (so they kill those columns)."""
    walk = {}
    for r, st in enumerate(states):
        for (l, g), v in st.items():
            walk.setdefault((l, S.pos(g)), []).append((r, g))
    out = []
    for (l, g), v in target.items():
        for (r, h) in walk.get((l, S.pos(g)), []):
            p = S.mul(g, S.inv(h))
            for (l2, g2), v2 in states[r].items():
                if (l2, g2) == (l, h): continue
                row = (l2, S.mul(p, g2)); X = S.pos(row[1])
                nprod = len(walk.get((l2, X), []))
                out.append(dict(col=(r, p), row=row, pos=X, nprod=nprod, level=S.level(row[1])))
    return out

# -------------------------------------------------------------------- the psi screen (left P-coinvariants)
def psi_vec(S, c):
    """psi(c): the position profile, dict (l, position) -> sum of coefficients"""
    out = {}
    for (l, g), v in c.items():
        k = (l, S.pos(g)); out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}
def psi_screen(S, states, l):
    """necessary condition: (l, L1) in the C-span of the psi(S_r) (psi(beta x) = aug(beta) psi(x)).  Exact."""
    target = {(l, S.pos(S.K)): Fr(1)}
    vecs = [psi_vec(S, st) for st in states]
    rows = sorted({k for v in vecs for k in v} | set(target), key=str)
    ridx = {k: i for i, k in enumerate(rows)}
    colent = [{ridx[k]: v for k, v in vec.items()} for vec in vecs]
    ok, sol = solve_exact(list(range(len(rows))), list(range(len(vecs))), colent, {ridx[k]: v for k, v in target.items()})
    return ok

# -------------------------------------------------------------------- polynomial coefficients (injection scalars as unknowns)
class Poly:
    """multilinear-ish polynomial in symbols c_i with Fraction coefficients: dict frozenset-of-(i,exp) -> Fraction"""
    __slots__ = ("t",)
    def __init__(s, t=None):
        s.t = {k: v for k, v in (t or {}).items() if v != 0}
    @staticmethod
    def const(x): return Poly({(): Fr(x)}) if x != 0 else Poly()
    @staticmethod
    def var(i): return Poly({((i, 1),): Fr(1)})
    def _co(s, o): return o if isinstance(o, Poly) else Poly.const(o)
    def __add__(s, o):
        o = s._co(o); t = dict(s.t)
        for k, v in o.t.items(): t[k] = t.get(k, 0) + v
        return Poly(t)
    __radd__ = __add__
    def __neg__(s): return Poly({k: -v for k, v in s.t.items()})
    def __sub__(s, o): return s + (-s._co(o))
    def __mul__(s, o):
        o = s._co(o); t = {}
        for k1, v1 in s.t.items():
            for k2, v2 in o.t.items():
                d = dict(k1)
                for i, e in k2: d[i] = d.get(i, 0) + e
                k = tuple(sorted(d.items())); t[k] = t.get(k, 0) + v1 * v2
        return Poly(t)
    __rmul__ = __mul__
    def __eq__(s, o):
        o = s._co(o); return s.t == o.t
    def __ne__(s, o): return not s.__eq__(o)
    def __hash__(s): return hash(tuple(sorted(s.t.items())))
    def __repr__(s): return "Poly(%r)" % s.t

def psi_lin_screen(S, states, l):
    """psi screen with every (state, monomial in the injection scalars) as its own column: if this fails, no choice of
    injection scalars and readout scalars passes psi."""
    target = {(l, S.pos(S.K)): Fr(1)}
    cols = []
    for st in states:
        bymono = {}
        for (ll, g), v in st.items():
            pv = v if isinstance(v, Poly) else Poly.const(v)
            key = (ll, S.pos(g))
            for mono, co in pv.t.items():
                d = bymono.setdefault(mono, {}); d[key] = d.get(key, 0) + co
        for mono, d in bymono.items():
            d = {k: v for k, v in d.items() if v != 0}
            if d: cols.append(d)
    rows = sorted({k for v in cols for k in v} | set(target), key=str)
    ridx = {k: i for i, k in enumerate(rows)}
    colent = [{ridx[k]: v for k, v in c.items()} for c in cols]
    ok, sol = solve_exact(list(range(len(rows))), list(range(len(cols))), colent, {ridx[k]: v for k, v in target.items()})
    return ok

# -------------------------------------------------------------------- constructive transitivity of P on each sphere
def normal_P(S, Ub):
    """p in P with p(U) = U_d = <(e_1|0)..(e_d|0), (0|e_d+1)..(0|e_m)>, for a Lagrangian U given by m basis vectors"""
    m = S.m; mask = (1 << m) - 1
    # x-parts
    xs = [v & mask for v in Ub]
    # basis of A_U = span of x-parts (row reduce)
    basis = []
    for v in xs:
        for b in basis:
            if v & (b & -b): v ^= b
        if v:
            # keep reduced: eliminate the new pivot from others
            lb = v & -v
            basis = [b ^ v if b & lb else b for b in basis]
            basis.append(v)
    d = len(basis)
    # extend to a basis of F_2^m
    full = list(basis)
    for i in range(m):
        cand = 1 << i
        test = full + [cand]
        if S.rank_rows(test) == len(test): full.append(cand)
    # M maps full[k] -> e_k: M = (column matrix of full)^-1 ; Bmat has columns full[k]
    Brows = [sum(((full[k] >> i) & 1) << k for k in range(m)) for i in range(m)]
    Mrows = S._inv_m(Brows)
    cols_M = [sum(((Mrows[i] >> k) & 1) << i for i in range(m)) for k in range(m)]   # images of e_k
    g = S.addr(cols_M)
    U1 = [S.apply(g, v) for v in Ub]
    # reduce U1: vectors with x-part e_i (i < d) and pure (0|b') with b' in span(e_d..e_m-1)
    R = [0] * m          # R[i] = b-part (restricted to coords < d) of the vector with x-part e_i
    rows = list(U1)
    red = []
    for i in range(d):
        piv = next(v for v in rows if (v >> i) & 1)
        rows = [v ^ piv if ((v >> i) & 1) and v != piv else v for v in rows]
        rows.remove(piv)
        rows = [v for v in rows]
        red.append(piv)
    # make red[i] have x-part exactly e_i
    for i in range(d):
        for j in range(d):
            if j != i and (red[i] >> j) & 1: red[i] ^= red[j]
    Ssym = [0] * m
    for i in range(d):
        bpart = (red[i] >> m) & ((1 << d) - 1)
        for j in range(d):
            if (bpart >> j) & 1: Ssym[j] |= 1 << i      # S e_i has bit j  -> row j gets bit i
    # symmetrise check
    for i in range(d):
        for j in range(d):
            assert ((Ssym[i] >> j) & 1) == ((Ssym[j] >> i) & 1), "R block not symmetric"
    p = S.mul(S.ph(Ssym), g)
    return p, d
def chart_to(S, U, V):
    """mu in P with mu(U) = V (U, V Lagrangians as basis-vector tuples, same distance from L0), or None"""
    pU, dU = normal_P(S, list(U)); pV, dV = normal_P(S, list(V))
    if dU != dV: return None
    return S.mul(S.inv(pV), pU)
def find_chart_c(S, w, g, X):
    """constructive find_chart: mu in P with pos(w mu g) = X"""
    U = S._rref([S.apply(g, v) for v in X])
    V = S._rref([S.apply(S.inv(w), 1 << (S.m + j)) for j in range(S.m)])
    mu = chart_to(S, U, V)
    if mu is None: return None
    assert S.pos(S.mul(w, S.mul(mu, g))) == X
    return mu
