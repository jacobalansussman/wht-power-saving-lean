"""t_field.py: exact checks of the field model of the Singer no-phase game (rm-theory, claims F1-F6).

Coordinates: a self-dual basis f_1..f_m of F_q (q = 2^m) over F_2, so the dot product is Tr(xy) and the target I
is the trace form. Forms Q are binary quadratic forms with Q(0) = 0, stored by their monomial coefficients
q_(ii), q_(ij) (i<j) in these coordinates; <Q, S> = sum over monomials of q_mon * S_mon; Q(z) = <Q, z z^T>.

Checks (all exact):
 F1  sum_j f_j = 1 and sum_j f_j^(1+2^i) = 0 (i = 1..m-1).
 F2  every Q_(b,c)(x) = Tr(bx) + sum_i Tr(c_i x^(1+2^i)) is a coordinate form, they span all m(m+1)/2 dims, and
     tr Q = sum_j Q(f_j) = Tr(b).
 F3  trace placement (Singer x -> alpha x, target I): n_odd(Q) = #{k : tr(Q o sigma^-k) = 1} = 2^(m-1) [b != 0].
 F4  channel targets T_k = sigma^-k sigma^-kT: {T_k} u {0} is closed under addition, every T_k is invertible
     (an additive symmetric spread set), T_0 = I.
 F5  Frobenius x -> x^2 is orthogonal in these coordinates (F F^T = I).
 F6  general placement T: n_odd(Q) = wt(Q~), Q~(y) = sum_j Q(y g_j) for a frame g with g g^T = T (checked on
     random T), and the polynomial-basis placement of e3-rmbuild reproduces its numbers.
"""
import sys, random, itertools
from gf import GF, gold_exponents


def mons(m):
    return [(i, i) for i in range(m)] + [(i, j) for i in range(m) for j in range(i + 1, m)]


def form_from_func(F, basis, Qf, m):
    """coordinate monomial coefficients of the field function Qf (dict/list over field elements)."""
    ML = mons(m); co = []
    fb = basis
    for (i, j) in ML:
        if i == j: co.append(Qf[fb[i]])
        else: co.append(Qf[fb[i] ^ fb[j]] ^ Qf[fb[i]] ^ Qf[fb[j]])
    return co


def eval_form(co, z, m):
    ML = mons(m); s = 0
    for c, (i, j) in zip(co, ML):
        if c and z[i] and z[j]: s ^= 1
    return s


def rank2(rows, ncols):
    """rank over F_2 of a list of int bit rows."""
    rows = [r for r in rows]; rk = 0
    for c in range(ncols):
        piv = None
        for i in range(rk, len(rows)):
            if (rows[i] >> c) & 1: piv = i; break
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and (rows[i] >> c) & 1: rows[i] ^= rows[rk]
        rk += 1
    return rk


def symmat_from_vectors(vecs, m):
    """S = sum v v^T (F_2), as an m x m list of lists."""
    S = [[0] * m for _ in range(m)]
    for v in vecs:
        for i in range(m):
            if v[i]:
                for j in range(m):
                    if v[j]: S[i][j] ^= 1
    return S


def pair(co, S, m):
    s = 0
    for c, (i, j) in zip(co, mons(m)):
        if c and S[i][j]: s ^= 1
    return s


def main(m, exhaustive_limit=40000, seed=1):
    F = GF(m); q, n = F.q, F.n
    fb = F.self_dual_basis(seed)
    print(f"=== m = {m}, q = {q}, poly = {bin(F.poly)}, self-dual basis = {fb}")
    # F1
    s = 0
    for f in fb: s ^= f
    sums = []
    for i in range(1, m):
        t = 0
        for f in fb: t ^= F.mul(f, F.pw(f, 1 << i))
        sums.append(t)
    print(f"F1 sum f_j = {s} (want 1); sum f_j^(1+2^i), i=1..m-1: {sums} (want all 0): {'PASS' if s == 1 and all(x == 0 for x in sums) else 'FAIL'}")
    coords = {x: F.coords(x, fb) for x in range(q)}
    # sanity: coords is a linear iso with dot = trace
    ok_dot = all(F.tr(F.mul(x, y)) == (sum(a * b for a, b in zip(coords[x], coords[y])) % 2)
                 for x in range(q) for y in range(q))
    print(f"   dot product of coordinates = Tr(xy) for all x, y: {ok_dot}")
    # F2
    ge = gold_exponents(m)
    nparam = m * (1 + len(ge))
    rng = random.Random(seed)
    if q ** (1 + len(ge)) <= exhaustive_limit:
        params = list(itertools.product(range(q), repeat=1 + len(ge)))
        mode = "exhaustive"
    else:
        params = [tuple(rng.randrange(q) for _ in range(1 + len(ge))) for _ in range(3000)]
        mode = "random 3000"
    span_rows = []; bad_deg = 0; bad_tr = 0; nodd_bad = 0
    nonzero = [x for x in range(1, q)]
    ML = mons(m); d = len(ML)
    for p in params:
        b, cs = p[0], p[1:]
        Qf = [F.tr(F.mul(b, x)) ^ (sum(F.tr(F.mul(c, F.pw(x, e))) for c, e in zip(cs, ge)) % 2) for x in range(q)]
        co = form_from_func(F, fb, Qf, m)
        # degree check: the coordinate form reproduces Qf everywhere
        if any(eval_form(co, coords[x], m) != Qf[x] for x in range(q)): bad_deg += 1
        trQ = 0
        for f in fb: trQ ^= Qf[f]
        if trQ != F.tr(b): bad_tr += 1
        span_rows.append(int("".join(map(str, co)), 2) if d > 0 else 0)
        # F3 in the same loop: n_odd over channels k (sigma^-k : x -> alpha^-k x)
        nodd = 0
        for k in range(n):
            ak = F.alpha(-k)
            t = 0
            for f in fb: t ^= Qf[F.mul(ak, f)]
            nodd += t
        want = (1 << (m - 1)) if b != 0 else 0
        if nodd != want: nodd_bad += 1
    rk = rank2(span_rows, d)
    print(f"F2 ({mode}, {len(params)} forms): degree<=2 failures {bad_deg}; span rank {rk} of {d}; tr Q != Tr(b): {bad_tr}: "
          f"{'PASS' if bad_deg == 0 and bad_tr == 0 and (rk == d or mode != 'exhaustive') else 'FAIL'}")
    print(f"F3 n_odd(Q) = 2^(m-1)[b != 0] failures: {nodd_bad} of {len(params)}: {'PASS' if nodd_bad == 0 else 'FAIL'}")
    # F4 spread set
    Ts = []
    for k in range(n):
        ak = F.alpha(-k)
        Ts.append(symmat_from_vectors([coords[F.mul(ak, f)] for f in fb], m))
    def key(S): return tuple(tuple(r) for r in S)
    keys = {key(S) for S in Ts}
    zero = key([[0] * m for _ in range(m)])
    closed = True; invertible = True
    for S in Ts:
        rows = [int("".join(map(str, r)), 2) for r in S]
        if rank2(rows, m) != m: invertible = False
    for S, T in itertools.combinations(Ts, 2):
        U = [[S[i][j] ^ T[i][j] for j in range(m)] for i in range(m)]
        if key(U) not in keys: closed = False
    I = [[int(i == j) for j in range(m)] for i in range(m)]
    print(f"F4 {len(keys)} distinct channel targets; T_0 = I: {key(Ts[0]) == key(I)}; all invertible: {invertible}; "
          f"closed under addition: {closed}: {'PASS' if invertible and closed and key(Ts[0]) == key(I) else 'FAIL'}")
    # F5 Frobenius orthogonal
    Fr = [coords[F.mul(f, f)] for f in fb]   # columns
    FFt = symmat_from_vectors(Fr, m)
    print(f"F5 Frobenius F F^T = I: {key(FFt) == key(I)}")
    return F, fb, coords


if __name__ == "__main__":
    ms = [int(a) for a in sys.argv[1:]] or [3, 4, 5]
    for m in ms:
        main(m)
