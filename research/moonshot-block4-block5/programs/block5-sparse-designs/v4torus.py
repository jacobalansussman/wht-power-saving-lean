"""v4torus.py (w8-b5sparse): exact V4 for d = 4 closures whose rows are all BINOMIAL, any number of columns.

Row r: v1 * x^{m1} + v2 * x^{m2} = 0 with monomials in the unknowns x = (alpha_1..alpha_n, w_c (c != seed)), w_seed = 1.
Only the seed's connected component (columns linked by rows) constrains alpha: every other component may be set to 0.
On the torus (all unknowns of the component nonzero) the rows say x^{E_r} = rho_r (E_r = m1 - m2, rho_r = -v2/v1).
Smith normal form E = U D V^{-1}...: we compute unimodular P (rows), Q (cols) with P E Q = D = diag(d_1..d_k, 0..).
With x = y^{Q^T} (x_j = prod_i y_i^{Q[j][i]}): x^{E_r} = y^{(E Q)_r}; P E Q = D gives y^{D_i} = prod_r rho_r^{P[i][r]}
=: sigma_i for i < k, and consistency sigma_i = 1 for i >= k (rows of P beyond the rank).  So y_i (i < k) is any d_i-th
root of sigma_i, y_i (i >= k, up to N) is free.  Everything is done exactly; roots of unity live in F_p,
p = 2305842973415982001 = 1 + k * lcm(1..24) (we need d_i | 2 * lcm(1..24) -> checked).
API: solve(stash, nal) -> V (dict) or (None, reason); points(V, rng) -> iterator of alpha vectors (mod p), one per
choice of roots, free parameters random; parts(V, choice) -> for every alpha_k: (constant mod p, exponent vector over the
free parameters) for the relaxed all-parameters test."""
from fractions import Fraction
import itertools, random

P2 = 2305842973415982001
ZL = 1465896768020624744
import os as _os
if _os.environ.get('PRIME2'):                       # second prime, same root-of-unity property (cross-check)
    P2 = 2305837608478644241; ZL = 558701126638716937
LL = 5354228880


def snf(E):
    """integer Smith form: returns (P, D, Q) with P E Q = D (lists), P (R x R), Q (N x N) unimodular."""
    R = len(E); N = len(E[0]) if R else 0
    A = [row[:] for row in E]
    P = [[int(i == j) for j in range(R)] for i in range(R)]
    Q = [[int(i == j) for j in range(N)] for i in range(N)]

    def swap_rows(i, j):
        A[i], A[j] = A[j], A[i]; P[i], P[j] = P[j], P[i]

    def swap_cols(i, j):
        for row in A:
            row[i], row[j] = row[j], row[i]
        for row in Q:
            row[i], row[j] = row[j], row[i]

    def add_row(src, dst, f):          # row dst += f * row src
        A[dst] = [a + f * b for a, b in zip(A[dst], A[src])]; P[dst] = [a + f * b for a, b in zip(P[dst], P[src])]

    def add_col(src, dst, f):          # col dst += f * col src
        for row in A:
            row[dst] += f * row[src]
        for row in Q:
            row[dst] += f * row[src]

    t = 0
    while t < min(R, N):
        # pivot: smallest nonzero |entry| in the remaining block
        best = None
        for i in range(t, R):
            for j in range(t, N):
                if A[i][j] and (best is None or abs(A[i][j]) < abs(A[best[0]][best[1]])):
                    best = (i, j)
        if best is None:
            break
        swap_rows(t, best[0]); swap_cols(t, best[1])
        done = False
        while not done:
            done = True
            for i in range(t + 1, R):
                if A[i][t]:
                    f = A[i][t] // A[t][t]; add_row(t, i, -f)
                    if A[i][t]:
                        swap_rows(t, i); done = False
            for j in range(t + 1, N):
                if A[t][j]:
                    f = A[t][j] // A[t][t]; add_col(t, j, -f)
                    if A[t][j]:
                        swap_cols(t, j); done = False
            if done:
                # divisibility condition
                for i in range(t + 1, R):
                    for j in range(t + 1, N):
                        if A[i][j] % A[t][t]:
                            add_row(i, t, 1); done = False; break
                    if not done:
                        break
        if A[t][t] < 0:
            A[t] = [-a for a in A[t]]; P[t] = [-a for a in P[t]]
        t += 1
    return P, A, Q


def fmod(fr):
    return (fr.numerator % P2) * pow(fr.denominator % P2, P2 - 2, P2) % P2


def solve(stash, nal):
    eqs = []; varidx = {}
    for l in (0, 1):
        s = stash[l]
        rows, cols, ent, kl = s['cl']
        if len(s['seedcols']) != 1:
            return None, 'out %d: %d seed columns' % (l + 1, len(s['seedcols']))
        seed = next(iter(s['seedcols']))
        byrow = {}
        for (r, c, cf) in ent:
            byrow.setdefault(r, []).append((c, cf))
        # components of the column graph
        adj = {}
        for r, lst in byrow.items():
            cs = [c for c, _ in lst]
            for c in cs:
                adj.setdefault(c, set()).update(cs)
        comp = {seed}; st = [seed]
        while st:
            c = st.pop()
            for c2 in adj.get(c, ()):
                if c2 not in comp:
                    comp.add(c2); st.append(c2)
        for r, lst in byrow.items():
            if not any(c in comp for c, _ in lst):
                continue
            terms = []
            for (c, cf) in lst:
                if len(cf) != 1:
                    return None, 'out %d: entry with %d coefficients' % (l + 1, len(cf))
                (k, v), = cf.items()
                mono = {}
                if k:
                    mono[('a', k)] = 1
                if c != seed:
                    mono[('w', l, c)] = mono.get(('w', l, c), 0) + 1
                terms.append((Fraction(v).limit_denominator(10 ** 6), mono))
            # merge equal monomials
            mm = {}
            for v, mono in terms:
                key = tuple(sorted(mono.items()))
                mm[key] = mm.get(key, 0) + v
            mm = {k_: v for k_, v in mm.items() if v != 0}
            if len(mm) == 1:
                return None, 'out %d: a monomial row (forces a zero)' % (l + 1)
            if len(mm) != 2:
                return None, 'out %d: a row with %d terms' % (l + 1, len(mm))
            (m1, v1), (m2, v2) = mm.items()
            for key in list(m1) + list(m2):
                if key[0] not in varidx:
                    varidx[key[0]] = len(varidx)
            eqs.append((dict(m1), dict(m2), -v2 / v1))
    if not eqs:
        return None, 'no equations'
    N = len(varidx)
    E = []; rho = []
    for (m1, m2, r) in eqs:
        row = [0] * N
        for k_, e in m1.items():
            row[varidx[k_]] += e
        for k_, e in m2.items():
            row[varidx[k_]] -= e
        E.append(row); rho.append(r)
    Pm, D, Q = snf(E)
    R = len(E)
    ds = [D[i][i] for i in range(min(R, N)) if D[i][i]]
    k = len(ds)
    # sigma_i = prod_r rho_r^{P[i][r]}  (rho are rationals; for the consistency rows they must give 1)
    sig = []
    for i in range(R):
        x = Fraction(1)
        for r in range(R):
            e = Pm[i][r]
            if e:
                x *= rho[r] ** e
        sig.append(x)
    for i in range(k, R):
        if sig[i] != 1:
            return None, 'inconsistent torus system (sigma_%d = %s): V4 has no point with all unknowns nonzero' % (i, sig[i])
    for i in range(k):
        if sig[i] not in (1, -1):
            return None, 'sigma_%d = %s not +-1 (roots not handled)' % (i, sig[i])
        if (2 * LL) % (2 * ds[i]) and (LL % ds[i]):
            return None, 'd_%d = %d: roots not in F_p' % (i, ds[i])
    return dict(varidx=varidx, Q=Q, ds=ds, sig=sig[:k], N=N, k=k, nal=nal, neq=R), None


def root_choices(V):
    """for each i < k, the list of y_i in F_p with y_i^{d_i} = sigma_i"""
    out = []
    for d, s in zip(V['ds'], V['sig']):
        if s == 1:
            assert LL % d == 0
            z = pow(ZL, LL // d, P2)
            out.append([pow(z, j, P2) for j in range(d)])
        else:
            # y^d = -1: y = zeta_{2d}^{odd}
            assert LL % (2 * d) == 0, 'need 2d | lcm(1..24)'
            z = pow(ZL, LL // (2 * d), P2)
            out.append([pow(z, j, P2) for j in range(1, 2 * d, 2)])
    return out


def alpha_from(V, roots, free):
    """alpha_k (k = 1..nal) mod p; alpha not in the system -> None (free).  x_j = prod_i y_i^{Q[j][i]}."""
    y = list(roots) + list(free)
    al = [None] * V['nal']
    for key, j in V['varidx'].items():
        if key[0] != 'a':
            continue
        x = 1
        for i in range(V['N']):
            e = V['Q'][j][i]
            if e:
                x = x * pow(y[i], e % (P2 - 1), P2) % P2
        al[key[1] - 1] = x
    return al


def alpha_parts(V, roots):
    """for the relaxed test: alpha_k = const_k * prod_{i >= k} y_i^{Q[j][i]}: returns {k: (const, exponent tuple)}"""
    out = {}
    for key, j in V['varidx'].items():
        if key[0] != 'a':
            continue
        c = 1
        for i in range(V['k']):
            e = V['Q'][j][i]
            if e:
                c = c * pow(roots[i], e % (P2 - 1), P2) % P2
        ex = tuple(V['Q'][j][i] for i in range(V['k'], V['N']))
        out[key[1]] = (c, ex)
    return out


if __name__ == "__main__":
    # self-test of snf on random small integer matrices
    rng = random.Random(1)
    for t in range(300):
        R, N = rng.randrange(1, 6), rng.randrange(1, 6)
        E = [[rng.choice([0, 0, 1, -1, 2]) for _ in range(N)] for _ in range(R)]
        Pm, D, Q = snf(E)
        PE = [[sum(Pm[i][r] * E[r][j] for r in range(R)) for j in range(N)] for i in range(R)]
        PEQ = [[sum(PE[i][j] * Q[j][c] for j in range(N)) for c in range(N)] for i in range(R)]
        assert PEQ == D, (E, D)
        for i in range(R):
            for j in range(N):
                if i != j:
                    assert D[i][j] == 0
        dg = [D[i][i] for i in range(min(R, N))]
        nz = [x for x in dg if x]
        assert all(nz[i + 1] % nz[i] == 0 for i in range(len(nz) - 1)), dg
    print("snf self-test: 300 random matrices OK (P E Q = D, diagonal, divisibility)")
