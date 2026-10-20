"""frames.py -- sl-algebra (15 Oct 2026): exact frame bookkeeping for the engine game at small m.

Frame of a group element g (a Clifford matrix up to scalars): lam(g) = g^-1 L_0, the image of the Z-type Lagrangian
L_0 = {(0|z)} under conjugation by g^-1, as a set of vectors (a|b) of F_2^2m (encoded a + (b << m)).
Free steps p (Paulis, affine maps, Clifford Z4 phases) satisfy lam(p g) = lam(g); a paid kernel moves lam by one step
of the dual polar graph DSp(2m,2) (distance d(A,B) = m - dim(A cap B)).

Checks (all exact up to snapping dyadic/Gaussian entries; matrices at m <= 3 or 4):
 (a) lam(free generator) = L_0; lam(one directional kernel C_j) is adjacent to L_0; lam(C^(x)m) = L_Y = {(z|z)};
     lam(H^(x)m) = L_X = {(x|0)}.
 (b) subspace frames L_U = {(x | x + w): x in U, w in U^perp} are Lagrangian for every subspace U (m <= 4);
     lam(Phi_U) = L_U for every nondegenerate U (Phi_U = J diag(i^{|P_U x|}) J, OpenAI's frame), m <= 3;
     d(L_U, L_V) = dim(U+V) - dim(U cap V) for all pairs, m <= 4.
 (c) OpenAI's pair in coordinate-free form: for every nondegenerate U, L_U and L_{U^perp} are opposite; the overlaps
     dim(L_U cap L_whole) = dim U and dim(L_0 cap L_{U^perp}) = dim U; for degenerate U, L_U and L_{U^perp} are NOT
     opposite (d < m).
 (d) the route operator Phi_{U^perp} Phi_U^-1 is C^(x)m times a Pauli (a free step) for every nondegenerate U, m <= 3.
 (e) cut ratios of DSp(2m,2): for the point cuts S_v = {L : v in L}, P[separate an opposite pair]/P[separate an edge]
     = 2 - 2^(1-m); and the spectral (eigenvector) ceiling over ALL cuts, m = 2, 3.
"""
import itertools, sys, json, os
PYLIB = os.environ.get('PYLIB', '')   # optional folder holding numpy
if os.path.isdir(PYLIB):
    sys.path.insert(0, PYLIB)
import numpy as np

OUT = {}

def popc(x):
    return bin(x).count("1")

def dot(a, b):
    return popc(a & b) & 1

# ---------------------------------------------------------------- matrices
def pauli(a, b, m):
    n = 1 << m
    M = np.zeros((n, n), dtype=complex)
    for x in range(n):
        M[x ^ a, x] = (-1) ** dot(b, x)
    return M

def ident_pauli(M, m):
    """M = c X(a) Z(b) -> (a, b); asserts it is a Pauli up to phase"""
    n = 1 << m
    a = int(np.argmax(np.abs(M[:, 0])))
    c = M[a, 0]
    assert abs(abs(c) - 1) < 1e-9, "not unitary-Pauli"
    b = 0
    for j in range(m):
        x = 1 << j
        v = M[a ^ x, x] / c
        if abs(v + 1) < 1e-9:
            b |= x
        else:
            assert abs(v - 1) < 1e-9
    P = c * pauli(a, b, m)
    assert np.max(np.abs(P - M)) < 1e-9, "not a Pauli"
    return a, b

def span(vecs):
    S = {0}
    for v in vecs:
        S |= {s ^ v for s in S}
    return frozenset(S)

def frame(g, m):
    """lam(g) = g^-1 L_0 as a frozenset of encoded vectors a + (b << m)"""
    gi = np.linalg.inv(g)
    gens = []
    for j in range(m):
        a, b = ident_pauli(gi @ pauli(0, 1 << j, m) @ g, m)
        gens.append(a | (b << m))
    L = span(gens)
    assert len(L) == 1 << m
    return L

def symp(u, v, m):
    msk = (1 << m) - 1
    return (dot(u & msk, v >> m) + dot(v & msk, u >> m)) & 1

def is_lagrangian(L, m):
    return len(L) == 1 << m and all(symp(u, v, m) == 0 for u in L for v in L)

def dist(A, B, m):
    k = len(A & B).bit_length() - 1
    return m - k

def kron(*Ms):
    out = np.array([[1.0 + 0j]])
    for M in Ms:
        out = np.kron(M, out)   # bit j = j-th factor from the right: factor 0 acts on bit 0
    return out

I2 = np.eye(2, dtype=complex)
Cm = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])
Hm = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
Sm = np.diag([1, 1j])

def on_bit(M2, j, m):
    return kron(*[(M2 if k == j else I2) for k in range(m)])

def perm_matrix(f, m):
    n = 1 << m
    M = np.zeros((n, n), dtype=complex)
    for x in range(n):
        M[f(x), x] = 1
    return M

def cnot(c, t, m):
    return perm_matrix(lambda x: x ^ (((x >> c) & 1) << t), m)

def cz(c, t, m):
    n = 1 << m
    return np.diag([(-1) ** (((x >> c) & (x >> t)) & 1) for x in range(n)]).astype(complex)

# ---------------------------------------------------------------- subspaces of F_2^m
def all_subspaces(m):
    seen = set(); out = []
    vecs = list(range(1, 1 << m))
    for k in range(m + 1):
        for basis in itertools.combinations(vecs, k):
            S = span(basis)
            if len(S) == 1 << k and S not in seen:
                seen.add(S); out.append(S)
    return out

def perp(U, m):
    return frozenset(w for w in range(1 << m) if all(dot(w, u) == 0 for u in U))

def L_of(U, m):
    Up = perp(U, m)
    return frozenset(x | ((x ^ w) << m) for x in U for w in Up)

def dimv(S):
    return len(S).bit_length() - 1

def nondeg(U, m):
    return len(U & perp(U, m)) == 1

def proj_matrix(U, m):
    """P_U over F_2 for nondegenerate U: P_U x = the U-component of x in F_2^m = U + U^perp"""
    Up = perp(U, m)
    P = {}
    for u in U:
        for w in Up:
            P[u ^ w] = u
    return P

def Phi(U, m):
    n = 1 << m
    J = kron(*[Hm] * m)
    P = proj_matrix(U, m)
    D = np.diag([1j ** (popc(P[x]) % 4) for x in range(n)])
    return J @ D @ J

def is_pauli_times(M, m):
    try:
        ident_pauli(M, m); return True
    except AssertionError:
        return False

# ---------------------------------------------------------------- (a) generators
def check_a():
    res = {}
    for m in (1, 2, 3):
        n = 1 << m
        L0 = span([1 << (j + m) for j in range(m)])
        LX = span([1 << j for j in range(m)])
        LY = span([(1 << j) | (1 << (j + m)) for j in range(m)])
        free = [pauli(1 << j, 0, m) for j in range(m)] + [pauli(0, 1 << j, m) for j in range(m)] + \
               [on_bit(Sm, j, m) for j in range(m)] + \
               [cnot(c, t, m) for c in range(m) for t in range(m) if c != t] + \
               [cz(c, t, m) for c in range(m) for t in range(m) if c < t]
        ok_free = all(frame(p, m) == L0 for p in free)
        ok_adj = all(dist(frame(on_bit(Cm, j, m), m), L0, m) == 1 for j in range(m))
        # pg: a free step on the left does not change the frame; a paid move changes it by exactly 1
        rng = np.random.default_rng(m)
        ok_left = True; ok_move = True
        for _ in range(30):
            g = np.eye(n, dtype=complex)
            for _ in range(6):
                if rng.random() < 0.5:
                    g = free[rng.integers(len(free))] @ g
                else:
                    g = on_bit(Cm, int(rng.integers(m)), m) @ g
            p = free[rng.integers(len(free))]
            ok_left &= frame(p @ g, m) == frame(g, m)
            h = on_bit(Cm, int(rng.integers(m)), m)
            ok_move &= dist(frame(h @ g, m), frame(g, m), m) == 1
        Cfull = kron(*[Cm] * m); Hfull = kron(*[Hm] * m)
        res[m] = dict(free_fix_L0=bool(ok_free), kernel_adjacent=bool(ok_adj), left_free_invariant=bool(ok_left),
                      paid_moves_one_step=bool(ok_move), C_full_is_LY=bool(frame(Cfull, m) == LY),
                      H_full_is_LX=bool(frame(Hfull, m) == LX), L0_LY_opposite=dist(L0, LY, m) == m)
    return res

# ---------------------------------------------------------------- (b), (c), (d)
def check_bcd():
    res = {}
    for m in (1, 2, 3, 4):
        subs = all_subspaces(m)
        Ls = {U: L_of(U, m) for U in subs}
        lag = all(is_lagrangian(Ls[U], m) for U in subs)
        dist_ok = all(dist(Ls[U], Ls[V], m) == dimv(span(list(U | V))) - dimv(U & V) for U in subs for V in subs)
        whole = frozenset(range(1 << m)); zero = frozenset([0])
        pair_ok = True; degen_not_opp = True; n_nd = 0; n_deg = 0
        for U in subs:
            Up = perp(U, m)
            if nondeg(U, m):
                n_nd += 1
                pair_ok &= dist(Ls[U], Ls[Up], m) == m
                pair_ok &= dimv(Ls[U] & Ls[whole]) == dimv(U)
                pair_ok &= dimv(Ls[zero] & Ls[Up]) == dimv(U)
            else:
                n_deg += 1
                degen_not_opp &= dist(Ls[U], Ls[Up], m) < m
        r = dict(n_subspaces=len(subs), n_nondegenerate=n_nd, all_L_U_lagrangian=bool(lag),
                 subspace_distance_formula=bool(dist_ok), pair_overlaps_and_opposite=bool(pair_ok),
                 degenerate_U_not_opposite=bool(degen_not_opp))
        if m <= 3:
            phi_ok = True; route_ok = True
            Cfull = kron(*[Cm] * m)
            for U in subs:
                if not nondeg(U, m):
                    continue
                phi_ok &= frame(Phi(U, m), m) == Ls[U]
                Up = perp(U, m)
                R = np.linalg.inv(Cfull) @ Phi(Up, m) @ np.linalg.inv(Phi(U, m))
                # R must be a Pauli up to a global phase
                route_ok &= is_pauli_times(R, m)
            r.update(lam_Phi_U_equals_L_U=bool(phi_ok), route_is_C_times_pauli=bool(route_ok))
        res[m] = r
    return res

# ---------------------------------------------------------------- (e) cut ratios
def all_lagrangians(m):
    """all Lagrangians of F_2^2m as frozensets (m <= 3)"""
    N = 1 << (2 * m)
    out = set()
    # grow isotropic spans greedily by DFS over sorted bases
    def rec(S, start, k):
        if k == m:
            out.add(S); return
        for v in range(start, N):
            if v in S:
                continue
            if all(symp(v, s, m) == 0 for s in S):
                rec(frozenset(S | {s ^ v for s in S}), v + 1, k + 1)
    rec(frozenset([0]), 1, 0)
    return sorted(out, key=lambda S: sorted(S))

def check_e():
    res = {}
    for m in (2, 3):
        Ls = all_lagrangians(m)
        n = len(Ls)
        idx = {L: i for i, L in enumerate(Ls)}
        D = np.array([[dist(A, B, m) for B in Ls] for A in Ls])
        A1 = (D == 1).astype(float); Am = (D == m).astype(float)
        k1 = A1.sum(1)[0]; km = Am.sum(1)[0]
        # point cuts
        N = 1 << (2 * m)
        ratios = []
        for v in range(1, N):
            s = np.array([1.0 if v in L else 0.0 for L in Ls])
            e_cut = s @ A1 @ (1 - s); a_cut = s @ Am @ (1 - s)
            ratios.append((a_cut / km) / (e_cut / k1))
        # spectral ceiling: A1 and Am commute (association scheme); joint eigenvectors
        w1, V = np.linalg.eigh(A1 + 1e-3 * Am)  # perturb to split
        th1 = np.round(np.einsum('ij,jk,ki->i', V.T, A1, V), 6)
        thm = np.round(np.einsum('ij,jk,ki->i', V.T, Am, V), 6)
        pairs = sorted(set(zip(th1, thm)))
        ceil = max((1 - tm / km) / (1 - t1 / k1) for (t1, tm) in pairs if abs(t1 - k1) > 1e-6)
        res[m] = dict(n_frames=n, valency=int(k1), opposite_per_frame=int(km),
                      point_cut_ratio=float(np.round(ratios[0], 6)), point_cut_ratio_all_equal=bool(np.allclose(ratios, ratios[0])),
                      formula_2_minus_2pow=2 - 2 ** (1 - m), eigen_pairs=[(float(a), float(b)) for a, b in pairs],
                      spectral_cut_ceiling=float(np.round(ceil, 6)), plain_m=m)
    return res

if __name__ == "__main__":
    OUT["a_generators"] = check_a()
    print("(a)", json.dumps(OUT["a_generators"]), flush=True)
    OUT["bcd_subspace_frames"] = check_bcd()
    print("(b-d)", json.dumps(OUT["bcd_subspace_frames"]), flush=True)
    OUT["e_cuts"] = check_e()
    print("(e)", json.dumps(OUT["e_cuts"]), flush=True)
    json.dump(OUT, open(__file__.replace(".py", ".json"), "w"), indent=1, default=str)
