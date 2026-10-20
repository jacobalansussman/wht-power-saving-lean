#!/usr/bin/env python3
"""vfield.py (v-rm-singer): independent rebuild of the field side of rm-singer's Theorem S2, Props S3, S4, S5 and
Remark S9.  Written from the statements only (rm-singer's sfield.py was NOT read or imported).

Field GF(2^m) with my own primitive polynomials; elements are ints (bit i = coefficient of a^i).
Coordinates: 'sd' = a self-dual basis found by backtracking (Tr(w_i w_j) = delta_ij); 'poly' = 1, a, ..., a^(m-1).
Forms: Q(c) = sum_{i<=j} q_ij c_i c_j over F_2 in the chosen coordinates.
Target channel k: <Q, S_k> with S_k = s^-k s^-kT (s = matrix of multiplication by a), computed from the MATRIX
(<Q,S> = sum_i q_ii S_ii + sum_{i<j} q_ij S_ij), never from the trace formula.  Move channel k: Q(coords(a^-k)).
"""
import sys, itertools, random
import numpy as np

PRIM = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011}


class GF:
    def __init__(s, m):
        s.m = m; s.q = 1 << m; s.n = s.q - 1
        p = PRIM[m]
        s.exp = [0] * (2 * s.n); s.log = [None] * s.q
        x = 1
        for i in range(s.n):
            s.exp[i] = x; s.log[x] = i
            x <<= 1
            if x & s.q: x ^= p
        assert x == 1 and len(set(s.exp[:s.n])) == s.n, "polynomial not primitive"
        for i in range(s.n, 2 * s.n): s.exp[i] = s.exp[i - s.n]

    def mul(s, x, y):
        if x == 0 or y == 0: return 0
        return s.exp[s.log[x] + s.log[y]]

    def inv(s, x): return s.exp[(s.n - s.log[x]) % s.n]
    def pw(s, k): return s.exp[k % s.n]   # a^k

    def tr(s, x):
        t = 0; y = x
        for _ in range(s.m):
            t ^= y; y = s.mul(y, y)
        assert t in (0, 1)
        return t


def bits(x, m): return np.array([(x >> i) & 1 for i in range(m)], dtype=np.int64)


def rank2(M):
    M = [list(r) for r in (np.array(M) % 2)]
    r = 0; rows = len(M); cols = len(M[0]) if rows else 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i][c]: piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(rows):
            if i != r and M[i][c]: M[i] = [(u ^ v) for u, v in zip(M[i], M[r])]
        r += 1
    return r


def selfdual(F):
    m = F.m
    def rec(chosen):
        if len(chosen) == m: return chosen
        for x in range(1, F.q):
            if F.tr(F.mul(x, x)) != 1: continue
            if any(F.tr(F.mul(x, w)) for w in chosen): continue
            if rank2([bits(w, m) for w in chosen + [x]]) < len(chosen) + 1: continue
            r = rec(chosen + [x])
            if r: return r
        return None
    return rec([])


class Coords:
    """coordinate system: basis b_1..b_m of F over F_2; coords(x) solves x = sum c_i b_i."""
    def __init__(s, F, basis):
        s.F = F; s.basis = basis; m = F.m
        Bm = np.array([bits(b, m) for b in basis]).T % 2   # columns = basis in poly bits
        # brute-force coordinate table
        s.tab = {}
        for c in itertools.product([0, 1], repeat=m):
            x = 0
            for ci, b in zip(c, basis):
                if ci: x ^= b
            s.tab[x] = np.array(c, dtype=np.int64)
        assert len(s.tab) == F.q
        s.mat_of = lambda f: np.array([s.tab[f(b)] for b in basis]).T % 2   # matrix of an F_2-linear map

    def c(s, x): return s.tab[x]


def forms(m):
    pairs = [(i, j) for i in range(m) for j in range(i, m)]
    return pairs


def Qval(q, pairs, c):
    v = 0
    for (i, j), qij in zip(pairs, q):
        if qij: v ^= int(c[i] & c[j])
    return v


def Qpair(q, pairs, S):
    v = 0
    for (i, j), qij in zip(pairs, q):
        if qij: v ^= int(S[i, j] & 1)
    return v


def matpow2(A, k):
    m = A.shape[0]; R = np.eye(m, dtype=np.int64); B = A.copy() % 2
    while k:
        if k & 1: R = (R @ B) % 2
        B = (B @ B) % 2; k >>= 1
    return R


def setup(m, coords):
    F = GF(m)
    if coords == 'sd':
        basis = selfdual(F)
    else:
        basis = [F.pw(i) for i in range(m)]
    C = Coords(F, basis)
    sig = C.mat_of(lambda x: F.mul(F.pw(1), x))
    n = F.n
    siginv = matpow2(sig, n - 1)
    assert ((sig @ siginv) % 2 == np.eye(m, dtype=np.int64)).all()
    Sk = [(matpow2(siginv, k) @ matpow2(siginv, k).T) % 2 for k in range(n)]
    pairs = forms(m)
    return F, basis, C, sig, Sk, pairs


def channel_data(F, C, Sk, pairs, q):
    n = F.n
    K = tuple(Qval(q, pairs, C.c(F.pw(-k))) for k in range(n))     # move channel k: Q(a^-k)
    J = tuple(Qpair(q, pairs, Sk[k]) for k in range(n))           # target channel k from the matrix
    return K, J


def aQ(F, C, pairs, q):
    a = 0
    for y in range(1, F.q):
        if Qval(q, pairs, C.c(y)): a ^= F.inv(y)
    return a


def main():
    ms = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [2, 3, 4, 5]
    full_upto = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    rng = random.Random(12345)
    for m in ms:
        for coords in ('sd', 'poly'):
            F, basis, C, sig, Sk, pairs = setup(m, coords)
            n = F.n; N = len(pairs)
            print(f"\n== m={m} coords={coords} n={n} N={N} basis={basis}")
            sym = bool((sig == sig.T).all())
            print(f"  sigma symmetric: {sym}")
            if coords == 'sd':
                G = np.array([[F.tr(F.mul(u, v)) for v in basis] for u in basis])
                assert (G == np.eye(m, dtype=np.int64)).all()
                Fr = C.mat_of(lambda x: F.mul(x, x))
                orth = bool(((Fr @ Fr.T) % 2 == np.eye(m, dtype=np.int64)).all())
                conj = bool(((Fr @ sig) % 2 == (matpow2(sig, 2) @ Fr) % 2).all())
                print(f"  S9 Frobenius orthogonal: {orth};  phi sigma phi^-1 = sigma^2: {conj}")
            # forms to examine
            if m <= full_upto:
                Qs = list(itertools.product([0, 1], repeat=N))
                mode = 'all'
            else:
                Qs = [tuple(int(i == j) for i in range(N)) for j in range(N)]
                Qs += [tuple(rng.randint(0, 1) for _ in range(N)) for _ in range(400)]
                mode = 'basis+400 random'
            viol_S2b = 0; viol_lin = 0; odd_stab_bad = 0; nodd = 0; Jbad = 0
            seen = {}
            orbit_rows = []
            for q in Qs:
                if not any(q): continue
                K, J = channel_data(F, C, Sk, pairs, q)
                if coords == 'sd':
                    a = aQ(F, C, pairs, q)
                    # S2(b): sum_i Q(d w_i) = Tr(a_Q d) for every d in F*
                    for d in range(1, F.q):
                        lhs = 0
                        for w in basis: lhs ^= Qval(q, pairs, C.c(F.mul(d, w)))
                        if lhs != F.tr(F.mul(a, d)): viol_S2b += 1
                    # S2(c): target channel k (from matrix) = Tr(a_Q a^-k)
                    for k in range(n):
                        if J[k] != F.tr(F.mul(a, F.pw(-k))): viol_lin += 1
                    if a != 0:
                        nodd += 1
                        if sum(J) != 2 ** (m - 1): Jbad += 1
                    else:
                        if sum(J) != 0: Jbad += 1
                # orbit of K under rotation (K as the value sequence of Q along the orbit)
                if mode == 'all':
                    key = min(K[r:] + K[:r] for r in range(n))
                    if key not in seen:
                        size = len(set(K[r:] + K[:r] for r in range(n)))
                        seen[key] = (size, sum(K), sum(J), q)
            if coords == 'sd':
                print(f"  forms checked: {mode}; S2(b) violations {viol_S2b}; target-from-matrix vs Tr(a_Q a^-k)"
                      f" violations {viol_lin}; |J| != 2^(m-1) (odd) or != 0 (even): {Jbad}")
            if mode == 'all':
                # orbits: stabiliser, counts, deletion
                sizes = {}
                for key, (size, wK, wJ, q) in seen.items():
                    sizes.setdefault((size, wJ > 0), []).append(wK)
                print(f"  orbits of nonzero forms: {len(seen)}")
                for (size, odd), ws in sorted(sizes.items()):
                    print(f"    size {size:3d} {'odd ' if odd else 'even'}: count {len(ws)};"
                          f" |K| multiset {sorted(set(ws)) if size == n else sorted(ws)}")
                ratios = [(wJ / wK if wK else float('inf'), wK, wJ) for (size, wK, wJ, q) in seen.values() if wJ > 0]
                best = max(ratios)
                print(f"  deletion bound max |J|/|K| over odd fibres: {best[0]:.6f} (|K|={best[1]}, |J|={best[2]})")
                if coords == 'sd':
                    nodd_orb = sum(1 for (size, wK, wJ, q) in seen.values() if wJ > 0)
                    print(f"  S4: odd orbits {nodd_orb} (claim 2^(N-m) = {2 ** (N - m)}); odd orbits all free:"
                          f" {all(size == n for (size, wK, wJ, q) in seen.values() if wJ > 0)};"
                          f" min odd |K| = {min(wK for (size, wK, wJ, q) in seen.values() if wJ > 0)}"
                          f" (claim 2^(m-2) = {2 ** (m - 2)})")
                    # each odd orbit contains exactly one form with a_Q = 1
                    cnt1 = 0
                    for q in Qs:
                        if any(q) and aQ(F, C, pairs, q) == 1: cnt1 += 1
                    print(f"  S4: forms with a_Q = 1: {cnt1} (= number of odd orbits if one per orbit)")


if __name__ == '__main__':
    main()
