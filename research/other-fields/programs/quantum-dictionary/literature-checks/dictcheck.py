#!/usr/bin/env python3
"""x-lit dictionary checks (small m), to back the literature statements with computation.

Conventions: a vector of F_2^{2m} is (x, z) packed as an int: bits 0..m-1 = x, bits m..2m-1 = z.
Symplectic form w((x,z),(x',z')) = x.z' + z.x'  (mod 2)   [Pauli commutation form].
L0 = {(0,z)} = Z-type Lagrangian = stabilizer group of |0...0> = computational-basis frame.

Checks
 (a) number of Lagrangians = prod_{i=1..m} (2^i+1)                       [Aaronson-Gottesman Prop 2 / Gross Thm 20, mod signs]
 (b) split by r = m - dim(L cap L0): [m r]_2 * 2^{r(r+1)/2}               [Pllaha-Tirkkonen-Calderbank Lemma 1 / Cor 1]
     (r = distance from L0 in the dual polar graph = Hadamard count of the Bruhat cell)
 (c) neighbours at distance 1 of any Lagrangian: 2(2^m - 1)                [Garcia-Markov-Cross Thm 13 counts 4(2^m-1) states]
 (d) subspace frames L_U = {(u, u+v): u in U, v in U^perp}: count = number of subspaces (Galois numbers 2,5,16,67),
     and d(L_U, L_V) = dim U + dim V - 2 dim(U cap V)  (Koetter-Kschischang subspace metric) for all pairs
 (e) real case: maximal totally singular subspaces of Q(x,z) = x.z: count prod_{i=0..m-1}(2^i+1)  (real stabilizer states mod signs)
 (f) 0-1 BFS over Sp(2m,2) with free {CNOT, S} and paid {H_j}: min #H = rank of the B block = d(L0, M L0)
     (Bravyi-Maslov eq. (4) / Maslov-Roetteler Thm 13: Hadamard count of the Bruhat cell; optimality claim of M-R abstract)
"""
import itertools, sys
from collections import deque

def popc(v): return bin(v).count('1')

def symp(a, b, m):
    mask = (1 << m) - 1
    ax, az = a & mask, a >> m
    bx, bz = b & mask, b >> m
    return (popc(ax & bz) + popc(az & bx)) & 1

def span(vecs):
    s = {0}
    for v in vecs:
        s |= {w ^ v for w in s}
    return frozenset(s)

def rank(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)

def gauss_binom(m, r):
    num = den = 1
    for i in range(r):
        num *= (2 ** (m - i) - 1)
        den *= (2 ** (i + 1) - 1)
    return num // den

def lagrangians(m, form):
    """All maximal (dim m) subspaces totally isotropic for `form` (a predicate on pairs, incl. a==b for quadratic)."""
    N = 1 << (2 * m)
    found = set()
    # DFS over isotropic bases in increasing canonical order
    def rec(basis, S):
        if len(basis) == m:
            found.add(S)
            return
        for v in range(1, N):
            if v in S:
                continue
            if not form(v, v):
                continue
            if all(form(v, b) for b in basis):
                S2 = span(list(basis) + [v])
                # canonical: only extend if v is the smallest new element (avoid blowup); cheap dedupe via set
                if S2 in seen:
                    continue
                seen.add(S2)
                rec(basis + [v], S2)
    seen = set()
    rec([], frozenset({0}))
    return found

def main():
    out = []
    for m in (1, 2, 3):
        mask = (1 << m) - 1
        w = lambda a, b: symp(a, b, m) == 0 if a != b else True
        Ls = lagrangians(m, w)
        L0 = span([1 << (m + j) for j in range(m)])
        prod = 1
        for i in range(1, m + 1):
            prod *= (2 ** i + 1)
        dist = lambda A, B: m - (len(A & B).bit_length() - 1)
        hist = {}
        for L in Ls:
            r = dist(L, L0)
            hist[r] = hist.get(r, 0) + 1
        pred = {r: gauss_binom(m, r) * 2 ** (r * (r + 1) // 2) for r in range(m + 1)}
        nb = {L: sum(1 for K in Ls if dist(L, K) == 1) for L in Ls}
        # subspace frames
        subs = set()
        for k in range(m + 1):
            for vecs in itertools.combinations(range(1, 1 << m), k):
                if rank(list(vecs)) == k:
                    subs.add(span(list(vecs)))
        def perp(U):
            return frozenset(v for v in range(1 << m) if all(popc(v & u) % 2 == 0 for u in U))
        def LU(U):
            Up = perp(U)
            return frozenset((u | ((u ^ v) << m)) for u in U for v in Up)
        Lsub = {U: LU(U) for U in subs}
        assert all(L in Ls for L in Lsub.values())
        iso_ok = all(dist(Lsub[U], Lsub[V]) == (len(U).bit_length() - 1) + (len(V).bit_length() - 1) - 2 * (len(U & V).bit_length() - 1)
                     for U in subs for V in subs)
        # real case: Q(x,z) = x.z ; polar form = symplectic form
        Q = lambda v: popc((v & mask) & (v >> m)) & 1
        wr = lambda a, b: (Q(a) == 0) if a == b else (symp(a, b, m) == 0)
        Lr = lagrangians(m, wr)
        prodr = 1
        for i in range(m):
            prodr *= (2 ** i + 1)
        out.append(f"m={m}: #Lagrangians={len(Ls)} (formula {prod}); cells r->count {dict(sorted(hist.items()))} "
                   f"(formula {pred}); neighbours per Lagrangian {sorted(set(nb.values()))} (formula {2*(2**m-1)}); "
                   f"#subspace frames={len(Lsub)}; subspace-metric isometry {iso_ok}; "
                   f"#real (Q=x.z) maximal singular={len(Lr)} (formula {prodr}); "
                   f"subspace frames that are Q-singular={sum(1 for L in Lsub.values() if L in Lr)}")
    # (f) BFS on Sp(2m,2)
    for m in (2, 3):
        n2 = 2 * m
        def mat_mul_left(G, M):
            # matrices as tuple of column images (ints); (G M) e_j = G (M e_j)
            return tuple(apply(G, c) for c in M)
        def apply(G, v):
            r = 0
            j = 0
            while v:
                if v & 1:
                    r ^= G[j]
                v >>= 1
                j += 1
            return r
        I = tuple(1 << j for j in range(n2))
        def gen_cnot(c, t):
            cols = list(I)
            # x_t += x_c : column of x_c gets +x_t ; z_c += z_t : column of z_t gets + z_c
            cols[c] ^= (1 << t)
            cols[m + t] ^= (1 << (m + c))
            return tuple(cols)
        def gen_s(j):
            cols = list(I)
            cols[j] ^= (1 << (m + j))   # z_j += x_j
            return tuple(cols)
        def gen_h(j):
            cols = list(I)
            cols[j], cols[m + j] = cols[m + j], cols[j]
            return tuple(cols)
        free = [gen_cnot(c, t) for c in range(m) for t in range(m) if c != t] + [gen_s(j) for j in range(m)]
        paid = [gen_h(j) for j in range(m)]
        cost = {I: 0}
        dq = deque([I])
        while dq:
            M = dq.popleft()
            cM = cost[M]
            for G in free:
                N_ = mat_mul_left(G, M)
                if N_ not in cost or cost[N_] > cM:
                    cost[N_] = cM
                    dq.appendleft(N_)
            for G in paid:
                N_ = mat_mul_left(G, M)
                if N_ not in cost or cost[N_] > cM + 1:
                    cost[N_] = cM + 1
                    dq.append(N_)
        maskx = (1 << m) - 1
        bad = 0
        hist = {}
        for M, c in cost.items():
            rB = rank([M[m + j] & maskx for j in range(m)])
            if rB != c:
                bad += 1
            hist[c] = hist.get(c, 0) + 1
        order = 2 ** (m * m)
        for i in range(1, m + 1):
            order *= (4 ** i - 1)
        out.append(f"m={m}: |group reached|={len(cost)} (|Sp(2m,2)|={order}); min-#H histogram {dict(sorted(hist.items()))}; "
                   f"mismatches min#H vs rank(B) = {bad}")
    print("\n".join(out))

if __name__ == "__main__":
    main()
