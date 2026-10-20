#!/usr/bin/env python3
"""x-dict: exact checks of the dictionary  engine frames <-> Lagrangians <-> stabilizer states (13 Oct 2026).

Standard library only.  Run:  python3 -I dictcheck.py [sections]   (sections: A B C D E F G; default all)

Conventions (same as ../lib/frames.py).  m qubits = m address bits.  A Pauli vector is an int
v = x | (z << m), standing for X^x Z^z (bit i of x acts on qubit/address bit i); omega(u,v) = x.z' + x'.z.
Engine matrices act on C^(2^m) with (shift(z) f)(x) = f(x+z) (= X^z), dir(z) = ca*1 + cb*shift(z),
ca = (1+i)/2, cb = (1-i)/2, kernel = littleC^(tensor m) (BinaryFrames.lean:69,92,137,152).
The frame label of an invertible S is L(S) = {P : S P S^-1 diagonal} (frames.py; hub general-clifford-frames.tex).

Exactness: Gaussian integers are pairs (re, im) of Python ints; every matrix is an integer matrix times a power of 2
that we track or that the test does not need (all label tests are scale-free).
"""
import sys
import itertools
from math import prod

OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def popc(x):
    return bin(x).count("1")


# ---------------------------------------------------------------------------------------------------------------
# F_2 symplectic helpers
# ---------------------------------------------------------------------------------------------------------------
def omega(u, v, m):
    M = (1 << m) - 1
    return (popc((u & M) & (v >> m)) + popc((v & M) & (u >> m))) & 1


def span(vs):
    s = {0}
    for v in vs:
        if v not in s:
            s |= {x ^ v for x in s}
    return frozenset(s)


def dim2(U):
    return len(U).bit_length() - 1


def basis_of(U):
    b, sp = [], {0}
    for v in sorted(U):
        if v not in sp:
            b.append(v)
            sp |= {x ^ v for x in sp}
    return b


def L0(m):
    return frozenset(z << m for z in range(1 << m))


def LF(m):
    return frozenset(z | (z << m) for z in range(1 << m))


def perp_dot(U, m):
    return frozenset(v for v in range(1 << m) if all(popc(u & v) % 2 == 0 for u in U))


def L_sub(U, m):
    """the hub's / frames.py subspace frame L_U = {(u | u+v) : u in U, v in U^perp}"""
    return frozenset(u | ((u ^ v) << m) for u in U for v in perp_dot(U, m))


def L_css(U, m):
    """the CSS Lagrangian U (+) U^perp = {(u | v)} : X-part in U, Z-part in U^perp"""
    return frozenset(u | (v << m) for u in U for v in perp_dot(U, m))


def all_subspaces(m):
    found, fr = {frozenset([0])}, [frozenset([0])]
    while fr:
        nx = []
        for U in fr:
            for v in range(1, 1 << m):
                if v not in U:
                    V = span(list(U) + [v])
                    if V not in found:
                        found.add(V)
                        nx.append(V)
        fr = nx
    return found


def all_lagrangians_direct(m):
    """every totally isotropic subspace of F_2^(2m), built by isotropic extension; keep the m-dimensional ones"""
    N = 1 << (2 * m)
    level = {frozenset([0])}
    for _ in range(m):
        nxt = set()
        for U in level:
            for v in range(1, N):
                if v not in U and all(omega(v, u, m) == 0 for u in U):
                    nxt.add(span(list(U) + [v]))
        level = nxt
    return level


def gl_count(m):
    return prod((1 << m) - (1 << i) for i in range(m))


# symplectic action of gates on Pauli vectors (Heisenberg picture, signs dropped)
def sym_H(v, j, m):
    x, z = v & ((1 << m) - 1), v >> m
    bx, bz = (x >> j) & 1, (z >> j) & 1
    x = (x & ~(1 << j)) | (bz << j)
    z = (z & ~(1 << j)) | (bx << j)
    return x | (z << m)


def sym_S(v, j, m):       # S X S^-1 = Y : (1|0) -> (1|1)
    x, z = v & ((1 << m) - 1), v >> m
    z ^= ((x >> j) & 1) << j
    return x | (z << m)


def sym_CX(v, c, t, m):   # CNOT c->t : X_c -> X_c X_t, Z_t -> Z_c Z_t
    x, z = v & ((1 << m) - 1), v >> m
    x ^= ((x >> c) & 1) << t
    z ^= ((z >> t) & 1) << c
    return x | (z << m)


# ---------------------------------------------------------------------------------------------------------------
# Gaussian-integer vectors and matrices
# ---------------------------------------------------------------------------------------------------------------
def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def gconj(a):
    return (a[0], -a[1])


Z0 = (0, 0)
ONE = (1, 0)
IU = (0, 1)
IPOW = [(1, 0), (0, 1), (-1, 0), (0, -1)]


def gdiv_exact(a, b):
    n = b[0] * b[0] + b[1] * b[1]
    p = gmul(a, gconj(b))
    assert p[0] % n == 0 and p[1] % n == 0, ("inexact", a, b)
    return (p[0] // n, p[1] // n)


def mmul(A, B):
    n = len(A)
    out = []
    for i in range(n):
        Ai = A[i]
        row = []
        for j in range(n):
            re = im = 0
            for t in range(n):
                a = Ai[t]
                if a[0] or a[1]:
                    b = B[t][j]
                    re += a[0] * b[0] - a[1] * b[1]
                    im += a[0] * b[1] + a[1] * b[0]
            row.append((re, im))
        out.append(row)
    return out


def mdag(A):
    n = len(A)
    return [[gconj(A[j][i]) for j in range(n)] for i in range(n)]


def mscale(A, c):
    return [[gmul(c, e) for e in r] for r in A]


def ident(n):
    return [[ONE if i == j else Z0 for j in range(n)] for i in range(n)]


def kron(A, B):
    return [[gmul(A[i // len(B)][j // len(B)], B[i % len(B)][j % len(B)]) for j in range(len(A) * len(B))]
            for i in range(len(A) * len(B))]


def shiftM(z, m):            # (shift z) x y = 1 iff x + z = y   (BinaryFrames.lean:137)
    n = 1 << m
    return [[ONE if (x ^ z) == y else Z0 for y in range(n)] for x in range(n)]


def dir2M(z, m):             # 2 * dir(z) = (1+i) 1 + (1-i) shift(z)
    n = 1 << m
    return [[(1, 1) if x == y else ((1, -1) if (x ^ z) == y else Z0) for y in range(n)] for x in range(n)]


def ZM(a, m):                # Z^a = diag (-1)^{a.x}
    n = 1 << m
    return [[((-1, 0) if popc(a & x) & 1 else ONE) if x == y else Z0 for y in range(n)] for x in range(n)]


def XZM(v, m):               # X^x Z^z as a matrix: (X^x Z^z f)(y) = (-1)^{z.(y+x)} f(y+x)
    x, z = v & ((1 << m) - 1), v >> m
    return mmul(shiftM(x, m), ZM(z, m))


def permM(Gimg, m):          # Letters.permMat: (permMat G) x y = 1 iff G x = y ; G given by images of basis vectors
    n = 1 << m

    def G(x):
        r = 0
        for i in range(m):
            if (x >> i) & 1:
                r ^= Gimg[i]
        return r
    return [[ONE if G(x) == y else Z0 for y in range(n)] for x in range(n)]


def phaseM(z, c, m):         # Letters.phaseMat: diag i^[z.x + c]  ([.] in {0,1})
    n = 1 << m
    return [[(IU if ((popc(z & x) + c) & 1) else ONE) if x == y else Z0 for y in range(n)] for x in range(n)]


def littleC2():              # 2 * littleC = [[1+i, 1-i], [1-i, 1+i]]
    return [[(1, 1), (1, -1)], [(1, -1), (1, 1)]]


def littleH():
    return [[ONE, ONE], [ONE, (-1, 0)]]


def kernel_scaled(m):        # 2^m * kernel
    K = [[ONE]]
    for _ in range(m):
        K = kron(littleC2(), K)
    return K


def tensor_power(A, m):
    K = [[ONE]]
    for _ in range(m):
        K = kron(A, K)
    return K


def on_qubit(A2, j, m):      # 2x2 matrix A2 acting on address bit j (bit j of x)
    n = 1 << m
    out = [[Z0] * n for _ in range(n)]
    for x in range(n):
        for y in range(n):
            if (x ^ y) & ~(1 << j) == 0:
                out[x][y] = A2[(x >> j) & 1][(y >> j) & 1]
    return out


def pauli_of(Q, m):
    """Q = c * X^x Z^z (any nonzero c)  ->  v = x | z << m ; else None"""
    n = 1 << m
    nz = [j for j in range(n) if Q[0][j] != Z0]
    if len(nz) != 1:
        return None
    x = nz[0]
    c = Q[0][x]
    z = 0
    for i in range(n):
        for j in range(n):
            if j != (i ^ x):
                if Q[i][j] != Z0:
                    return None
        e = Q[i][i ^ x]
        if e == c:
            s = 0
        elif e == (-c[0], -c[1]):
            s = 1
        else:
            return None
        # sign at row i is (-1)^{z.(i+x)} relative to row 0 sign (-1)^{z.x}
        if i and (i & (i - 1)) == 0:
            z |= s << (i.bit_length() - 1)
    # verify
    for i in range(n):
        want = c if (popc(z & i) & 1) == 0 else (-c[0], -c[1])
        if Q[i][i ^ x] != want:
            return None
    return x | (z << m)


def label_of(S, m):
    """L(S) = { P : S P S^-1 diagonal } = span of S^-1 Z^a S ; for S unitary up to scale, S^-1 ~ S^dagger"""
    Sd = mdag(S)
    vs = []
    for j in range(m):
        P = pauli_of(mmul(mmul(Sd, ZM(1 << j, m)), S), m)
        assert P is not None, "not Clifford"
        vs.append(P)
    L = span(vs)
    assert len(L) == 1 << m and all(omega(a, b, m) == 0 for a in vs for b in vs)
    return L


def is_clifford_mat(S, m):
    Sd = mdag(S)
    for j in range(2 * m):
        if pauli_of(mmul(mmul(S, XZM(1 << j, m)), Sd), m) is None:
            return False
    return True


def to_monomial(M):
    """M = t * (monomial matrix with entries i^e), t > 0 common  ->  (perm tuple, exps tuple) else None"""
    pi, ex, t = [], [], None
    for row in M:
        nz = [j for j in range(len(row)) if row[j] != Z0]
        if len(nz) != 1:
            return None
        re, im = row[nz[0]]
        if re and im:
            return None
        a = abs(re) + abs(im)
        if t is None:
            t = a
        if a != t:
            return None
        pi.append(nz[0])
        ex.append(0 if re > 0 else 2 if re < 0 else 1 if im > 0 else 3)
    return tuple(pi), tuple(ex)


def monomial_is_free(mono, m):
    """the closed form of the engine's Free group: perm x -> pi(x) affine over F_2 and exponents
    e(x) = c + sum a_j x_j + 2 sum_{j<k} b_jk x_j x_k (mod 4)  (row x has its entry i^e(x) in column pi(x))"""
    pi, ex = mono
    n = len(pi)
    b0 = pi[0]
    lin = [pi[1 << j] ^ b0 for j in range(m)]
    for x in range(n):
        y = b0
        for j in range(m):
            if (x >> j) & 1:
                y ^= lin[j]
        if pi[x] != y:
            return False
    if len(set(lin)) != m or 0 in lin or len(span(lin)) != n:
        return False
    c = ex[0]
    a = [(ex[1 << j] - c) % 4 for j in range(m)]
    b = {}
    for j in range(m):
        for k in range(j + 1, m):
            d = (ex[(1 << j) | (1 << k)] - c - a[j] - a[k]) % 4
            if d % 2:
                return False
            b[(j, k)] = d // 2
    for x in range(n):
        e = c + sum(a[j] for j in range(m) if (x >> j) & 1) + 2 * sum(b[(j, k)] for (j, k) in b
                                                                       if (x >> j) & 1 and (x >> k) & 1)
        if e % 4 != ex[x]:
            return False
    return True


# ---------------------------------------------------------------------------------------------------------------
# A. Lagrangians and stabilizer states
# ---------------------------------------------------------------------------------------------------------------
def stabilizer_states(m):
    """BFS from |0...0> under H_j, S_j, CNOT_jk.  A state is a tuple of exponents (None = 0 amplitude, e = i^e),
    normalised so that the first non-zero amplitude is 1.  Each state carries the Lagrangian obtained by applying
    the same gates' symplectic actions to L_0; returns {state: Lagrangian}."""
    n = 1 << m

    def norm(vec):  # vec: list of Gaussian ints
        c = next(e for e in vec if e != Z0)
        out = []
        for e in vec:
            if e == Z0:
                out.append(None)
            else:
                q = gdiv_exact(e, c)
                assert q in IPOW, ("amplitudes not flat", vec)
                out.append(IPOW.index(q))
        return tuple(out)

    def to_vec(st):
        return [Z0 if e is None else IPOW[e] for e in st]

    gates = []
    for j in range(m):
        def gH(vec, j=j):
            out = list(vec)
            for x in range(n):
                if not (x >> j) & 1:
                    a, b = vec[x], vec[x | (1 << j)]
                    out[x], out[x | (1 << j)] = gadd(a, b), (a[0] - b[0], a[1] - b[1])
            return out
        gates.append((gH, lambda v, j=j: sym_H(v, j, m)))

        def gS(vec, j=j):
            return [gmul(IU, e) if (x >> j) & 1 else e for x, e in enumerate(vec)]
        gates.append((gS, lambda v, j=j: sym_S(v, j, m)))
    for c in range(m):
        for t in range(m):
            if c != t:
                def gC(vec, c=c, t=t):
                    return [vec[x ^ (((x >> c) & 1) << t)] for x in range(n)]
                gates.append((gC, lambda v, c=c, t=t: sym_CX(v, c, t, m)))
    s0 = norm([ONE] + [Z0] * (n - 1))
    lab = {s0: L0(m)}
    fr = [s0]
    while fr:
        nx = []
        for s in fr:
            vec = to_vec(s)
            for g, sg in gates:
                t = norm(g(vec))
                Lt = frozenset(sg(v) for v in lab[s])
                if t in lab:
                    assert lab[t] == Lt, "Lagrangian of a state depends on the path"
                else:
                    lab[t] = Lt
                    nx.append(t)
        fr = nx
    return lab


def stabilizes(st, v, m):
    """is the Hermitian Pauli i^{x.z} X^x Z^z a +-1 eigen-operator of the state?"""
    n = 1 << m
    x, z = v & ((1 << m) - 1), v >> m
    ph = popc(x & z) % 4
    sgn = None
    for y in range(n):
        a = st[y ^ x]
        b = st[y]
        if (a is None) != (b is None):
            return False
        if a is None:
            continue
        e = (a + ph + 2 * (popc(z & (y ^ x)) & 1)) % 4    # amplitude exponent of P psi at y
        d = (e - b) % 4
        if d % 2:
            return False
        if sgn is None:
            sgn = d
        elif d != sgn:
            return False
    return True


def section_A():
    say("== A. Lagrangians of F_2^(2m)  <->  stabilizer states of m qubits (mod phase and Pauli)")
    res = {}
    for m in (1, 2, 3, 4):
        want = prod((1 << i) + 1 for i in range(1, m + 1))
        Ls = all_lagrangians_direct(m)
        lab = stabilizer_states(m)
        nstates = len(lab)
        fib = {}
        for s, L in lab.items():
            fib.setdefault(L, []).append(s)
        ok_fib = all(len(v) == 1 << m for v in fib.values())
        onto = set(fib) == set(Ls)
        # direct check: the generators of the tracked Lagrangian stabilise the state up to sign
        bad = 0
        for s, L in lab.items():
            for v in basis_of(L):
                if not stabilizes(s, v, m):
                    bad += 1
        # a fibre is an orthonormal basis (m <= 3: every pair; m = 4: the fibre of L_0 and 30 others)
        orth_bad = 0
        fibs = list(fib.values()) if m <= 3 else list(fib.values())[:31] + [fib[L0(m)]]
        for F in fibs:
            vecs = [[Z0 if e is None else IPOW[e] for e in s] for s in F]
            for a, b in itertools.combinations(vecs, 2):
                ip = (0, 0)
                for p, q in zip(a, b):
                    ip = gadd(ip, gmul(gconj(p), q))
                if ip != Z0:
                    orth_bad += 1
        AG = (1 << m) * want
        say("  m=%d: Lagrangians (direct isotropic enumeration) %d, prod(2^i+1) = %d; stabilizer states by BFS %d,"
            " Aaronson-Gottesman 2^m*prod = %d; every Lagrangian has exactly 2^m states: %s; onto: %s;"
            " generator-stabilises failures %d; non-orthogonal pairs inside fibres %d (%d fibres checked)"
            % (m, len(Ls), want, nstates, AG, ok_fib, onto, bad, orth_bad, len(fibs)))
        assert len(Ls) == want and nstates == AG and ok_fib and onto and bad == 0 and orth_bad == 0
        res[m] = (len(Ls), nstates)
    return res


# ---------------------------------------------------------------------------------------------------------------
# B. CSS / subspace frames
# ---------------------------------------------------------------------------------------------------------------
def section_B():
    say("== B. subspace frames L_U  <->  CSS Lagrangians U (+) U^perp, twisted by S on every qubit")
    for m in (1, 2, 3, 4):
        subs = all_subspaces(m)
        Lsub = {L_sub(U, m) for U in subs}
        Lcss = {L_css(U, m) for U in subs}
        Ls = all_lagrangians_direct(m)
        # CSS-type Lagrangians = those spanned by their pure-X and pure-Z elements
        X = lambda v: (v >> m) == 0
        Zt = lambda v: (v & ((1 << m) - 1)) == 0
        css_type = {L for L in Ls if len({a ^ b for a in L if X(a) for b in L if Zt(b)}) == len(L)}
        # (Y,Z)-type: spanned by elements (u|u) and (0|v)
        Yt = lambda v: (v & ((1 << m) - 1)) == (v >> m)
        yz_type = {L for L in Ls if len({a ^ b for a in L if Yt(a) for b in L if Zt(b)}) == len(L)}
        Sall = lambda v: v ^ ((v & ((1 << m) - 1)) << m)      # symplectic action of S on every qubit
        twist_ok = all(frozenset(Sall(v) for v in L_css(U, m)) == L_sub(U, m) for U in subs)
        Fk = lambda v: ((v & ((1 << m) - 1)) ^ (v >> m)) | ((v >> m) << m)   # kernel (x|z) -> (x+z|z)
        kern_ok = all(frozenset(Fk(v) for v in L_sub(U, m)) == L_sub(perp_dot(U, m), m) for U in subs)
        Hall = lambda v: (v >> m) | ((v & ((1 << m) - 1)) << m)            # H on every qubit swaps x, z
        hk_ok = all(frozenset(Hall(v) for v in L_css(U, m)) == L_css(perp_dot(U, m), m) for U in subs)
        all_lag = all(L in Ls for L in Lsub | Lcss)
        say("  m=%d: subspaces %d; L_U distinct %d, CSS U+U^perp distinct %d; CSS-type Lagrangians (X/Z-spanned) %d;"
            " (Y,Z)-spanned %d = L_U set: %s; L_U = S^(x m)(U + U^perp): %s; kernel F: L_U -> L_(U^perp): %s;"
            " H^(x m): U+U^perp -> U^perp+U: %s; all Lagrangian: %s"
            % (m, len(subs), len(Lsub), len(Lcss), len(css_type), len(yz_type), yz_type == Lsub, twist_ok, kern_ok,
               hk_ok, all_lag))
        assert len(Lsub) == len(subs) == len(Lcss) == len(css_type) and css_type == Lcss and yz_type == Lsub
        assert twist_ok and kern_ok and hk_ok and all_lag
    say("  start L_0 = L_{0} (all Z: |0..0>), end L_F = L_{F_2^m} (all Y-type: S^(x m)|+..+>), and the third m=1 frame"
        " <X> (the Hadamard frame, |+>) is not a subspace frame")


# ---------------------------------------------------------------------------------------------------------------
# C. free steps = Hadamard-free Clifford group ; moves = Hadamards
# ---------------------------------------------------------------------------------------------------------------
def free_closure(m):
    """closure of the engine's Free generators (Free.lean:27-32) as monomial (perm, exps) pairs"""
    n = 1 << m
    gens = []
    # permMat over GL(m,2) generators: transvections x_i += x_j and swaps
    for i in range(m):
        for j in range(m):
            if i != j:
                img = [(1 << k) | ((1 << j) if k == i else 0) for k in range(m)]
                gens.append(permM(img, m))
    for z in range(1 << m):
        for c in (0, 1):
            if z or c:
                gens.append(phaseM(z, c, m))
    for z in range(1, 1 << m):
        gens.append(shiftM(z, m))
    cg = [to_monomial(g) for g in gens]
    mul = lambda a, b: (tuple(b[0][a[0][x]] for x in range(n)), tuple((a[1][x] + b[1][a[0][x]]) % 4 for x in range(n)))
    e = (tuple(range(n)), (0,) * n)
    seen, fr = {e}, [e]
    while fr:
        nx = []
        for a in fr:
            for g in cg:
                b = mul(a, g)
                if b not in seen:
                    seen.add(b)
                    nx.append(b)
        fr = nx
    return seen


def mono_matrix(mono):
    pi, ex = mono
    n = len(pi)
    return [[IPOW[ex[x]] if pi[x] == y else Z0 for y in range(n)] for x in range(n)]


def sp_order_generated(gen_maps, m):
    """order of the subgroup of Sp(2m,2) generated by maps given as tables on basis vectors (tuple of 2m images)"""
    N = 2 * m

    def comp(a, b):     # (a o b)(e_k) = a(b(e_k))
        out = []
        for k in range(N):
            v = b[k]
            r = 0
            for t in range(N):
                if (v >> t) & 1:
                    r ^= a[t]
            out.append(r)
        return tuple(out)
    e = tuple(1 << k for k in range(N))
    seen, fr = {e}, [e]
    while fr:
        nx = []
        for a in fr:
            for g in gen_maps:
                b = comp(g, a)
                if b not in seen:
                    seen.add(b)
                    nx.append(b)
        fr = nx
    return seen


def section_C():
    say("== C. free steps = the Hadamard-free Clifford group (stabiliser of L_0) ; paid move = one Hadamard")
    for m in (1, 2):
        cl = free_closure(m)
        want = 4 * (4 ** m) * gl_count(m) * (2 ** (m * (m + 1) // 2))
        n = 1 << m
        # every monomial {+-1, +-i} matrix with an affine permutation (all perms at m <= 2): which are Clifford?
        cliff_mono = set()
        for pi in itertools.permutations(range(n)):
            for ex in itertools.product(range(4), repeat=n):
                M = mono_matrix((pi, ex))
                if is_clifford_mat(M, m):
                    cliff_mono.add((pi, ex))
        closed_ok = all(monomial_is_free(a, m) for a in cl)
        say("  m=%d: closure of Free generators has %d elements, formula 4*4^m*|GL(m,2)|*2^(m(m+1)/2) = %d;"
            " monomial i-power matrices that are Clifford: %d; equal to the closure: %s; closed form holds: %s"
            % (m, len(cl), want, len(cliff_mono), cliff_mono == cl, closed_ok))
        assert len(cl) == want and cliff_mono == cl and closed_ok
    # m = 3: phase subgroup and AGL separately (the full group is their semidirect product, order 2,752,512)
    m = 3
    n = 8
    ph = set()
    gens = [tuple((1 if ((popc(z & x) + c) & 1) else 0) for x in range(n)) for z in range(n) for c in (0, 1) if z or c]
    e = (0,) * n
    ph, fr = {e}, [e]
    while fr:
        nx = []
        for a in fr:
            for g in gens:
                b = tuple((p + q) % 4 for p, q in zip(a, g))
                if b not in ph:
                    ph.add(b)
                    nx.append(b)
        fr = nx
    say("  m=3: diagonal phases generated by phaseMat: %d = 4*4^3*2^3 = %d; every one of the closed form: %s"
        % (len(ph), 4 * 64 * 8, all(monomial_is_free((tuple(range(n)), p), 3) for p in ph)))
    assert len(ph) == 2048
    # symplectic images: the Free generators generate exactly Stab_Sp(L_0); with one direction, all of Sp(2m,2)
    for m in (1, 2, 3):
        N = 2 * m
        tab = lambda f: tuple(f(1 << k) for k in range(N))
        fg = []
        for i in range(m):
            for j in range(m):
                if i != j:
                    fg.append(tab(lambda v, i=i, j=j: sym_CX(v, j, i, m)))
            fg.append(tab(lambda v, i=i: sym_S(v, i, m)))
        for i in range(m):
            for j in range(i + 1, m):     # CZ: X_i -> X_i Z_j, X_j -> X_j Z_i
                def cz(v, i=i, j=j):
                    x, z = v & ((1 << m) - 1), v >> m
                    z ^= (((x >> i) & 1) << j) ^ (((x >> j) & 1) << i)
                    return x | (z << m)
                fg.append(tab(cz))
        G = sp_order_generated(fg, m)
        stab_ok = all(frozenset(_apply(g, v, N) for v in L0(m)) == L0(m) for g in G)
        want = gl_count(m) * (2 ** (m * (m + 1) // 2))
        # add one H (= one direction up to S): the whole Sp(2m,2)
        hg = fg + [tab(lambda v: sym_H(v, 0, m))]
        spo = None
        if m <= 3:
            spo = len(sp_order_generated(hg, m))
        sp_want = (2 ** (m * m)) * prod((4 ** i) - 1 for i in range(1, m + 1))
        say("  m=%d: symplectic image of Free = %d elements = |GL(m,2)|*2^(m(m+1)/2) = %d, fixes L_0: %s;"
            " with one Hadamard added: %s = |Sp(2m,2)| = %d" % (m, len(G), want, stab_ok, spo, sp_want))
        assert len(G) == want and stab_ok and spo == sp_want
    # moves versus Hadamards, explicit matrices
    for m in (1, 2, 3):
        n = 1 << m
        ok1 = ok2 = ok3 = ok4 = True
        Sq = [[ONE, Z0], [Z0, IU]]
        for j in range(m):
            # dir(e_j) = littleC on bit j = H S H = SX (exactly, scale 2 tracked)
            D = dir2M(1 << j, m)
            ok1 &= D == on_qubit(littleC2(), j, m)
            # S_j * dir(e_j) * S_j = ((1+i)/2) littleH on bit j, i.e. e^{i pi/4} H_j  (hub: (1+i)/2 H = S C S)
            Sj = on_qubit(Sq, j, m)
            ok2 &= mmul(mmul(Sj, D), Sj) == on_qubit([[(1, 1), (1, 1)], [(1, 1), (-1, -1)]], j, m)
        # dir(z) = P_G dir(e_0) P_G^-1 with G e_0 = z (an address map, i.e. a CNOT circuit)
        for z in range(1, n):
            low = (z & -z).bit_length() - 1
            # G: e_low -> z, other basis vectors fixed (invertible since z has bit low)
            img = [(z if k == low else (1 << k)) for k in range(m)]
            P = permM(img, m)
            Pinv = mdag(P)
            ok3 &= dir2M(z, m) == mmul(mmul(Pinv, dir2M(1 << low, m)), P) or \
                dir2M(z, m) == mmul(mmul(P, dir2M(1 << low, m)), Pinv)
        # a block of k independent directions = CNOT-conjugate of SX on k qubits
        for k in range(1, m + 1):
            for zs in itertools.combinations(range(1, n), k):
                if len(span(zs)) != 1 << k:
                    continue
                B = ident(n)
                for z in zs:
                    B = mmul(dir2M(z, m), B)
                # find G with G e_t = zs[t] for t < k, completed to a basis
                basis = list(zs)
                for t in range(m):
                    if len(span(basis + [1 << t])) > len(span(basis)):
                        basis.append(1 << t)
                P = permM(basis, m)
                SXk = ident(n)
                for t in range(k):
                    SXk = mmul(dir2M(1 << t, m), SXk)
                ok4 &= B == mmul(mmul(mdag(P), SXk), P) or B == mmul(mmul(P, SXk), mdag(P))
        say("  m=%d: dir(e_j) = SX_j = H S H exactly: %s; S_j dir(e_j) S_j = e^(i pi/4) H_j exactly: %s;"
            " dir(z) = (address map) dir(e_low) (address map)^-1: %s; every block of k independent directions ="
            " (address map) SX^(x k) (address map)^-1: %s" % (m, ok1, ok2, ok3, ok4))
        assert ok1 and ok2 and ok3 and ok4


def _apply(g, v, N):
    r = 0
    for k in range(N):
        if (v >> k) & 1:
            r ^= g[k]
    return r


# ---------------------------------------------------------------------------------------------------------------
# D. distance
# ---------------------------------------------------------------------------------------------------------------
def lag_graph_bfs(m, sources):
    Ls = sorted(all_lagrangians_direct(m), key=lambda s: sorted(s))
    idx = {L: i for i, L in enumerate(Ls)}
    n = len(Ls)
    # neighbours by one transvection tau_w(v) = v + omega(v,w) w  (the symplectic action of a pi/4 Pauli rotation)
    adj = [set() for _ in range(n)]
    N = 1 << (2 * m)
    for i, L in enumerate(Ls):
        for w in range(1, N):
            if w in L:
                continue
            Lw = frozenset(v ^ (omega(v, w, m) * w) for v in L)
            adj[i].add(idx[Lw])
    out = {}
    for s in sources:
        si = idx[s]
        d = {si: 0}
        fr = [si]
        while fr:
            nx = []
            for u in fr:
                for v in adj[u]:
                    if v not in d:
                        d[v] = d[u] + 1
                        nx.append(v)
            fr = nx
        out[s] = {Ls[j]: dd for j, dd in d.items()}
    return Ls, adj, out


def section_D():
    say("== D. frame distance d(A,B) = m - dim(A cap B)")
    for m in (1, 2, 3, 4):
        Ls = all_lagrangians_direct(m)
        srcs = list(Ls) if m <= 3 else [L0(m), LF(m)]
        Lsorted, adj, dist = lag_graph_bfs(m, srcs)
        ok = all(dist[s][t] == m - dim2(s & t) for s in srcs for t in Ls)
        nb = {len(a) for a in adj}
        edges_ok = all(dim2(Lsorted[i] & Lsorted[j]) == m - 1 for i in range(len(Lsorted)) for j in adj[i])
        # every Lagrangian at intersection dim m-1 is reached by one transvection
        cnt_nb = {sum(1 for t in Ls if dim2(s & t) == m - 1) for s in srcs}
        subs = all_subspaces(m)
        sub_ok = all(m - dim2(L_sub(U, m) & L_sub(V, m)) == dim2(U) + dim2(V) - 2 * dim2(U & V)
                     for U in subs for V in subs)
        say("  m=%d: graph 'one pi/4 Pauli rotation' (transvection) has valency %s (dual polar graph Sp(2m,2):"
            " 2(2^m-1) = %d), edges have dim(A cap B) = m-1: %s, Lagrangians with dim(A cap B)=m-1 per vertex %s;"
            " BFS distance = m - dim(A cap B) for %d sources x all: %s; on subspace frames"
            " = dim U + dim V - 2 dim(U cap V) (Koetter-Kschischang subspace distance): %s"
            % (m, nb, 2 * ((1 << m) - 1), edges_ok, cnt_nb, len(srcs), ok, sub_ok))
        assert ok and edges_ok and sub_ok and nb == {2 * ((1 << m) - 1)}


def min_H_count(m):
    """0-1 BFS over Sp(2m,2): H_j costs 1; S_j, CNOT, (and CZ) cost 0.  Returns {symplectic map: min #H}."""
    N = 2 * m
    tab = lambda f: tuple(f(1 << k) for k in range(N))
    free, paid = [], []
    for i in range(m):
        free.append(tab(lambda v, i=i: sym_S(v, i, m)))
        paid.append(tab(lambda v, i=i: sym_H(v, i, m)))
        for j in range(m):
            if i != j:
                free.append(tab(lambda v, i=i, j=j: sym_CX(v, i, j, m)))

    def comp(a, b):
        return tuple(_apply(a, b[k], N) for k in range(N))
    from collections import deque
    e = tuple(1 << k for k in range(N))
    dist = {e: 0}
    dq = deque([e])
    while dq:
        a = dq.popleft()
        da = dist[a]
        for g in free:
            b = comp(g, a)
            if b not in dist or dist[b] > da:
                dist[b] = da
                dq.appendleft(b)
        for g in paid:
            b = comp(g, a)
            if b not in dist:
                dist[b] = da + 1
                dq.append(b)
    return dist


def rank2(rows):
    p = {}
    r = 0
    for x in rows:
        while x:
            b = x.bit_length() - 1
            if b not in p:
                p[b] = x
                r += 1
                break
            x ^= p[b]
    return r


def section_D2():
    say("== D2. minimal Hadamard count of a Clifford operator = rank of its Z->X block = d(L_0, L(T)) = Fourier rank")
    for m in (1, 2, 3):
        dist = min_H_count(m)
        M = (1 << m) - 1
        ok = True
        hist = {}
        for g, h in dist.items():
            # B block: X-parts of the images of Z_k  (g[m+k] is the image of Z_k)
            B = [g[m + k] & M for k in range(m)]
            rb = rank2(B)
            # label of T: T^-1 L_0 ; frame distance from L_0 = m - dim(L_0 cap g(L_0)) as well
            img = frozenset(_apply(g, v, 2 * m) for v in L0(m))
            d = m - dim2(img & L0(m))
            ok &= (rb == h == d)
            hist[h] = hist.get(h, 0) + 1
        sp_want = (2 ** (m * m)) * prod((4 ** i) - 1 for i in range(1, m + 1))
        say("  m=%d: |Sp(2m,2)| reached %d (= %d); min #H = rank B = m - dim(L_0 cap T L_0) on every element: %s;"
            " histogram of min #H %s" % (m, len(dist), sp_want, ok, dict(sorted(hist.items()))))
        assert ok and len(dist) == sp_want


def section_D3():
    say("== D3. frame distance = overlap exponent of stabilizer states:  max |<psi|phi>|^2 = 2^-d(A,B)")
    for m in (1, 2, 3, 4):
        lab = stabilizer_states(m)
        states = list(lab.items())
        n = 1 << m
        vec = {s: [Z0 if e is None else IPOW[e] for e in s] for s in lab}
        nn = {s: sum(1 for e in s if e is not None) for s in lab}
        pairs = itertools.combinations(states, 2) if m <= 3 else ((a, b) for a in states[:1] for b in states)
        ok = True
        cnt = 0
        best = {}
        for (s, A), (t, B) in pairs:
            ip = (0, 0)
            for p, q in zip(vec[s], vec[t]):
                if p != Z0 and q != Z0:
                    ip = gadd(ip, gmul(gconj(p), q))
            num = ip[0] ** 2 + ip[1] ** 2            # |<s|t>|^2 * nn[s] * nn[t]
            d = m - dim2(A & B)
            # |<psi|phi>|^2 in {0, 2^-d}
            if num:
                ok &= num * (1 << d) == nn[s] * nn[t]
            if A != B:
                key = frozenset((A, B))
                best[key] = best.get(key, False) or bool(num)
            cnt += 1
        some_nonzero = all(best.values())
        say("  m=%d: %d pairs; every overlap |<psi|phi>|^2 is 0 or 2^-(m - dim(A cap B)): %s;"
            " for every pair of distinct Lagrangians some sign choice gives the non-zero value: %s"
            " (%d Lagrangian pairs; equal Lagrangians: distinct states are orthogonal, section A)"
            % (m, cnt, ok, some_nonzero, len(best)))
        assert ok and some_nonzero
    say("  caveat (state picture): S maps |+> (frame <X>) to |+i> (frame <Y>) with no Hadamard, but d(<X>,<Y>) = 1;"
        " so d is NOT the H-count of a state-to-state circuit whose free gates act on the state; it is the H-count of"
        " the transition operator S' S^-1 (D2) and the overlap exponent (D3)")


# ---------------------------------------------------------------------------------------------------------------
# E. kernel, WHT, DFT pieces, several arrays
# ---------------------------------------------------------------------------------------------------------------
def section_E():
    say("== E. kernel = SX^(x m); WHT = H^(x m) up to free S layers; DFT pieces")
    for m in (1, 2, 3):
        n = 1 << m
        K2 = kernel_scaled(m)                       # 2^m K
        W = tensor_power(littleH(), m)              # unnormalised WHT, W^2 = 2^m
        # W K W^-1 = diag(i^{wt x}) :  W K2 W = 2^m * 2^m diag(i^wt)   (W^-1 = W / 2^m)
        lhs = mmul(mmul(W, K2), W)
        D = [[gmul((1 << (2 * m), 0), IPOW[popc(x) % 4]) if x == y else Z0 for y in range(n)] for x in range(n)]
        ok1 = lhs == D
        # S^(x m) K S^(x m) = e^{i m pi/4} H^(x m):  scaled: Sm K2 Sm = (1+i)^m W
        Sm = tensor_power([[ONE, Z0], [Z0, IU]], m)
        c = (1, 0)
        for _ in range(m):
            c = gmul(c, (1, 1))
        ok2 = mmul(mmul(Sm, K2), Sm) == mscale(W, c)
        # the kernel's label is L_F = {(z|z)}; the WHT's label is the X frame {(x|0)}
        Kl = label_of(K2, m)
        Wl = label_of(W, m)
        ok3 = Kl == LF(m) and Wl == frozenset(range(1 << m))
        # is the kernel Clifford, is the WHT Clifford
        say("  m=%d: W K W^-1 = diag(i^wt(x)) = S^(x m): %s; S^(x m) K S^(x m) = e^(i m pi/4) H^(x m): %s;"
            " L(K) = L_F = {(z|z)}, L(WHT) = {(x|0)} (X frame): %s; d(L_0, L(WHT)) = %d"
            % (m, ok1, ok2, ok3, m - dim2(Wl & L0(m))))
        assert ok1 and ok2 and ok3
    # DFT_n with root w = e^{2 pi i/n}: Clifford only for n <= 2 among n = 2, 4 (n = 8 leaves Q(i))
    F2 = [[ONE, ONE], [ONE, (-1, 0)]]
    F4 = [[IPOW[(j * k) % 4] for k in range(4)] for j in range(4)]    # w = i
    say("  DFT_2 is Clifford: %s; DFT_4 (root i) is Clifford: %s; DFT_8 needs e^(i pi/4), outside Q(i) (not checked)"
        % (is_clifford_mat(F2, 1), is_clifford_mat(F4, 2)))
    # DFT_4 = (H on bit 1) CS (H on bit 0) SWAP-type: the twiddle diag(1,1,1,i) = controlled-S is not Clifford
    CS = [[ONE if (x == y and x != 3) else (IU if x == y == 3 else Z0) for y in range(4)] for x in range(4)]
    say("  the radix-2 twiddle diag(1,1,1,i) (= controlled-S) is Clifford: %s; CZ = diag(1,1,1,-1) is Clifford: %s"
        % (is_clifford_mat(CS, 2), is_clifford_mat([[ONE if x == y and x != 3 else ((-1, 0) if x == y == 3 else Z0)
                                                       for y in range(4)] for x in range(4)], 2)))
    # several arrays: a sum of two Clifford terms is not a multiple of a unitary
    I2 = ident(2)
    C2 = littleC2()
    Ssum = [[gadd((2 * I2[i][j][0], 2 * I2[i][j][1]), C2[i][j]) for j in range(2)] for i in range(2)]   # 2(I + C)
    G = mmul(mdag(Ssum), Ssum)
    say("  free addition of two arrays at m=1 holding v and Cv gives (I + C) v; (I+C)^dag (I+C) = %s / 4 is not a"
        " multiple of 1 (not unitary up to scale): %s" % (G, not (G[0][1] == Z0 and G[0][0] == G[1][1])))


# ---------------------------------------------------------------------------------------------------------------
# F. E8 / E7 line systems  <->  4-qubit / 3-qubit Paulis ; extended Hamming code ; octads
# ---------------------------------------------------------------------------------------------------------------
def e8_roots():
    from fractions import Fraction as Fr
    R = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [0] * 8
                    v[i], v[j] = 2 * si, 2 * sj          # doubled coordinates
                    R.append(tuple(v))
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            R.append(tuple(signs))                       # doubled: (+-1/2)^8 -> (+-1)^8
    return R                                              # all coordinates doubled; inner product /4


def ip4(a, b):
    return sum(x * y for x, y in zip(a, b))              # = 4 * true inner product


def b_q(a, c, qbasis, b):
    va = vc = 0
    for k in range(len(qbasis)):
        if (a >> k) & 1:
            va ^= qbasis[k]
        if (c >> k) & 1:
            vc ^= qbasis[k]
    return b(va, vc)


def section_F():
    say("== F. E8 / E7 lines, the [8,4,4] code, octads  <->  multi-qubit Paulis and stabilizer codes")
    R = e8_roots()
    assert len(R) == 240
    lines = []
    seen = set()
    for r in R:
        if r not in seen:
            neg = tuple(-x for x in r)
            seen |= {r, neg}
            lines.append(r)
    assert len(lines) == 120
    # E8 / 2E8 as F_2^8: coordinates of a vector in a Z-basis of E8 (simple roots), mod 2.
    simple = [(2, -2, 0, 0, 0, 0, 0, 0), (0, 2, -2, 0, 0, 0, 0, 0), (0, 0, 2, -2, 0, 0, 0, 0), (0, 0, 0, 2, -2, 0, 0, 0),
              (0, 0, 0, 0, 2, -2, 0, 0), (0, 0, 0, 0, 0, 2, -2, 0), (0, 0, 0, 0, 0, 2, 2, 0),
              (-1, -1, -1, -1, -1, -1, -1, -1)]
    # solve integer coordinates by exact Gaussian elimination over Q
    from fractions import Fraction as Fr
    A = [[Fr(simple[j][i]) for j in range(8)] for i in range(8)]       # columns = simple roots

    def coords(v):
        M = [row[:] + [Fr(v[i])] for i, row in enumerate(A)]
        n = 8
        for c in range(n):
            p = next(r for r in range(c, n) if M[r][c] != 0)
            M[c], M[p] = M[p], M[c]
            pv = M[c][c]
            M[c] = [x / pv for x in M[c]]
            for r in range(n):
                if r != c and M[r][c] != 0:
                    f = M[r][c]
                    M[r] = [x - f * y for x, y in zip(M[r], M[c])]
        out = [M[i][n] for i in range(n)]
        assert all(x.denominator == 1 for x in out)
        return [int(x) for x in out]
    cls = {}
    for r in lines:
        c = coords(r)
        cls[r] = sum((x % 2) << i for i, x in enumerate(c))
    assert len(set(cls.values())) == 120                # distinct, non-zero classes
    # quadratic form q(v) = (v.v)/2 mod 2 on E8/2E8 and its polar form b(u,v) = u.v mod 2
    Bm = [[(ip4(simple[i], simple[j]) // 4) % 2 for j in range(8)] for i in range(8)]
    qb = [(ip4(simple[i], simple[i]) // 8) % 2 for i in range(8)]

    def b(u, v):
        return sum(Bm[i][j] for i in range(8) if (u >> i) & 1 for j in range(8) if (v >> j) & 1) % 2

    def q(u):
        return (sum(qb[i] for i in range(8) if (u >> i) & 1) + sum(Bm[i][j] for i in range(8) for j in range(i + 1, 8)
                                                                    if (u >> i) & 1 and (u >> j) & 1)) % 2
    # build a hyperbolic basis e_k, f_k (q = 0, b(e_k, f_k) = 1, others 0)
    space = list(range(1, 256))
    E, Fv = [], []
    cur = list(range(256))
    for k in range(4):
        e = next(v for v in cur if v and q(v) == 0)
        f = next(v for v in cur if b(e, v) == 1)
        if q(f):
            f ^= e
        assert q(f) == 0 and b(e, f) == 1
        E.append(e)
        Fv.append(f)
        cur = [v for v in cur if b(v, e) == 0 and b(v, f) == 0]
    # coordinates: v = sum alpha_k e_k + beta_k f_k with alpha_k = b(v, f_k), beta_k = b(v, e_k)
    def to_pauli(v):
        x = sum(b(v, Fv[k]) << k for k in range(4))
        z = sum(b(v, E[k]) << k for k in range(4))
        return x | (z << 4)
    P = {r: to_pauli(cls[r]) for r in lines}
    odd_y = all(popc((P[r] & 15) & (P[r] >> 4)) % 2 == 1 for r in lines)
    allodd = {v for v in range(1, 256) if popc((v & 15) & (v >> 4)) % 2 == 1}
    bij = set(P.values()) == allodd
    comm = all(((ip4(r, s) // 4) % 2 == 0) == (omega(P[r], P[s], 4) == 0) for r in lines for s in lines)
    say("  E8: 120 root lines -> 120 distinct classes of E8/2E8 (quadric O+(8,2), q = norm/2 mod 2); a hyperbolic basis"
        " sends them to 4-qubit Paulis: all have an odd number of Y: %s; they are exactly the 120 odd-Y Paulis: %s;"
        " root lines orthogonal <-> Paulis commute (all 120^2 pairs): %s" % (odd_y, bij, comm))
    assert odd_y and bij and comm
    # E7 = lines orthogonal to one root line r0 -> r0^perp / <r0> = 6-dim symplectic space = 3-qubit Paulis
    r0 = lines[0]
    e7 = [r for r in lines if ip4(r, r0) == 0]
    c0 = cls[r0]
    quo = {}
    for r in e7:
        v = cls[r]
        quo.setdefault(frozenset([v, v ^ c0]), []).append(r)
    comm7 = all(((ip4(r, s) // 4) % 2 == 0) == (b(cls[r], cls[s]) == 0) for r in e7 for s in e7)
    say("  E7: %d lines orthogonal to a root line; their classes give %d distinct non-zero vectors of r0^perp/<r0>"
        " (6-dim symplectic space = 3-qubit Paulis, 63 of them); orthogonal <-> symplectically orthogonal"
        " (commuting 3-qubit Paulis): %s" % (len(e7), len(quo), comm7))
    assert len(e7) == 63 and len(quo) == 63 and comm7
    # octads: sets of 8 mutually orthogonal root lines
    idx = {r: i for i, r in enumerate(lines)}
    orth = [set(j for j, s in enumerate(lines) if ip4(r, s) == 0 and s != r) for r in lines]
    octads = set()

    def grow(cl, cand):
        if len(cl) == 8:
            octads.add(frozenset(cl))
            return
        for c in sorted(cand):
            if c > max(cl):
                grow(cl + [c], cand & orth[c])
    for i in range(120):
        grow([i], orth[i])
    ranks = {}
    paulis_comm = True
    for o in octads:
        vs = [P[lines[i]] for i in o]
        ranks[dim2(span(vs))] = ranks.get(dim2(span(vs)), 0) + 1
        paulis_comm &= all(omega(a, c, 4) == 0 for a in vs for c in vs)
    say("  octads (8 mutually orthogonal root lines): %d; F_2-span dimension of their 8 Paulis: %s;"
        " each octad = 8 mutually commuting odd-Y 4-qubit Paulis: %s" % (len(octads), ranks, paulis_comm))
    # octad -> the Lagrangian its Paulis span; compare with the 2,295 Lagrangian frames at m = 4
    L4 = all_lagrangians_direct(4)
    oddY = lambda v: popc((v & 15) & (v >> 4)) % 2 == 1
    real = {L for L in L4 if not any(oddY(v) for v in L)}
    olag = {}
    for o in octads:
        Lo = span([P[lines[i]] for i in o])
        olag[Lo] = o
    exact8 = all(sorted(P[lines[i]] for i in o) == sorted(v for v in Lo if oddY(v)) for Lo, o in olag.items())
    say("  m=4 Lagrangian frames: %d = %d with no odd-Y element (totally singular: real stabilizer groups) + %d others;"
        " octad -> span is injective (%d distinct) and onto the others: %s; an octad is exactly the odd-Y half of its"
        " Lagrangian: %s" % (len(L4), len(real), len(L4) - len(real), len(olag),
                            set(olag) == (L4 - real), exact8))
    assert set(olag) == (L4 - real) and len(olag) == len(octads) and exact8 and len(real) == 270
    # E7: maximal sets of mutually orthogonal lines (heptads) <-> the 135 Lagrangian frames at m = 3
    e7i = [idx[r] for r in e7]
    hept = set()

    def grow7(cl, cand):
        if len(cl) == 7:
            hept.add(frozenset(cl))
            return
        for c in sorted(cand):
            if c > max(cl):
                grow7(cl + [c], cand & orth[c])
    e7set = set(e7i)
    for i in e7i:
        grow7([i], orth[i] & e7set)
    # images in r0^perp/<r0>: map each E7 line to its quotient class, then to F_2^6 coordinates
    qbasis = []
    for r in e7:
        v = cls[r]
        cand = qbasis + [v]
        if dim2(span(cand + [c0])) > dim2(span(qbasis + [c0])):
            qbasis.append(v)
    assert len(qbasis) == 6
    def qcoord(v):     # coordinates of v modulo c0 in qbasis
        for mask in range(64):
            w = 0
            for k in range(6):
                if (mask >> k) & 1:
                    w ^= qbasis[k]
            if w == v or w == v ^ c0:
                return mask
        raise AssertionError
    hl = set()
    for h in hept:
        hl.add(frozenset(qcoord(cls[lines[i]]) for i in h) | {0})
    iso = all(len(H) == 8 and all(b_q(a, c, qbasis, b) == 0 for a in H for c in H) for H in hl)
    say("  E7: maximal sets of mutually orthogonal lines: %d of size 7; each maps to the 7 non-zero vectors of a"
        " totally isotropic 3-space of r0^perp/<r0> (a Lagrangian of 3-qubit Paulis): %s; distinct images %d"
        " (= the 135 Lagrangian frames at m = 3)" % (len(hept), iso, len(hl)))
    assert len(hept) == 135 and len(hl) == 135 and iso
    # the extended Hamming [8,4,4] code: doubly even, self-dual; C (+) C is a Lagrangian at m = 8, fixed by the kernel
    gens = [0b11110000, 0b11001100, 0b10101010, 0b11111111]
    C = span(gens)
    de = all(popc(c) % 4 == 0 for c in C)
    sd = perp_dot(C, 8) == C
    m = 8
    Lc = L_sub(C, m)
    Fk = lambda v: ((v & 255) ^ (v >> 8)) | ((v >> 8) << 8)
    fixed = frozenset(Fk(v) for v in Lc) == Lc
    css = Lc == L_css(C, m)
    say("  extended Hamming [8,4,4]: doubly even %s, self-dual %s; its subspace frame L_C equals the CSS Lagrangian"
        " C (+) C ([[8,0,4]] stabilizer state): %s; fixed by the kernel F (L_U -> L_U^perp): %s;"
        " d(L_0, L_C) = %d, d(L_C, L_F) = %d" % (de, sd, css, fixed, m - dim2(Lc & L0(m)), m - dim2(Lc & LF(m))))
    assert de and sd and fixed and css
    # Construction A: E8 (doubled coordinates scaled by sqrt 2) = {x in Z^8 : x mod 2 in C} / sqrt 2 :
    # norm-4 vectors of that lattice: +-2e_i (16) and +-1 on the 14 weight-4 codewords (224) = 240; the 8 lines
    # +-2e_i form one octad.
    n4 = 16 + 16 * sum(1 for c in C if popc(c) == 4)
    say("  Construction A from [8,4,4]: minimal vectors 16 + 16*%d = %d = 240; the coordinate frame {+-2e_i} is one"
        " octad, so each octad presents E8 as Construction A of an [8,4,4] code (stated, not enumerated)"
        % (sum(1 for c in C if popc(c) == 4), n4))
    return len(octads)


def section_G():
    say("== G. engine frames at matrix level (explicit Gaussian matrices): equal labels <-> equal up to a free matrix")
    for m in (1, 2, 3):
        n = 1 << m
        reps = {label_of(ident(n), m): ident(n)}
        order = [label_of(ident(n), m)]
        coinc = 0
        prods = 0
        for L in order:
            S = reps[L]
            for z in range(1, n):
                # F in {1, phaseMat(c,0)} with c.z = 1 reaches both neighbours along z
                cz = next(c for c in range(1, n) if popc(c & z) & 1)
                for F in (ident(n), phaseM(cz, 0, m)):
                    N = mmul(dir2M(z, m), mmul(F, S))
                    prods += 1
                    Lb = label_of(N, m)
                    assert dim2(Lb & L) == m - 1, "a move did not go to an adjacent frame"
                    if Lb not in reps:
                        reps[Lb] = N
                        order.append(Lb)
                    else:
                        R = to_monomial(mmul(N, mdag(reps[Lb])))
                        assert R is not None and monomial_is_free(R, m), "same label, ratio not free"
                        coinc += 1
        want = prod((1 << i) + 1 for i in range(1, m + 1))
        # different labels: the ratio is never free (a free matrix fixes L_0, so it cannot change the label)
        diff_ok = True
        labs = list(reps)
        for a, bL in itertools.islice(itertools.combinations(labs, 2), 400):
            R = to_monomial(mmul(reps[a], mdag(reps[bL])))
            if R is not None and monomial_is_free(R, m):
                diff_ok = False
        say("  m=%d: frames reached from 1 by single directions (BFS): %d = prod(2^i+1) = %d; %d products, %d label"
            " coincidences, each ratio a free matrix (closed form); sampled distinct-label pairs never differ by a"
            " free matrix: %s" % (m, len(reps), want, prods, coinc, diff_ok))
        assert len(reps) == want and diff_ok


if __name__ == "__main__":
    import time
    secs = sys.argv[1:] or list("ABCGDE") + ["D2", "D3", "F"]
    t0 = time.time()
    for s in secs:
        t = time.time()
        {"A": section_A, "B": section_B, "C": section_C, "D": section_D, "D2": section_D2, "D3": section_D3,
         "E": section_E, "F": section_F, "G": section_G}[s]()
        say("  (%s: %.1fs)" % (s, time.time() - t))
    say("ALL CHECKS PASSED in %.1fs" % (time.time() - t0))
