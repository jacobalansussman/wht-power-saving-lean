#!/usr/bin/env python3
"""v-dict-compute: an independent re-check of the [computed] entries of the frames <-> stabilizer dictionary.
Written from the Lean definitions (BinaryFrames.lean, Letters.lean, Free.lean, Split.lean), not from x-dict's code.
Run: python3 -I recheck.py (needs numpy) [sections]
Convention: address y (int), bit i = qubit i.  Pauli (x|z) = X^x Z^z, X^x|y> = |y^x>, Z^z|y> = (-1)^{z.y}|y>.
Symplectic vector v = x | (z << m).  M_C: C P(a) C^-1 ~ P(M_C a).  Label L(S) = {a : S P(a) S^-1 diagonal}.
"""
import sys, itertools, random, time
from collections import deque, Counter
import numpy as np

I_ = 1j
TOL = 1e-9
pc = lambda v: bin(v).count("1")
dot = lambda a, b: pc(a & b) & 1
OUT = []


def report(name, ok, extra=""):
    s = "%s  %s  %s" % ("PASS" if ok else "FAIL", name, extra)
    print(s, flush=True)
    OUT.append((name, ok))


# ------------------------------------------------------------------ engine matrices (from the Lean source)
ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
littleC = np.array([[ca, cb], [cb, ca]])
littleH = np.array([[1, 1], [1, -1]], dtype=complex)


def kron_all(mats):
    out = np.array([[1]], dtype=complex)
    for M in reversed(mats):  # qubit 0 = least significant bit -> rightmost kron factor
        out = np.kron(out, M)
    return out


def digitProd(M2, m):
    return kron_all([M2] * m)


def shift(z, m):
    d = 1 << m
    S = np.zeros((d, d), dtype=complex)
    for i in range(d):
        S[i, i ^ z] = 1  # Lean: shift z i j = if i+z = j then 1
    return S


def dirm(z, m):
    return ca * np.eye(1 << m) + cb * shift(z, m)


def permMat(G, m):  # G: list of images of basis vectors (linear map on F_2^m); Lean: permMat x y = [G x = y]
    d = 1 << m
    P = np.zeros((d, d), dtype=complex)
    for x in range(d):
        P[x, applyG(G, x)] = 1
    return P


def applyG(G, x):
    r = 0
    for i, g in enumerate(G):
        if x >> i & 1:
            r ^= g
    return r


def phaseMat(z, c, m):
    return np.diag([1j if (dot(z, x) ^ c) else 1 for x in range(1 << m)])


_PC = {}


def pauli(v, m):
    if (v, m) in _PC:
        return _PC[(v, m)]
    _PC[(v, m)] = P = _pauli(v, m)
    return P


def _pauli(v, m):
    x, z = v & ((1 << m) - 1), v >> m
    d = 1 << m
    P = np.zeros((d, d), dtype=complex)
    for y in range(d):
        P[y ^ x, y] = (-1) ** dot(z, y)
    return P


def is_diag(A):
    return np.allclose(A - np.diag(np.diag(A)), 0, atol=TOL)


def label(S, m, Sinv=None):
    Sinv = np.linalg.inv(S) if Sinv is None else Sinv
    return frozenset(v for v in range(1 << 2 * m) if is_diag(S @ pauli(v, m) @ Sinv))


def symp_of(C, m, Cinv=None):
    """M_C as list of images of the 2m basis vectors, or None if C is not Clifford."""
    Cinv = np.linalg.inv(C) if Cinv is None else Cinv
    cols = []
    for b in range(2 * m):
        A = C @ pauli(1 << b, m) @ Cinv
        found = None
        for w in range(1 << 2 * m):
            P = pauli(w, m)
            k = np.vdot(P.ravel(), A.ravel()) / (1 << m)
            if abs(abs(k) - 1) < 1e-7 and np.allclose(A, k * P, atol=1e-7):
                found = w
                break
        if found is None:
            return None
        cols.append(found)
    return cols


def is_clifford(C, m):
    return symp_of(C, m) is not None


def sapply(M, v):
    r = 0
    b = 0
    while v:
        if v & 1:
            r ^= M[b]
        v >>= 1
        b += 1
    return r


def scompose(A, B):  # A o B
    return tuple(sapply(A, c) for c in B)


def omega(u, v, m):
    M = (1 << m) - 1
    return (pc((u & M) & (v >> m)) + pc((v & M) & (u >> m))) & 1


def span(vs):
    s = {0}
    for v in vs:
        if v not in s:
            s |= {x ^ v for x in s}
    return frozenset(s)


def dimof(S):
    return len(S).bit_length() - 1


def all_lagrangians(m):
    """Exhaustive: BFS over isotropic subspaces by adding one vector."""
    N = 1 << 2 * m
    level = {frozenset([0])}
    for _ in range(m):
        nxt = set()
        for L in level:
            for v in range(1, N):
                if v not in L and all(omega(v, u, m) == 0 for u in L):
                    nxt.add(span(list(L) + [v]))
        level = nxt
    return level


def perp1(U, m):
    return frozenset(v for v in range(1 << m) if all(dot(u, v) == 0 for u in U))


def subspaces(m):
    found = {frozenset([0])}
    fr = [frozenset([0])]
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


def L_U(U, m):
    return frozenset(u | ((u ^ v) << m) for u in U for v in perp1(U, m))


def CSS(U, m):
    return frozenset(u | (v << m) for u in U for v in perp1(U, m))


L0 = lambda m: frozenset(z << m for z in range(1 << m))
LF = lambda m: frozenset(z | (z << m) for z in range(1 << m))


# ------------------------------------------------------------------ section A: conventions, paid moves, kernel, WHT
def secA():
    for m in (1, 2, 3):
        d = 1 << m
        ok = all(np.allclose(shift(z, m), pauli(z, m)) for z in range(d))
        ok &= all(np.allclose(np.diag([(-1) ** dot(a, x) for x in range(d)]), pauli(a << m, m)) for a in range(d))
        report("A1 m=%d shift(z) = X^z, diag (-1)^(a.x) = Z^a" % m, ok)
        # dir(e_j) = SX_j ; S_j dir S_j = ((1+i)/2) littleH_j ; dir(z) = e^{i pi/4} exp(-i pi/4 X^z)
        SX = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])
        Sg = np.diag([1, 1j])
        ok = True
        for j in range(m):
            one = [np.eye(2)] * m
            Mj = list(one); Mj[j] = SX
            ok &= np.allclose(dirm(1 << j, m), kron_all(Mj))
            Sj = list(one); Sj[j] = Sg; Sj = kron_all(Sj)
            Hj = list(one); Hj[j] = littleH; Hj = kron_all(Hj)
            ok &= np.allclose(Sj @ dirm(1 << j, m) @ Sj, ca * Hj)
        for z in range(1, d):
            X = pauli(z, m)
            ok &= np.allclose(dirm(z, m), np.exp(1j * np.pi / 4) * (np.eye(d) - 1j * X) / np.sqrt(2))
        report("A2 m=%d dir(e_j)=SX_j, S_j dir(e_j) S_j=((1+i)/2)littleH_j, dir(z)=e^{i pi/4}exp(-i pi/4 X^z)" % m, ok)
        # dir(z) = P_G dir(e_low) P_G^-1 for an address map with G^-1 e_low = z  (P_G X^e P_G^-1 = X^{G^-1 e})
        ok = True
        for z in range(1, d):
            # build G^-1 with column 0 = z, completed to a basis; G = inverse
            low = (z & -z).bit_length() - 1
            Ginv = [z] + [1 << i for i in range(m) if i != low]
            # order columns: e_0 -> z, others e_i
            cols = [0] * m
            cols[0] = z
            rest = [1 << i for i in range(m) if i != low]
            for k, i in enumerate(range(1, m)):
                cols[i] = rest[k]
            Ginvmat = cols
            # invert over F_2 by brute force
            G = None
            for cand in itertools.product(range(d), repeat=m):
                if all(applyG(cand, applyG(Ginvmat, 1 << i)) == 1 << i for i in range(m)):
                    G = list(cand); break
            P = permMat(G, m)
            ok &= np.allclose(dirm(z, m), P @ dirm(1, m) @ np.linalg.inv(P))
        report("A3 m=%d dir(z) = P_G dir(e_0) P_G^-1 for every z != 0" % m, ok)
        # blocks: product of k independent directions = P (SX^(x k) (x) 1) P^-1, and Clifford
        ok = True
        for k in range(1, m + 1):
            for zs in itertools.combinations(range(1, d), k):
                if dimof(span(zs)) < k:
                    continue
                B = np.eye(d, dtype=complex)
                for z in zs:
                    B = B @ dirm(z, m)
                # address map with G^-1 e_i = z_i (i < k), completed
                cols = list(zs)
                for i in range(m):
                    if len(cols) == m:
                        break
                    if dimof(span(cols + [1 << i])) > len(cols):
                        cols.append(1 << i)
                G = None
                for cand in itertools.product(range(d), repeat=m):
                    if all(applyG(cand, applyG(cols, 1 << i)) == 1 << i for i in range(m)):
                        G = list(cand); break
                P = permMat(G, m)
                target = kron_all([SX] * k + [np.eye(2)] * (m - k)) if True else None
                # qubit i = bit i; kron_all puts list[0] on bit 0
                ok &= np.allclose(B, P @ target @ np.linalg.inv(P))
        report("A4 m=%d blockMat of k independent directions = P_G (SX^(x k) (x) 1) P_G^-1, all families" % m, ok)
        # kernel, WHT
        K = digitProd(littleC, m)
        W = digitProd(littleH, m)
        Hn = W / np.sqrt(d)
        Sm = digitProd(Sg, m)
        ok = np.allclose(K, digitProd(SX, m))
        ok &= label(K, m) == LF(m)
        ok &= label(W, m) == frozenset(range(d))  # {(x|0)}
        ok &= np.allclose(W @ K @ np.linalg.inv(W), np.diag([1j ** pc(x) for x in range(d)]))
        ok &= np.allclose(Sm @ K @ Sm, np.exp(1j * m * np.pi / 4) * Hn)
        report("A5 m=%d K=SX^(x m), L(K)=L_F (Y frame), L(W)=X frame, W K W^-1=diag(i^wt), S K S=e^{i m pi/4}H" % m, ok)


# ------------------------------------------------------------------ section B: Lagrangians, stabilizer states, fibres
def stab_states(m):
    d = 1 << m
    one = [np.eye(2)] * m
    gens = []
    Hh = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    for j in range(m):
        a = list(one); a[j] = Hh; gens.append(kron_all(a))
        a = list(one); a[j] = np.diag([1, 1j]); gens.append(kron_all(a))
    for c in range(m):
        for t in range(m):
            if c != t:
                P = np.zeros((d, d))
                for y in range(d):
                    P[y ^ ((y >> c & 1) << t), y] = 1
                gens.append(P)

    def key(v):
        k = np.flatnonzero(np.abs(v) > 1e-9)[0]
        v = v * (abs(v[k]) / v[k])
        return tuple(np.round(np.concatenate([v.real, v.imag]), 6)), v

    v0 = np.zeros(d, dtype=complex); v0[0] = 1
    k0, v0 = key(v0)
    seen = {k0: v0}
    fr = [v0]
    while fr:
        nx = []
        for v in fr:
            for g in gens:
                k, w = key(g @ v)
                if k not in seen:
                    seen[k] = w
                    nx.append(w)
        fr = nx
    return np.array(list(seen.values()))


def stab_group(Psi, m):
    """for each state, frozenset of Paulis v with |<psi|P|psi>| = 1."""
    d = 1 << m
    n = len(Psi)
    groups = [[] for _ in range(n)]
    for v in range(1 << 2 * m):
        P = pauli(v, m)
        ev = (Psi.conj() * (Psi @ P.T)).sum(axis=1)
        for i in np.flatnonzero(np.abs(np.abs(ev) - 1) < 1e-7):
            groups[i].append(v)
        # also: every expectation is 0 or modulus 1
        bad = np.flatnonzero((np.abs(ev) > 1e-7) & (np.abs(np.abs(ev) - 1) > 1e-7))
        assert len(bad) == 0
    return [frozenset(g) for g in groups]


LAG, STATES, SGRP = {}, {}, {}


def secB(ms=(1, 2, 3, 4)):
    for m in ms:
        t = time.time()
        Ls = all_lagrangians(m)
        LAG[m] = Ls
        exp = 1
        for i in range(1, m + 1):
            exp *= (1 << i) + 1
        report("B1 m=%d Lagrangians by isotropic enumeration" % m, len(Ls) == exp, "%d (prod (2^i+1) = %d)" % (len(Ls), exp))
        Psi = stab_states(m)
        STATES[m] = Psi
        report("B2 m=%d stabilizer states by BFS from |0> under H,S,CNOT" % m, len(Psi) == (1 << m) * exp,
               "%d" % len(Psi))
        G = stab_group(Psi, m)
        SGRP[m] = G
        fib = {}
        for i, g in enumerate(G):
            fib.setdefault(g, []).append(i)
        ok = set(fib) == set(Ls) and all(len(f) == 1 << m for f in fib.values())
        report("B3 m=%d stabilizer group (signs dropped) of each state is a Lagrangian; map onto; fibres 2^m" % m, ok,
               "%d fibres" % len(fib))
        ok = True
        for f in fib.values():
            A = Psi[f]
            ok &= np.allclose(A.conj() @ A.T, np.eye(1 << m), atol=1e-7)
        report("B4 m=%d every fibre is an orthonormal basis (ALL %d fibres)" % (m, len(fib)), ok, "%.1fs" % (time.time() - t))


# ------------------------------------------------------------------ section C: free group
def free_gens(m):
    d = 1 << m
    gens = []
    # all invertible linear maps (permMat) -- use elementary transvections + swaps, enough to generate GL
    for i in range(m):
        for j in range(m):
            if i != j:
                G = [1 << k for k in range(m)]
                G[j] ^= 1 << i
                gens.append(permMat(G, m))
    for z in range(d):
        for c in (0, 1):
            gens.append(phaseMat(z, c, m))
    for z in range(1, d):
        gens.append(shift(z, m))
    return gens


def mkey(A):
    return tuple(np.round(np.concatenate([A.real.ravel(), A.imag.ravel()]), 6))


def closure(gens):
    d = gens[0].shape[0]
    e = np.eye(d, dtype=complex)
    seen = {mkey(e): e}
    fr = [e]
    while fr:
        nx = []
        for A in fr:
            for g in gens:
                B = g @ A
                k = mkey(B)
                if k not in seen:
                    seen[k] = B
                    nx.append(B)
        fr = nx
    return seen


def sp_close(gens, m):
    e = tuple(1 << b for b in range(2 * m))
    seen = {e}
    fr = [e]
    while fr:
        nx = []
        for A in fr:
            for g in gens:
                B = scompose(g, A)
                if B not in seen:
                    seen.add(B)
                    nx.append(B)
        fr = nx
    return seen


def secC():
    for m in (1, 2):
        d = 1 << m
        Fr = closure(free_gens(m))
        # all monomial matrices with entries in {1,i,-1,-i}
        cands = 0
        cliff = set()
        for perm in itertools.permutations(range(d)):
            for ph in itertools.product(range(4), repeat=d):
                A = np.zeros((d, d), dtype=complex)
                for x in range(d):
                    A[perm[x], x] = 1j ** ph[x]
                cands += 1
                if is_clifford(A, m):
                    cliff.add(mkey(A))
        report("C1 m=%d Free closure = all Clifford monomial i-power matrices" % m,
               set(Fr) == cliff, "closure %d, Clifford monomials %d of %d candidates" % (len(Fr), len(cliff), cands))
        # Bravyi-Maslov h-free generators {X, CNOT, CZ, S} plus global i
        one = [np.eye(2)] * m
        gens = [np.eye(d) * 1j]
        for j in range(m):
            a = list(one); a[j] = np.array([[0, 1], [1, 0]]); gens.append(kron_all(a))
            a = list(one); a[j] = np.diag([1, 1j]); gens.append(kron_all(a))
        for c in range(m):
            for t in range(m):
                if c != t:
                    P = np.zeros((d, d))
                    for y in range(d):
                        P[y ^ ((y >> c & 1) << t), y] = 1
                    gens.append(P)
                    if c < t:
                        gens.append(np.diag([(-1) ** ((y >> c & 1) & (y >> t & 1)) for y in range(d)]))
        BM = closure(gens)
        report("C2 m=%d Free closure = <X, CNOT, CZ, S, iI>" % m, set(BM) == set(Fr), "%d" % len(BM))
    # m = 3 diagonal subgroup
    m = 3
    D = closure([phaseMat(z, c, m) for z in range(8) for c in (0, 1)])
    report("C3 m=3 group generated by phaseMat (diagonal part)", len(D) == 2048, "%d" % len(D))
    # phase letter i^[x_j + x_k] = S_j S_k CZ_jk
    m = 2
    S0 = np.diag([1j ** (x & 1) for x in range(4)]); S1 = np.diag([1j ** (x >> 1 & 1) for x in range(4)])
    ok = np.allclose(phaseMat(3, 0, 2), S0 @ S1 @ np.diag([1, 1, 1, -1]))
    report("C4 phaseMat(e0+e1,0) = S_0 S_1 CZ_01", ok)
    # symplectic images
    gl = {1: 1, 2: 6, 3: 168}
    spo = {1: 6, 2: 720, 3: 1451520}
    for m in (1, 2, 3):
        fg = [symp_of(g, m) for g in free_gens(m)]
        assert all(f is not None for f in fg)
        img = sp_close([tuple(f) for f in fg], m)
        exp = gl[m] * 2 ** (m * (m + 1) // 2)
        fixes = all(frozenset(sapply(M, v) for v in L0(m)) == L0(m) for M in img)
        report("C5 m=%d symplectic image of Free" % m, len(img) == exp and fixes, "%d (= %d), all fix L_0: %s" % (len(img), exp, fixes))
        Hq = np.kron(np.eye(1 << (m - 1)), littleH) if m > 1 else littleH
        h = symp_of(Hq, m)
        if m < 3:
            full = sp_close([tuple(f) for f in fg] + [tuple(h)], m)
            report("C6 m=%d Free + one H generates Sp(2m,2)" % m, len(full) == spo[m], "%d" % len(full))
        else:
            FG3 = [tuple(f) for f in fg]
            globals()["_FG3"] = (FG3, tuple(h))


# ------------------------------------------------------------------ section D: frames from explicit engine matrices
def is_free_matrix(T, m, Dset):
    d = 1 << m
    # monomial with i-power entries, affine permutation, diagonal in the phaseMat group
    nz = np.abs(T) > 1e-9
    if not (nz.sum(axis=0) == 1).all() or not (nz.sum(axis=1) == 1).all():
        return False
    perm = [int(np.flatnonzero(nz[:, x])[0]) for x in range(d)]  # T e_x ~ e_perm[x]
    b = perm[0]
    lin = [perm[x] ^ b for x in range(d)]
    if any(lin[x ^ y] != lin[x] ^ lin[y] for x in range(d) for y in range(d)):
        return False
    Pm = np.zeros((d, d), dtype=complex)
    for x in range(d):
        Pm[perm[x], x] = 1
    Dg = np.linalg.inv(Pm) @ T  # T = Pm Dg
    if not is_diag(Dg):
        return False
    return mkey(Dg) in Dset


def secD():
    rng = random.Random(7)
    for m, nwords in ((1, 2000), (2, 6000), (3, 12000)):
        d = 1 << m
        Dset = set(closure([phaseMat(z, c, m) for z in range(d) for c in (0, 1)]))
        F = free_gens(m)
        dirs = [dirm(z, m) for z in range(1, d)]
        reps = {}
        coinc = 0
        bad = 0
        nonlag = 0
        trans_bad = 0
        Lset = LAG.get(m) or all_lagrangians(m)
        for w in range(nwords):
            S = np.eye(d, dtype=complex)
            ln = rng.randint(0, 6 * m)
            for _ in range(ln):
                if rng.random() < 0.35:
                    S = dirs[rng.randrange(len(dirs))] @ S
                else:
                    S = F[rng.randrange(len(F))] @ S
            Sinv = np.linalg.inv(S)
            L = label(S, m, Sinv)
            if L not in Lset:
                nonlag += 1
            # label = M_S^-1 L_0
            MS = symp_of(S, m, Sinv)
            if frozenset(v for v in range(1 << 2 * m) if sapply(MS, v) in L0(m)) != L:
                bad += 1
            # one more direction acts as the transvection along P = M_S^-1 (z|0)
            if w % 5 == 0:
                z = rng.randrange(1, d)
                L2 = label(dirm(z, m) @ S, m)
                MSinv = None
                P = next(a for a in range(1 << 2 * m) if sapply(MS, a) == z)
                tau = frozenset(a ^ (P if omega(a, P, m) else 0) for a in L)
                if tau != L2 or P in L:
                    trans_bad += 1
            if L in reps:
                coinc += 1
                if not is_free_matrix(S @ np.linalg.inv(reps[L]), m, Dset):
                    bad += 1
            else:
                reps[L] = S
        exp = {1: 3, 2: 15, 3: 135}[m]
        report("D1 m=%d frames reached by random engine words (free letters + dir)" % m, len(reps) == exp,
               "%d distinct labels, all Lagrangian: %s" % (len(reps), nonlag == 0))
        report("D2 m=%d equal labels => S' S^-1 is a free matrix (affine perm x phaseMat group); label = M_S^-1 L_0" % m,
               bad == 0, "%d coincidences checked, %d failures" % (coinc, bad))
        report("D3 m=%d dir(z) acts on labels as transvection tau_P, P = M_S^-1(z|0) outside L" % m, trans_bad == 0)
        globals().setdefault("_REPS", {})[m] = reps


# ------------------------------------------------------------------ section E: subspace frames, CSS, DDM forms, Phi_U
def secE():
    for m in (1, 2, 3, 4):
        Ls = LAG.get(m) or all_lagrangians(m)
        Us = subspaces(m)
        subs = {L_U(U, m) for U in Us}
        css = {CSS(U, m) for U in Us}
        yz = set()
        xz = set()
        M = (1 << m) - 1
        for L in Ls:
            Y = [v for v in L if (v & M) == (v >> m) or (v & M) == 0]
            if span(Y) == L:
                yz.add(L)
            Xz = [v for v in L if (v >> m) == 0 or (v & M) == 0]
            if span(Xz) == L:
                xz.add(L)
        ok = subs <= Ls and css <= Ls and subs == yz and css == xz
        report("E1 m=%d subspace frames L_U = (Y,Z)-spanned Lagrangians; CSS = (X,Z)-spanned" % m, ok,
               "#subspaces %d, #L_U %d, #CSS %d, #(Y,Z) %d, sets equal L_U==CSS: %s" % (len(Us), len(subs), len(css), len(yz), subs == css))
        # twist S^(x m): (x|z) -> (x|z+x) maps CSS_U to L_U; kernel map (x|z)->(x+z|z) maps L_U -> L_{U^perp};
        # H^(x m): swap maps CSS_U -> CSS_{U^perp}
        tw = lambda v: v ^ ((v & M) << m)
        ker = lambda v: v ^ (v >> m)
        hh = lambda v: (v >> m) | ((v & M) << m)
        ok = all(frozenset(map(tw, CSS(U, m))) == L_U(U, m) for U in Us)
        ok2 = all(frozenset(map(ker, L_U(U, m))) == L_U(perp1(U, m), m) for U in Us)
        ok3 = all(frozenset(map(hh, CSS(U, m))) == CSS(perp1(U, m), m) for U in Us)
        report("E2 m=%d twist S^m CSS_U = L_U; kernel L_U -> L_Uperp; H^m CSS_U -> CSS_Uperp" % m, ok and ok2 and ok3)
        if m <= 3:
            # the symplectic maps used above are those of the actual matrices
            Sm = digitProd(np.diag([1, 1j]), m)
            K = digitProd(littleC, m)
            okm = all(sapply(symp_of(Sm, m), 1 << b) == tw(1 << b) for b in range(2 * m))
            okm &= all(sapply(symp_of(K, m), 1 << b) == ker(1 << b) for b in range(2 * m))
            # the twisted CSS state sum_{y in U} i^{|y|}|y> is stabilized (up to sign) by L_U
            d = 1 << m
            okst = True
            for U in Us:
                psi = np.array([1j ** pc(y) if y in U else 0 for y in range(d)])
                psi = psi / np.linalg.norm(psi)
                g = stab_group(psi[None, :], m)[0]
                okst &= g == L_U(U, m)
            report("E3 m=%d symplectic maps of S^(x m), K match; sum_{y in U} i^|y| |y> has stabilizer L_U" % m, okm and okst)
        if m == 1:
            non = [L for L in Ls if L not in subs]
            report("E4 m=1 the one non-subspace frame is <X>", non == [frozenset({0, 1})], str([sorted(x) for x in non]))
        report("E5 m=%d L_0 = L_{0}, L_F = L_{F^m}" % m, L0(m) == L_U(frozenset([0]), m) and LF(m) == L_U(frozenset(range(1 << m)), m))
        # Dehaene-De Moor: X-projection U and symmetric beta
        byU = {}
        okb = True
        for L in Ls:
            U = frozenset(v & M for v in L)
            Z0 = frozenset(v >> m for v in L if (v & M) == 0)
            okb &= Z0 == perp1(U, m)
            zof = {}
            for v in L:
                zof.setdefault(v & M, v >> m)
            beta = {(u, w): dot(zof[u], w) for u in U for w in U}
            okb &= all(beta[(u, w)] == beta[(w, u)] for u in U for w in U)
            isdot = all(beta[(u, w)] == dot(u, w) for u in U for w in U)
            iszero = all(beta[(u, w)] == 0 for u in U for w in U)
            okb &= (isdot == (L == L_U(U, m))) and (iszero == (L == CSS(U, m)))
            byU.setdefault(U, []).append((isdot, iszero))
        okc = all(len(v) == 2 ** (dimof(U) * (dimof(U) + 1) // 2) and sum(a for a, _ in v) == 1 and sum(b for _, b in v) == 1
                  for U, v in byU.items())
        report("E6 m=%d DDM: each U carries 2^(k(k+1)/2) Lagrangians, beta symmetric, exactly one beta=dot (=L_U), one beta=0 (=CSS)" % m,
               okb and okc and len(byU) == len(Us), "non-subspace frames %d" % (len(Ls) - len(subs)))
    # Phi_U (family 130), nondegenerate U
    for m in (1, 2, 3):
        d = 1 << m
        W = digitProd(littleH, m)
        n = 0
        ok = True
        for U in subspaces(m):
            Up = perp1(U, m)
            if len(U & Up) != 1:
                continue
            n += 1
            # projection onto U along U^perp
            proj = {}
            for u in U:
                for v in Up:
                    proj[u ^ v] = u
            Phi = W @ np.diag([1j ** pc(proj[x]) for x in range(d)]) @ W
            ok &= label(Phi, m) == L_U(U, m)
        report("E7 m=%d L(Phi_U) = L_U for every nondegenerate U" % m, ok, "%d nondegenerate U" % n)


# ------------------------------------------------------------------ section F: graph, H-count, overlaps, KK distance
def secF():
    for m in (1, 2, 3, 4):
        Ls = sorted(LAG.get(m) or all_lagrangians(m), key=sorted)
        idx = {L: i for i, L in enumerate(Ls)}
        n = len(Ls)
        masks = [sum(1 << v for v in L) for L in Ls]
        half = 1 << (m - 1)
        adj = [[j for j in range(n) if pc(masks[i] & masks[j]) == half] for i in range(n)]
        val = set(len(a) for a in adj)
        report("F1 m=%d dual polar valency" % m, val == {2 * ((1 << m) - 1)}, str(val))
        # engine one-direction neighbours of L_0: orbit of {M_dir(z)^-1 L_0} under the free symplectic image
        Mm = (1 << m) - 1
        def dirL(z):  # label of dir(z): {a : M_dir a in L_0}, M_dir = transvection along (z|0) (involution)
            P = z
            return frozenset(a for a in range(1 << 2 * m) if (a ^ (P if omega(a, P, m) else 0)) in L0(m))
        if m <= 3:
            for z in range(1, 1 << m):
                assert label(dirm(z, m), m) == dirL(z)
        # simpler: build free generators from matrices for m<=3, analytic for m=4
        def gl(i, j):
            def f(v):
                x, z = v & Mm, v >> m
                x ^= ((x >> j) & 1) << i
                z ^= ((z >> i) & 1) << j
                return x | (z << m)
            return f
        def ph(i, j):
            def f(v):
                x, z = v & Mm, v >> m
                if i == j:
                    z ^= ((x >> i) & 1) << i
                else:
                    z ^= (((x >> j) & 1) << i) ^ (((x >> i) & 1) << j)
                return x | (z << m)
            return f
        fg = [gl(i, j) for i in range(m) for j in range(m) if i != j] + [ph(i, j) for i in range(m) for j in range(i, m)]
        if m <= 3:  # check these are the symplectic maps of CNOT-type permMat and S / CZ
            F = free_gens(m)
            imgs = {tuple(symp_of(g, m)) for g in F}
            # each analytic generator is in the image group (closure check would be heavy at m=3; check membership in
            # the closure of the matrix images)
            cl = sp_close(list(imgs), m)
            assert all(tuple(f(1 << b) for b in range(2 * m)) in cl for f in fg)
        N0 = {dirL(z) for z in range(1, 1 << m)}
        fr = list(N0)
        while fr:
            nx = []
            for L in fr:
                for f in fg:
                    L2 = frozenset(map(f, L))
                    if L2 not in N0:
                        N0.add(L2); nx.append(L2)
            fr = nx
        nb = {Ls[j] for j in adj[idx[L0(m)]]}
        report("F2 m=%d one engine direction (with free gates) from L_0 reaches exactly the dual-polar neighbours" % m, N0 == nb,
               "%d vs %d" % (len(N0), len(nb)))
        srcs = range(n) if m <= 3 else [idx[L0(m)], idx[LF(m)]]
        ok = True
        for s in srcs:
            dist = {s: 0}
            q = deque([s])
            while q:
                u = q.popleft()
                for v in adj[u]:
                    if v not in dist:
                        dist[v] = dist[u] + 1; q.append(v)
            ok &= all(dist[t] == m - dimof(Ls[s] & Ls[t]) for t in range(n))
        report("F3 m=%d BFS distance = m - dim(A cap B) from %s" % (m, "every source" if m <= 3 else "L_0 and L_F"), ok)
        # KK subspace distance on subspace frames
        Us = list(subspaces(m))
        ok = all(m - dimof(L_U(U, m) & L_U(V, m)) == dimof(U) + dimof(V) - 2 * dimof(U & V) for U in Us for V in Us)
        report("F4 m=%d on subspace frames d = dim U + dim V - 2 dim(U cap V) (all pairs)" % m, ok)
    # overlaps (vectorised): lag index per state, distance matrix between Lagrangians via bitmask intersections
    for m in (1, 2, 3, 4):
        Psi, G = STATES[m], SGRP[m]
        n = len(Psi)
        lags = sorted(set(G), key=sorted)
        li = {L: k for k, L in enumerate(lags)}
        lab = np.array([li[g] for g in G])
        masks = [sum(1 << v for v in L) for L in lags]
        nl = len(lags)
        Dl = np.array([[m - (pc(masks[a] & masks[b]).bit_length() - 1) for b in range(nl)] for a in range(nl)])
        if m <= 3:
            rows = np.arange(n)
        else:
            rng = np.random.default_rng(3)
            rows = np.concatenate([[0], rng.choice(n, 400, replace=False)])
        Gm = np.abs(Psi[rows].conj() @ Psi.T) ** 2
        Dm = Dl[lab[rows][:, None], lab[None, :]]
        expv = 2.0 ** (-Dm)
        is0 = np.abs(Gm) < 1e-9
        isd = np.abs(Gm - expv) < 1e-9
        ok = bool((is0 | isd).all())
        # attained: for every (lag_a, lag_b) pair present, some pair of states reaches 2^-d
        pa = lab[rows][:, None] * nl + lab[None, :]
        present = set(np.unique(pa).tolist())
        hit = set(np.unique(pa[isd]).tolist())
        report("F5 m=%d |<a|b>|^2 in {0, 2^-d(A,B)}; 2^-d attained for every Lagrangian pair seen" % m, ok and present == hit,
               "%d x %d states, %d Lagrangian pairs" % (len(rows), n, len(present)))
    # min #H over Sp(2m,2) by 0-1 BFS vs rank B vs m - dim(L_0 cap T L_0)
    for m in (1, 2, 3):
        Mm = (1 << m) - 1
        def gl(i, j):
            def f(v):
                x, z = v & Mm, v >> m
                x ^= ((x >> j) & 1) << i
                z ^= ((z >> i) & 1) << j
                return x | (z << m)
            return f
        def ph(i, j):
            def f(v):
                x, z = v & Mm, v >> m
                if i == j:
                    z ^= ((x >> i) & 1) << i
                else:
                    z ^= (((x >> j) & 1) << i) ^ (((x >> i) & 1) << j)
                return x | (z << m)
            return f
        def hd(i):
            def f(v):
                x, z = v & Mm, v >> m
                xi, zi = (x >> i) & 1, (z >> i) & 1
                x = (x & ~(1 << i)) | (zi << i)
                z = (z & ~(1 << i)) | (xi << i)
                return x | (z << m)
            return f
        tab = lambda f: tuple(f(1 << b) for b in range(2 * m))
        free = [tab(gl(i, j)) for i in range(m) for j in range(m) if i != j] + [tab(ph(i, j)) for i in range(m) for j in range(i, m)]
        hs = [tab(hd(0))]   # one H on qubit 0 only; free swaps move it
        T64 = lambda g: [sapply(g, v) for v in range(1 << 2 * m)]
        free = [T64(g) for g in free]
        hs = [T64(g) for g in hs]
        scompose = lambda tg, A: tuple(tg[c] for c in A)
        e = tuple(1 << b for b in range(2 * m))
        dist = {e: 0}
        dq = deque([e])
        t0 = time.time()
        while dq:
            A = dq.popleft()
            dA = dist[A]
            for g in free:
                B = scompose(g, A)
                if B not in dist or dist[B] > dA:
                    dist[B] = dA
                    dq.appendleft(B)
            for g in hs:
                B = scompose(g, A)
                if B not in dist:
                    dist[B] = dA + 1
                    dq.append(B)
        # 0-1 BFS with "B not in dist" for weight-1 edges can be wrong if later improved; redo a correctness pass:
        changed = True
        while changed:
            changed = False
            for A, dA in list(dist.items()):
                for g in free:
                    B = scompose(g, A)
                    if dist[B] > dA:
                        dist[B] = dA; changed = True
                for g in hs:
                    B = scompose(g, A)
                    if dist[B] > dA + 1:
                        dist[B] = dA + 1; changed = True
        hist = Counter()
        ok = True
        for T, dT in dist.items():
            Bcols = [T[m + j] & Mm for j in range(m)]  # x-part of T(0|e_j)
            rB = dimof(span(Bcols))
            inter = dimof(frozenset(sapply(T, z << m) for z in range(1 << m)) & L0(m))
            ok &= (dT == rB == m - inter)
            hist[dT] += 1
        report("F6 m=%d min #H (0-1 BFS over Sp(2m,2)) = rank B = m - dim(L_0 cap T L_0) on every element" % m, ok,
               "|Sp|=%d hist %s (%.0fs)" % (len(dist), dict(sorted(hist.items())), time.time() - t0))
    # the caveat: S |+> = |+i>, d(<X>,<Y>) = 1
    plus = np.array([1, 1]) / np.sqrt(2)
    ok = np.allclose(np.diag([1, 1j]) @ plus, np.array([1, 1j]) / np.sqrt(2))
    ok &= 1 - dimof(frozenset({0, 1}) & frozenset({0, 3})) == 1
    report("F7 caveat: S|+> = |+i> with no H, yet d(<X>,<Y>) = 1", ok)


# ------------------------------------------------------------------ section G: DFT pieces, multi-array
def secG():
    F2 = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    F4 = np.array([[1j ** (j * k) for k in range(4)] for j in range(4)]) / 2
    F4b = np.array([[(-1j) ** (j * k) for k in range(4)] for j in range(4)]) / 2
    cs = np.diag([1, 1, 1, 1j])
    cz = np.diag([1, 1, 1, -1])
    ok = is_clifford(F2, 1) and not is_clifford(F4, 2) and not is_clifford(F4b, 2) and not is_clifford(cs, 2) and is_clifford(cz, 2)
    # bit-reversed DFT_4 also not Clifford
    R = np.zeros((4, 4)); R[[0, 2, 1, 3], [0, 1, 2, 3]] = 1
    ok &= not is_clifford(R @ F4, 2)
    report("G1 DFT_2 Clifford; DFT_4 (either sign, with/without bit reversal) not; controlled-S not; CZ is", ok)
    A = np.eye(2) + littleC
    G = A.conj().T @ A
    report("G2 (I+C)^dag (I+C) = [[12,4],[4,12]]/4", np.allclose(G, np.array([[12, 4], [4, 12]]) / 4), str(np.round(G * 4, 9).real.tolist()))


# ------------------------------------------------------------------ section H: E8, E7, octads, heptads, Hamming
def hnf_basis(gens):
    """Z-basis of the lattice spanned by integer row vectors (simple row reduction)."""
    rows = [list(g) for g in gens]
    n = len(rows[0])
    basis = []
    col = 0
    while rows and col < n:
        rows = [r for r in rows if any(r)]
        piv = [r for r in rows if r[col] != 0]
        if not piv:
            col += 1
            continue
        while len([r for r in rows if r[col] != 0]) > 1:
            piv = sorted([r for r in rows if r[col] != 0], key=lambda r: abs(r[col]))
            p = piv[0]
            for r in piv[1:]:
                q = r[col] // p[col]
                for k in range(n):
                    r[k] -= q * p[k]
            rows = [r for r in rows if any(r)]
        p = [r for r in rows if r[col] != 0][0]
        basis.append(p)
        rows = [r for r in rows if r is not p]
        col += 1
    return basis


def secH():
    Hgen = [0b00001111, 0b00110011, 0b01010101, 0b11111111]
    C = sorted(span(Hgen))
    vec = lambda w: [w >> i & 1 for i in range(8)]
    # extended Hamming: doubly even, self-dual
    ok = len(C) == 16 and all(pc(c) % 4 == 0 for c in C) and all(dot(a, b) == 0 for a in C for b in C)
    # Construction A (scaled by sqrt 2): Lambda = {x in Z^8 : x mod 2 in C}; norm-4 vectors = roots
    roots = [x for x in itertools.product(range(-2, 3), repeat=8)
             if sum(t * t for t in x) == 4 and sum(((t % 2) << i) for i, t in enumerate(x)) in C]
    ok &= len(roots) == 240 and sum(1 for r in roots if sorted(map(abs, r)) == [0] * 7 + [2]) == 16
    report("H1 extended Hamming doubly even self-dual; Construction A has 240 = 16 + 224 minimal vectors", ok)
    # L_C = C (+) C, d's at m = 8 and kernel fixed
    m = 8
    Cs = frozenset(C)
    LC = L_U(Cs, m)
    ok = LC == frozenset(a | (b << m) for a in C for b in C)
    ok &= m - dimof(LC & L0(m)) == 4 and m - dimof(LC & LF(m)) == 4
    ok &= frozenset(v ^ (v >> m) for v in LC) == LC
    report("H2 m=8: L_C = C(+)C, d(L_0,L_C) = d(L_C,L_F) = 4, kernel fixes L_C", ok)
    # E8/2E8 coordinates
    gens = [vec(g) for g in Hgen] + [[2 if i == j else 0 for i in range(8)] for j in range(8)]
    B = hnf_basis([list(g) for g in gens])
    Bm = np.array(B, dtype=float)
    assert len(B) == 8 and abs(abs(np.linalg.det(Bm)) - 16) < 1e-6  # covolume 2^8/16 = 16
    Binv = np.linalg.inv(Bm)
    def coords(x):
        c = np.array(x, dtype=float) @ Binv
        ci = np.round(c).astype(int)
        assert np.allclose(c, ci)
        return sum((int(t) % 2) << i for i, t in enumerate(ci))
    Bint = [list(map(int, b)) for b in B]
    def lift(v):
        return [sum(((v >> i) & 1) * Bint[i][k] for i in range(8)) for k in range(8)]
    Bf = lambda u, v: (sum(a * b for a, b in zip(lift(u), lift(v))) // 2) % 2   # true inner product mod 2
    qf = lambda v: (sum(a * a for a in lift(v)) // 4) % 2                          # true norm / 2 mod 2
    # hyperbolic basis
    space = list(range(1, 256))
    W = set(range(256))
    pairs = []
    for _ in range(4):
        e = next(v for v in sorted(W) if v and qf(v) == 0)
        f = next(v for v in sorted(W) if Bf(e, v) == 1)
        if qf(f):
            f ^= e
        assert qf(f) == 0 and Bf(e, f) == 1
        pairs.append((e, f))
        W = {v for v in W if Bf(v, e) == 0 and Bf(v, f) == 0}
    def tox(v):
        x = sum(Bf(v, f) << i for i, (e, f) in enumerate(pairs))
        z = sum(Bf(v, e) << i for i, (e, f) in enumerate(pairs))
        return x | (z << 4)
    allv = {v: tox(v) for v in range(256)}
    ok = len(set(allv.values())) == 256
    ok &= all(qf(v) == pc((allv[v] & 15) & (allv[v] >> 4)) % 2 for v in range(256))
    ok &= all(Bf(u, v) == omega(allv[u], allv[v], 4) for u in range(256) for v in range(0, 256, 7))
    report("H3 (E8/2E8, norm/2 mod 2) isometric to (F_2^8, x.z), polar form = omega", ok)
    # lines
    lines = []
    seen = set()
    for r in roots:
        k = tuple(r) if tuple(r) > tuple(-t for t in r) else tuple(-t for t in r)
        if k not in seen:
            seen.add(k); lines.append(list(k))
    P = [allv[coords(r)] for r in lines]
    oddY = [v for v in range(1, 256) if pc((v & 15) & (v >> 4)) % 2 == 1]
    ok = len(lines) == 120 and sorted(P) == sorted(oddY)
    orth = lambda a, b: sum(s * t for s, t in zip(a, b)) == 0
    ok &= all(orth(lines[i], lines[j]) == (omega(P[i], P[j], 4) == 0) for i in range(120) for j in range(120) if i != j)
    report("H4 120 root lines <-> 120 odd-Y 4-qubit Paulis; orthogonal <-> commuting (all 120^2 pairs)", ok)
    # hub's line list (copyfloor e8_lines construction): check its knows-graph is this graph
    hubL = []
    for a in range(8):
        hubL.append([2 if i == a else 0 for i in range(8)])
    sups = sorted(tuple(i for i in range(8) if w >> i & 1) for w in C if pc(w) == 4)
    for sup in sups:
        for signs in itertools.product((1, -1), repeat=3):
            v = [0] * 8; v[sup[0]] = 1
            for s, i in zip(signs, sup[1:]):
                v[i] = s
            hubL.append(v)
    Ph = [allv[coords(r)] for r in hubL]
    ok = sorted(Ph) == sorted(oddY) and all(orth(hubL[i], hubL[j]) == (omega(Ph[i], Ph[j], 4) == 0) for i in range(120) for j in range(120) if i != j)
    report("H5 the hub/copyfloor P120 lines: knows-graph (N = 0, i != j) = commutation graph of odd-Y Paulis", ok)
    # E7: lines orthogonal to hub line 0 (2e_0)
    r0 = hubL[0]
    e7 = [i for i in range(1, 120) if orth(hubL[i], r0)]
    rb = Ph[0]
    # quotient r^perp / <r>: pick representative odd-Y element; build symplectic basis of the quotient
    perpv = [v for v in range(256) if omega(v, rb, 4) == 0]
    # symplectic basis of perp/<rb>
    def reduce(v):
        return min(v, v ^ rb)
    Q = sorted({reduce(v) for v in perpv if v not in (0, rb)})
    Wq = set(Q)
    qp = []
    for _ in range(3):
        e = min(Wq)
        f = next(v for v in sorted(Wq) if omega(e, v, 4) == 1)
        qp.append((e, f))
        Wq = {v for v in Wq if omega(v, e, 4) == 0 and omega(v, f, 4) == 0}
    def to3(v):
        x = sum(omega(v, f, 4) << i for i, (e, f) in enumerate(qp))
        z = sum(omega(v, e, 4) << i for i, (e, f) in enumerate(qp))
        return x | (z << 3)
    P7 = [to3(Ph[i]) for i in e7]
    ok = len(e7) == 63 and sorted(P7) == list(range(1, 64))
    ok &= all(orth(hubL[e7[a]], hubL[e7[b]]) == (omega(P7[a], P7[b], 3) == 0) for a in range(63) for b in range(63) if a != b)
    report("H6 E7: 63 lines orthogonal to a line <-> 63 nontrivial 3-qubit Paulis; orthogonal <-> commuting (all 63^2)", ok)
    # octads: 8 mutually orthogonal lines; heptads in E7
    def cliques(verts, adjf, k):
        res = []
        def rec(cur, cand):
            if len(cur) == k:
                res.append(tuple(cur)); return
            for i, v in enumerate(cand):
                rec(cur + [v], [w for w in cand[i + 1:] if adjf(v, w)])
        rec([], list(verts))
        return res
    O = cliques(range(120), lambda a, b: orth(hubL[a], hubL[b]), 8)
    L4 = LAG.get(4) or all_lagrangians(4)
    ts = [L for L in L4 if all(pc((v & 15) & (v >> 4)) % 2 == 0 for v in L)]
    imgs = set()
    ok = True
    for o in O:
        S = span([Ph[i] for i in o])
        ok &= S in L4 and frozenset(v for v in S if pc((v & 15) & (v >> 4)) % 2) == frozenset(Ph[i] for i in o)
        imgs.add(S)
    ok &= len(O) == 2025 and len(imgs) == 2025 and imgs == set(L4) - set(ts) and len(ts) == 270
    report("H7 2025 octads <-> 4-qubit Lagrangians not totally singular (odd-Y half); 2295 = 2025 + 270", ok,
           "octads %d, totally singular %d" % (len(O), len(ts)))
    Hp = cliques(range(63), lambda a, b: orth(hubL[e7[a]], hubL[e7[b]]), 7)
    L3 = LAG.get(3) or all_lagrangians(3)
    im3 = {span([P7[i] for i in h]) for h in Hp}
    ok = len(Hp) == 135 and im3 == set(L3) and all(len(span([P7[i] for i in h])) == 8 for h in Hp)
    report("H8 135 E7 heptads <-> 135 three-qubit Lagrangians", ok, "heptads %d" % len(Hp))
    # real (totally singular) stabilizer groups 2, 6, 30, 270
    cnt = []
    for m in (1, 2, 3, 4):
        Mm = (1 << m) - 1
        cnt.append(sum(1 for L in (LAG.get(m) or all_lagrangians(m)) if all(pc((v & Mm) & (v >> m)) % 2 == 0 for v in L)))
    report("H9 Lagrangians with no odd-Y element (real stabilizer groups) m=1..4", cnt == [2, 6, 30, 270], str(cnt))


if __name__ == "__main__":
    secs = sys.argv[1:] or ["A", "B", "C", "D", "E", "F", "G", "H"]
    t = time.time()
    if any(s in secs for s in "DEFH") and "B" not in secs:
        secs = ["B"] + secs
    for s in secs:
        print("== section", s, flush=True)
        globals()["sec" + s]()
    print("TOTAL %d checks, %d FAIL, %.0fs" % (len(OUT), sum(1 for _, ok in OUT if not ok), time.time() - t))
