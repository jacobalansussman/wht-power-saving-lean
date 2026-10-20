"""vgroup.py (v-fp-build, 15 Oct 2026): independent check of fp-build Prop 1.3 and Theorem 2.2's lift, block m = 2.
Own code (numpy only; none of fp-build's or ledger.py's code).
(1) Engine group G from actual 4x4 matrices: |G|, scalars in G, symplectic images, fibre sizes mod mu4,
    linear independence of each fibre.
(2) Lambda_1(lam) = sum_g a_g [g] with sum_g a_g g = I over the fibre of lam (block f = 1 mod 4), and
    Lambda_3(lam) with sum_g a_g conj(g) = I (block f = 3 mod 4, not in the claim). Homomorphism test on raw
    random complex arrays at f = 1..7 and 9: Lambda(l)Lambda(u) v == Lambda(lu) v.
(3) Lambda(T_z) has exactly two terms, [dir z] and [X^z dir z]; Lambda(p) for p in P is supported on free letters.
"""
import itertools, sys
import numpy as np
rng = np.random.default_rng(20261015)
m = 2; n = 1 << m
ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
def Xz(z):
    M = np.zeros((n, n), complex)
    for x in range(n): M[x ^ z, x] = 1
    return M
def Zz(b):
    return np.diag([(-1) ** bin(b & x).count("1") for x in range(n)]).astype(complex)
def dirz(z): return ca * np.eye(n) + cb * Xz(z)
def par(x): return bin(x).count("1") & 1
def phase(z, c): return np.diag([1j ** ((par(z & x) + c) % 2) for x in range(n)])
def addr(h):  # h: list of images of basis vectors (ints); |x> -> |h x>
    M = np.zeros((n, n), complex)
    for x in range(n):
        y = 0
        for k in range(m):
            if (x >> k) & 1: y ^= h[k]
        M[y, x] = 1
    return M
GLgens = [[2, 1], [1, 3]]  # swap, and e1->e1, e2->e1+e2
free_gens = [Xz(z) for z in range(1, n)] + [phase(z, c) for z in range(n) for c in (0, 1)] + [addr(h) for h in GLgens]
paid_gens = [dirz(z) for z in range(1, n)]
def key(M):
    return tuple(np.round(M.real * 1e6).astype(np.int64).ravel()) + tuple(np.round(M.imag * 1e6).astype(np.int64).ravel())
def closure(gens):
    I = np.eye(n, dtype=complex)
    seen = {key(I): I}; fr = [I]
    while fr:
        nf = []
        for g in fr:
            for h in gens:
                k = h @ g
                kk = key(k)
                if kk not in seen: seen[kk] = k; nf.append(k)
        fr = nf
    return seen
G = closure(free_gens + paid_gens)
B = closure(free_gens)
print(f"|G| = {len(G)}, |B| = {len(B)}")
scal = [g for g in G.values() if np.allclose(g, g[0, 0] * np.eye(n))]
print("scalars in G:", sorted(np.round(np.angle([g[0, 0] for g in scal]) / (np.pi / 4)).astype(int) % 8), "(units of pi/4)")
# symplectic image: label (a|b) as int a | b << m ; Pauli P(a|b) = X^a Z^b
def pauli(v): return Xz(v & (n - 1)) @ Zz(v >> m)
PA = {v: pauli(v) for v in range(1 << (2 * m))}
def which(M):
    for v, P in PA.items():
        c = np.vdot(P, M) / n
        if abs(abs(c) - 1) < 1e-9 and np.allclose(M, c * P): return v
    raise ValueError
def sym(g):
    gi = np.linalg.inv(g)
    return tuple(which(g @ PA[1 << i] @ gi) for i in range(2 * m))  # images of basis labels
fib = {}
for k, g in G.items():
    fib.setdefault(sym(g), []).append(g)
print("distinct symplectic images:", len(fib), " (|Sp(4,2)| = 720)")
def modmu4(gs):
    reps = {}
    for g in gs:
        # canonical: multiply by i^k making first nonzero entry have angle in [0, pi/2)
        idx = np.flatnonzero(np.abs(g.ravel()) > 1e-9)[0]
        a = np.angle(g.ravel()[idx]) % (2 * np.pi)
        k = int(np.floor(a / (np.pi / 2) + 1e-9)) % 4
        h = g * (1j) ** (-k)
        reps.setdefault(key(h), h)
    return list(reps.values())
fibr = {lam: modmu4(gs) for lam, gs in fib.items()}
print("fibre sizes mod mu4:", sorted(set(len(v) for v in fibr.values())), " |G/mu4| =", sum(len(v) for v in fibr.values()))
bad = sum(1 for v in fibr.values() if np.linalg.matrix_rank(np.array([g.ravel() for g in v])) != len(v))
print("fibres with linearly dependent matrices:", bad)
# which images are in P (stabiliser of L0 = z-labels {(0|b)}): images of (0|e_i) have zero x-part
def inP(lam): return all((lam[m + i] & (n - 1)) == 0 for i in range(m))
Bimgs = set(sym(g) for g in B.values())
Pimgs = set(l for l in fib if inP(l))
print("|image of B| =", len(Bimgs), " |P| (images fixing L0) =", len(Pimgs), " equal:", Bimgs == Pimgs)
# free part over P: are all fibre elements over P in B (i.e., preimage of P = B)?
print("preimage of P equals B:", sum(len(fib[l]) for l in Pimgs) == len(B))
# ---------------------------------------------------------------- the lifts
def coeffs(lam, conj=False):
    gs = fibr[lam]
    A = np.array([(np.conj(g) if conj else g).ravel() for g in gs]).T
    a = np.linalg.solve(A.T.conj() @ A, A.T.conj() @ np.eye(n).ravel())
    assert np.allclose(A @ a, np.eye(n).ravel())
    return list(zip(a, gs))
def compose(l1, l2):  # l1 o l2 on labels; image as tuple of images of basis labels
    def app(l, v):
        out = 0
        for i in range(2 * m):
            if (v >> i) & 1: out ^= l[i]
        return out
    return tuple(app(l1, l2[i]) for i in range(2 * m))
def act(terms, v, f):
    """sum_g a_g g^(x)f applied to tensor v of shape (n,)*f"""
    out = np.zeros_like(v)
    for a, g in terms:
        w = v
        for ax in range(f):
            w = np.moveaxis(np.tensordot(g, w, axes=([1], [ax])), 0, ax)
        out += a * w
    return out
lams = list(fib.keys())
Sid = tuple(1 << i for i in range(2 * m))
for conj, name in ((False, "Lambda_1"), (True, "Lambda_3")):
    res = {}
    for f in (1, 2, 3, 4, 5, 6, 7, 9):
        worst = 0.0
        for trial in range(3 if f >= 7 else 6):
            l1 = lams[rng.integers(len(lams))]; l2 = lams[rng.integers(len(lams))]
            v = rng.standard_normal((n,) * f) + 1j * rng.standard_normal((n,) * f)
            lhs = act(coeffs(l1, conj), act(coeffs(l2, conj), v, f), f)
            rhs = act(coeffs(compose(l1, l2), conj), v, f)
            worst = max(worst, np.linalg.norm(lhs - rhs) / np.linalg.norm(rhs))
        res[f] = worst
    print(name, "homomorphism, worst relative error by f:", {f: f"{e:.1e}" for f, e in res.items()})
# Lambda_1(T_z): support
for z in range(1, n):
    lam = sym(dirz(z))
    t = [(a, g) for a, g in coeffs(lam) if abs(a) > 1e-9]
    names = []
    for a, g in t:
        for cand, nm in ((dirz(z), "dir z"), (Xz(z) @ dirz(z), "X^z dir z")):
            for k in range(4):
                if np.allclose(g, (1j) ** k * cand): names.append((nm, np.round(a * (1j) ** k, 6)))
    print(f"Lambda_1(T_{z}) terms: {len(t)}:", names)
# Lambda_1(p) for p in P is supported on B
Bkeys = set(B.keys())
ok = all(all(key(g) in Bkeys or any(key(g * (1j) ** k) in Bkeys for k in range(4)) for a, g in coeffs(l) if abs(a) > 1e-12) for l in Pimgs)
print("Lambda_1(P) supported on free letters:", ok)
