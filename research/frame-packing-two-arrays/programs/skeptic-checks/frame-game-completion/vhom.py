"""vhom.py (v-pk2-complete): structure checks of pk2-complete Theorem 1, FORMAL in C[G] with exact coefficients.
 (a) homomorphism + lift independence: for random words w = g_r ... g_1 in letters (paid dir z, phase letters,
     address maps, shifts, scalar i), prod_t Lambda(gbar_t) (each from its own letter) == [g n] D(g n) for g = the
     matrix product and a random n in N;   and pi(Lambda(gbar)) == gbar (pi groups coefficients by symplectic image)
 (b) D(dir z) == c_z (e_1 + e_2) + cbar_z (e_0 + e_3) exactly, and its number of group elements, every z
 (c) Remark 3 sizes: phase letter, elementary address map, shift, scalar, random address map (and even-block part)
 (d) character orbits of the letters in block 0 (w -> gbar w) and block 2 (w -> label of g P_w g^T)
usage: vhom.py m nwords wordlen"""
import sys, random, itertools
import numpy as np
from fractions import Fraction as Fr
from veng import *
m = int(sys.argv[1]); nw = int(sys.argv[2]); L = int(sys.argv[3]); rng = random.Random(1000 + m)
E = VEng(m); n = E.n
def Lam(g):
    """Lambda(gbar) as an element of C[G]: {key: GQ}"""
    out = {}
    for (u, j), c in E.D(g).items():
        k = E.key(g @ E.nmat(u, j)); out[k] = out.get(k, GQ(0)) + c
    return {k: v for k, v in out.items() if not v.iszero()}
def cmul(X, Y):
    out = {}
    for kx, vx in X.items():
        Mx = E.store[kx]
        for ky, vy in Y.items():
            k = E.key(Mx @ E.store[ky]); out[k] = out.get(k, GQ(0)) + vx * vy
    return {k: v for k, v in out.items() if not v.iszero()}
def eq(X, Y):
    return {k: (v.r, v.i) for k, v in X.items()} == {k: (v.r, v.i) for k, v in Y.items()}
def rand_letter():
    t = rng.randrange(5)
    if t <= 1: return E.dirm(rng.randrange(1, n)), 'dir'
    if t == 2: return E.phase_letter(rng.randrange(0, n), rng.randrange(2)), 'ph'
    if t == 3:
        while True:
            cols = [rng.randrange(1, n) for _ in range(m)]
            span = {0}
            for c in cols: span |= {x ^ c for x in span}
            if len(span) == n: return E.addr_mat(cols), 'addr'
    return E._X(rng.randrange(n)) * (1j ** rng.randrange(4)), 'shift'
# (a)
okh = 0; okpi = 0
for t in range(nw):
    lets = [rand_letter() for _ in range(L)]
    X = {E.key(E.I): GQ(1)}; g = E.I.copy()
    for M, _ in lets: X = cmul(Lam(M), X); g = M @ g
    nrand = E.nmat(rng.randrange(E.NL), rng.randrange(4))
    Y = Lam(g @ nrand)
    okh += eq(X, Y)
    # pi o Lambda = id
    byS = {}
    for k, v in Y.items():
        sk_ = E.skey(E.gbar(E.store[k])); byS[sk_] = byS.get(sk_, GQ(0)) + v
    byS = {k: (v.r, v.i) for k, v in byS.items() if not v.iszero()}
    okpi += byS == {E.skey(E.gbar(g)): (1, 0)}
print(f"m={m} (a) homomorphism/lift-independence: {okh}/{nw} words of {L} letters equal; pi o Lambda = id: {okpi}/{nw}; "
      f"letter types {[t for _, t in lets]}", flush=True)
# (b)
def e_k(k):
    return {(0, j): GQ.of((1j) ** ((-j * k) % 4)).scale(Fr(1, 4)) for j in range(4)}
def nmul(A, B):   # product in C[N] on (u, j) keys: i^j P_u i^j' P_u' = i^(j+j') (-1)^(b.a') P_(u+u')
    out = {}
    for (u, j), a in A.items():
        for (v, k), b in B.items():
            sgn = pc((u >> m) & (v & E.mask)) & 1
            key = (u ^ v, (j + k + 2 * sgn) % 4); out[key] = out.get(key, GQ(0)) + a * b
    return {k: v for k, v in out.items() if not v.iszero()}
def nadd(A, B):
    out = dict(A)
    for k, v in B.items(): out[k] = out.get(k, GQ(0)) + v
    return {k: v for k, v in out.items() if not v.iszero()}
def norm(A):   # collapse to group elements: (u, j) with i^j scalar is a distinct group element; already canonical
    return {k: (v.r, v.i) for k, v in A.items()}
okb = 0; sizes = set()
for z in range(1, n):
    cz = {(0, 0): GQ(Fr(1, 2), Fr(-1, 2)), (z, 0): GQ(Fr(1, 2), Fr(1, 2))}
    cbz = {(0, 0): GQ(Fr(1, 2), Fr(1, 2)), (z, 0): GQ(Fr(1, 2), Fr(-1, 2))}
    want = nadd(nmul(cz, nadd(e_k(1), e_k(2))), nmul(cbz, nadd(e_k(0), e_k(3))))
    got = E.D(E.dirm(z))
    okb += norm(got) == norm(want); sizes.add(len(got))
print(f"m={m} (b) D(dir z) == closed form for {okb}/{n-1} z; group-element counts {sorted(sizes)}", flush=True)
# (c)
def sz(M):
    D = E.D(M); Ds = E.Dk(M)
    even = [int(sum(abs(Ds[k]) > 1e-12)) for k in (0, 2)]; odd = [int(sum(abs(Ds[k]) > 1e-12)) for k in (1, 3)]
    return len(D), even, odd
print(f"m={m} (c) sizes (total group elements, Pauli terms in blocks 0,2, in blocks 1,3):")
print("   phase letters z!=0:", sorted({sz(E.phase_letter(z, c)) [0] for z in range(1, n) for c in (0, 1)}))
print("   phase letter z=0,c=1 (scalar i):", sz(E.phase_letter(0, 1)), "  shift X^1:", sz(E._X(1)))
elem = [i for i in range(m)]
cols = [1 << k for k in range(m)]; cols[0] |= 2          # x_1 += x_0 style elementary map
print("   elementary address map:", sz(E.addr_mat(cols)))
rs = []
for _ in range(6):
    while True:
        cols = [rng.randrange(1, n) for _ in range(m)]
        span = {0}
        for c in cols: span |= {x ^ c for x in span}
        if len(span) == n: break
    rs.append(sz(E.addr_mat(cols)))
print("   random address maps:", rs, flush=True)
# (d)
gens = [E.dirm(z) for z in range(1, n)] + [E.phase_letter(1 << i, 0) for i in range(m)] + \
       [E.phase_letter((1 << i) | (1 << j), 0) for i in range(m) for j in range(i + 1, m)]
for i in range(m):
    for j in range(m):
        if i != j:
            cols = [1 << k for k in range(m)]; cols[i] |= 1 << j; gens.append(E.addr_mat(cols))
def orbits(act):
    seen = set(); out = []
    for w0 in range(E.NL):
        if w0 in seen: continue
        orb = {w0}; fr = [w0]
        while fr:
            w = fr.pop()
            for g in gens:
                v = act(g, w)
                if v not in orb: orb.add(v); fr.append(v)
        seen |= orb; out.append(len(orb))
    return sorted(out)
o0 = orbits(lambda g, w: E.ident(g @ E.P[w] @ g.conj().T)[1])
o2 = orbits(lambda g, w: E.ident(g @ E.P[w] @ g.T)[1])
Np, Nm = 2 ** (m - 1) * (2 ** m + 1), 2 ** (m - 1) * (2 ** m - 1)
print(f"m={m} (d) block-0 orbits {o0} (claim [1, {4**m-1}]); block-2 orbits {o2} (claim [{Nm}, {Np}])", flush=True)
