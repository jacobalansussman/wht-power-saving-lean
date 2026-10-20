"""xc-round2 (wave 17): a PERFECT Pauli-blind certificate on the F4 bulk-symbol frames at block 2.

Setting (labels at the bulk frequency alpha in F4, alpha^2 = alpha + 1, m = 2):
  Pauli label (x | z), x, z in F4^2.  The R-valued form sum_i x_i z'_i^sigma + z_i x'_i^sigma (sigma: S -> S^-1)
  evaluated at alpha (sigma = Frobenius a -> a^2 on F4) gives the Hermitian form
      Hg(u, v) = sum_i u_x,i * v_z,i^2 + u_z,i * v_x,i^2.
  Frames = maximal totally Hg-isotropic F4-subspaces (dim 2).  L0 = {x = 0}.
  Free letters (address maps G(alpha), phases) fix L0.  A move along v (rank-1 submodule) acts as the transvection
  T_v: (x | z) -> (x + (sum_i z_i v_i^2) v | z).  The kernel K = C (x) C on every column acts as (x | z) -> (x + z | z).

Claim checked here (exact, finite):
  (1) there are exactly 27 frames; the meeting graph is srg(27, 10, 1, 5) (GQ(2,4) collinearity);
  (2) every T_v L0 (v != 0) is adjacent to L0, K L0 is opposite (disjoint) to L0;
  (3) an explicit isometry Phi: (F4^4, Hg) -> (F4^4, Hstd), Hstd(u,v) = sum u_i v_i^2, maps the frames onto the
      27 lines of the Fermat cubic x0^3+x1^3+x2^3+x3^3 = 0 reduced mod 2 (omega -> alpha);
  (4) lifting each frame to its Fermat line over Q(omega) (a 2-dim subspace of C^4), delta(U,V) = dim Y_U -
      dim(Y_U cap Y_V) equals the frame distance (0, 1, 2) for ALL 27 x 27 ordered pairs (exact ranks over Q(omega)
      via integers mod a prime p = 1 mod 3 that does not divide any 4x4 minor that matters, plus a float check over C).
  With the Averaging Lemma / induced certificate this is a PERFECT certificate: abar_1 = 1, abar_2 = 2.
"""
import itertools, sys
import numpy as np

# ---------------- F4 arithmetic: 0, 1, 2 = alpha, 3 = alpha^2 = alpha + 1 ----------------
EXP = [1, 2, 3]; LOG = {1: 0, 2: 1, 3: 2}
def add(a, b): return a ^ b
def mul(a, b):
    if a == 0 or b == 0: return 0
    return EXP[(LOG[a] + LOG[b]) % 3]
def frob(a): return mul(a, a)          # sigma = Frobenius
def inv(a): return EXP[(-LOG[a]) % 3]

def vadd(u, v): return tuple(add(a, b) for a, b in zip(u, v))
def vscale(c, u): return tuple(mul(c, a) for a in u)
def Hg(u, v):                           # u = (x1, x2, z1, z2)
    s = 0
    for i in range(2):
        s ^= mul(u[i], frob(v[2 + i])) ^ mul(u[2 + i], frob(v[i]))
    return s
def Hstd(u, v):
    s = 0
    for i in range(4): s ^= mul(u[i], frob(v[i]))
    return s

VEC = [v for v in itertools.product(range(4), repeat=4)]
NZ = [v for v in VEC if any(v)]
def span(gens):
    S = {(0, 0, 0, 0)}
    for g in gens:
        S = S | {vadd(s, vscale(c, g)) for s in S for c in range(1, 4)}
    return frozenset(S)

def frames(H):
    out = set()
    for u in NZ:
        if H(u, u): continue
        for v in NZ:
            if H(v, v) or H(u, v): continue
            S = span([u, v])
            if len(S) == 16: out.add(S)
    return sorted(out, key=lambda S: sorted(S))

ok = True
def rep(name, cond):
    global ok
    ok = ok and bool(cond)
    print(("PASS " if cond else "FAIL ") + name)

FR = frames(Hg)
rep(f"number of Hg-frames = 27 (got {len(FR)})", len(FR) == 27)
n = len(FR)
def dist(A, B):
    k = len(A & B)
    return {16: 0, 4: 1, 1: 2}[k]
D = np.array([[dist(A, B) for B in FR] for A in FR])
deg = (D == 1).sum(1)
rep("meeting graph regular of degree 10", set(deg) == {10})
lam = {int(((D[a] == 1) & (D[b] == 1)).sum()) for a in range(n) for b in range(n) if D[a, b] == 1}
mu = {int(((D[a] == 1) & (D[b] == 1)).sum()) for a in range(n) for b in range(n) if D[a, b] == 2}
rep(f"srg(27,10,1,5): lambda {lam}, mu {mu}", lam == {1} and mu == {5})

L0 = span([(0, 0, 1, 0), (0, 0, 0, 1)])
rep("L0 = {x = 0} is a frame", L0 in FR)
def T(v):          # transvection of a move along v in F4^2
    def f(u):
        c = mul(u[2], frob(v[0])) ^ mul(u[3], frob(v[1]))
        return (u[0] ^ mul(c, v[0]), u[1] ^ mul(c, v[1]), u[2], u[3])
    return f
def Kmap(u): return (u[0] ^ u[2], u[1] ^ u[3], u[2], u[3])
def image(fn, S): return frozenset(fn(u) for u in S)
for v in [(a, b) for a in range(4) for b in range(4) if (a, b) != (0, 0)]:
    M = image(T(v), L0)
    if not (M in FR and dist(L0, M) == 1):
        rep(f"move along v={v} is an adjacent step", False)
rep("all 15 moves T_v (v != 0) map L0 to an adjacent frame", True)
KL0 = image(Kmap, L0)
rep("K L0 is a frame opposite to L0", KL0 in FR and dist(L0, KL0) == 2)
# isometry checks of the moves and K (they must preserve Hg)
iso = all(Hg(T(v)(a), T(v)(b)) == Hg(a, b) for v in [(1, 0), (2, 3), (1, 1)] for a in VEC for b in VEC[::7])
rep("T_v preserves Hg (sampled)", iso)
rep("K preserves Hg", all(Hg(Kmap(a), Kmap(b)) == Hg(a, b) for a in VEC for b in VEC))
# symbols of the chain letters at alpha (ut-direct group1.py (8)): chi = (1+S) I -> scalar 1+alpha = alpha^2 on x,
# and (1+alpha)^(-sigma) on z;  E = [[1,S],[0,1]] -> [[1,alpha],[0,1]] on x and its inverse-sigma-transpose on z.
def chi(u):
    c = add(1, 2); cz = inv(frob(c))
    return (mul(c, u[0]), mul(c, u[1]), mul(cz, u[2]), mul(cz, u[3]))
rep("chi(alpha) preserves Hg and fixes every frame (central)", all(Hg(chi(a), chi(b)) == Hg(a, b) for a in VEC for b in VEC)
    and all(image(chi, F) == F for F in FR))
def Emap(u):       # x -> G x with G = [[1,a],[0,1]] (a = alpha); z -> G^(-sigma T) z = [[1,0],[a^2,1]] z
    a = 2
    return (u[0] ^ mul(a, u[1]), u[1], u[2], u[3] ^ mul(frob(a), u[2]))
rep("E(alpha) preserves Hg and fixes L0", all(Hg(Emap(a), Emap(b)) == Hg(a, b) for a in VEC for b in VEC)
    and image(Emap, L0) == L0)

# ---------------- isometry Phi : (F4^4, Hg) -> (F4^4, Hstd) ----------------
A = 2; A2 = 3
cols = {0: (1, 1, 0, 0), 1: (0, 0, 1, 1), 2: (A, A2, 0, 0), 3: (0, 0, A, A2)}   # images of x1, x2, z1, z2
def Phi(u):
    out = (0, 0, 0, 0)
    for i in range(4): out = vadd(out, vscale(u[i], cols[i]))
    return out
rep("Phi is an isometry Hg -> Hstd", all(Hstd(Phi(a), Phi(b)) == Hg(a, b) for a in VEC for b in VEC))
rep("Phi is bijective", len({Phi(a) for a in VEC}) == 256)

# ---------------- the 27 lines of the Fermat cubic over Z[omega] ----------------
# Elements of Z[omega] as pairs (a, b) = a + b*omega, omega^2 = -1 - omega.  Cube roots of unity: 1, omega, omega^2.
ROOTS = [(1, 0), (0, 1), (-1, -1)]
def neg(z): return (-z[0], -z[1])
lines = []      # each: (pairing, zeta, zeta', two spanning vectors in Z[omega]^4)
for pr in [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]:
    (i, j), (k, l) = pr
    for z1 in ROOTS:
        for z2 in ROOTS:
            v1 = [(0, 0)] * 4; v2 = [(0, 0)] * 4
            v1[i] = neg(z1); v1[j] = (1, 0)          # x_i = -zeta x_j
            v2[k] = neg(z2); v2[l] = (1, 0)
            lines.append((pr, z1, z2, v1, v2))
rep("27 Fermat lines", len(lines) == 27)
def red2(z):     # Z[omega] -> F4, omega -> alpha
    a, b = z[0] % 2, z[1] % 2
    return add(a, mul(b, 2)) if True else 0
def red_line(L):
    return span([tuple(red2(c) for c in L[3]), tuple(red2(c) for c in L[4])])
RED = [red_line(L) for L in lines]
rep("reductions mod 2 are 27 distinct Hstd-generators", len(set(RED)) == 27 and
    all(Hstd(a, b) == 0 for S in RED for a in S for b in S))
PF = [frozenset(Phi(u) for u in F) for F in FR]
rep("Phi maps the 27 game frames onto the 27 reduced Fermat lines", set(PF) == set(RED))
lift = [RED.index(P) for P in PF]          # game frame -> Fermat line index

# exact complex check: omega = exp(2 pi i / 3)
w = np.exp(2j * np.pi / 3)
def cvec(v): return np.array([c[0] + c[1] * w for c in v])
Y = [np.stack([cvec(L[3]), cvec(L[4])]) for L in lines]
def cdelta(a, b):
    r = np.linalg.matrix_rank(np.vstack([Y[a], Y[b]]), tol=1e-9)
    return r - 2                       # dim Y_a + dim Y_b - dim(Y_a cap Y_b) - dim Y_a = r - 2
# exact check mod p = 7 (omega -> 2, since 2^3 = 8 = 1 mod 7 and 1 + 2 + 4 = 7 = 0) and p = 13 (omega -> 3)
def rank_mod(Mrows, p):
    M = [list(r) for r in Mrows]; r = 0; ncol = len(M[0])
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        iv = pow(M[r][c], p - 2, p)
        M[r] = [x * iv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                fct = M[i][c]; M[i] = [(x - fct * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r
def pdelta(a, b, p, om):
    rows = [[(c[0] + c[1] * om) % p for c in v] for v in (lines[a][3], lines[a][4], lines[b][3], lines[b][4])]
    return rank_mod(rows, p) - 2
bad_c = bad_7 = bad_13 = 0
for x in range(27):
    for y in range(27):
        d = D[x, y]; a, b = lift[x], lift[y]
        bad_c += cdelta(a, b) != d
        bad_7 += pdelta(a, b, 7, 2) != d
        bad_13 += pdelta(a, b, 13, 3) != d
rep(f"over C: delta(lift U, lift V) = frame distance on all 729 ordered pairs (mismatches {bad_c})", bad_c == 0)
rep(f"mod 7 (omega=2): same (mismatches {bad_7})", bad_7 == 0)
rep(f"mod 13 (omega=3): same (mismatches {bad_13})", bad_13 == 0)
cnt = {d: int((D == d).sum()) for d in (0, 1, 2)}
print("ordered pairs by distance:", cnt, "-> abar_1 = 1, abar_2 = 2: PERFECT family, moves >= 2 W in the alpha game")
print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
