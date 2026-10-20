"""placements.py (v-t2-pair): how many inequivalent 'Singer games' are there at block m?
The no-phase game with address group H = <sigma> (Singer) and target [I] is the group algebra C[Sym_m x| H] with
moves [zz^T] and target [T], T = I in the chosen coordinates.  Changing coordinates by g in GL(m,2) replaces H by
g H g^-1; equivalently keep H and use the target T = g^-1 g^-T.  Games with targets T, T' are isomorphic (as
presented systems) if T' = n T n^T with n in the normaliser N(H) = <sigma, Frobenius>.  So placements <-> N(H)-orbits
on {g g^T : g in GL(m,2)} = nondegenerate NON-alternating symmetric matrices.  Independent code (no sfield import)."""
import sys, itertools
m = int(sys.argv[1])
POLY = {3: 0b1011, 4: 0b10011, 5: 0b100101}[m]
q = 1 << m
def mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> m: a ^= POLY
    return r
# matrices over F_2 in the polynomial basis 1, a, ..., a^{m-1}; column j = image of a^j
def matof(f):
    cols = [f(1 << j) for j in range(m)]
    return tuple(tuple((cols[j] >> i) & 1 for j in range(m)) for i in range(m))
def mm(A, B):
    return tuple(tuple(sum(A[i][t] & B[t][j] for t in range(m)) & 1 for j in range(m)) for i in range(m))
def tp(A): return tuple(zip(*A))
sig = matof(lambda x: mul(x, 2)); fro = matof(lambda x: mul(x, x))
# group N = <sig, fro>
I = tuple(tuple(int(i == j) for j in range(m)) for i in range(m))
N = {I}; frontier = [I]
while frontier:
    new = []
    for g in frontier:
        for s in (sig, fro):
            h = mm(s, g)
            if h not in N: N.add(h); new.append(h)
    frontier = new
print("m", m, "|N(H)|", len(N))
def det2(A):
    A = [list(r) for r in A]; r = 0
    for c in range(m):
        p = next((i for i in range(r, m) if A[i][c]), None)
        if p is None: return 0
        A[r], A[p] = A[p], A[r]
        for i in range(m):
            if i != r and A[i][c]: A[i] = [x ^ y for x, y in zip(A[i], A[r])]
        r += 1
    return 1
syms = []
idx = [(i, j) for i in range(m) for j in range(i, m)]
for bits in itertools.product((0, 1), repeat=len(idx)):
    S = [[0] * m for _ in range(m)]
    for b, (i, j) in zip(bits, idx): S[i][j] = S[j][i] = b
    S = tuple(tuple(r) for r in S)
    if det2(S) and any(S[i][i] for i in range(m)): syms.append(S)
print("nondegenerate non-alternating symmetric:", len(syms))
left = set(syms); orbits = []
while left:
    T = next(iter(left)); orb = {mm(mm(n, T), tp(n)) for n in N}
    orbits.append(orb); left -= orb
print("N(H)-orbits:", len(orbits), "sizes", sorted(len(o) for o in orbits))
# where do the two placements of t2-pair sit?  poly placement: T = I in poly coordinates.
# sd placement: coordinates c_i = Tr(x w_i) for a self-dual basis; in poly coordinates the target is
# T_sd = sum_i w_i w_i^T (w_i as poly-coordinate column vectors).
def tr(x):
    t = 0; y = x
    for _ in range(m): t ^= y; y = mul(y, y)
    return t
def vec(x): return [(x >> i) & 1 for i in range(m)]
# find a self-dual basis
def sdb():
    for ws in itertools.combinations(range(1, q), m):
        if all(tr(mul(a, b)) == (1 if a == b else 0) for a in ws for b in ws): return ws
w = sdb(); print("self-dual basis (poly ints):", w)
Tsd = tuple(tuple(sum(vec(wi)[i] & vec(wi)[j] for wi in w) & 1 for j in range(m)) for i in range(m))
for k, o in enumerate(orbits):
    tags = []
    if I in o: tags.append("poly placement (T=I)")
    if Tsd in o: tags.append("sd placement")
    if tags: print(" orbit", k, "size", len(o), ":", ", ".join(tags))
