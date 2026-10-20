# geo.py -- skeptic's own F_4 / PG(2,4) helpers (v-meta2-hermitian-compute). Shares no code with c2-hermitian.
# F_4 = {0,1,2,3}: 2 = alpha (root of x^2+x+1), 3 = alpha^2 = alpha+1. Addition = XOR. Multiplication via logs.
import itertools
EXP = [1, 2, 3]
LOG = {1: 0, 2: 1, 3: 2}
def fmul(a, b):
    if a == 0 or b == 0: return 0
    return EXP[(LOG[a] + LOG[b]) % 3]
def finv(a): return EXP[(-LOG[a]) % 3]
def fconj(a): return fmul(a, a)          # Frobenius x -> x^2 (the involution of F_4/F_2)
def vadd(u, v): return tuple(x ^ y for x, y in zip(u, v))
def vscale(c, u): return tuple(fmul(c, x) for x in u)
VEC3 = list(itertools.product(range(4), repeat=3))
def span(gens, n=3):
    S = {tuple([0] * n)}
    for g in gens:
        new = set()
        for v in S:
            for c in range(4):
                new.add(vadd(v, vscale(c, g)))
        S = new
    return frozenset(S)
def dim_of(S):
    n = len(S); d = 0
    while 4 ** d < n: d += 1
    assert 4 ** d == n
    return d
# all subspaces of F_4^3
SUBS = sorted({span(list(g)) for k in range(4) for g in itertools.combinations(VEC3, k)}, key=lambda S: (dim_of(S), sorted(S)))
assert len(SUBS) == 1 + 21 + 21 + 1
IDX = {S: i for i, S in enumerate(SUBS)}
DIM = [dim_of(S) for S in SUBS]
def meet(i, j): return IDX[SUBS[i] & SUBS[j]]
def join(i, j): return IDX[span(list(SUBS[i] | SUBS[j]))]
D = [[DIM[join(i, j)] - DIM[meet(i, j)] for j in range(44)] for i in range(44)]
COVER = [(i, j) for i in range(44) for j in range(44) if SUBS[i] < SUBS[j] and DIM[j] == DIM[i] + 1]
def herm_perp(i):     # Hermitian: {z : sum conj(x_k) z_k = 0 for all x in s}
    S = SUBS[i]
    out = []
    for z in VEC3:
        ok = True
        for x in S:
            t = 0
            for k in range(3): t ^= fmul(fconj(x[k]), z[k])
            if t: ok = False; break
        if ok: out.append(z)
    return IDX[frozenset(out)]
def bil_perp(i):
    S = SUBS[i]
    out = [z for z in VEC3 if all((fmul(x[0], z[0]) ^ fmul(x[1], z[1]) ^ fmul(x[2], z[2])) == 0 for x in S)]
    return IDX[frozenset(out)]
def frame_of_basis(basis):
    return IDX[span([tuple(v) for v in basis])]
