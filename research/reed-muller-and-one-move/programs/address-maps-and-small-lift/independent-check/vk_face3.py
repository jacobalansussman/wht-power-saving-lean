"""v-rm-theory: independent recomputation of the block-3 face certificate (Theorem C3 inputs).
Own code. Exact ranks over Q(i) via realification (integer matrices, exact Fraction elimination),
plus checks: N kept by all 8 shifts; caps per direction; demand of the kernel in every unordered basis
(all 28) and every ordered basis (168); which address maps keep N.
dir(w) = (1+i)/2 I + (1-i)/2 X^w  ->  2 dir(w) = (1+i) I + (1-i) X^w (Gaussian integers).
Third exterior power: Lambda^3(g)[S,T] = det g[S,T] (S,T sorted 3-subsets of addresses).
"""
import itertools
from fractions import Fraction

m = 3
n = 1 << m
addrs = list(range(n))
triples = list(itertools.combinations(addrs, 3))
def bit(x, j): return (x >> j) & 1
def is_face(T):
    return any(len({bit(x, j) for x in T}) == 1 for j in range(m))
F = [T for T in triples if is_face(T)]
NF = [T for T in triples if not is_face(T)]
assert len(F) == 24 and len(NF) == 32

# Gaussian integers as Python complex with int parts -> use tuples (a,b)
def gmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gadd(x, y): return (x[0]+y[0], x[1]+y[1])
def gsub(x, y): return (x[0]-y[0], x[1]-y[1])
Z = (0, 0)

def matmul(A, B):
    return [[ (lambda r,c: (sum(gmul(A[r][k],B[k][c])[0] for k in range(n)), sum(gmul(A[r][k],B[k][c])[1] for k in range(n))))(r,c)
              for c in range(n)] for r in range(n)]

def shift(a):
    M = [[Z]*n for _ in range(n)]
    for x in range(n): M[x ^ a][x] = (1, 0)
    return M
def dir2(w):  # 2*dir(w)
    M = [[Z]*n for _ in range(n)]
    for x in range(n):
        M[x][x] = gadd(M[x][x], (1, 1))
        M[x ^ w][x] = gadd(M[x ^ w][x], (1, -1))
    return M
def addrmap(h):  # h: 3x3 over F2 as list of column images h(e_j); P e_x = e_{hx}
    def app(x):
        y = 0
        for j in range(m):
            if bit(x, j): y ^= h[j]
        return y
    M = [[Z]*n for _ in range(n)]
    for x in range(n): M[app(x)][x] = (1, 0)
    return M

def det3(A, S, T):
    a = [[A[s][t] for t in T] for s in S]
    tot = Z
    for p in itertools.permutations(range(3)):
        sgn = 1
        for i in range(3):
            for j in range(i+1, 3):
                if p[i] > p[j]: sgn = -sgn
        prod = (1, 0)
        for i in range(3): prod = gmul(prod, a[i][p[i]])
        tot = gadd(tot, prod) if sgn == 1 else gsub(tot, prod)
    return tot

def rank_gauss(rows):
    """exact rank over Q(i) of a matrix with Gaussian-integer entries: realify -> rank/2."""
    R = []
    for r in rows:
        R.append([Fraction(e[0]) for e in r] + [Fraction(-e[1]) for e in r])
        R.append([Fraction(e[1]) for e in r] + [Fraction(e[0]) for e in r])
    rk = 0
    ncol = len(R[0]) if R else 0
    for c in range(ncol):
        piv = None
        for r in range(rk, len(R)):
            if R[r][c] != 0: piv = r; break
        if piv is None: continue
        R[rk], R[piv] = R[piv], R[rk]
        for r in range(len(R)):
            if r != rk and R[r][c] != 0:
                f = R[r][c] / R[rk][c]
                R[r] = [x - f*y for x, y in zip(R[r], R[rk])]
        rk += 1
    assert rk % 2 == 0
    return rk // 2

def leak(A, rowset=NF, colset=F):
    """dim(N + gN) - dim N for coordinate subspace N = span(F): rank of Lambda^3(g)[NF, F]."""
    return rank_gauss([[det3(A, S, T) for T in colset] for S in rowset])

out = []
# kept by shifts
for a in range(n):
    out.append(f"shift {a}: leak {leak(shift(a))}")
caps = {}
for w in range(1, n):
    caps[w] = leak(dir2(w))
out.append("caps per direction: " + str(caps))
# bases
vecs = list(range(1, n))
def indep(b):
    b1, b2, b3 = b
    span = {0, b1, b2, b1 ^ b2}
    return b3 not in span and b1 != b2
bases = [b for b in itertools.combinations(vecs, 3) if indep(b)]
dems = {}
for b in bases:
    K = matmul(matmul(dir2(b[0]), dir2(b[1])), dir2(b[2]))
    dems[b] = leak(K)
out.append(f"unordered bases: {len(bases)}; demands histogram: " + str({v: list(dems.values()).count(v) for v in set(dems.values())}))
# order independence (A abelian) spot check on 6 orderings of one non-standard basis
b = (3, 5, 6) if (3,5,6) in dems else bases[5]
vals = []
for p in itertools.permutations(b):
    K = matmul(matmul(dir2(p[0]), dir2(p[1])), dir2(p[2]))
    vals.append(leak(K))
out.append(f"orderings of {b}: {vals}")
# address maps keeping N
keep = 0; leaks = {}
for c in itertools.product(vecs, repeat=3):
    if not indep(c): continue
    L = leak(addrmap(list(c)))
    leaks[L] = leaks.get(L, 0) + 1
    if L == 0: keep += 1
out.append(f"address maps keeping N: {keep}; leak histogram over 168 maps: {dict(sorted(leaks.items()))}")
# rank of N itself sanity, and that dir caps are <= 8 also when two moves combine (rank-2 block)
c2 = {}
for w1 in range(1, n):
    for w2 in range(w1+1, n):
        c2[(w1, w2)] = leak(matmul(dir2(w1), dir2(w2)))
out.append("rank-2 block caps histogram: " + str({v: list(c2.values()).count(v) for v in set(c2.values())}))
print("\n".join(out))
