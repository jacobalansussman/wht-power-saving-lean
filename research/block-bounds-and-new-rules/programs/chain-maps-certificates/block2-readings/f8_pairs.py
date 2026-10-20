"""xc3-s1s2: is the group of letter symbols at lambda of order 7 (F8) transitive on ordered pairs of planes at each
distance?  (Needed in PROOF.md 7.4/7.5 to pass from an averaged frame family to a uniform one.)
Method: BFS transversal from L0 over the 4745 planes; Schreier generators t_(gU)^-1 g t_U of Stab(L0) (a random sample);
orbits of the subgroup they generate on the planes.  A subgroup with one orbit at distance 1 and one at distance 2
proves transitivity of the whole stabiliser, hence of the group on ordered pairs at each distance."""
import itertools, random
from collections import deque
random.seed(7)
q = 8; MOD = 0b1011
EXP = [0] * 14; LOG = [0] * 8
x = 1
for i in range(7):
    EXP[i] = x; LOG[x] = i; x <<= 1
    if x & 8: x ^= MOD
for i in range(7, 14): EXP[i] = EXP[i - 7]
def mul(a, b): return 0 if a == 0 or b == 0 else EXP[LOG[a] + LOG[b]]
def inv(a): return EXP[(7 - LOG[a]) % 7]
def pw(a, e): return 1 if e == 0 else EXP[(LOG[a] * e) % 7]
lam = EXP[1]
def ev(poly, sign=1):
    s = 0
    for d, c in poly.items():
        if c: s ^= pw(lam, (sign * d) % 7)
    return s
def mmul(A, B):
    return [[(lambda i, j: __import__('functools').reduce(lambda u, v: u ^ v, (mul(A[i][k], B[k][j]) for k in range(len(B))), 0))(i, j)
             for j in range(len(B[0]))] for i in range(len(A))]
def minv(A):
    n = len(A); M = [A[i][:] + [int(i == j) for j in range(n)] for i in range(n)]
    r = 0
    for c in range(n):
        p = next(i for i in range(r, n) if M[i][c]); M[r], M[p] = M[p], M[r]
        iv = inv(M[r][c]); M[r] = [mul(iv, t) for t in M[r]]
        for i in range(n):
            if i != r and M[i][c]:
                fc = M[i][c]; M[i] = [t ^ mul(fc, y) for t, y in zip(M[i], M[r])]
        r += 1
    return [row[n:] for row in M]
def addr(G):
    Gl = [[ev(G[i][j]) for j in range(2)] for i in range(2)]
    Gi = [[ev(G[i][j], -1) for j in range(2)] for i in range(2)]
    Git = minv(Gi); Gz = [[Git[j][i] for j in range(2)] for i in range(2)]
    s = [[0] * 4 for _ in range(4)]
    for i in range(2):
        for j in range(2): s[i][j] = Gl[i][j]; s[2 + i][2 + j] = Gz[i][j]
    return s
def phase(B):
    s = [[int(i == j) for j in range(4)] for i in range(4)]
    for i in range(2):
        for j in range(2): s[2 + i][j] = ev(B[i][j])
    return s
one = {0: 1}; zero = {}
GEN = [addr([[{0: 1, 1: 1}, zero], [zero, {0: 1, 1: 1}]]), addr([[one, {1: 1}], [zero, one]]),
       addr([[one, {1: 1}], [{1: 1}, {0: 1, 1: 1, 2: 1}]]), addr([[{0: 1, 1: 1, 2: 1}, zero], [zero, {0: 1, 1: 1, 2: 1}]]),
       addr([[{0: 1, 1: 1}, {1: 1}], [{1: 1}, one]]), phase([[zero, zero], [zero, one]]), phase([[zero, one], [one, zero]]),
       phase([[zero, zero], [zero, {1: 1, -1: 1}]]),
       [[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 1]], [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
       [[1, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], [[1, 0, 0, 0], [0, 1, 0, 1], [0, 0, 1, 0], [0, 0, 0, 1]],
       [[1, 0, 1, 1], [0, 1, 1, 1], [0, 0, 1, 0], [0, 0, 0, 1]], [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]]
GEN = GEN + [minv(g) for g in GEN]
def app(s, u): return tuple(__import__('functools').reduce(lambda a, b: a ^ b, (mul(s[i][k], u[k]) for k in range(4)), 0) for i in range(4))
def rref(rows):
    rows = [list(r) for r in rows]; r = 0; n = len(rows)
    for c in range(4):
        p = next((i for i in range(r, n) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        iv = inv(rows[r][c]); rows[r] = [mul(iv, t) for t in rows[r]]
        for i in range(n):
            if i != r and rows[i][c]:
                fc = rows[i][c]; rows[i] = [t ^ mul(fc, y) for t, y in zip(rows[i], rows[r])]
        r += 1
    return tuple(map(tuple, rows[:r]))
vecs = [v for v in itertools.product(range(q), repeat=4) if any(v)]
lead1 = [v for v in vecs if next(c for c in v if c) == 1]
planes = set()
for i, u in enumerate(lead1):
    for v in lead1[i + 1:]:
        r = rref([u, v])
        if len(r) == 2: planes.add(r)
planes = sorted(planes); idx = {p: i for i, p in enumerate(planes)}
print("planes:", len(planes), flush=True)
def act(s, i): return idx[rref([app(s, planes[i][0]), app(s, planes[i][1])])]
def dist(i, j): return len(rref(list(planes[i]) + list(planes[j]))) - 2
L0 = idx[rref([(0, 0, 1, 0), (0, 0, 0, 1)])]
I4 = [[int(i == j) for j in range(4)] for i in range(4)]
trans = {L0: I4}; dq = deque([L0])
while dq:
    a = dq.popleft()
    for g in GEN:
        b = act(g, a)
        if b not in trans:
            trans[b] = mmul(g, trans[a]); dq.append(b)
print("orbit of L0:", len(trans), flush=True)
# Schreier generators
keys = list(trans)
sch = []
for _ in range(60):
    U = random.choice(keys); g = random.choice(GEN)
    gU = act(g, U)
    h = mmul(minv(trans[gU]), mmul(g, trans[U]))
    assert act(h, L0) == L0
    sch.append(h)
perms = [[act(h, i) for i in range(len(planes))] for h in sch]
seen = [False] * len(planes); orbits = []
for st in range(len(planes)):
    if seen[st]: continue
    orb = [st]; seen[st] = True; k = 0
    while k < len(orb):
        a = orb[k]; k += 1
        for p in perms:
            b = p[a]
            if not seen[b]: seen[b] = True; orb.append(b)
    orbits.append(orb)
desc = sorted((len(o), sorted({dist(L0, i) for i in o})) for o in orbits)
print("orbits of the sampled stabiliser subgroup on planes (size, distances):", desc)
ok = len(orbits) == 3 and all(len(d) == 1 for _, d in desc)
print(("PASS" if ok else "FAIL") + " Stab(L0) transitive on planes at distance 1 and at distance 2 -> the symbol group is"
      " transitive on ordered pairs at each distance")
