# v-xc3-cert3-compute: independent recomputation of xc3-cert3's alpha-game lift certificate (section 3.3).
# Own code.  Frames = generators (totally isotropic 3-spaces) of the Hermitian form h(x,y) = sum x_i y_i^2 on F_4^6
# (Fermat coordinates), enumerated directly from all 3x6 RREF matrices over F_4.
# Lift: Y_U in Q(omega)^6, 3-dim, reducing mod 2 (omega -> alpha) to U.
#   (a) Fermat lift: the 405 planes x_a = -zeta x_b of the Fermat cubic for the frames they reduce to,
#       Teichmueller lift of the F_4-RREF for the others;  (b) all-Teichmueller-RREF lift.
# Exact ranks over Q(omega) (Fractions) on every edge (distance-1 pair): excess = dim(Y_U + Y_V) - 3 - 1.
import itertools, sys
from fractions import Fraction as Fr

# F4: 0, 1, 2 = alpha, 3 = alpha^2 = alpha + 1; addition = xor
EXP = [1, 2, 3]; LOG = {1: 0, 2: 1, 3: 2}
def fm(a, b): return 0 if a == 0 or b == 0 else EXP[(LOG[a] + LOG[b]) % 3]
def fconj(a): return fm(a, a)
def finv(a): return EXP[(-LOG[a]) % 3]

def herm(x, y):  # sum x_i * conj(y_i)
    s = 0
    for a, b in zip(x, y): s ^= fm(a, fconj(b))
    return s

def enc(v): return sum(c << (2 * i) for i, c in enumerate(v))

def span_set(rows):
    S = set()
    for cs in itertools.product(range(4), repeat=len(rows)):
        v = [0] * 6
        for c, r in zip(cs, rows):
            if c:
                for i in range(6): v[i] ^= fm(c, r[i])
        S.add(enc(v))
    return frozenset(S)

def rref_f4(rows):
    M = [list(r) for r in rows]; piv = []; r = 0
    for c in range(6):
        pr = next((i for i in range(r, len(M)) if M[i][c]), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        inv = finv(M[r][c]); M[r] = [fm(inv, v) for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [M[i][k] ^ fm(f, M[r][k]) for k in range(6)]
        piv.append(c); r += 1
    return tuple(tuple(x) for x in M[:r]), tuple(piv)

# ---- enumerate all generators via RREF matrices ----
frames = {}   # span_set -> rref rows
for piv in itertools.combinations(range(6), 3):
    free = [(i, c) for i in range(3) for c in range(6) if c > piv[i] and c not in piv]
    for vals in itertools.product(range(4), repeat=len(free)):
        M = [[0] * 6 for _ in range(3)]
        for i in range(3): M[i][piv[i]] = 1
        for (i, c), v in zip(free, vals): M[i][c] = v
        ok = all(herm(M[i], M[j]) == 0 for i in range(3) for j in range(i, 3))
        if ok:
            frames[span_set(M)] = tuple(tuple(r) for r in M)
F = list(frames.keys()); idx = {U: k for k, U in enumerate(F)}
out = []
def say(s): print(s, flush=True); out.append(s)
say(f'generators of H(5,4) (Fermat form) by direct RREF enumeration: {len(F)} (want 891)')

def dist(U, V):
    n = len(U & V); d = {1: 0, 4: 1, 16: 2, 64: 3}[n]
    return 3 - d

# neighbour lists and distance spheres
nbr = [[] for _ in F]; sph_ok = True; sphere_hist = {}
for a in range(len(F)):
    cnt = [0, 0, 0, 0]
    for b in range(len(F)):
        d = dist(F[a], F[b]); cnt[d] += 1
        if d == 1: nbr[a].append(b)
    sphere_hist[tuple(cnt)] = sphere_hist.get(tuple(cnt), 0) + 1
say(f'distance spheres (0,1,2,3) -> number of frames: {sphere_hist}')
edges = [(a, b) for a in range(len(F)) for b in nbr[a] if a < b]
say(f'edges (unordered distance-1 pairs): {len(edges)}')

# ---- Q(omega) arithmetic: (a, b) = a + b w, w^2 = -1 - w ----
def qm(x, y):
    a, b = x; c, d = y
    return (a * c - b * d, a * d + b * c - b * d)
def qs(x, y): return (x[0] - y[0], x[1] - y[1])
def qinv(x):
    a, b = x; n = a * a - a * b + b * b
    return (Fr(a - b) / n, Fr(-b) / n)
Z0 = (Fr(0), Fr(0))
def qrank(rows):
    M = [[(Fr(e[0]), Fr(e[1])) for e in r] for r in rows]; r = 0
    for c in range(len(M[0])):
        pr = next((i for i in range(r, len(M)) if M[i][c] != Z0), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        inv = qinv(M[r][c]); M[r] = [qm(inv, v) for v in M[r]]
        for i in range(r + 1, len(M)):
            if M[i][c] != Z0:
                f = M[i][c]; M[i] = [qs(M[i][k], qm(f, M[r][k])) for k in range(len(M[0]))]
        r += 1
    return r

TEICH = {0: (0, 0), 1: (1, 0), 2: (0, 1), 3: (-1, -1)}   # alpha -> w, alpha^2 -> w^2 = -1 - w
def red(q):  # reduce an element of Z[w] mod 2 to F4 (q integral)
    a, b = q; a %= 2; b %= 2
    return {(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 3}[(a, b)]

# ---- Fermat planes ----
def matchings(s):
    if not s: yield []; return
    a = s[0]
    for i in range(1, len(s)):
        b = s[i]; rest = s[1:i] + s[i + 1:]
        for m in matchings(rest): yield [(a, b)] + m
ZETA = [(1, 0), (0, 1), (-1, -1)]
fermat = {}
for m in matchings(list(range(6))):
    for zs in itertools.product(ZETA, repeat=3):
        rows = []
        for (a, b), z in zip(m, zs):
            v = [(0, 0)] * 6; v[b] = (1, 0); v[a] = (-z[0], -z[1])   # x_a = -z x_b  <=>  e_b - z e_a
            rows.append(v)
        U = span_set([[red(e) for e in r] for r in rows])
        assert U in frames, 'Fermat plane does not reduce to a generator'
        fermat.setdefault(U, []).append(rows)
say(f'Fermat planes: {sum(len(v) for v in fermat.values())}, distinct frames: {len(fermat)}')
isF = [F[k] in fermat for k in range(len(F))]
fn = {}
for k in range(len(F)):
    c = sum(1 for b in nbr[k] if isF[b]); fn.setdefault((isF[k], c), 0); fn[(isF[k], c)] += 1
say(f'(is Fermat, number of Fermat neighbours) -> count: {fn}')

def teich_rows(U): return [[TEICH[e] for e in r] for r in frames[U]]

def run(lift, name):
    hist = {}
    for a, b in edges:
        r = qrank(lift[a] + lift[b])
        assert qrank(lift[a]) == 3
        ex = r - 3 - 1
        hist[ex] = hist.get(ex, 0) + 1
    E = sum(k * v for k, v in hist.items())
    abar1 = 1 + Fr(E, len(edges))
    say(f'{name}: exact excess histogram {dict(sorted(hist.items()))}, E = {E}, abar_1 = {abar1}, bound 3/abar_1 = {Fr(3) / abar1} = {float(Fr(3) / abar1):.6f}')

liftA = [fermat[F[k]][0] if isF[k] else teich_rows(F[k]) for k in range(len(F))]
run(liftA, 'Fermat (405) + Teichmueller-RREF (486)')
liftB = [teich_rows(F[k]) for k in range(len(F))]
run(liftB, 'all Teichmueller-RREF (891)')
# perfection of the Fermat family on Fermat-Fermat pairs (all distances), exact
hist = {}
FF = [k for k in range(len(F)) if isF[k]]
for i, a in enumerate(FF):
    for b in FF[i + 1:]:
        d = dist(F[a], F[b]); dl = qrank(liftA[a] + liftA[b]) - 3
        hist.setdefault(d, {}); hist[d][dl] = hist[d].get(dl, 0) + 1
say(f'Fermat-Fermat unordered pairs: delta by distance (exact) {dict(sorted(hist.items()))}')
with open('out_lift891_check.txt', 'w') as fh: fh.write('\n'.join(out) + '\n')
