# vhook.py -- v-w7-mixed3: own rebuild of the hook single-letter defect ranks and stable ranks (Lemma Q3's input).
# X = V (x) Lambda^4 V, N = span{e_x (x) e_T : x in T}, Y' = span{e_x (x) e_T : x not in T}; rho(g) = g (x) Lambda^4 g.
# D_z = ((1+i) I + (1-i) X_z)/2.  E_T = N cap rho(T)^-1 N.  A_b = pi_{Y'} rho(D_1 b)|_{E_T}.
# Letters: engine = affine permutation x -> Ax + c times phase i^(c0 + a.x + 2 q(x)) (q quadratic, no squares);
#          widened = any permutation times Z4 phases.  Classes: D_1 b up to a left monomial (rows up to scalar, sorted).
import itertools, sys, time, collections, numpy as np
t0 = time.time()
n = 8; S4 = list(itertools.combinations(range(n), 4)); BAS = [(x, T) for x in range(n) for T in S4]
NI = np.array([i for i, (x, T) in enumerate(BAS) if x in T]); OI = np.array([i for i, (x, T) in enumerate(BAS) if x not in T])
S4a = np.array(S4)
def lam4(g):
    # 70x70 matrix of 4x4 minors: (Lambda^4 g)[I, J] = det g[I, J]
    sub = g[S4a[:, None, :, None], S4a[None, :, None, :]]      # 70 x 70 x 4 x 4
    return np.linalg.det(sub)
def rho(g): return np.kron(g, lam4(g))
def Xz(z):
    M = np.zeros((8, 8)); 
    for x in range(8): M[x ^ z, x] = 1
    return M
def Dz(z): return ((1 + 1j) * np.eye(8) + (1 - 1j) * Xz(z)) / 2
D1, D2, D4 = Dz(1), Dz(2), Dz(4)
TG = {'K': D1 @ D2 @ D4, 'R': D1 @ D2}
def nullsp(M, tol=1e-9):
    u, s, vh = np.linalg.svd(M); r = int((s > tol * max(1, s[0])).sum()); return vh[r:].conj().T
E = {}
for t, T in TG.items():
    RT = rho(T); E[t] = nullsp(RT[np.ix_(OI, NI)])       # coefficient vectors on N coordinates
cap = np.linalg.matrix_rank(rho(D1)[np.ix_(OI, NI)], tol=1e-8)
print('dim N', len(NI), 'cap', cap, 'dim E_K', E['K'].shape[1], 'dim E_R', E['R'].shape[1], ' (expect 280, 80, 56, 136)', flush=True)
RD1 = None
def mono(perm, ph):
    M = np.zeros((8, 8), complex)
    for x in range(8): M[perm[x], x] = 1j ** (ph[x] % 4)
    return M
def rand_widened(rng): return mono(list(rng.permutation(8)), list(rng.integers(0, 4, 8)))
GL = []
for cols in itertools.product(range(1, 8), repeat=3):
    if len({cols[0], cols[1], cols[2], cols[0] ^ cols[1], cols[0] ^ cols[2], cols[1] ^ cols[2], cols[0] ^ cols[1] ^ cols[2]}) == 7 and 0 not in {cols[0] ^ cols[1] ^ cols[2]}:
        GL.append(cols)
assert len(GL) == 168, len(GL)
def rand_engine(rng):
    A = GL[rng.integers(168)]; c = int(rng.integers(8))
    perm = []
    for x in range(8):
        y = c
        for j in range(3):
            if x >> j & 1: y ^= A[j]
        perm.append(y)
    a = rng.integers(0, 4, 3); q = rng.integers(0, 2, 3); c0 = int(rng.integers(4))
    ph = []
    for x in range(8):
        b = [(x >> j) & 1 for j in range(3)]
        ph.append(int(c0 + sum(int(a[j]) * b[j] for j in range(3)) + 2 * (q[0] * b[0] * b[1] + q[1] * b[0] * b[2] + q[2] * b[1] * b[2])) % 4)
    return mono(perm, ph)
def canon(M):
    rows = []
    for r in M:
        nz = np.nonzero(np.abs(r) > 1e-9)[0]; r = r / r[nz[0]]
        rows.append(tuple(np.round(r.real, 6) + 0.0) + tuple(np.round(r.imag, 6) + 0.0))
    return tuple(sorted(rows))
game = sys.argv[1]; NS = int(sys.argv[2]); seed = int(sys.argv[3])
def rand_cont(rng):
    M = np.zeros((8, 8), complex); perm = rng.permutation(8); ph = np.exp(2j * np.pi * rng.random(8))
    for x in range(8): M[perm[x], x] = ph[x]
    return M
rng = np.random.default_rng(seed); gen = {'engine': rand_engine, 'widened': rand_widened, 'cont': rand_cont}[game]
reps = {}
for _ in range(NS):
    b = gen(rng); k = canon(D1 @ b)
    if k not in reps: reps[k] = b
print(game, 'samples', NS, 'distinct left classes of D1 b:', len(reps), f'({time.time()-t0:.0f}s)', flush=True)
tab = collections.Counter(); worst = {'R': 99.0, 'K': 99.0}; worst_non8 = 99.0; gaps = [1e9, 0.0]
for k, b in reps.items():
    Rb = rho(D1 @ b); row = []
    for t in ('R', 'K'):
        A = Rb[np.ix_(OI, NI)] @ E[t]; s = np.linalg.svd(A, compute_uv=False)
        r = int((s > 1e-7 * max(1, s[0])).sum())
        if r: gaps[0] = min(gaps[0], s[r - 1] / s[0])
        if r < len(s): gaps[1] = max(gaps[1], s[r] / max(1, s[0]))
        sr = float((s ** 2).sum() / s[0] ** 2) if r else 0.0
        row.append((r, round(sr, 3)))
        if r: worst[t] = min(worst[t], sr)
        if t == 'R' and r and r != 8: worst_non8 = min(worst_non8, sr)
    tab[tuple(row)] += 1
print('((rank R, srank R), (rank K, srank K)) -> count')
for k, v in sorted(tab.items()): print('   ', k, v)
print('least srank R', worst['R'], '; least srank R among rank != 8:', worst_non8, '; least srank K (nonzero):', worst['K'])
print(f'min kept relative sv {gaps[0]:.3e}, max dropped relative sv {gaps[1]:.3e}  ({time.time()-t0:.0f}s)')
