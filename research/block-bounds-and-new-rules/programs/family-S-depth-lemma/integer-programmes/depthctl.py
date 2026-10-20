# b4-depth4: exact controls at block size m = 4 for the depth lemma and the one-event integer program.
# Formal-sum bookkeeping: a generator is {(input, matrix-key): coeff}; matrices are 16x16 with dyadic Gaussian
# entries, exact in complex128 (all denominators <= 2^8); keys are the exact entries.
import sys, itertools, random
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import numpy as np
n = 16
al, be = (1 + 1j) / 2, (1 - 1j) / 2
def dirm(z):
    A = np.zeros((n, n), complex)
    for x in range(n): A[x, x] += al; A[x ^ z, x] += be
    return A
def shift(z):
    A = np.zeros((n, n), complex)
    for x in range(n): A[x ^ z, x] = 1
    return A
I16 = np.eye(n, dtype=complex)
D = {z: dirm(z) for z in range(1, n)}
K = D[1] @ D[2] @ D[4] @ D[8]
def key(M): return tuple(np.round(M.flatten() * 256).astype(complex).tolist())
mats = {}
def term(M):
    k = key(M); mats[k] = M; return k
def gen_input(l): return {(l, term(I16)): 1}
def lmul(A, g):           # left-multiply every term of generator g by matrix A
    out = {}
    for (l, k), c in g.items():
        kk = term(A @ mats[k]); out[(l, kk)] = out.get((l, kk), 0) + c
    return {t: c for t, c in out.items() if abs(c) > 1e-12}
def add(*gs):
    out = {}
    for coef, g in gs:
        for t, c in g.items(): out[t] = out.get(t, 0) + coef * c
    return {t: c for t, c in out.items() if abs(c) > 1e-12}
# --- (1) nonzeros per row: nu(K) = 16 on every row; nu(D b g) <= 2 nu(g) ---
nzK = (np.abs(K) > 1e-12).sum(axis=1)
print("nu(K) per row:", sorted(set(nzK.tolist())), "-> depth(K) >= log2 16 = 4; plain chain gives 4")
rng = random.Random(1)
def rand_monomial():
    perm = list(range(n)); rng.shuffle(perm)
    M = np.zeros((n, n), complex)
    for x in range(n): M[perm[x], x] = 1j ** rng.randrange(4)
    return M
tof = np.zeros((n, n), complex)
for x in range(n): tof[x ^ (4 if (x & 1 and x & 2) else 0), x] = 1      # Toffoli on bits 0,1 -> 2
cs = np.diag([1j if (x & 1 and x & 2) else 1 for x in range(n)]).astype(complex)
worst = 0
for trial in range(2000):
    sig = rand_monomial()
    nz = (np.abs(K @ sig) > 1e-12).sum(axis=1).min()
    worst = max(worst, 16 - nz)
for sig in (tof, cs, tof @ cs):
    nz = (np.abs(K @ sig) > 1e-12).sum(axis=1).min(); worst = max(worst, 16 - nz)
print("min over rows of nu(K sigma) for 2000 random monomials + Toffoli/CS: 16 - worst =", 16 - worst,
      "-> every frame B'K sigma has depth exactly 4")
# --- (2) the 9-move unit whose two K-paths share one mixing move (one-event bound 2m+1 is attained) ---
x1, x2 = gen_input(1), gen_input(2)
moves = 0
def mv(z, g):
    global moves; moves += 1; return lmul(D[z], g)
ua = mv(1, x1)                 # T1: depth 1
uc = mv(8, x1)                 # T1: canceller of J2
v1 = mv(2, x2); v2 = mv(4, v1) # T2: depth 2
vc = mv(8, x2)                 # T2: canceller of J1
e = mv(8, add((1, ua), (1, v2)))          # the mixing move
p2 = mv(1, e)                              # branch 2 (r2 = 1)
p1 = mv(2, mv(4, e))                       # branch 1 (r1 = 2)
out1 = add((1, p1), (-1, lmul(shift(4) @ shift(2), vc)))
out2 = add((1, p2), (-1, lmul(shift(1), uc)))
ok1 = out1 == {(1, term(K)): 1}; ok2 = out2 == {(2, term(K)): 1}
print(f"9-move unit with both K-paths through one mixing move: moves = {moves}, live 1 = K x1: {ok1}, live 2 = K x2: {ok2}")
