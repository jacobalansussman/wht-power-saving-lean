# v-xc3-cert3-compute: independent check of xc3-cert3 Lemma 1 (the block-3 hook certificate is a frame family).
# Own code, nothing imported from xc3-cert3 or d-tof-cert.
# V = C^8 (basis e_x, x in F_2^3, bit j = qubit j), W = V (x) Lambda^4 V (dim 560), N = span{e_x (x) e_T : x in T} (dim 280).
# For a Clifford g: Y(g) = (g (x) Lambda^4 g) N.  Frame of g = { v : g Z_b g^-1 = c P_v } (8 labels, a Lagrangian of F_2^6).
# Part 1 (EXACT over Q, Bareiss on integer matrices): delta(L0, g L0) for g = H0, H0H1, H0H1H2 (scaled Hadamards, integer).
# Part 2 (mod two primes p = 1 mod 4): BFS over H_j, S_j, CNOT_jk; all 135 frames; Y well defined (several g per frame);
#   delta(L0, V) for all 134 V; delta on random pairs; monomial letters keep N; the direction D_z = (1+i) + (1-i) X_z.
import sys, itertools, random
import numpy as np

T4 = list(itertools.combinations(range(8), 4))           # 70 subsets
TIDX = {t: i for i, t in enumerate(T4)}
DIM = 8 * 70
NCOLS = [x * 70 + TIDX[t] for x in range(8) for t in T4 if x in t]
NSET = set(NCOLS)
NROWS_C = [r for r in range(DIM) if r not in NSET]
assert len(NCOLS) == 280

def det4(m):  # exact integer (or mod-reduced later) determinant of a 4x4 list-of-lists, Laplace
    a = m
    def det3(b):
        return (b[0][0]*(b[1][1]*b[2][2]-b[1][2]*b[2][1]) - b[0][1]*(b[1][0]*b[2][2]-b[1][2]*b[2][0])
                + b[0][2]*(b[1][0]*b[2][1]-b[1][1]*b[2][0]))
    s = 0
    for j in range(4):
        minor = [[a[i][k] for k in range(4) if k != j] for i in range(1, 4)]
        s += (-1)**j * a[0][j] * det3(minor)
    return s

def wedge4(g, p=None):  # (Lambda^4 g)[S,T] = det g[S,T]; g as python int lists
    L = [[0]*70 for _ in range(70)]
    for si, S in enumerate(T4):
        for ti, T in enumerate(T4):
            d = det4([[g[r][c] for c in T] for r in S])
            L[si][ti] = d % p if p else d
    return L

def kron_cols(g, L, cols, p=None):
    # columns (x,T) of g (x) L: entry at row (y,S) = g[y][x] * L[S][T]
    out = []
    for c in cols:
        x, t = divmod(c, 70)
        col = [0]*DIM
        for y in range(8):
            gyx = g[y][x]
            if gyx == 0: continue
            for s in range(70):
                v = gyx * L[s][t]
                col[y*70+s] = v % p if p else v
        out.append(col)
    return out  # list of columns

def bareiss_rank(M):  # exact rank of integer matrix (list of rows), fraction-free
    A = [row[:] for row in M]
    nr, nc = len(A), len(A[0])
    r = 0; prev = 1
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if A[i][c] != 0: piv = i; break
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        prc = A[r][c]; rowr = A[r]
        for i in range(r+1, nr):
            ai = A[i]; aic = ai[c]
            if aic == 0:
                # still need the Bareiss scaling: ai = (prc*ai - aic*rowr)/prev = prc*ai/prev
                A[i] = [(prc*v)//prev for v in ai]
            else:
                A[i] = [(prc*ai[k] - aic*rowr[k])//prev for k in range(nc)]
        prev = prc; r += 1
        if r == nr: break
    return r

def rank_mod(M, p):  # M numpy int64 array, entries in [0,p)
    A = M.copy() % p
    nr, nc = A.shape; r = 0
    for c in range(nc):
        if r == nr: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        piv = r + nz[0]
        if piv != r: A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p-2, p)
        A[r] = (A[r] * inv) % p
        rows = r + 1 + np.nonzero(A[r+1:, c])[0]
        if len(rows):
            f = A[rows, c].reshape(-1, 1)
            A[rows] = (A[rows] - (f * A[r]) % p) % p
        r += 1
    return r

# ---------- basic gates (integer / Gaussian-integer; scaled Hadamard) ----------
def Hgate(j):
    g = [[0]*8 for _ in range(8)]
    for x in range(8):
        for y in range(8):
            if (x ^ y) & ~(1 << j) & 7: continue
            g[y][x] = -1 if ((x >> j) & 1) and ((y >> j) & 1) else 1
    return g

def matmul_int(a, b, p=None):
    n = len(a)
    r = [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    if p: r = [[v % p for v in row] for row in r]
    return r

ID8 = [[int(i == j) for j in range(8)] for i in range(8)]

def delta_L0_exact(g):
    L = wedge4(g)
    cols = kron_cols(g, L, NCOLS)
    sub = [[cols[k][r] for k in range(280)] for r in NROWS_C]   # complement rows x 280 columns
    return bareiss_rank(sub)

out = []
def say(s):
    print(s, flush=True); out.append(s)

if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if part in ('exact', 'all'):
        H0, H1, H2 = Hgate(0), Hgate(1), Hgate(2)
        for name, g in [('H0 (distance 1)', H0), ('H0H1 (distance 2)', matmul_int(H0, H1)),
                        ('H0H1H2 (distance 3)', matmul_int(matmul_int(H0, H1), H2))]:
            d = delta_L0_exact(g)
            say(f'EXACT over Q: delta(L0, g L0) for g = {name}: {d}')

# ---------------------------------------------------------------- mod-p part ----------------------------------------
def is_prime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True

def primes_1mod4(count, below=2**31):
    ps = []; n = below - 1
    while len(ps) < count:
        if n % 4 == 1 and is_prime(n): ps.append(n)
        n -= 1
    return ps

def sqrt_m1(p):
    for a in range(2, 1000):
        r = pow(a, (p - 1) // 4, p)
        if r * r % p == p - 1: return r
    raise ValueError

def run_modp(p, seed, npairs, lines):
    I = sqrt_m1(p)
    def arr(g): return np.array(g, dtype=np.int64) % p
    def Hg(j): return arr(Hgate(j))
    def Sg(j): return arr([[ (pow(I, (x >> j) & 1, p) if x == y else 0) for x in range(8)] for y in range(8)])
    def perm(fn):
        g = np.zeros((8, 8), dtype=np.int64)
        for x in range(8): g[fn(x), x] = 1
        return g
    def CNOT(c, t): return perm(lambda x: x ^ (((x >> c) & 1) << t))
    def Xg(j): return perm(lambda x: x ^ (1 << j))
    def TOF(): return perm(lambda x: x ^ ((((x >> 0) & (x >> 1)) & 1) << 2))
    def CZ(a, b): return arr([[ ((p - 1) if (((x >> a) & (x >> b)) & 1) else 1) if x == y else 0 for x in range(8)] for y in range(8)])
    def pauli(v):  # v = (a | b), a = low 3 bits (X part), b = high 3 bits (Z part):  X^a Z^b e_x = (-1)^{b.x} e_{x^a}
        a, b = v & 7, v >> 3
        g = np.zeros((8, 8), dtype=np.int64)
        for x in range(8):
            g[x ^ a, x] = (p - 1) if bin(b & x).count('1') % 2 else 1
        return g
    PAULI = [pauli(v) for v in range(64)]
    def mm(a, b): return (a @ b) % p
    def proportional(A, B):
        idx = np.argwhere(B != 0)
        if len(idx) == 0: return False
        r, c = idx[0]
        if A[r, c] == 0: return False
        cst = A[r, c] * pow(int(B[r, c]), p - 2, p) % p
        return np.array_equal(A % p, (cst * B) % p)
    def frame(g):
        labs = []
        for j in range(3):
            A = mm(g, PAULI[1 << (3 + j)])        # g Z_j
            hit = [v for v in range(64) if proportional(A, mm(PAULI[v], g))]
            assert len(hit) == 1, hit
            labs.append(hit[0])
        span = set()
        for c in range(8):
            v = 0
            for j in range(3):
                if (c >> j) & 1: v ^= labs[j]
            span.add(v)
        assert len(span) == 8
        return frozenset(span)
    def sympl(u, v):
        return (bin((u & 7) & (v >> 3)).count('1') + bin((u >> 3) & (v & 7)).count('1')) % 2
    def dist(U, V): return 3 - (len(U & V).bit_length() - 1)
    def Y(g):
        L = np.array(wedge4([[int(v) for v in row] for row in g], p), dtype=np.int64)
        G = np.kron(g, L) % p
        return G[:, NCOLS]
    ID = arr(ID8)
    L0 = frame(ID)
    assert all(sympl(u, v) == 0 for u in L0 for v in L0)
    gens = [('H%d' % j, Hg(j)) for j in range(3)] + [('S%d' % j, Sg(j)) for j in range(3)] + \
           [('CX%d%d' % (c, t), CNOT(c, t)) for c in range(3) for t in range(3) if c != t]
    reps = {L0: [ID]}; order = [L0]; q = [L0]
    while q:
        U = q.pop(0)
        for nm, h in gens:
            g2 = mm(h, reps[U][0])
            V = frame(g2)
            if V not in reps:
                reps[V] = [g2]; order.append(V); q.append(V)
            elif len(reps[V]) < 3:
                reps[V].append(g2)
    lines.append(f'[p={p}] frames reached by BFS over H_j, S_j, CNOT: {len(reps)} (want 135)')
    assert all(all(sympl(u, v) == 0 for u in U for v in U) for U in reps)
    # well-definedness of Y_U: every extra representative gives the same subspace
    Ycache = {U: Y(reps[U][0]) for U in order}
    wd_checked = 0; wd_bad = 0
    for U in order:
        r0 = rank_mod(Ycache[U], p)
        if r0 != 280: wd_bad += 1
        for g2 in reps[U][1:]:
            r = rank_mod(np.hstack([Ycache[U], Y(g2)]), p)
            wd_checked += 1
            if r != 280: wd_bad += 1
    lines.append(f'[p={p}] Y_U well defined: dim Y_U = 280 for all 135 frames; {wd_checked} second/third representatives span the same Y_U; failures {wd_bad}')
    # monomial letters keep N
    compl = np.array(NROWS_C)
    mono = [('S0', Sg(0)), ('X1', Xg(1)), ('CX01', CNOT(0, 1)), ('CZ12', CZ(1, 2)), ('Toffoli', TOF())]
    res = []
    for nm, g in mono:
        res.append((nm, rank_mod(Y(g)[compl, :], p)))
    lines.append(f'[p={p}] monomial letters: rank of (g N) outside N = {res} (0 means g keeps N)')
    # delta(L0, V) for all V
    hist = {}
    for V in order:
        d = dist(L0, V)
        dl = rank_mod(Ycache[V][compl, :], p)
        hist.setdefault(d, {}); hist[d][dl] = hist[d].get(dl, 0) + 1
    lines.append(f'[p={p}] delta(L0, V) by distance (mod p, a lower bound for the char-0 value): {dict(sorted(hist.items()))}')
    # direction D_z = (1+i) I + (1-i) X_z, z = e_0 and z = e_0 + e_2
    for z in (1, 5):
        Xz = perm(lambda x, z=z: x ^ z)
        D = ((1 + I) % p * ID + (1 - I) % p * Xz) % p
        V = frame(D)
        lines.append(f'[p={p}] direction D_z, z={z:03b}: distance(L0, D L0) = {dist(L0, V)}, delta(L0, D L0) = {rank_mod(Y(D)[compl, :], p)}')
    # random ordered pairs (U, V)
    rng = random.Random(seed)
    ph = {}
    for _ in range(npairs):
        U, V = rng.sample(order, 2)
        r = rank_mod(np.hstack([Ycache[U], Ycache[V]]), p) - 280
        d = dist(U, V)
        ph.setdefault(d, {}); ph[d][r] = ph[d].get(r, 0) + 1
    lines.append(f'[p={p}] random ordered pairs (U, V): delta by distance {dict(sorted(ph.items()))}')

if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if part in ('modp', 'all', 'modp1', 'modp2'):
        ps = primes_1mod4(2)
        if part == 'modp1': ps = ps[:1]
        if part == 'modp2': ps = ps[1:]
        npairs = int(sys.argv[2]) if len(sys.argv) > 2 else 60
        for k, p in enumerate(ps):
            lines = []
            run_modp(p, 1000 + k + (1 if part == 'modp2' else 0), npairs, lines)
            for s in lines: say(s)
    with open('out_hook_frames_%s.txt' % part, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
