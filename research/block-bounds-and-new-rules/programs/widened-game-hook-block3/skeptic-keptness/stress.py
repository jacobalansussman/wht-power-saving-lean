# v28-kept: randomized stress test of Theorem R with the hook certificate, EXACT mod p (no floats).
# Module W = V (x) Lambda^4 V over F_p (p = 1 mod 4, i = sqrt(-1) mod p), rho(g) = g (x) Lambda^4 g computed from
# 4x4 minors (my code; cross-checked against g^(x)5 on V^(x)5 for 3 matrices).  N = {x in T} coordinates.
# Registers: W_live live + S scratch; state M[r][l] : N -> W (register r as a function of input y_l in N).
# Potential Phi = dim(A + N^R) - R dim N = rank of the out-of-N projection of the stacked state.
# Steps: free letters = random monomials (random permutation of the 8 addresses x random F_p^* diagonal: covers
# Toffoli, CS, arbitrary diagonals), named Toffoli / CCZ / CS, engine letters; mixing = random integer R x R matrix
# (often singular); paid = b dir(z) b^-1 (rank 1) or b dir(z1) dir(z2) b^-1 (rank 2), b random monomial.
# Checks: free and mixing steps never raise Phi; rank-1 raises <= 80; rank-2 raises <= 160 (claim: f(B2) = 144).
# Controls: a 'free' Hadamard on bit 0 (not monomial) must break keptness; plain schedule must end with Phi >= 224 W.
# Usage: python3 -I stress.py SEED STEPS
import sys, itertools, time
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import numpy as np
seed, STEPS = int(sys.argv[1]), int(sys.argv[2])
rng = np.random.default_rng(seed)
n = 8
p = 1000033            # prime, 1 mod 4
assert all(p % d for d in range(2, 1001)) and p % 4 == 1
I = next(pow(g, (p - 1) // 4, p) for g in range(2, 100) if pow(g, (p - 1) // 2, p) == p - 1)
inv2 = pow(2, p - 2, p)
S4 = list(itertools.combinations(range(n), 4)); sidx = {T: i for i, T in enumerate(S4)}
basis = [(x, T) for x in range(n) for T in S4]
Nidx = np.array([i for i, (x, T) in enumerate(basis) if x in T]); Oidx = np.array([i for i, (x, T) in enumerate(basis) if x not in T])
DW, DN = len(basis), len(Nidx)
def det4(A):
    s = 0
    for pi in itertools.permutations(range(4)):
        sg = 1
        for a in range(4):
            for b in range(a + 1, 4):
                if pi[a] > pi[b]: sg = -sg
        t = sg
        for a in range(4): t = t * int(A[a, pi[a]]) % p
        s += t
    return s % p
def rhoW(g):
    g = np.array(g, dtype=object) % p
    L = np.zeros((70, 70), dtype=np.int64)
    for i, U in enumerate(S4):
        for j, T in enumerate(S4):
            L[i, j] = det4(g[np.ix_(U, T)])
    return (np.kron(np.array(g, dtype=np.int64), L)) % p
def dir_modp(z):
    M = np.zeros((n, n), dtype=np.int64)
    for x in range(n):
        M[x, x] = (1 + I) * inv2 % p; M[x, x ^ z] = (M[x, x ^ z] + (1 - I) * inv2) % p
    return M
def mono(perm, d):
    M = np.zeros((n, n), dtype=np.int64)
    for x in range(n): M[x, perm[x]] = d[x] % p
    return M
def mono_inv(perm, d):
    M = np.zeros((n, n), dtype=np.int64)
    for x in range(n): M[perm[x], x] = pow(int(d[x]), p - 2, p)
    return M
def mm(A, B): return (A.astype(object) @ B.astype(object) % p).astype(np.int64)
def rand_mono():
    perm = rng.permutation(n); d = rng.integers(1, p, n)
    return mono(perm, d), mono_inv(perm, d)
def rank_modp(M):
    M = M.copy() % p; r = 0; rows, cols = M.shape
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        piv = r + nz[0]
        if piv != r: M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): M[nzr] = (M[nzr] - (col[nzr, None] * M[r][None, :]) % p) % p
        r += 1
    return r
def matmul_p(A, B):    # A (560x560), B (560 x k) mod p, exact via splitting into 10-bit pieces
    lo = B & 1023; hi = B >> 10
    return ((A @ lo) % p + ((A @ hi) % p) * 1024) % p
def keeps_N(R):        # R = rho(g): image of N inside N?
    return not (R[np.ix_(Oidx, Nidx)] % p).any()
def fW(R):
    return rank_modp(R[np.ix_(Oidx, Nidx)])

t0 = time.time()
# cross-check rhoW against g^(x)5 on V^(x)5 at the sorted positions (independent path)
def via_tensor(g, cols):
    out = np.zeros((DW, len(cols)), dtype=np.int64)
    for c, j in enumerate(cols):
        x, T = basis[j]
        v = np.zeros((n,) * 5, dtype=object)
        for pi in itertools.permutations(range(4)):
            sg = 1
            for a in range(4):
                for b in range(a + 1, 4):
                    if pi[a] > pi[b]: sg = -sg
            v[(x,) + tuple(T[k] for k in pi)] = sg
        go = np.array(g, dtype=object)
        for ax in range(5):
            v = np.moveaxis(np.tensordot(go, v, axes=([1], [ax])), 0, ax)
        for i, (y, U) in enumerate(basis): out[i, c] = int(v[(y,) + U]) % p
    return out
for name, g in [("dir e1", dir_modp(1)), ("random monomial", rand_mono()[0]), ("random dense", rng.integers(0, p, (n, n)))]:
    cols = list(rng.choice(DW, 3, replace=False))
    ok = np.array_equal(rhoW(g)[:, cols], via_tensor(g, cols))
    print(f"rhoW == g^(x)5 restricted to W (3 random basis columns), {name}: {ok}", flush=True)
    assert ok
# homomorphism check
a, b = rand_mono()[0], dir_modp(3)
print("rhoW(a b) == rhoW(a) rhoW(b):", np.array_equal(rhoW(mm(a, b)), mm(rhoW(a), rhoW(b))), flush=True)

# control: Hadamard on bit 0 is not monomial and must NOT keep N; dir must not keep N
H0 = np.zeros((n, n), dtype=np.int64)
for x in range(n):
    H0[x, x] = (p - 1) if (x & 1) else 1; H0[x, x ^ 1] = 1
print("control: Hadamard on bit 0 keeps N?", keeps_N(rhoW(H0)), "  dir e1 keeps N?", keeps_N(rhoW(dir_modp(1))), flush=True)
TOF = mono([x ^ (4 if (x & 3) == 3 else 0) for x in range(n)], [1] * n)
CS = mono(list(range(n)), [I if (x & 3) == 3 else 1 for x in range(n)])
CCZ = mono(list(range(n)), [p - 1 if x == 7 else 1 for x in range(n)])
named = {"Toffoli": TOF, "CS": CS, "CCZ": CCZ}
for k, g in named.items(): print(f"{k} keeps N: {keeps_N(rhoW(g))}", flush=True)
RD = [rhoW(dir_modp(z)) for z in range(1, 8)]
print("f(dir z), z=1..7:", [fW(R) for R in RD], f"[{time.time()-t0:.0f}s]", flush=True)

def run(Wl, S, steps, plain=False):
    Rn = Wl + S
    M = [[(np.eye(DW, dtype=np.int64)[:, Nidx] if r == l else np.zeros((DW, DN), dtype=np.int64)) for l in range(Wl)] for r in range(Rn)]
    def Phi():
        big = np.vstack([np.hstack([M[r][l][Oidx] for l in range(Wl)]) for r in range(Rn)])
        return rank_modp(big)
    phi = Phi(); viol = []; incs = {1: [], 2: []}; nfree = 0
    sched = []
    if plain:
        for l in range(Wl):
            for z in (1, 2, 4): sched.append(('paid1', l, z))
    else:
        for s in range(steps):
            kind = rng.choice(['free', 'free', 'mix', 'paid1', 'paid1', 'paid2', 'named'])
            sched.append((kind, int(rng.integers(Rn)), None))
    for (kind, r, z) in sched:
        if kind == 'free' or kind == 'named':
            g = rand_mono()[0] if kind == 'free' else named[rng.choice(list(named))]
            R = rhoW(g); assert keeps_N(R)
            M[r] = [matmul_p(R, X) for X in M[r]]
        elif kind == 'mix':
            C = rng.integers(-2, 3, (Rn, Rn))
            M = [[sum(int(C[r2, s]) * M[s][l] for s in range(Rn)) % p for l in range(Wl)] for r2 in range(Rn)]
        else:
            b, bi = rand_mono()
            if kind == 'paid1':
                zz = z if z is not None else int(rng.integers(1, 8))
                g = dir_modp(zz) if plain else mm(mm(b, dir_modp(zz)), bi)
            else:
                z1 = int(rng.integers(1, 8)); z2 = int(rng.integers(1, 8))
                while z2 == z1: z2 = int(rng.integers(1, 8))
                g = mm(mm(b, mm(dir_modp(z1), dir_modp(z2))), bi)
            R = rhoW(g)
            M[r] = [matmul_p(R, X) for X in M[r]]
        new = Phi()
        if kind in ('free', 'named', 'mix'):
            nfree += 1
            if new > phi: viol.append((kind, phi, new))
        else:
            rk = 1 if kind == 'paid1' else 2
            incs[rk].append(new - phi)
            if new - phi > 80 * rk: viol.append((kind, phi, new))
        phi = new
    return phi, viol, incs, nfree
phi, viol, incs, nfree = run(2, 0, 0, plain=True)
print(f"plain schedule W=2: Phi_end = {phi} (Theorem R end count needs >= 448); increments {incs[1]}; violations {viol}", flush=True)
tot = {1: [], 2: []}; allviol = []; tf = 0
for trial in range(3):
    phi, viol, incs, nfree = run(2, 2, STEPS)
    tot[1] += incs[1]; tot[2] += incs[2]; allviol += viol; tf += nfree
    print(f"trial {trial}: Phi_end {phi}; rank-1 incs {incs[1]}; rank-2 incs {incs[2]}; free/mix steps {nfree}; violations {viol}  [{time.time()-t0:.0f}s]", flush=True)
print(f"SUMMARY seed {seed}: paid rank-1 steps {len(tot[1])} (max raise {max(tot[1]) if tot[1] else None}), "
      f"rank-2 steps {len(tot[2])} (max raise {max(tot[2]) if tot[2] else None}), free/mix steps {tf}, violations {len(allviol)}")
