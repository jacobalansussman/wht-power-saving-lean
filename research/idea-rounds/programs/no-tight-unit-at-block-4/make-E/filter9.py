# tb4-tight: the f = 9 TIGHTNESS FILTER for formal words (wave 16, block 4).
# In a tight 2-in-7 (blocks of rank 2,2,2,1) every content (beta_1, beta_2) that ever appears satisfies
#   rho_9(beta_l) E  inside N   (E = N cap K^-1 N, dim 16368)          [N-side, Theorem R with Y = N]
# This file applies rho_9 = g (x) Lambda^8 g to random vectors of E, letter by letter, on arrays Y[x, S]
# (16 x 12870, S = 8-subset as bitmask index), exact mod p.  A word is a list of letters applied RIGHT TO LEFT
# (word[-1] acts first): ('dir', z) = dir z = alpha I + beta S_z; ('dirinv', z); ('perm', sigma) = e_x -> e_sigma(x);
# ('phase', q) = diag(i^q[x]).  A random e in E is built from E_0 (chi = 0 block, blk9.py 'E' output) and
# E_chi = Z_c E_0 (Z_c keeps E: K Z_c K^-1 is a Pauli times a scalar) [proved, NOTE].
# One random e with rho(word) e not in N disproves the condition for that word (certainly); passing for two
# independent random e's means the condition holds except with probability <= 2/p.
import sys, os, time, itertools
import numpy as np
import scipy.sparse as sp
n = 16; k = 8
p = int(os.environ.get('PRIME', '33554393'))
I = next(pow(g, (p - 1) // 4, p) for g in range(2, 100) if pow(g, (p - 1) // 2, p) == p - 1)
inv2 = (p + 1) // 2; alpha = (1 + I) * inv2 % p; beta = (1 - I) * inv2 % p
SK = sorted([sum(1 << t for t in c) for c in itertools.combinations(range(n), k)])
IDX = {b: i for i, b in enumerate(SK)}
NS = len(SK)
ELS = [[t for t in range(n) if S >> t & 1] for S in SK]
# N mask: entry (x, S) with x in S
NMASK = np.array([[(S >> x) & 1 for S in SK] for x in range(n)], dtype=bool)

def parity(seq):
    return sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq)) if seq[i] > seq[j]) & 1

_lam_cache = {}
def lam_dir(z, inv=False):
    """Lambda^8 of dir z (or its inverse) as a sparse matrix, column T -> sum_S det D[S,T] e_S"""
    key = (z, inv)
    if key in _lam_cache: return _lam_cache[key]
    a, b = (beta, alpha) if inv else (alpha, beta)
    det2 = (a * a - b * b) % p
    rows, cols, vals = [], [], []
    for ci, T in enumerate(SK):
        els = ELS[ci]
        single = [t for t in els if not (T >> (t ^ z)) & 1]
        nfull = (len(els) - len(single)) // 2
        base = pow(det2, nfull, p)
        for choice in itertools.product((0, 1), repeat=len(single)):
            phi = {}
            coef = base
            for t, c in zip(single, choice):
                if c: phi[t] = t ^ z; coef = coef * b % p
                else: phi[t] = t; coef = coef * a % p
            img = [phi.get(t, t) for t in els]
            Smask = sum(1 << s for s in img)
            # rows of D[S,T] taken in phi-order; sign = parity of sorting img
            sgn = -1 if parity(img) else 1
            rows.append(IDX[Smask]); cols.append(ci); vals.append((sgn * coef) % p)
    M = sp.csr_matrix((np.array(vals, dtype=np.int64), (np.array(rows), np.array(cols))), shape=(NS, NS))
    _lam_cache[key] = M
    return M

def dirmat(z, inv=False):
    a, b = (beta, alpha) if inv else (alpha, beta)
    A = np.zeros((n, n), dtype=np.int64)
    for x in range(n):
        A[x, x] = (A[x, x] + a) % p; A[x ^ z, x] = (A[x ^ z, x] + b) % p
    return A

def apply_dir(Y, z, inv=False):
    A = dirmat(z, inv)
    Y1 = (A @ Y) % p                                  # first factor
    L = lam_dir(z, inv)
    return np.asarray((L @ Y1.T) % p).T % p           # (Lambda^8 A) on the second factor

def perm_tables(sigma):
    tgt = np.zeros(NS, dtype=np.int64); sg = np.zeros(NS, dtype=np.int64)
    for i, S in enumerate(SK):
        img = [sigma[t] for t in ELS[i]]
        tgt[i] = IDX[sum(1 << s for s in img)]; sg[i] = -1 if parity(img) else 1
    return tgt, sg

def apply_perm(Y, sigma):
    tgt, sg = perm_tables(sigma)
    out = np.zeros_like(Y)
    for x in range(n):
        out[sigma[x], tgt] = (sg * Y[x]) % p
    return out

def apply_phase(Y, q):
    ip = [pow(I, j, p) for j in range(4)]
    # weight of (x, S): q[x] + sum_{s in S} q[s]
    qs = np.array([sum(q[t] for t in ELS[i]) for i in range(NS)], dtype=np.int64)
    out = np.zeros_like(Y)
    for x in range(n):
        ph = np.array([ip[(q[x] + v) % 4] for v in qs], dtype=np.int64)
        out[x] = (Y[x] * ph) % p
    return out

def apply_word(word, Y):
    for let in reversed(word):
        kind = let[0]
        if kind == 'dir': Y = apply_dir(Y, let[1])
        elif kind == 'dirinv': Y = apply_dir(Y, let[1], inv=True)
        elif kind == 'perm': Y = apply_perm(Y, let[1])
        elif kind == 'phase': Y = apply_phase(Y, let[1])
        else: raise ValueError(kind)
    return Y

def out_part(Y):
    return Y[~NMASK]

def in_N(Y): return not out_part(Y).any()

# ---- random vectors of E from the chi = 0 block basis ----
def load_E0(r=1):
    here = os.environ.get('E_DIR', os.path.dirname(os.path.abspath(__file__)))
    return np.load(os.path.join(here, f"E_r{r}_p{p}.npy"))

def expand_chi0(v):
    """chi = 0 block vector over N-reps (T containing 0, in SK order restricted) -> full array"""
    Ncols = [i for i, S in enumerate(SK) if S & 1]
    Y = np.zeros((n, NS), dtype=np.int64)
    for j, ci in enumerate(Ncols):
        if v[j] == 0: continue
        T = SK[ci]
        for d in range(n):
            img = [t ^ d for t in ELS[ci]]
            s = -1 if parity(img) else 1
            Y[d, IDX[sum(1 << u for u in img)]] = (Y[d, IDX[sum(1 << u for u in img)]] + s * v[j]) % p
    return Y

_exp_cache = {}
def random_E(E0, rng, nchi=16):
    Y = np.zeros((n, NS), dtype=np.int64)
    for c in range(nchi):
        coeff = rng.integers(0, p, size=len(E0))
        v = (coeff @ (E0 % p)) % p if False else np.zeros(E0.shape[1], dtype=np.int64)
        # coeff @ E0 can overflow int64 (1023 * p^2 > 2^63): accumulate in chunks
        for s in range(0, len(E0), 128):
            v = (v + (coeff[s:s + 128] @ E0[s:s + 128]) % p) % p
        Yc = expand_chi0(v)
        Zq = [2 * (bin(c & x).count('1') & 1) for x in range(n)]     # Z_c = diag((-1)^{c.x}) = phase with q = 2 c.x
        Y = (Y + apply_phase(Yc, Zq)) % p
    return Y

K_WORD = [('dir', 1), ('dir', 2), ('dir', 4), ('dir', 8)]

if __name__ == '__main__':
    t0 = time.time()
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    # self-test 1: Lambda^8 dir against brute-force minors on a few columns
    D = dirmat(1)
    L = lam_dir(1).tocsc()
    for ci in rng.choice(NS, 5, replace=False):
        T = ELS[ci]
        col = np.zeros(NS, dtype=np.int64)
        for ri, S in enumerate(SK):
            sub = D[np.ix_(ELS[ri], T)] % p
            # det mod p by elimination
            Mx = [list(map(int, r_)) for r_ in sub]; det = 1
            for c_ in range(k):
                pv = next((r_ for r_ in range(c_, k) if Mx[r_][c_] % p), None)
                if pv is None: det = 0; break
                if pv != c_: Mx[c_], Mx[pv] = Mx[pv], Mx[c_]; det = -det
                det = det * Mx[c_][c_] % p; iv = pow(Mx[c_][c_], p - 2, p)
                for r_ in range(c_ + 1, k):
                    fct = Mx[r_][c_] * iv % p
                    if fct: Mx[r_] = [(Mx[r_][j] - fct * Mx[c_][j]) % p for j in range(k)]
            col[ri] = det % p
        assert np.array_equal(col % p, np.asarray(L[:, ci].todense()).ravel() % p), "Lambda^8 dir mismatch"
    print(f"self-test Lambda^8 dir: 5 columns match brute-force minors  ({time.time()-t0:.0f}s)", flush=True)
    E0 = load_E0(1)
    e = random_E(E0, rng)
    print(f"random e in E built; e in N: {in_N(e)}; nonzeros {np.count_nonzero(e)}  ({time.time()-t0:.0f}s)", flush=True)
    Ke = apply_word(K_WORD, e)
    print(f"rho(K) e in N: {in_N(Ke)}   (checks E_0 and the Z_c claim)  ({time.time()-t0:.0f}s)", flush=True)
    d1e = apply_word([('dir', 1)], e)
    print(f"rho(d1) e in N: {in_N(d1e)}; out-part nonzeros {np.count_nonzero(out_part(d1e))}  ({time.time()-t0:.0f}s)", flush=True)

# ---- the filter for a beam search ----
# beta = list of (lam, word): lam an integer mod p (use F.I for i), word as above.  Returns True when
# rho_9(beta) e is in N for every test vector e (necessary for any content of a tight 2-in-7, ranks 2,2,2,1).
def e_condition(beta, es):
    for e in es:
        acc = np.zeros_like(e)
        for lam, word in beta:
            acc = (acc + (lam % p) * apply_word(word, e)) % p
        if out_part(acc).any():
            return False
    return True

# P-side: E_P (+) L = null of the N x out block of K at chi = 0 (blk9.py 'EP'); a content must map E_P into P (+) L,
# i.e. rho_9(beta) u has zero N-part for u in E_P (L is G-invariant with zero N-part, so it never matters).
def load_EP0(r=1):
    here = os.environ.get('E_DIR', os.path.dirname(os.path.abspath(__file__)))
    return np.load(os.path.join(here, f"EP_r{r}_p{p}.npy"))

def expand_chi0_out(v):
    Ocols = [i for i, S in enumerate(SK) if not S & 1]
    Y = np.zeros((n, NS), dtype=np.int64)
    for j, ci in enumerate(Ocols):
        if v[j] == 0: continue
        for d in range(n):
            img = [t ^ d for t in ELS[ci]]
            s = -1 if parity(img) else 1
            ii = IDX[sum(1 << u for u in img)]
            Y[d, ii] = (Y[d, ii] + s * v[j]) % p
    return Y

def random_EP(EP0, rng, nchi=16):
    Y = np.zeros((n, NS), dtype=np.int64)
    for c in range(nchi):
        coeff = rng.integers(0, p, size=len(EP0))
        v = np.zeros(EP0.shape[1], dtype=np.int64)
        for s in range(0, len(EP0), 128):
            v = (v + (coeff[s:s + 128] @ EP0[s:s + 128]) % p) % p
        Zq = [2 * (bin(c & x).count('1') & 1) for x in range(n)]
        Y = (Y + apply_phase(expand_chi0_out(v), Zq)) % p
    return Y

def p_condition(beta, us):
    for u in us:
        acc = np.zeros_like(u)
        for lam, word in beta:
            acc = (acc + (lam % p) * apply_word(word, u)) % p
        if acc[NMASK].any():
            return False
    return True
