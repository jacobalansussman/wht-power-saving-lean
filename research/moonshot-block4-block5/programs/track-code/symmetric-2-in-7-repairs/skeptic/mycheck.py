"""mycheck.py (v-mw2-absorb4, skeptic): independent closure + least-squares test of a dumped program (dump.py output).
Own code: exact free-group membership for the engine game E (Z4 phases i^q with q a Clifford quadratic form, affine
address maps, shifts, global i^k), and, for comparison, the looser 'any Z4 monomial' test the claimant's library uses.
Formal sums: a cell is an exact matrix (no identification of g with i g).  Columns are (generator, free letter on the
left); cells whose matrix is free are dropped (cancelled by the input generator of that array).
Coefficients are recomputed from the move maps for any alpha.
Usage: mycheck.py dump.pkl mode [starts]   mode: exact | z4"""
import sys, pickle, itertools
import numpy as np

SC = 4096.0
def key(M):
    return (np.round(M.real * SC).astype(np.int64).tobytes() + np.round(M.imag * SC).astype(np.int64).tobytes())

Z4 = [1, 1j, -1, -1j]
def mono_parts(M, tol=1e-9):
    nz = np.abs(M) > tol
    if not ((nz.sum(axis=0) == 1).all() and (nz.sum(axis=1) == 1).all()):
        return None
    n = M.shape[0]
    sig = [0] * n; q = [0] * n
    for y in range(n):
        x = int(np.argmax(np.abs(M[:, y])))
        v = M[x, y]
        ks = [k for k in range(4) if abs(v - Z4[k]) < 1e-9]
        if not ks:
            return None
        sig[y] = x; q[y] = ks[0]
    return sig, q

def is_affine(sig, m):
    t = sig[0]
    P = [sig[1 << j] ^ t for j in range(m)]
    for x in range(1 << m):
        y = t
        for j in range(m):
            if (x >> j) & 1:
                y ^= P[j]
        if sig[x] != y:
            return False
    return True

def is_clifford_phase(q, m):
    """q: Z4 exponent as function of the (input) address y.  In Q iff q(y) = a0 + sum a_j y_j + 2 sum b_jk y_j y_k mod 4."""
    a0 = q[0]
    a = [(q[1 << j] - a0) % 4 for j in range(m)]
    b = {}
    for j in range(m):
        for k in range(j + 1, m):
            r = (q[(1 << j) | (1 << k)] - a0 - a[j] - a[k]) % 4
            if r % 2:
                return False
            b[(j, k)] = r // 2
    for y in range(1 << m):
        bits = [(y >> j) & 1 for j in range(m)]
        v = a0 + sum(a[j] * bits[j] for j in range(m)) + 2 * sum(b[(j, k)] * bits[j] * bits[k] for (j, k) in b)
        if (v - q[y]) % 4:
            return False
    return True

def make_free_test(mode, m):
    cache = {}
    def f(M):
        k = key(M)
        v = cache.get(k)
        if v is None:
            mp = mono_parts(M)
            if mp is None:
                v = False
            elif mode == 'z4':
                v = True
            else:
                sig, q = mp
                v = is_affine(sig, m) and is_clifford_phase(q, m)
            cache[k] = v
        return v
    return f

def coefs(dmp, alpha):
    out = []
    for i, cf in enumerate(dmp['coef_of']):
        kind = cf[0]
        if kind == 'inp':
            out.append(np.ones(1, dtype=complex))
        elif kind == 'move':
            mv = dmp['moves'][cf[1]]
            v = np.zeros(len(dmp['gens'][i]), dtype=complex)
            for (a, src, si, ri, ch) in mv['maps']:
                np.add.at(v, np.array(ri, dtype=int), ch * alpha[a] * out[src][np.array(si, dtype=int)])
            out.append(v)
        else:
            out.append(out[dmp['moves'][cf[1]]['rep']])
    return out

def closure(dmp, l0, isfree, cap=200000):
    gens = dmp['gens']; K = dmp['K']
    rows = {}; rowmat = []; cols = {}; ent = []
    def addrow(l, M):
        k = (l, key(M))
        r = rows.get(k)
        if r is None:
            r = len(rows); rows[k] = r; rowmat.append((l, M)); queue.append(r)
        return r
    queue = []
    addrow(l0, K)
    # index generator terms per array
    terms = [(i, j, l, h, np.conj(h).T) for i, T in enumerate(gens) for j, (l, h) in enumerate(T)]
    while queue:
        r = queue.pop()
        l, C = rowmat[r]
        for (i, j, lt, h, hinv) in terms:
            if lt != l:
                continue
            b = C @ hinv
            if not isfree(b):
                continue
            ck = (i, key(b))
            if ck in cols:
                continue
            c = len(cols); cols[ck] = c
            if c > cap:
                return None
            for jj, (l2, h2) in enumerate(gens[i]):
                M = b @ h2
                if isfree(M):
                    continue
                r2 = addrow(l2, M)
                ent.append((r2, c, i, jj))
    return rows, cols, ent

def residual(nr, nc, ent, cv):
    A = np.zeros((nr, nc), dtype=complex)
    for (r, c, i, jj) in ent:
        A[r, c] += cv[i][jj]
    t = np.zeros(nr, dtype=complex); t[0] = 1
    b = np.linalg.lstsq(A, t, rcond=None)[0]
    res = A @ b - t
    return float(np.sum(np.abs(res) ** 2)), A

if __name__ == '__main__':
    dmp = pickle.load(open(sys.argv[1], 'rb'))
    mode = sys.argv[2]; starts = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    m, W = dmp['m'], dmp['W']
    isfree = make_free_test(mode, m)
    # legitimacy of the content letters, the symmetry elements: report which are exact-E free
    fe = make_free_test('exact', m)
    for t, mv in enumerate(dmp['moves']):
        pass
    print("Hc (c, c') exact-free:", [(fe(c), fe(cp)) for (c, pi, cp) in dmp['Hc']])
    cls = []
    for l in range(W):
        cl = closure(dmp, l, isfree)
        rows, cols, ent = cl
        print(f"array {l}: closure {len(rows)} x {len(cols)} (mode {mode})")
        cls.append(cl)
    free_idx = dmp['free_idx']; N = len(dmp['alpha'])
    def total(alpha):
        cv = coefs(dmp, alpha)
        return sum(residual(len(r), len(c), e, cv)[0] for (r, c, e) in cls)
    al = np.array(dmp['alpha'], dtype=complex)
    print("residual^2 at the dumped alpha (sum over arrays):", total(al))
    cvd = coefs(dmp, al)
    print("coef check vs dump:", max(np.max(np.abs(a - b)) for a, b in zip(cvd, dmp['coef_alpha'])))
    from scipy.optimize import minimize
    rng = np.random.default_rng(11)
    nf = len(free_idx)
    def f(x):
        v = al.copy(); v[free_idx] = x[:nf] + 1j * x[nf:]
        return total(v)
    best = []
    for st in range(starts):
        x0 = rng.normal(size=2 * nf) * (1 if st % 2 else 2)
        r = minimize(f, x0, method='BFGS', options={'maxiter': 400})
        v = r.x[:nf] + 1j * r.x[nf:]
        best.append((r.fun, v))
        print(f"start {st}: total res^2 {r.fun:.8f}  alpha_free {np.round(v, 5)}", flush=True)
    best.sort(key=lambda t: t[0])
    print("BEST", best[0][0], np.round(best[0][1], 6))
    pickle.dump({'best': best[:5], 'closure_dims': [(len(r), len(c)) for (r, c, e) in cls]},
                open(sys.argv[1].replace('.pkl', f'_my_{mode}.pkl'), 'wb'))
