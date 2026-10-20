"""paid.py (copy in vex2lim_rb4, menu cliff2s added 23:21; 23:15 restart; menus cliff2 and h1small added) -- v-ex2-limits: exact Grams of words with paid letters (independent of e1_paid.py).
Words L.D<z>.R (one paid letter) or L.D<z>.M.D<z'>.R (two).  With a = (1+i)/2, b = (1-i)/2,
a^(f-k) b^k = 2^(-f/2) w^f i^(-k), so every word with k paid letters (no inverses) is 2^(-kf/2) w^(kf) times a matrix
with GAUSSIAN-INTEGER entries; columns are aggregated exactly over the paths s in F_2^f.
N*G_f(i,j) = 2^(-kf) sum_(v,y) conj(c_i(v,y)) c_j(v,y)  (same k for all words in a menu, so the factor is common).
Dedupe: |G(i,j)| = sqrt(G(i,i) G(j,j)) at the dedupe f's.
Usage: python3 paid.py <menu> <f> ...      menu in {one, two}"""
import sys, os, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, scipy.sparse as sp, myxw
SIDE = ["1", "CS", "A", "Ai", "Pr0", "Pr1", "lCS", "fCS", "B", "a0", "S1", "C01", "lX1", "fX1", "E01", "P"]
def menu_words(menu):
    if menu == 'one':
        return ["%s.D%d.%s" % (L, z, R) for L in SIDE for z in (1, 2, 3) for R in SIDE]
    if menu == 'two':   # H1 screen: CS combined with two paid letters
        LS = ["1", "CS", "A"]; MS = ["1", "CS", "A", "Pr0", "CS.A", "A.CS", "lCS", "P"]; RS = ["1", "CS", "Ai"]
        return ["%s.D%d.%s.D%d.%s" % (L, z, M, y, R) for L in LS for z in (1, 2, 3) for M in MS for y in (1, 2, 3)
                for R in RS]
    if menu == 'cliff2':   # Theorem 9 screen: CLIFFORD letters only (no CS), two paid letters
        LS = ["1", "A", "Pr0", "S1"]; MS = ["1", "A", "Pr0", "P", "CZ", "S1", "Ai", "B", "lS1", "A.Pr0", "C01", "E01"]
        RS = ["1", "Ai", "Pr1"]
        return ["%s.D%d.%s.D%d.%s" % (L, z, M, y, R) for L in LS for z in (1, 2, 3) for M in MS for y in (1, 2, 3)
                for R in RS]
    if menu == 'cliff2s':  # Theorem 9 screen to f = 5: Clifford letters only, two paid letters, reduced menu
        MS = ["1", "A", "Pr0", "P", "CZ", "S1", "Ai", "B", "lS1", "A.Pr0", "C01", "E01"]
        return ["%s.D%d.%s.D%d.1" % (L, z, M, y) for L in ("1", "A") for z in (1, 2, 3) for M in MS for y in (1, 2, 3)]
    if menu == 'h1small':  # H1 screen to f = 5: CS with two paid letters, small menu
        MS = ["CS", "A.CS", "CS.A", "lCS", "CS.Pr0", "Pr0.CS", "A", "Pr0", "1", "CS.CS"]
        return ["1.D%d.%s.D%d.1" % (z, M, y) for z in (1, 2, 3) for M in MS for y in (1, 2, 3)] + \
               ["CS.D%d.%s.D%d.1" % (z, M, y) for z in (1, 2, 3) for M in ("CS", "A", "1") for y in (1, 2, 3)]
    raise SystemExit('menu?')
def popcount(x):
    c = np.zeros_like(x)
    while x.any(): c += x & 1; x = x >> 1
    return c
def spread(s, z, f):            # s in F_2^f -> address with z in every column j where s_j = 1
    out = np.zeros_like(s)
    for j in range(f): out |= ((s >> j) & 1) * (z << (2 * j))
    return out
def column_data(w, f):
    """returns (keys, vals): keys = v*N + target, vals = Gaussian-integer amplitude (complex), aggregated."""
    N = 4 ** f; parts = w.split('.'); segs = [[]]; zs = []
    for p in parts:
        if p.startswith('D') and p[1:].isdigit(): zs.append(int(p[1:])); segs.append([])
        else: segs[-1].append(p)
    mono = [myxw.word('.'.join(s) if s else '1', f) for s in segs]
    S = np.arange(2 ** f, dtype=np.int64); sh = {z: spread(S, z, f) for z in set(zs)}; ws = popcount(S)
    v = np.arange(N, dtype=np.int64)
    t, ph = mono[0][0][v], mono[0][1][v]                     # after the first free segment
    vv = v
    for k, z in enumerate(zs):
        t = (t[:, None] ^ sh[z][None, :]).ravel()
        ph = ((ph[:, None] - ws[None, :]) % 4).ravel()
        vv = np.repeat(vv, 2 ** f)
        tn, pn = mono[k + 1]
        ph = (ph + pn[t]) % 4; t = tn[t]
    key = vv * N + t
    u, inv = np.unique(key, return_inverse=True)
    val = np.zeros(u.size, complex)
    np.add.at(val, inv, np.array([1, 1j, -1, -1j])[ph])
    nz = val != 0
    return u[nz], val[nz], len(zs)
def gram(words, f):
    cols = []; ks = []
    for w in words:
        k, v, npaid = column_data(w, f); cols.append((k, v)); ks.append(npaid)
    assert len(set(ks)) == 1
    allk = np.unique(np.concatenate([c[0] for c in cols]))
    rows = np.concatenate([np.full(c[0].size, i) for i, c in enumerate(cols)])
    cidx = np.concatenate([np.searchsorted(allk, c[0]) for c in cols])
    data = np.concatenate([c[1] for c in cols])
    M = sp.csr_matrix((data, (rows, cidx)), shape=(len(words), allk.size))
    GN = (M.conj() @ M.T).toarray() / 2.0 ** (ks[0] * f)
    return GN / 4.0 ** f
def control(words, f):
    """compare my column data with a direct dense computation of R D M D L at f (small f only)."""
    N = 4 ** f; a, b = (1 + 1j) / 2, (1 - 1j) / 2
    def dmat(z):
        d = np.zeros((4, 4), complex)
        for x in range(4): d[x, x] += a; d[x ^ z, x] += b
        out = np.ones((1, 1))
        for j in range(f): out = np.kron(d, out)
        return out
    def mmat(seg):
        t, p = myxw.word(seg if seg else '1', f); M = np.zeros((N, N), complex)
        M[t, np.arange(N)] = np.array([1, 1j, -1, -1j])[p]; return M
    err = 0
    for w in words:
        parts = w.split('.'); segs = [[]]; zs = []
        for p in parts:
            if p.startswith('D') and p[1:].isdigit(): zs.append(int(p[1:])); segs.append([])
            else: segs[-1].append(p)
        W = mmat('.'.join(segs[0]))
        for k, z in enumerate(zs): W = mmat('.'.join(segs[k + 1])) @ dmat(z) @ W
        k, v, npaid = column_data(w, f)
        fac = 2.0 ** (-npaid * f / 2) * np.exp(1j * np.pi / 4) ** (npaid * f)
        D = np.zeros((N, N), complex); D[k % N, k // N] = v * fac
        err = max(err, np.abs(D - W).max())
    return err
if __name__ == '__main__':
    menu = sys.argv[1]; FS = [int(a) for a in sys.argv[2:]]; t0 = time.time()
    words = menu_words(menu)
    e = control(words[::37], 2); e3 = control(words[::97], 3)
    print('control (dense R D M D L vs my exact columns): max err %.2e at f=2, %.2e at f=3' % (e, e3), flush=True)
    DF = (3, 4)
    same = None
    for f in DF:
        G = gram(words, f); d = np.sqrt(np.real(np.diag(G)))
        s = np.abs(np.abs(G) - np.outer(d, d)) < 1e-10
        same = s if same is None else (same & s)
    keep = []
    for i in range(len(words)):
        if not any(same[i, k] for k in keep): keep.append(i)
    kw = [words[i] for i in keep]
    open('my_words_paid_%s.txt' % menu, 'w').write('\n'.join(kw) + '\n')
    print('menu %s: %d words -> %d distinct up to a global phase (f = 3, 4); %.1f s' % (menu, len(words), len(kw), time.time() - t0), flush=True)
    for f in FS:
        t1 = time.time(); G = gram(kw, f)
        np.save('my_G_paid_%s_f%d.npy' % (menu, f), G)
        ev = np.linalg.eigvalsh((G + G.conj().T) / 2); nz = ev[ev > 1e-9]
        print('f=%d: rank %d of %d; bottom nonzero %s; %.1f s' % (f, nz.size, len(kw), ' '.join('%.4f' % x for x in nz[:5]), time.time() - t1), flush=True)
