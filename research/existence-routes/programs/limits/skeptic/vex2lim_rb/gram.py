"""gram.py -- v-ex2-limits independent Gram computation for free (monomial) words.
Exact: N*G_f(i,j) = sum_v [t_i(v) = t_j(v)] i^(p_j(v) - p_i(v)), accumulated as INTEGER counts per Z4 difference.
Two algorithms: 'dense' (chunked pairwise comparison) and 'sort' (group (v, target) pairs by sorting).
Projective dedupe: words i, j merged if |G_f(i,j)| = 1 (= G_f(i,i)) at f = 3, 4, 5 (equal up to an f-dependent phase).
Usage: python3 gram.py <menu> <alg> <f> ...   menu in {main, chain3}"""
import sys, os, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, myxw
FREE = ["1", "P", "C01", "C10", "X1", "X2", "X3", "Z1", "Z2", "Z3", "S1", "S2", "S3", "CZ", "CS", "CSi",
        "A", "Ai", "a0", "a1", "a0i", "a1i", "E01", "E10", "B", "Bi", "Pr0", "Pr1",
        "lX1", "lX2", "lZ1", "lZ2", "lS1", "lS2", "lCZ", "lCS", "fX1", "fX2", "fZ1", "fZ2", "fS1", "fS2", "fCZ", "fCS",
        "Ai.lS1.A", "Ai.lS2.A", "Ai.lCS.A", "A.CS.Ai", "Ai.CS.A", "B.CS.Bi", "A.B", "A.CS.B", "B.A", "Ai.fX1.A",
        "lS1.lS2", "Pr0.Pr1", "fS1.lS1"]
GEN = ["P", "C01", "X1", "Z1", "Z2", "S1", "S3", "CS", "CZ", "A", "Ai", "a0", "a1i", "E01", "B", "Bi",
       "Pr0", "Pr1", "lX1", "lZ1", "lS1", "lCS", "lC01", "fX1", "fZ2", "fCS"]
G3 = ["A", "Ai", "a0", "E01", "B", "Pr0", "Pr1", "lX1", "lS1", "lCS", "fX1", "fCS", "CS"]
def menu_words(menu):
    if menu == 'main':
        return list(FREE) + list(GEN) + ["%s.%s" % (a, b) for a in GEN for b in GEN]
    if menu == 'chain3':
        return ["1"] + G3 + ["%s.%s" % (a, b) for a in G3 for b in G3] + \
               ["%s.%s.%s" % (a, b, c) for a in G3 for b in G3 for c in G3]
    raise SystemExit('menu?')
PH = np.array([1, 1j, -1, -1j])
def gram_dense(T, P, chunk=8):
    n, N = T.shape; cnt = np.zeros((4, n, n), np.int64)
    for s in range(0, N, chunk):
        t = T[:, s:s + chunk]; p = P[:, s:s + chunk]
        eq = t[:, None, :] == t[None, :, :]
        d = (p[None, :, :] - p[:, None, :]) % 4
        for k in range(4):
            cnt[k] += (eq & (d == k)).sum(axis=2)
    return cnt
def gram_sort(T, P):
    import scipy.sparse as sp
    n, N = T.shape
    key = (np.arange(N, dtype=np.int64)[None, :] * N + T).ravel()
    u, inv = np.unique(key, return_inverse=True)
    rows = np.repeat(np.arange(n), N)
    cnt = np.zeros((4, n, n), np.int64)
    # integer one-hot per phase value, then counts of pairs per phase difference
    Ms = [sp.csr_matrix(((P.ravel() == k).astype(np.int64), (rows, inv.ravel())), shape=(n, u.size)) for k in range(4)]
    for a in range(4):
        for b in range(4):
            cnt[(b - a) % 4] += (Ms[a] @ Ms[b].T).toarray()
    return cnt
def to_complex(cnt):
    return cnt[0] - cnt[2] + 1j * (cnt[1] - cnt[3])
if __name__ == '__main__':
    menu, alg = sys.argv[1], sys.argv[2]; FS = [int(a) for a in sys.argv[3:]]
    words = menu_words(menu); t0 = time.time()
    # dedupe at f = 3, 4, 5 by |G| = 1 (computed with the sort algorithm on all words)
    same = None
    for f in (3, 4, 5):
        T = np.empty((len(words), 4 ** f), np.int64); P = np.empty_like(T)
        for i, w in enumerate(words): T[i], P[i] = myxw.word(w, f)
        Gc = to_complex(gram_sort(T, P)) / 4 ** f
        s = np.abs(np.abs(Gc) - 1) < 1e-12
        same = s if same is None else (same & s)
    keep = []; rep = {}
    for i, w in enumerate(words):
        j = next((k for k in keep if same[i, k]), None)
        if j is None: keep.append(i)
    kw = [words[i] for i in keep]
    open('my_words_%s.txt' % menu, 'w').write('\n'.join(kw) + '\n')
    print('menu %s: %d words -> %d projective families (my dedupe, |G|=1 at f=3,4,5); %.1f s' % (menu, len(words), len(kw), time.time() - t0), flush=True)
    for f in FS:
        t1 = time.time(); N = 4 ** f
        T = np.empty((len(kw), N), np.int64); P = np.empty_like(T)
        for i, w in enumerate(kw): T[i], P[i] = myxw.word(w, f)
        cnt = gram_dense(T, P) if alg == 'dense' else gram_sort(T, P)
        np.save('my_cnt_%s_f%d.npy' % (menu, f), cnt)
        Gf = to_complex(cnt) / N
        ev = np.linalg.eigvalsh(Gf)
        nz = ev[ev > 1e-9]
        print('f=%d: rank %d of %d; smallest nonzero %.5f; %.1f s' % (f, nz.size, len(kw), nz.min(), time.time() - t1), flush=True)
