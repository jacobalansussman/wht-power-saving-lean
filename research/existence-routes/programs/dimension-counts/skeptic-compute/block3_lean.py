# (lean version of block3_full.py, which was killed by the memory guard)
# Block 3 "blindness" check, exact mod p (two primes), independent code.  Rep V(x)V (d = 64), inputs restricted to
# U = span{e_x (x) e_x} (8-dim).  For a shape (W, n rank-1 steps, targets), pick SPECIFIC step coefficients (random
# combinations of random words in the free letters), compute the module R = A.{inputs, step outputs} (closure under
# the free letters), and test whether every target rho(W_k) iota_U (slot t) lies in R.  If so, the shape is exactly
# solvable in this relaxation (final coefficients = a linear solve), i.e. the relaxation is blind to it.
import sys, os, json, time, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import letters, rep, mm, Echelon
p = int(sys.argv[1]); shapes = json.loads(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
rng = np.random.default_rng(seed); m = 3; n0 = 8
free, Wr = letters(m, p); F = [rep('VV', g, p) for g in free]; Wv = [rep('VV', g, p) for g in Wr]; d = 64
U = np.zeros((d, n0), np.int64)
for x in range(n0): U[x * n0 + x, x] = 1
for g in F: assert not ((g @ U) % p)[[i for i in range(d) if i not in [x * n0 + x for x in range(n0)]]].any(), 'U not B-invariant'
def closure(gens, X0s, c):
    E = Echelon(d * c, p); fr = E.add(np.stack([X.reshape(-1) for X in X0s]))
    while fr:
        nf = []
        for g in gens:                                   # one generator at a time (memory)
            nf += E.add(np.stack([mm(g, v.reshape(d, c), p).reshape(-1) for v in fr]))
        fr = nf
    return E
# random words in the free letters
words = []
for _ in range(400):
    M = np.eye(d, dtype=np.int64)
    for g in rng.integers(0, len(F), size=50): M = mm(F[g], M, p)
    words.append(M.reshape(-1))
words = np.stack(words)
def randA():
    c = rng.integers(0, p, size=(1, len(words)))
    acc = np.zeros(d * d, np.int64)
    for a in range(0, len(words), 200): acc = (acc + (c[:, a:a + 200] @ words[a:a + 200]) % p) % p
    return acc[0].reshape(d, d)
Ecg = closure(F + [Wv[1]], [U], n0); print('dim C_G iota_U =', Ecg.dim(), flush=True)
for sh in shapes:
    W, nst, tg = sh['W'], sh['n'], sh['targets']; c = n0 * W; t0 = time.time()
    src = []
    for l in range(W):
        X = np.zeros((d, c), np.int64); X[:, l * n0:(l + 1) * n0] = U; src.append(X)
    growth = [closure(F, src, c).dim()]
    for j in range(nst):
        cj = np.zeros((d, c), np.int64)
        for s in src: cj = (cj + mm(randA(), s, p)) % p
        src.append(mm(Wv[1], cj, p))
        growth.append(closure(F, src, c).dim() if (j == nst - 1 or sh.get('growth')) else None)
    E = closure(F, src, c)
    ins = []
    for t in range(W):
        T = np.zeros((d, c), np.int64); T[:, t * n0:(t + 1) * n0] = mm(Wv[tg[t]], U, p)
        ins.append(not E.reduce(T.reshape(1, -1)).any())
    print(json.dumps({'p': p, 'shape': sh['label'], 'dimR_growth': growth, 'max': W * Ecg.dim(), 'targets_in_R': ins,
                      'exactly_solvable_in_relaxation': all(ins), 'sec': round(time.time() - t0, 1)}), flush=True)
