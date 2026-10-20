"""v-w3-blueprint: numeric test of Lemma R1/R2's core identity with own code.
Shared term t (start: t = b in B; middle: t = l D(z) a), branches A delta c_h t c'_h on array h, weight psi(h).
Array-0 chi-slice: branch (A, h) -> cell u_A . q_h  (u_A = B A delta t, q_h = t^-1 c_h t), weight w_A psi(h) chi(h).
R2's core: the slice v satisfies v(u . q_g) = phi(g) v(u) for every cell u and g (phi = psi chi) -- so the kernel's
coefficient and each junk cell's coefficient [K q_g] have equal modulus.  We compute v from the ACTUAL cells
B (A delta c_h t c'_h) c'_(-h) (not from the formula) for random paths and test the identity on every touched cell."""
import sys, random, numpy as np
np.seterr(all='ignore')
sys.path.insert(0, __import__('os').path.dirname(__file__))
from groups import G, z7, z2cube
def fkey(M):
    rows = []
    for r in M:
        nz = np.nonzero(np.abs(r) > 1e-9)[0]; v = r[nz[0]]
        for z in (1, 1j, -1, -1j):
            w = v * z; a = np.angle(w)
            if -1e-9 <= a < np.pi / 2 - 1e-9: rr = r * z; break
        rows.append(tuple(np.round(rr.real * 4096).astype(int)) + tuple(np.round(rr.imag * 4096).astype(int)))
    return tuple(sorted(rows))
grp, npath, seed, kind = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
els = z7()[0] if grp == 'Z7' else z2cube()[0]
W = len(els); rng = random.Random(seed); nrng = np.random.default_rng(seed)
if grp == 'Z7':
    chars = [[np.exp(2j * np.pi * r * h / 7) for h in range(7)] for r in range(7)]; neg = lambda a: (-a) % 7
else:
    chars = [[(-1) ** bin(r & h).count('1') for h in range(8)] for r in range(8)]; neg = lambda a: a
cp = [G.mul(G.Ki, G.inv(c), G.K) for c in els]
worst = 0; ntest = 0
for trial in range(3):
    if kind == 'start': t = G.rand_engine(rng)
    else: t = G.mul(G.rand_engine(rng), G.D[rng.randrange(1, 16)], G.rand_engine(rng))
    delta = G.D[rng.randrange(1, 16)]
    psi = chars[rng.randrange(W)]
    q = [G.mul(G.inv(t), c, t) for c in els]
    reps = {}; paths = []
    for p in range(npath):
        A = np.eye(16, dtype=complex)
        for k in range(rng.randrange(2, 4)):
            A = G.mul(G.rand_engine(rng), G.D[rng.randrange(1, 16)], A)
        A = G.mul(G.rand_engine(rng), A)
        w = complex(nrng.normal(), nrng.normal())
        cells = []
        for h in range(W):
            M = G.mul(A, delta, els[h], t, cp[h], cp[neg(h)])
            k = fkey(M); reps.setdefault(k, M); cells.append(k)
        paths.append((cells, w))
    for r in range(W):
        chi = chars[r]
        v = {}
        for cells, w in paths:
            for h in range(W):
                v[cells[h]] = v.get(cells[h], 0) + w * psi[h] * chi[h]
        phi = [psi[g] * chi[g] for g in range(W)]
        for k, M in list(reps.items()):
            for g in range(W):
                k2 = fkey(G.mul(M, q[g]))
                err = abs(v.get(k2, 0) - phi[g] * v.get(k, 0))
                worst = max(worst, err); ntest += 1
    print(f'{grp} {kind} trial {trial}: {len(reps)} distinct cells touched by {npath} paths x {W} branches', flush=True)
print(f'identity v(u.q_g) = phi(g) v(u): {ntest} checks, max error {worst:.2e}')
