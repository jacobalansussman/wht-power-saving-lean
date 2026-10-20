"""t_formal.py: the strongest check available at block 2 -- the lifted engine words are FORMAL identities in the group
algebra C[G] (RULES section 0: right at every f = formal identity), not only right at f = 1..8.
G (block 2) is enumerated exactly (46,080 matrices, entries in Z[i]/2); each register holds, for each input, a vector of
coefficients over G (modulo a prime p = 1 mod 4, two primes); a letter acts by the left-multiplication permutation of
G. Output l must be exactly the basis vector [K] for input l and 0 for the other inputs.
Negative control: fp-build's f = 1 lift on the same programs must fail.
usage: t_formal.py"""
import sys, os, time, random, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lift import *
import t_units_cases as TC

m = 2; E = Eng(m); S = E.S; n = E.n
def key(M): return np.round(2 * np.asarray(M)).astype(np.int64).tobytes() if False else \
    (np.round(2 * M.real).astype(np.int8).tobytes() + np.round(2 * M.imag).astype(np.int8).tobytes())
t0 = time.time()
gens = [E.dirm(z) for z in range(1, n)] + [E.Xm[a] for a in range(1, n)] + \
       [np.diag([1j ** ((pc(z & x) + c) & 1) for x in range(n)]) for z in range(n) for c in (0, 1)]
for cols in ([2, 1], [1, 3], [3, 2]):
    gens.append(E.lift_P(S.addr(cols)))
mats = [np.eye(n, dtype=complex)]; idx = {key(mats[0]): 0}; fr = [0]
while fr:
    nf = []
    for i in fr:
        for g in gens:
            M = g @ mats[i]; k = key(M)
            assert np.abs(M - (np.round(2 * M.real) + 1j * np.round(2 * M.imag)) / 2).max() < 1e-9
            if k not in idx: idx[k] = len(mats); mats.append(M); nf.append(idx[k])
    fr = nf
NG = len(mats); GM = np.array(mats)
print(f"block 2: |G| = {NG} (expected 46080 = 4 x 11520), enumerated in {time.time() - t0:.1f} s", flush=True)
assert NG == 46080
def vkeys(A):
    B = np.concatenate([np.round(2 * A.real).astype(np.int8).reshape(len(A), -1),
                        np.round(2 * A.imag).astype(np.int8).reshape(len(A), -1)], axis=1)
    return np.ascontiguousarray(B).view(np.dtype((np.void, B.shape[1]))).ravel()
GK = vkeys(GM); order = np.argsort(GK); GKs = GK[order]
assert np.all(GKs[1:] != GKs[:-1])
permc = {}
def perm(M):
    k = key(M)
    if k not in permc:
        P = np.einsum('ij,njk->nik', M, GM)
        pos = np.searchsorted(GKs, vkeys(P)); assert np.all(GKs[pos] == vkeys(P))
        permc[k] = order[pos]
    return permc[k]
def run_formal(word, p):
    F = ModP(p)
    regs = []
    for st in word['steps']:
        if st[0] == 'in':
            v = np.zeros((word['W'], NG), dtype=np.int64); v[st[1], 0] = 1; regs.append(v)
        elif st[0] == 'paid':
            pm = perm(E.dirm(st[1])); src = regs[st[2]]; out = np.zeros_like(src); out[:, pm] = src; regs.append(out)
        else:
            out = np.zeros((word['W'], NG), dtype=np.int64)
            for (cf, x, c, h, u, j, s_) in st[1]:
                M = word['heads'][h] @ ((1j ** j) * E.P[u])
                cc = F.fr(x) * F.cx(complex(c)) % p
                pm = perm(M); src = regs[s_]
                tmp = np.zeros_like(src); tmp[:, pm] = src
                out = (out + cc * tmp) % p
            regs.append(out)
    kK = idx[key(E.K)]
    ok = True
    for l, o in enumerate(word['outs']):
        tgt = np.zeros((word['W'], NG), dtype=np.int64); tgt[l, kK] = 1
        ok &= np.array_equal(regs[o] % p, tgt)
    return ok
for name, P in TC.cases(E, random.Random(7 * m + 1)):
    if P.solve() != 'yes': print(name, "frame solve failed"); continue
    W1 = lift(E, P, 'full'); W0 = lift(E, P, 'f1')
    r1 = [run_formal(W1, p) for p in (167772161, 469762049)]
    r0 = [run_formal(W0, p) for p in (167772161, 469762049)]
    print(f"block 2 {name}: W = {P.W}, M = {P.moves()}: lifted word equals [K] e_l in C[G] (mod both primes): {r1};"
          f" f = 1 lift alone: {r0}  ({time.time() - t0:.0f} s)", flush=True)
