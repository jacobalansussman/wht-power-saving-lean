"""d4gen.py -- e3-bpexact: the exact distance-4 test of d4test.py for an irrep sigma whose row block Wsp is NOT full.
Wsp = {E_d1 R_sigma(A)} = sum_d' Alg . t_d' E_d',  Alg = sigma(C[B_d1]) (B_d1 = letters fixing the coset of D3),
t_d' = sigma(beta(d1, b')) for one letter b' with b'.d1 = d'.  Unknowns: the coordinates of w in this basis; the
conditions w M_a Pi_(S - d0) = 0 (a in req) are linear; null space by a Gram matrix; then the kernel-block ranks of a
generic null vector.  Representations: 'gl15' (linear part on the 15 nonzero addresses), 'gl14' (its 14-dim
irreducible part), 'aff16' (affine permutation on the 16 addresses), 'sq16' (b o b entrywise: the D16 block of V(x)V).
argv: seed psi nBm sigma_name"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4test import *

def make_sig(name):
    if name == 'gl15': return sig_gl15, 15
    if name == 'aff16': return sig_aff16, 16
    if name == 'sq16': return (lambda b: (b * b)), 16
    if name == 'gl14':
        v = np.ones(15) / np.sqrt(15)
        U, _, _ = np.linalg.svd(np.eye(15) - np.outer(v, v)); Q14 = U[:, :14].T
        return (lambda b: Q14 @ sig_gl15(b) @ Q14.T), 14
    raise ValueError(name)

def alg_basis(mats, rng, tol=1e-9, rounds=6):
    ds = mats[0].shape[0]
    def orth(vs):
        A = np.array([m.reshape(-1) for m in vs])
        U, s, Vh = np.linalg.svd(A, full_matrices=False); r = int((s > tol * s[0]).sum())
        return Vh[:r]
    B = orth(mats)
    for _ in range(rounds):
        prods = [mats[rng.randrange(len(mats))] @ mats[rng.randrange(len(mats))] for _ in range(4 * len(mats))]
        B2 = orth([b.reshape(ds, ds) for b in B] + prods)
        if B2.shape[0] == B.shape[0]: break
        B = B2
    return [b.reshape(ds, ds) for b in B]

def gen_test(D, sig, ds, ys, ks, req, letters, rng, nrng, tol=1e-9):
    L = D.L
    Bd1 = [b for b in letters if D.act(b)[D.d1][0] == D.d1]
    Alg = alg_basis([sig(D.act(b)[D.d1][1]) for b in Bd1], rng)
    t = [None] * 30
    for b in letters:
        dd, beta = D.act(b)[D.d1]
        if t[dd] is None: t[dd] = sig(beta)
    assert all(x is not None for x in t)
    Rk = sum(y * D.R(k, sig, ds) for y, k in zip(ys, ks))
    Rp = [D.R(p, sig, ds) for p in L.p]
    Ma = {a: sum(L.chars[a][j] * Rp[j] for j in range(7)) @ Rk for a in range(7)}
    colx = np.concatenate([np.arange(d * ds, (d + 1) * ds) for d in D.Sx])
    col0 = np.arange(D.d0 * ds, (D.d0 + 1) * ds)
    nA = len(Alg); N = 30 * nA
    Gm = np.zeros((N, N), dtype=complex)
    for a in req:
        T = Ma[a][:, colx]
        C = np.zeros((N, ds * len(colx)), dtype=complex)
        for dd in range(30):
            Td = t[dd] @ T[dd * ds:(dd + 1) * ds, :]
            for i, A in enumerate(Alg):
                C[dd * nA + i] = (A @ Td).reshape(-1)
        Gm += C @ C.conj().T
    ev, V = np.linalg.eigh(Gm)
    nul = V[:, ev < tol * max(1, ev.max())]
    res = {'dim_Alg': nA, 'dim_W': N, 'dim_null': nul.shape[1], 'Bd1_sample': len(Bd1)}
    if nul.shape[1] == 0:
        for a in range(7): res[a] = 0
        return res
    x = np.conj(nul @ (nrng.standard_normal(nul.shape[1]) + 1j * nrng.standard_normal(nul.shape[1])))   # row x with x C = 0
    w = np.zeros((ds, 30 * ds), dtype=complex)
    for dd in range(30):
        Z = sum(x[dd * nA + i] * A for i, A in enumerate(Alg))
        w[:, dd * ds:(dd + 1) * ds] = Z @ t[dd]
    wn = w / np.linalg.norm(w)
    for a in range(7):
        B0 = wn @ Ma[a][:, col0]
        sc = np.abs(Ma[a]).max()
        sv = np.linalg.svd(B0, compute_uv=False)
        res[a] = int((sv > 1e-8 * sc).sum())                     # absolute scale: w has unit norm
        res[('top', a)] = float(sv[0] / sc)
        res[('junk', a)] = float(np.abs(wn @ Ma[a][:, colx]).max() / sc)
        Ea = sum(L.chars[a][j] * sig(L.p[j]) for j in range(7)) / 7
        res[('m', a)] = int(round(np.trace(Ea).real))           # multiplicity of the a-eigenvalue of sigma(p)
    return res

if __name__ == '__main__':
    t0 = time.time()
    seed = int(sys.argv[1]); psi = int(sys.argv[2]); nBm = int(sys.argv[3]); names = sys.argv[4].split(',')
    L = Lift(seed=seed, psi=psi)
    rng = random.Random(seed + 1); nrng = np.random.default_rng(seed + 2)
    ys = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'ys_seed{seed}_psi{psi}.npy'))
    # regenerate run1's B_m letters (same rng order as run1.py)
    base = [b for (b, s) in L.Bm]
    Bm = list(base); seen = set(key(b) for b in Bm)
    while len(Bm) < nBm:
        w = mm(rng.choice(base), rng.choice(Bm))
        if key(w) not in seen: seen.add(key(w)); Bm.append(w)
    assert len(Bm) == len(ys)
    mi = minv(L.m); ks = [mm(L.m, b, mi) for b in Bm]
    D = D4(L)
    letters = [b for (b, s) in L.letters]
    # extra letters of B_d1 (products), to span its image
    Bd1 = [b for b in letters if D.act(b)[D.d1][0] == D.d1]
    rr = random.Random(5)
    while len(Bd1) < 600:
        w = mm(rr.choice(Bd1), rr.choice(Bd1)); Bd1.append(w)
    letters = letters + Bd1
    req = [0, 1, 2, 3, 4, 5, 6]   # all seven: F1's own character a = 0 too (the F2 criterion is read per character)
    for name in names:
        sig, ds = make_sig(name)
        # homomorphism check
        for _ in range(5):
            b1, b2 = random.choice(letters), random.choice(letters)
            assert np.allclose(sig(mm(b1, b2)), sig(b1) @ sig(b2)), name + ' not a homomorphism'
        res = gen_test(D, sig, ds, ys, ks, req, letters, rr, nrng)
        print(f'[d4gen] sigma = {name} (dim {ds}): dim Alg {res["dim_Alg"]}, dim W {res["dim_W"]}, null dim {res["dim_null"]}; '
              f'kernel-block rank (need {ds}) ' + ', '.join(f'a={a}:{res[a]}' for a in range(7)) +
              '; top sv/scale ' + ', '.join('%.0e' % res[('top', a)] for a in range(7) if ('top', a) in res) +
              '; junk/scale on required a ' + ', '.join('%.0e' % res[('junk', a)] for a in req if ('junk', a) in res) +
              '; F2 (twist p) needs r_a >= ds - m_a: m_a ' + ', '.join(str(res[('m', a)]) for a in range(7) if ('m', a) in res) +
              f'  {time.time() - t0:.0f}s', flush=True)
        # (characters decouple: alpha's a-components are independent, so the joint test equals the per-character tests)
    print('done %.0fs' % (time.time() - t0))
