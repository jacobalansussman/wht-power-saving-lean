"""wsp.py -- w8-bpgl14: is the STRUCTURED top row block as good as free rows?  The top's rows are w = E_d1 R_sigma(A),
A in C[B]: Wsp = sum_d Alg t_d E_d (d4gen.py), Alg = sigma(C[B_d1]).  What enters every test is w R(e_a) (rows in M e_a).
If dim(Wsp R(e_a)) = d_sigma x rank R(e_a) (all such matrices), the rows-free criterion is exact for the structured top.
argv: seed psi sigma"""
import sys, os, time, random
import numpy as np, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from decomp_lib import *
t0 = time.time()
seed, psi, sname = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
L, D, sigJ, DS, p, p2 = load_twists(seed, psi, [sname])
letters = [b for (b, s) in L.letters]
Bd1 = [b for b in letters if D.act(b)[D.d1][0] == D.d1]
rr = random.Random(5)
while len(Bd1) < 600:
    Bd1.append(mm(rr.choice(Bd1), rr.choice(Bd1)))
Alg = alg_basis([sigJ(D.act(b)[D.d1][1]) for b in Bd1 + letters if D.act(b)[D.d1][0] == D.d1], rr)
t = [None] * 30
for b in letters:
    dd, beta = D.act(b)[D.d1]
    if t[dd] is None: t[dd] = sigJ(beta)
n = 30 * DS
print(f'[wsp] {sname}: dim Alg {len(Alg)} (full would be {DS * DS}); all 30 coset blocks reached: {all(x is not None for x in t)}  {time.time() - t0:.0f}s', flush=True)
for name, pp in (('p', p), ('p2', p2)):
    for a in range(7):
        E = sum(L.chars[a][j] * D.R(pp[j], sigJ, DS) for j in range(7)) / 7
        rk_e = int(np.linalg.matrix_rank(E, tol=1e-8))
        vecs = []
        for dd in range(30):
            Eb = E[dd * DS:(dd + 1) * DS, :]
            for Z in Alg: vecs.append((Z @ t[dd] @ Eb).reshape(-1))
        rk = np.linalg.matrix_rank(np.array(vecs), tol=1e-8)
        print(f'[wsp] twist {name} a={a}: rank R(e_a) {rk_e}; dim Wsp R(e_a) = {rk} of {DS * rk_e} -> '
              f'{"FULL (structured = rows free)" if rk == DS * rk_e else "PROPER subspace"}  {time.time() - t0:.0f}s', flush=True)
        if name == 'p2' and a > 0: break
print('done %.0fs' % (time.time() - t0))
