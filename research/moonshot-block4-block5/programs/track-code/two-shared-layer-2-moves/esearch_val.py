"""esearch_val.py -- validates esearch.py's four condition types with EXACT matrices: for random private cores, pick
symplectic letters S_e' meeting each single condition, lift them to exact monomials e', compute the injection letter
the condition promises (c or kappa) and check (i) it is a monomial (a letter, up to a scalar) and (ii) the defining
identity holds exactly up to a scalar: (J,s->v) e' D c D b_j = z J; (J,x->u_l) e' D g_l D kappa = z J;
(k,s->v) D mu_e D c D b_k = z e' q_k; (k,x->u_l) D mu_e D g_l D kappa = z e' q_k."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from esearch import *
from exact import G
rng = random.Random(3)
Dm = G.D1; Dmi = np.conj(Dm).T
def ex(R, key): return R[key][1]
def prop(A, B):
    i = np.unravel_index(np.argmax(np.abs(B)), B.shape); z = A[i] / B[i]
    return abs(abs(z) - 1) < 1e-9 and np.allclose(A, z * B, atol=1e-9)
def mono_upto_scalar(M):
    nz = np.abs(M) > 1e-9
    if not ((nz.sum(0) == 1).all() and (nz.sum(1) == 1).all()): return False
    v = M[nz]; return np.allclose(np.abs(v), np.abs(v[0])) and np.allclose((v / v[0]) ** 4, 1)
checked = 0; bad = 0
for c in range(3):
    R0, R1, _, _, tg = core(rng, 'priv')
    R = [R0, R1]
    for r in R:
        r['qm'] = G.mul(Dm, ex(r, 'g'), Dm, ex(r, 'a'), Dm, ex(r, 'b'))
    for e in (0, 1):
        j = 1 - e
        y = sp.comp(sp.D, R[e]['mu'][0]); J = sp.comp(y, R[j]['q'])
        Jm = G.mul(Dm, ex(R[e], 'mu'), R[j]['qm'])
        X = sp.comp(sp.D, sp.comp(R[e]['mu'][0], sp.D))
        conds = []
        conds.append(('J s->v', basis(img_L0(sp.D)), img_L0(sp.comp(J, sp.inv(R[j]['s'])))))
        for l in (0, 1):
            conds.append((f'J x->u{l}', basis(img_L0(sp.comp(sp.D, sp.comp(R[l]['g'][0], sp.D)))), img_L0(J), l))
        for k in (0, 1):
            Yk = sp.comp(sp.D, sp.comp(R[k]['g'][0], sp.comp(sp.D, R[k]['a'][0])))
            conds.append((f'k{k} s->v', basis(img_L0(Yk)), img_L0(X), k))
            for l in (0, 1):
                W_ = sp.comp(sp.D, sp.comp(R[e]['mu'][0], sp.comp(sp.D, sp.comp(R[l]['g'][0], sp.D))))
                conds.append((f'k{k} x->u{l}', basis(img_L0(R[k]['q'])), img_L0(W_), k, l))
        for cd in conds:
            m = maps(cd[1], cd[2]); idx = np.argwhere(m)
            if len(idx) == 0: continue
            for (ip, ib) in idx[rng.sample(range(len(idx)), min(3, len(idx)))]:
                P = list(GL[ip]); Bm = SYM[ib]
                ep = sp.to_matrix(P, Bm, rng.randrange(16), rng.randrange(16))
                assert sp.symp_from(P, Bm) is not None
                name = cd[0]
                if name == 'J s->v':
                    cc = G.snap(Dmi @ np.conj(ep).T @ Jm @ np.conj(G.mul(Dm, ex(R[j], 'b'))).T)
                    ok = mono_upto_scalar(cc) and prop(G.snap(ep @ Dm @ cc @ Dm @ ex(R[j], 'b')), Jm)
                elif name.startswith('J x->u'):
                    l = cd[3]; DgD = G.mul(Dm, ex(R[l], 'g'), Dm)
                    kap = G.snap(np.conj(DgD).T @ np.conj(ep).T @ Jm)
                    ok = mono_upto_scalar(kap) and prop(G.snap(ep @ DgD @ kap), Jm)
                elif 's->v' in name:
                    k = cd[3]; Xm = G.mul(Dm, ex(R[e], 'mu'), Dm)
                    Ykm = G.mul(Dm, ex(R[k], 'g'), Dm, ex(R[k], 'a'))
                    cc = G.snap(np.conj(Xm).T @ ep @ Ykm)
                    ok = mono_upto_scalar(cc) and prop(G.snap(Xm @ cc @ Dm @ ex(R[k], 'b')), G.snap(ep @ R[k]['qm']))
                else:
                    k, l = cd[3], cd[4]; Wm = G.mul(Dm, ex(R[e], 'mu'), Dm, ex(R[l], 'g'), Dm)
                    kap = G.snap(np.conj(Wm).T @ ep @ R[k]['qm'])
                    ok = mono_upto_scalar(kap) and prop(G.snap(Wm @ kap), G.snap(ep @ R[k]['qm']))
                checked += 1; bad += 0 if ok else 1
                if not ok: print("FAIL", name, flush=True)
print("validated", checked, "condition instances with exact matrices; failures", bad)
