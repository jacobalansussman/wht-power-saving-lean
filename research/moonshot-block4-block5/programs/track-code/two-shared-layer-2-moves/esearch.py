"""esearch.py -- w5-twoshared: EXHAUSTIVE frame-level test of Lemma D (private variant of the shape: each array's own
two-move chain to layer 2, ONE shared layer-3 move v, W private tops; here W = 2).
Lemma D: if the cross junk J = y_e q_j (j != e) of top e has layer 2, output e needs a free letter e' in supp(eps_e)
and layer-2 terms of v (made by injections: s_k -> v, or x_k -> u_l) with
   (J)  e' w = zeta J                      (w on array j),
   (k)  D mu_e w_k = zeta_k e' q_k          (w_k on array k, for EVERY array k).
Each condition says: the symplectic image S of e' (an element of Stab(L0)) maps one Lagrangian A to another B:
   (J, s_j->v)   S(S_D L0)          = S_J S_(Db_j)^-1 L0
   (J, x_j->u_l) S(S_(D g_l D) L0)  = S_J L0
   (k, s_k->v)   S(S_(D g_k D a_k) L0) = S_(D mu_e D) L0
   (k, x_k->u_l) S(S_(q_k) L0)      = S_(D mu_e D g_l D) L0
This script enumerates ALL 20160 x 1024 symplectic engine letters (Stab(L0)) and counts the e' meeting the three
conditions (J; k = 0; k = 1) for every choice of injection types (3 x 3 x 3 = 27).  Zero for every choice => the output
cannot be completed (for any injections, coefficients, letters of the injections) => the core is excluded.
Usage: esearch.py ncores seed"""
import sys, os, random, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sp
from sweep import core

# ---- all of Stab(L0): P in GL(4,2) (columns), B symmetric
def all_gl4():
    out = []
    for cols in itertools.product(range(1, 16), repeat=4):
        if len(sp.span(cols)) == 16: out.append(cols)
    return out
GL = all_gl4(); nP = len(GL)
def mat_apply(cols, x):
    r = 0
    for j in range(4):
        if (x >> j) & 1: r ^= cols[j]
    return r
PX = np.array([[mat_apply(c, x) for x in range(16)] for c in GL], dtype=np.int64)     # P x
PTZ = np.zeros((nP, 16), dtype=np.int64)
for i, c in enumerate(GL):
    Pinv = {mat_apply(c, v): v for v in range(16)}
    PT = [sum(((Pinv[1 << j] >> r) & 1) << j for j in range(4)) for r in range(4)]   # P^-T e_r
    PTZ[i] = [mat_apply(PT, z) for z in range(16)]
SYM = []
for bits in range(1024):
    Bm = [[0] * 4 for _ in range(4)]; t = 0
    for i in range(4):
        for j in range(i, 4):
            Bm[i][j] = Bm[j][i] = (bits >> t) & 1; t += 1
    SYM.append(Bm)
BT = np.array([[sum((sum(Bm[i][j] * ((y >> j) & 1) for j in range(4)) & 1) << i for i in range(4)) for y in range(16)] for Bm in SYM], dtype=np.int64)  # (1024, 16)

def maps(Abasis, Bset):
    """boolean (nP, 1024): S_(P,B)(A) == Bset"""
    tab = np.zeros(256, dtype=bool); tab[list(Bset)] = True
    ok = np.ones((nP, 1024), dtype=bool)
    for a in Abasis:
        ax, az = a & 15, a >> 4
        px = PX[:, ax]                       # (nP,)
        pz = PTZ[:, az]                      # (nP,)
        z2 = pz[:, None] ^ BT[:, px].T       # (nP, 1024)
        img = px[:, None] | (z2 << 4)
        ok &= tab[img]
    return ok

def basis(L):
    b = []; s = {0}
    for v in sorted(L):
        if v not in s:
            b.append(v); s |= {w ^ v for w in s}
    return b

def img_L0(S):
    """S(L0)"""
    return sp.span([sp.app(S, 1 << (4 + j)) for j in range(4)])

if __name__ == '__main__':
    ncores, seed = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(seed)
    # sanity: the identity letter is (P = I, B = 0) and maps every A to itself
    iI = GL.index((1, 2, 4, 8))
    assert maps(basis(sp.LK), sp.LK)[iI, 0]
    D = sp.D; Di = sp.inv(D)
    tot = {'excluded_outputs': 0, 'open_outputs': 0}
    for c in range(ncores):
        t0 = time.time()
        cr = core(rng, 'priv')
        R0, R1, _, _, tg = cr
        R = [R0, R1]
        q = [r['q'][0] if isinstance(r['q'], tuple) and len(r['q']) == 8 else r['q'] for r in R]
        q = [R0['q'], R1['q']]
        g = [R0['g'][0], R1['g'][0]]; a = [R0['a'][0], R1['a'][0]]; s = [R0['s'], R1['s']]
        mu = [R0['mu'][0], R1['mu'][0]]
        line = [f"core {c}: tags {tg}"]
        for e in (0, 1):
            j = 1 - e
            y = sp.comp(D, mu[e]); J = sp.comp(y, q[j])
            lJ = sp.lam(J)
            if lJ != 2:
                line.append(f"top {e}: lambda(J) = {lJ} (not the layer-2 case)"); continue
            X = sp.comp(D, sp.comp(mu[e], D))
            # J conditions
            condJ = {}
            condJ['s->v'] = maps(basis(img_L0(D)), img_L0(sp.comp(J, sp.inv(s[j]))))
            for l in (0, 1):
                DgD = sp.comp(D, sp.comp(g[l], D))
                condJ[f'x->u{l}'] = maps(basis(img_L0(DgD)), img_L0(J))
            condK = {}
            for k in (0, 1):
                Yk = sp.comp(D, sp.comp(g[k], sp.comp(D, a[k])))
                condK[(k, 's->v')] = maps(basis(img_L0(Yk)), img_L0(X))
                for l in (0, 1):
                    W_ = sp.comp(D, sp.comp(mu[e], sp.comp(D, sp.comp(g[l], D))))
                    condK[(k, f'x->u{l}')] = maps(basis(img_L0(q[k])), img_L0(W_))
            types = ['s->v', 'x->u0', 'x->u1']
            counts = {}
            for tJ in types:
                for t0_ in types:
                    for t1_ in types:
                        m = condJ[tJ] & condK[(0, t0_)] & condK[(1, t1_)]
                        n = int(m.sum())
                        if n: counts[(tJ, t0_, t1_)] = n
            single = {k: int(v.sum()) for k, v in list(condJ.items()) + [(str(k), v) for k, v in condK.items()]}
            if counts:
                tot['open_outputs'] += 1
                line.append(f"top {e}: OPEN, combos with solutions {counts}")
            else:
                tot['excluded_outputs'] += 1
                line.append(f"top {e}: EXCLUDED (no e' for any of the 27 injection-type choices); single-condition counts {single}")
        print("\n   ".join(line), f"({time.time()-t0:.0f}s)", flush=True)
    print("TOTAL", tot, flush=True)
