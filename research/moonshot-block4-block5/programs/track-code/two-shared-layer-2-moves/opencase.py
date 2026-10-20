"""opencase.py -- w5-twoshared: decide the private cores where Lemma D's frame conditions HAVE solutions (esearch OPEN).
Rebuilds core c of esearch seed s, enumerates the symplectic e' of every solvable type combination for each top,
lifts e' to exact monomials (several Pauli parts), computes the three injection letters Lemma D prescribes
(J: x_j -> u_l or s_j -> v; k = 0, 1: x_k -> u_l or s_k -> v) exactly, adds them to the core and decides
(i) output e alone and (ii) both outputs (with the other top's injections) by closure + Singular, every coefficient free.
Usage: opencase.py seed core nlift cap"""
import sys, os, random, itertools, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from esearch import maps, basis, img_L0, GL, SYM
import sp
from sweep import core, moves_of
from exact import G, build, singular_script, decide
seed, cidx, nlift, cap = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = random.Random(seed)
for c in range(cidx + 1):
    cr = core(rng, 'priv')
R0, R1, _, _, tg = cr
R = [R0, R1]
Dm = G.D1; D = sp.D
ex = lambda r, k: r[k][1]
for r in R:
    r['qm'] = G.mul(Dm, ex(r, 'g'), Dm, ex(r, 'a'), Dm, ex(r, 'b'))
types = ['s->v', 'x->u0', 'x->u1']
rng2 = random.Random(1000 + seed)
os.makedirs('open_sing', exist_ok=True)
inj_by_top = {}
for e in (0, 1):
    j = 1 - e
    y = sp.comp(D, R[e]['mu'][0]); J = sp.comp(y, R[j]['q'])
    Jm = G.mul(Dm, ex(R[e], 'mu'), R[j]['qm'])
    X = sp.comp(D, sp.comp(R[e]['mu'][0], D))
    condJ = {'s->v': maps(basis(img_L0(D)), img_L0(sp.comp(J, sp.inv(R[j]['s']))))}
    for l in (0, 1):
        condJ[f'x->u{l}'] = maps(basis(img_L0(sp.comp(D, sp.comp(R[l]['g'][0], D)))), img_L0(J))
    condK = {}
    for k in (0, 1):
        Yk = sp.comp(D, sp.comp(R[k]['g'][0], sp.comp(D, R[k]['a'][0])))
        condK[(k, 's->v')] = maps(basis(img_L0(Yk)), img_L0(X))
        for l in (0, 1):
            W_ = sp.comp(D, sp.comp(R[e]['mu'][0], sp.comp(D, sp.comp(R[l]['g'][0], D))))
            condK[(k, f'x->u{l}')] = maps(basis(img_L0(R[k]['q'])), img_L0(W_))
    sols = []
    for tJ, t0_, t1_ in itertools.product(types, repeat=3):
        m = condJ[tJ] & condK[(0, t0_)] & condK[(1, t1_)]
        for (ip, ib) in np.argwhere(m):
            sols.append(((tJ, t0_, t1_), ip, ib))
    print(f"top {e}: {len(sols)} (type, e') solutions", flush=True)
    injs = []
    for (tt, ip, ib) in sols:
        for _ in range(nlift):
            ep = sp.to_matrix(list(GL[ip]), SYM[ib], rng2.randrange(16), rng2.randrange(16))
            new = []
            # J condition: e' w = J
            tJ = tt[0]
            if tJ == 's->v':
                cc = G.snap(np.conj(Dm).T @ np.conj(ep).T @ Jm @ np.conj(G.mul(Dm, ex(R[j], 'b'))).T)
                new.append((4, 2 + j, cc))                       # s_j into v (move 4); generator of s_j = 2 + j
            else:
                l = int(tJ[-1]); DgD = G.mul(Dm, ex(R[l], 'g'), Dm)
                kap = G.snap(np.conj(DgD).T @ np.conj(ep).T @ Jm)
                new.append((2 + l, j, kap))                      # x_j into u_l (move 2 + l)
            for k in (0, 1):
                tk = tt[1 + k]
                if tk == 's->v':
                    Xm = G.mul(Dm, ex(R[e], 'mu'), Dm); Ykm = G.mul(Dm, ex(R[k], 'g'), Dm, ex(R[k], 'a'))
                    cc = G.snap(np.conj(Xm).T @ ep @ Ykm)
                    new.append((4, 2 + k, cc))
                else:
                    l = int(tk[-1]); Wm = G.mul(Dm, ex(R[e], 'mu'), Dm, ex(R[l], 'g'), Dm)
                    kap = G.snap(np.conj(Wm).T @ ep @ R[k]['qm'])
                    new.append((2 + l, k, kap))
            ok = all(G.is_monomial(M) or np.allclose(np.abs(M[np.abs(M) > 1e-9]), np.abs(M[np.abs(M) > 1e-9]).max()) for (_, _, M) in new)
            # normalise scalar multiples to Z4 monomials
            norm = []
            for (mi, src, M) in new:
                nz = np.abs(M) > 1e-9; v = M[nz][0]
                norm.append((mi, src, G.snap(M / (v / abs(v)) / abs(v)) if not G.is_monomial(M) else M))
            if all(G.is_monomial(M) for (_, _, M) in norm):
                injs.append((tt, norm))
    inj_by_top[e] = injs
    print(f"top {e}: {len(injs)} exact injection sets", flush=True)
# (i) each output alone with its own top's injections
tally = {}
for e in (0, 1):
    for n, (tt, new) in enumerate(inj_by_top[e][:24]):
        # moves_of: move indices 0 s0, 1 s1, 2 u1(array0's merge), 3 u2, 4 v, 5 t0, 6 t1 ; extra = (mi, src, M)
        mvs = moves_of(R0, R1, None, None, [(mi, src, M) for (mi, src, M) in new])
        P = build(2, mvs)
        scr, sizes = singular_script(P, cap=cap, outputs=[e])
        if scr is None: v = 'CAP'
        else: v, out = decide(scr, f'open_sing/s{seed}_c{cidx}_t{e}.sing', timeout=120)
        tally[(e, v)] = tally.get((e, v), 0) + 1
        print(f"output {e} alone, types {tt}, set {n}: closure {sizes} {v}", flush=True)
        if v == 'CONSISTENT':
            pickle.dump({'moves': mvs, 'out': e}, open(f'open_sing/CONS_s{seed}_c{cidx}_t{e}_{n}.pkl', 'wb'))
print("TALLY", tally, flush=True)
