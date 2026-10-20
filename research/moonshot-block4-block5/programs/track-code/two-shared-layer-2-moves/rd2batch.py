"""rd2.py -- w5-twoshared (resumed): the SECOND ROUND of Lemma D in the private variant, over EVERY Lemma-D solution
of an open core (frame-exact, symplectic picture).
For a solution (types, e') of top e put Z = e'^-1 y_e (y_e = D mu_e).  Lemma D forces three layer-2 terms into v:
   w_J = Z q_j (array j, cancels J = y_e q_j through e' w_J),  w_k = Z^-1 q_k (k = 0, 1; y_e w_k cancels e' q_k).
Because t_e reads ALL of v through the single letter mu_e, t_e then contains
   N1 = y_e w_J = e' Z^2 q_j   (array j, layer <= 3).
Lemma D2 [proved in NOTE s3]: if lambda(N1) = 3, N1 can be cancelled in T_e = t_e + eps_e v only
   (ii)  by u_j injected into t_e with the letter D^-1 N1 p_j^-1 (needs Lag(D^-1 N1) = Lag(p_j)), or
   (iii) by a letter e' Z^2 of eps_e (needs Z^2 free, lambda(Z^2) = 0).
If neither holds the solution is dead.  Prints, per top, the counts of solutions by (lambda(N1), ii, iii) and the
layers of the other top's readings y_j w of the three forced terms (for the next round).
Batch version (rd2batch.py): Usage: rd2batch.py ncores seed [seed ...]; every private core of each seed is tested
(cores with no Lemma-D solution for some top are already excluded by Lemma D).
Lemma D2' (the other top j): a forced term w read by t_j at layer 3 on array k can be cancelled only by eps_j q_k
(y_j w in B'q_k) or by u_k injected into t_j (Lag(w) = Lag(p_k)); for w = w_J the latter is impossible (NOTE s3).
"""
import sys, os, random, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from esearch import maps, basis, img_L0, GL, SYM
import sp
from sweep import core
D = sp.D; comp = sp.comp; inv = sp.inv; lam = sp.lam; lag = sp.lag
types = ['s->v', 'x->u0', 'x->u1']
ncores = int(sys.argv[1]); TALLY = {}
for seed in [int(a) for a in sys.argv[2:]]:
  rng = random.Random(seed)
  for cidx in range(ncores):
    t0 = time.time()
    cr = core(rng, 'priv')
    if cr is None: print('core failed', flush=True); continue
    R0, R1, _, _, tg = cr
    R = [R0, R1]
    verdict = {}
    for e in (0, 1):
        j = 1 - e
        y = comp(D, R[e]['mu'][0]); J = comp(y, R[j]['q'])
        X = comp(D, comp(R[e]['mu'][0], D))
        condJ = {'s->v': maps(basis(img_L0(D)), img_L0(comp(J, inv(R[j]['s']))))}
        for l in (0, 1):
            condJ[f'x->u{l}'] = maps(basis(img_L0(comp(D, comp(R[l]['g'][0], D)))), img_L0(J))
        condK = {}
        for k in (0, 1):
            Yk = comp(D, comp(R[k]['g'][0], comp(D, R[k]['a'][0])))
            condK[(k, 's->v')] = maps(basis(img_L0(Yk)), img_L0(X))
            for l in (0, 1):
                W_ = comp(D, comp(R[e]['mu'][0], comp(D, comp(R[l]['g'][0], D))))
                condK[(k, f'x->u{l}')] = maps(basis(img_L0(R[k]['q'])), img_L0(W_))
        stats = {}; other = {}; both = {}; sanity_bad = 0; nsol = 0
        seenS = set()
        for tJ, t0_, t1_ in itertools.product(types, repeat=3):
            m = condJ[tJ] & condK[(0, t0_)] & condK[(1, t1_)]
            for (ip, ib) in np.argwhere(m):
                Se = sp.symp_from(list(GL[ip]), SYM[ib])
                nsol += 1
                Z = comp(inv(Se), y); Zi = inv(Z)
                wJ = comp(Z, R[j]['q']); wk = [comp(Zi, R[k]['q']) for k in (0, 1)]
                if lam(wJ) > 2 or any(lam(w) > 2 for w in wk) or lag(comp(Z, R[e]['q'])) != sp.LK:
                    sanity_bad += 1
                N1 = comp(Se, comp(Z, wJ))
                l1 = lam(N1)
                ii = lag(comp(inv(D), N1)) == lag(R[j]['p'])
                iii = lam(comp(Z, Z)) == 0
                uabs = lag(N1) == lag(R[j]['p'])          # for l1 = 2: N1 in the frame of u_j's main term
                key = (l1, 'ii' if ii else '-', 'iii' if iii else '-', 'u' if (l1 == 2 and uabs) else '-')
                stats[key] = stats.get(key, 0) + 1
                yj = comp(D, R[j]['mu'][0])
                ok = []
                for w, arr in ((wJ, j), (wk[0], 0), (wk[1], 1)):
                    g = comp(yj, w); L = lam(g)
                    if L == 4: ok.append('K' if lag(g) == sp.LK else '4')
                    elif L == 3: ok.append('3V' if lag(g) == lag(R[arr]['q']) else '3')
                    else: ok.append(str(L))
                other[tuple(ok)] = other.get(tuple(ok), 0) + 1
                # D2' for w_J read by the other top: layer 3 and not in B'q_j -> dead
                d2p = ok[0] == '3'
                d2 = (l1 == 3 and not ii and not iii)
                k2 = ('D2' if d2 else '-') + ('/D2p' if d2p else '/-')
                both[k2] = both.get(k2, 0) + 1
        dead = sum(v for k, v in both.items() if k != '-/-')
        verdict[e] = (nsol, dead)
        print(f"seed {seed} core {cidx} top {e}: {nsol} Lemma-D solutions; (lambda(N1), ii, iii, u-frame): {stats}; dead by Lemma D2: {dead}/{nsol}; sanity failures {sanity_bad}", flush=True)
        print(f"   other top's readings y_j (w_J, w_0, w_1): {other}; kills (D2 / D2'): {both}", flush=True)
    if any(verdict[e][0] == 0 for e in (0, 1)): v = 'excluded by Lemma D'
    elif any(verdict[e][0] == verdict[e][1] for e in (0, 1)): v = 'excluded by Lemma D2/D2p'
    else: v = 'OPEN after round 2'
    TALLY[v] = TALLY.get(v, 0) + 1
    print(f"seed {seed} core {cidx}: tags {tg}: {v} ({time.time()-t0:.0f}s)", flush=True)
print("TALLY", TALLY, flush=True)
