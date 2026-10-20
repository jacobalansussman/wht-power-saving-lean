"""esearchW.py -- w5-twoshared: Lemma D exhaustive test for the PRIVATE variant at W = 3 and 4 (3 in 10, 4 in 13:
W private two-move chains, ONE shared layer-3 move v, W tops).  Cores: a pool of independent routes
x_j -> s_j -> u_j -> v -> t_j is built, and W-cliques are taken in the compatibility graph
  e ~ j  iff  lambda(y_e q_j) = 2 = lambda(y_j q_e)  and  N(q_e) != N(q_j)
(Lemma K + Lemma C': in a unit every cross junk has layer 2 and the top frames differ; layer 3/4 junk or a shared top
frame is fatal).  For every top e and every cross junk J = y_e q_j, Lemma D needs ONE letter e' meeting the J-condition
and the k-condition for every array k (types: s -> v, or x -> u_l for each l); all of Stab(L0) is enumerated.
A core is excluded as soon as one junk has no e' for every choice of types.
Usage: esearchW.py W npool seed maxcores"""
import sys, os, random, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from esearch import maps, basis, img_L0, GL, SYM
import sp
from sweep import route

W, npool, seed, maxc = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = random.Random(seed)
D = sp.D
t0 = time.time()
pool = [route(rng) for _ in range(npool)]
for r in pool:
    r['y'] = sp.comp(D, r['mu'][0]); r['N'] = sp.lag(r['y'])
print(f"pool of {npool} routes built ({time.time()-t0:.0f}s)", flush=True)
ok = np.zeros((npool, npool), dtype=bool)
for e in range(npool):
    for j in range(npool):
        if e == j: continue
        ok[e, j] = sp.lam(sp.comp(pool[e]['y'], pool[j]['q'])) == 2
adj = ok & ok.T
for e in range(npool):
    for j in range(npool):
        if adj[e, j] and pool[e]['N'] == pool[j]['N']: adj[e, j] = False
print(f"edge density {adj.sum() / (npool * (npool - 1)):.4f}", flush=True)
cliques = []
def extend(cl, cand):
    if len(cliques) >= maxc: return
    if len(cl) == W: cliques.append(cl); return
    for v in cand:
        extend(cl + [v], [w for w in cand if w > v and adj[v, w]])
extend([], list(range(npool)))
print(f"{len(cliques)} {W}-cliques (cores) taken ({time.time()-t0:.0f}s)", flush=True)
tot = {'excluded': 0, 'open': 0}
for ci, cl in enumerate(cliques):
    R = [pool[i] for i in cl]
    verdict = 'excluded'; info = []
    open_tops = 0
    for e in range(W):
        X = sp.comp(D, sp.comp(R[e]['mu'][0], D))
        condK = {}
        for k in range(W):
            Yk = sp.comp(D, sp.comp(R[k]['g'][0], sp.comp(D, R[k]['a'][0])))
            condK[(k, 's')] = maps(basis(img_L0(Yk)), img_L0(X))
            for l in range(W):
                W_ = sp.comp(D, sp.comp(R[e]['mu'][0], sp.comp(D, sp.comp(R[l]['g'][0], D))))
                condK[(k, l)] = maps(basis(img_L0(R[k]['q'])), img_L0(W_))
        # OR over types for each k: e' must satisfy some type for every k -> AND_k OR_types
        allk = np.ones_like(condK[(0, 's')])
        for k in range(W):
            orK = condK[(k, 's')].copy()
            for l in range(W): orK |= condK[(k, l)]
            allk &= orK
        top_ok = True
        for j in range(W):
            if j == e: continue
            J = sp.comp(R[e]['y'], R[j]['q'])
            cj = maps(basis(img_L0(D)), img_L0(sp.comp(J, sp.inv(R[j]['s']))))
            for l in range(W):
                cj |= maps(basis(img_L0(sp.comp(D, sp.comp(R[l]['g'][0], D)))), img_L0(J))
            n = int((cj & allk).sum())
            info.append(f"t{e}/J{j}:{n}")
            if n == 0: top_ok = False
        if top_ok: open_tops += 1
    verdict = 'open' if open_tops == W else 'excluded'
    tot[verdict] += 1
    print(f"core {ci} {cl}: {verdict} (tops passing Lemma D: {open_tops}/{W}) {' '.join(info)}", flush=True)
print("TOTAL", tot, f"({time.time()-t0:.0f}s)", flush=True)
