# vgen.py -- v-nc-gain: generate my own independent layout samples for the Colab screens (my sampler, my seeds)
import sys, json, random, shutil
import numpy as np
sys.path.insert(0, '.')
import vlemmas
def load(m):
    E = np.load(f'gamma{m}_edges_dpgidx.npy'); D = np.load(f'gamma{m}_D_dpgidx.npy').astype(np.int16)
    class G: pass
    g = G(); g.m = m; g.n = D.shape[0]; g.D = D; g.edges = [tuple(map(int, e)) for e in E]
    g.adj = [np.nonzero(D[i] == 1)[0] for i in range(g.n)]
    return g
def hill(g, rng, iters=400):
    m = g.m; P = vlemmas.rand_layout(g, rng, 3); e = vlemmas.edef(g, P)[0]
    for it in range(iters):
        k = rng.randrange(3); end = rng.randrange(2); s, t = P[k]; x = (s, t)[end]; y = (t, s)[end]
        nb = int(rng.choice(list(g.adj[x])))
        if g.D[nb, y] != m or nb in {v for p in P for v in p}: continue
        Q = list(P); Q[k] = (nb, y) if end == 0 else (y, nb)
        e2 = vlemmas.edef(g, Q)[0]
        if e2 >= e: P, e = Q, e2
    return P, e
g3 = load(3); rng = random.Random(4242)
T = []; ncomp = 0; tot = 0
while len(T) < 15:
    P = vlemmas.rand_layout(g3, rng, 3); tot += 1
    if vlemmas.edef(g3, P)[0] == 0: ncomp += 1; continue
    T.append(dict(m=3, P=P, mode='undirected', claim='nc-gain: every block-3 three-pair layout found had 9'))
print(f"block 3: {tot} random layouts drawn, {ncomp} compatible, {len(T)} non-compatible sent")
k = 0
while k < 5:
    P, e = hill(g3, rng)
    if e >= 2: T.append(dict(m=3, P=P, mode='undirected', claim=f'hill-climbed e={e}')); k += 1
json.dump(T, open('../job_g/tasks.json', 'w'))
g4 = load(4); rng = random.Random(777)
T = []; ncomp = 0; tot = 0
while len(T) < 80:
    P = vlemmas.rand_layout(g4, rng, 3); tot += 1
    if vlemmas.edef(g4, P)[0] == 0: ncomp += 1; continue
    T.append(dict(m=4, P=P, mode='image', claim=f'e={vlemmas.edef(g4, P)[0]}'))
print(f"block 4: {tot} random layouts drawn, {ncomp} compatible ({100*ncomp/tot:.1f}%), {len(T)} non-compatible sent")
json.dump(T, open('../job_h/tasks.json', 'w'))
json.dump([dict(m=4, P=[(787,1231),(576,1200),(1677,534)], mode='directed', solver='pdlp', tl=5000, claim='NOT sent by nc-gain; e=3, image 10.5')], open('../job_e/tasks.json', 'w'))
json.dump([dict(m=4, P=[(954,702),(1849,1516),(1697,1353)], mode='directed', solver='pdlp', tl=5000, claim='NOT sent by nc-gain; e=3, image 10.57')], open('../job_f/tasks.json', 'w'))
for j in 'efgh':
    for f in ['vrun.py', 'vcut.py']: shutil.copy(f, f'../job_{j}/')
for j, m in (('e', 4), ('f', 4), ('g', 3), ('h', 4)):
    for f in [f'gamma{m}_edges_dpgidx.npy', f'gamma{m}_D_dpgidx.npy']: shutil.copy(f, f'../job_{j}/')
