# vtri.py -- v-nc-gain: the triangle + pendant W = 4 layout. Its distance pattern at block 3, and block-4 analogues
# (three mutually opposite frames a,b,c; pendant (x,y) with the same distance pattern scaled) for a full Gamma_4 LP.
import sys, json, random, itertools
import numpy as np
sys.path.insert(0, '.')
import vlemmas
from vgen import load
g3 = load(3)
P = [(0, 95), (0, 80), (80, 95), (84, 60)]
T = [0, 80, 95]
print("block 3 triangle", T, "pendant (84,60): d(84, T) =", [int(g3.D[84, t]) for t in T], " d(60, T) =", [int(g3.D[60, t]) for t in T])
g4 = load(4); rng = random.Random(11)
a = 0; opp = [int(v) for v in np.nonzero(g4.D[a] == 4)[0]]
sel = []
tries = 0
while len(sel) < 2 and tries < 200000:
    tries += 1
    b, c = rng.sample(opp, 2)
    if g4.D[b, c] != 4: continue
    x = int(rng.choice([int(v) for v in g4.adj[b]])); y = int(rng.choice([int(v) for v in g4.adj[c]]))
    if g4.D[x, y] != 4 or len({a, b, c, x, y}) < 5: continue
    Q = [(a, b), (a, c), (b, c), (x, y)]
    if vlemmas.two_class(g4, Q) is not None: continue
    prof = ([int(g4.D[x, t]) for t in (a, b, c)], [int(g4.D[y, t]) for t in (a, b, c)])
    if sel and prof == sel[-1][1]: continue
    sel.append((Q, prof))
print("block-4 analogues (x adjacent to b, y adjacent to c, x opposite y, not two-class stretchable):", sel)
json.dump([dict(m=4, P=q, mode='undirected', solver='pdlp', tl=5000, claim='W=4 block 4 triangle+pendant: cut gap (< 16)?') for q, p in sel], open('tri4_tasks.json', 'w'))
