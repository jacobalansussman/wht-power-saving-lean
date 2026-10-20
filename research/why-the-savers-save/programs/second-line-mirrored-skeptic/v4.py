#!/usr/bin/env python3
"""v-cx-mirror / v4.py (own code, stdlib).
 (a) The scatter.  The claim prices the V-centre as 'the merged scatter copy holds all 56 readers' lines and drops 8
     levels per total'.  A copy that has met y_t (frame contains v_t) can meet y_t' later without a drop only if
     v_t . u_t' = 0, i.e. N[t][t'] = 0 (y frames stay inside u^perp).  I compute, per total, the readers, the junk
     pairs among them and the largest drop-free scatter run (natural mirror and the best twist).
 (b) A random sample of twisted mirrors A = I + M (M = sum a_i b_i^T + b_i a_i^T over an isotropic set of even
     vectors), built by my own generator: orthogonality, the port condition, symmetry of N, and the largest drop-free
     run inside each support (exact, my own clique search) -- a spot check of the Colab job's 'never above 8'.
 (c) A longer annealing of the lazy-drop cleaning order for the best twist on two supports."""
import gzip, json, os, random, time, sys
T0 = time.time()
here = os.path.dirname(os.path.abspath(__file__))
d = json.load(gzip.open(os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))))
U = list(d['ports']); v = d['v']; J = (1 << 9) - 1
V1 = json.load(open(os.path.join(here, 'v1.json'))); supp = V1['supp']; readers = V1['readers']
def dot(a, b): return bin(a & b).count('1') & 1
def apply(A, x):
    y = 0
    for i in range(9):
        if x >> i & 1: y ^= A[i]
    return y
def bits(x):
    while x:
        l = x & -x; yield l.bit_length() - 1; x ^= l
def maxclique(cand, nb):
    best = [0]
    def bk(Rn, Pm, Xm):
        if not Pm and not Xm:
            best[0] = max(best[0], Rn); return
        if Rn + bin(Pm).count('1') <= best[0]: return
        piv = max(bits(Pm | Xm), key=lambda u: bin(Pm & nb[u]).count('1'))
        for u in list(bits(Pm & ~nb[piv])):
            bk(Rn + 1, Pm & nb[u], Xm & nb[u]); Pm &= ~(1 << u); Xm |= 1 << u
    bk(0, cand, 0); return best[0]
def Nmat(A):
    G = [J ^ apply(A, u) for u in U]
    return G, [[0 if a == b else dot(G[a], U[b]) for b in range(v)] for a in range(v)]
def free_nb(N):   # adjacency of 'can follow each other drop-free' (N = 0 both ways)
    return [sum(1 << b for b in range(v) if b != a and N[a][b] == 0 and N[b][a] == 0) for a in range(v)]
tw = json.load(open(os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../second-line-mirrored/all-twists/twistall.json'))))
I9 = [1 << i for i in range(9)]; Ab = [(1 << i) ^ m for i, m in enumerate(tw['best']['M'])]
out = {}
# (a)
for name, A in (('natural', I9), ('best_twist', Ab)):
    G, N = Nmat(A); nb = free_nb(N)
    same = [sorted(readers[k]) == sorted(supp[k]) for k in range(8)]
    jr = [sum(1 for i, a in enumerate(readers[k]) for b in readers[k][i + 1:] if dot(U[a], U[b]) == 0) for k in range(8)]
    runs = [maxclique(sum(1 << t for t in readers[k]), nb) for k in range(8)]
    # can ONE frame lie inside every reader's u^perp and contain every reader's v?  (needs v_t . u_t' = 0 for all pairs)
    allpairs = all(dot(S_ := (1 << 9) | G[t], U[t2]) == 0 for k in range(8) for t in readers[k] for t2 in readers[k])
    print('%-10s readers == support %s | orthogonal (junk) pairs among readers %s | longest drop-free scatter run %s | one frame for all readers: %s'
          % (name, same, jr, runs, allpairs))
    out['scatter_' + name] = dict(readers_eq_support=same, junk_among_readers=jr, scatter_run=runs)
# (b)
rnd = random.Random(7)
Ev = [e for e in range(1, 512) if bin(e).count('1') % 2 == 0]
def rand_twist():
    k = rnd.choice([1, 2])
    while True:
        vecs = []
        for _ in range(2 * k):
            c = [e for e in Ev if all(dot(e, w) == 0 for w in vecs)]
            # independence
            def indep(ws):
                b = []
                for x in ws:
                    for y in b: x = min(x, x ^ y)
                    if not x: return False
                    b.append(x)
                return True
            c = [e for e in c if indep(vecs + [e])]
            if not c: break
            vecs.append(rnd.choice(c))
        if len(vecs) == 2 * k: break
    M = [0] * 9
    for i in range(k):
        a, b = vecs[2 * i], vecs[2 * i + 1]
        for col in range(9):    # column col of a b^T + b a^T is a * b_col + b * a_col
            M[col] ^= (a if b >> col & 1 else 0) ^ (b if a >> col & 1 else 0)
    return [(1 << i) ^ M[i] for i in range(9)]
hist = {}; bad = 0; nsamp = 0
t1 = time.time()
while time.time() - t1 < 35:
    A = rand_twist(); nsamp += 1
    orth = all(dot(A[i], A[j]) == (1 if i == j else 0) for i in range(9) for j in range(9))
    cond = all(dot(u, apply(A, u)) == 1 for u in U)
    G, N = Nmat(A)
    sym = all(N[a][b] == N[b][a] for a in range(v) for b in range(v))
    if not (orth and cond and sym): bad += 1; continue
    nb = free_nb(N)
    runs = [maxclique(sum(1 << b for b in S), nb) for S in supp]
    key = (min(runs), max(runs)); hist[str(key)] = hist.get(str(key), 0) + 1
print('random twists sampled', nsamp, '| failing orthogonality / port condition / symmetry:', bad,
      '| (worst support, best support) longest run:', dict(sorted(hist.items())))
out['twist_sample'] = dict(n=nsamp, bad=bad, hist=hist)
# (c) longer anneal, best twist, supports 0 and 1
G = [J ^ apply(Ab, u) for u in U]
def red(B):
    o = []
    for x in B:
        for y in o: x = min(x, x ^ y)
        if x: o.append(x)
    return o
def sect(F, g):
    badv = [f for f in F if dot(f, g)]
    if not badv: return F, 0
    p = badv[0]
    return red([f ^ p if dot(f, g) else f for f in F if f != p]), 1
def cost(order):
    F = []; c = 0
    for b in order:
        F, x = sect(F, G[b]); c += x; F = red(F + [U[b]])
    return c
prev = json.load(open(os.path.join(here, 'v3order.json')))['best_twist']
res = []
for k in (0, 1):
    cur = list(prev[k]['order']); cc = cost(cur); best = cc
    t1 = time.time(); it = 0
    while time.time() - t1 < 18:
        it += 1
        nxt = list(cur); i, j = rnd.randrange(56), rnd.randrange(56)
        if rnd.random() < 0.5: nxt[i], nxt[j] = nxt[j], nxt[i]
        else: x = nxt.pop(i); nxt.insert(j, x)
        nc = cost(nxt)
        if nc <= cc or rnd.random() < 0.25 ** (nc - cc):
            cur, cc = nxt, nc; best = min(best, cc)
    res.append(dict(support=k, start=prev[k]['cost'], best=best, iters=it))
print('longer anneal, best twist:', res)
out['anneal'] = res
json.dump(out, open(os.path.join(here, 'v4.json'), 'w'), indent=1)
print('secs %.1f' % (time.time() - T0))
