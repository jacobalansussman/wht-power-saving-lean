"""check_caseB.py -- independent check of a K=10 case-B decoder JSON: fit, W0 not inside W, exact star cost per total."""
import sys, json, itertools
from fractions import Fraction
import numpy as np
gen = np.array([[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[0,1,0,1,0,1,0,1]])
code = {tuple((np.array(c) @ gen) % 2) for c in itertools.product([0,1], repeat=4)}
w4 = sorted(c for c in code if sum(c) == 4)
Ls = [tuple(2*(np.arange(8)==a).astype(int)) for a in range(8)]
for c in w4:
    idx = [i for i in range(8) if c[i]]
    for sg in itertools.product([1,-1], repeat=3):
        x = np.zeros(8, dtype=int); x[idx[0]] = 1
        for k in range(3): x[idx[k+1]] = sg[k]
        Ls.append(tuple(x))
L = np.array(Ls, dtype=np.int64); G = (L @ L.T) // 2
d = json.load(open(sys.argv[1])); live = d['live']; X = np.array(d['X']); Y = np.array(d['Y']); v = len(live)
B = (Y.T @ X) % 2
ok = all(B[s, s] == 1 for s in range(v)) and all(B[s, t] == 0 for s in range(v) for t in range(v) if s != t and G[live[s], live[t]] != 0)
def f2rank(rows):
    basis = []
    for x in rows:
        for b in basis: x = min(x, x ^ b)
        if x: basis.append(x); basis.sort(reverse=True)
    return len(basis)
toint = lambda r: int(''.join(str(int(b)) for b in r[::-1]), 2)
B0 = ((1 + G) % 2)[np.ix_(live, live)]
rW = f2rank([toint(r) for r in X]); rW0 = f2rank([toint(r) for r in B0]); rBoth = f2rank([toint(r) for r in X] + [toint(r) for r in B0])
def qrank(idx):
    M = [[Fraction(int(z)) for z in L[i]] for i in idx]; rk = 0
    for c in range(8):
        p = next((k for k in range(rk, len(M)) if M[k][c] != 0), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        for k in range(len(M)):
            if k != rk and M[k][c] != 0:
                f = M[k][c] / M[rk][c]; M[k] = [a - f * b for a, b in zip(M[k], M[rk])]
        rk += 1
    return rk
st = [qrank([live[t] for t in range(v) if X[j][t]]) for j in range(len(X))]
print(sys.argv[1].split('/')[-1], 'fits', ok, 'dim W', rW, 'dim W0', rW0, 'dim(W+W0)', rBoth, '-> W0 inside W:', rBoth == rW, ' star per total', st, 'sum', sum(st))
