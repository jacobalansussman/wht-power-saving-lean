"""aut.py (t-km): automorphism group of P112 inside GL(9,2): linear maps L of F_2^9 with L(LAB) = LAB and
par(L x & L y) = par(x & y) on labels.  Backtracking over images of 9 independent labels.  Saves generators and
the full element list (as images of the standard basis) to aut.npz.  usage: python -I aut.py"""
import sys, time, itertools
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from model import *
import numpy as np
LS = set(LAB); NL = len(LAB)
# a basis of 9 labels
basis = []; ech = []
for u in LAB:
    x = u
    for y in ech: x = min(x, x ^ y)
    if x: ech.append(x); ech.sort(reverse=True); basis.append(u)
    if len(basis) == 9: break
assert len(basis) == 9
# coordinates of every label in that basis: solve via span table
span = {0: 0}
for i, b in enumerate(basis):
    span.update({x ^ b: c | (1 << i) for x, c in list(span.items())})
coord = {u: span[u] for u in LAB}          # coefficient mask over basis
dot = lambda a, b: par(a & b)
G0 = [[dot(a, b) for b in basis] for a in basis]
# labels expressed in the first k basis vectors (for incremental checks)
bylev = [[] for _ in range(10)]
for u in LAB:
    c = coord[u]; bylev[c.bit_length()].append((u, c))
def img(c, ims):
    x = 0
    for i in range(len(ims)):
        if c >> i & 1: x ^= ims[i]
    return x
sols = []
t0 = time.time()
def bt(ims):
    k = len(ims)
    if k == 9:
        sols.append(tuple(ims)); return
    for cand in LAB:
        if any(dot(cand, ims[j]) != G0[k][j] for j in range(k)): continue
        ims2 = ims + [cand]
        ok = True
        for (u, c) in bylev[k + 1]:
            y = img(c, ims2)
            if y not in LS: ok = False; break
        if ok:
            # dot products among labels at this level vs earlier: ensured by basis dots (bilinear) -- check quickly
            bt(ims2)
bt([])
print('automorphisms (label-preserving, dot-preserving):', len(sols), '%.1fs' % (time.time() - t0), flush=True)
# convert to matrices acting on F_2^9 standard basis: L(e_i).  basis matrix Bm columns = basis labels
def solve_std(ims):
    # find L with L(basis_i) = ims_i; L(e) for unit e: write e in basis coords
    out = []
    for i in range(9):
        e = 1 << i; c = span[e]
        out.append(img(c, list(ims)))
    return tuple(out)
Ls = np.array([solve_std(s) for s in sols], dtype=np.int64)
def apply(L, x):
    y = 0
    for i in range(9):
        if x >> i & 1: y ^= int(L[i])
    return y
# port permutations
PI = {u: i for i, u in enumerate(LAB)}
perms = np.array([[PI[apply(L, u)] for u in LAB] for L in Ls], dtype=np.int16)
np.savez('aut.npz', L=Ls, perm=perms, basis=np.array(basis))
# element orders histogram
def order(p):
    p = np.asarray(p); q = p.copy(); k = 1
    ident = np.arange(len(p))
    while not np.array_equal(q, ident): q = p[q]; k += 1
    return k
from collections import Counter
print('order histogram:', sorted(Counter(order(p) for p in perms).items()))
