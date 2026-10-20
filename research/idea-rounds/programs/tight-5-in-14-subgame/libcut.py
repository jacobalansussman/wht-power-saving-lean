# libcut.py -- isotypic decomposition of a representation under a finite free group (all elements given as matrices)
# and exhaustive cut certificates: Y = sum of a subset of isotypic components; cap f(h), demand f(K), ratio.
import itertools, numpy as np
from lib514 import rank

def isotypic(mats, rng, ntest=60, tol=1e-7):
    """mats: list of orthogonal matrices rho(b) over the whole free group F. Returns list of components:
    dict(basis=orthonormal columns, d=irrep dim, m=multiplicity, chi=character on test elements)"""
    n = mats[0].shape[0]
    S = rng.standard_normal((n, n)); S = S + S.T
    R = sum(M @ S @ M.T for M in mats) / len(mats)
    w, V = np.linalg.eigh((R + R.T) / 2)
    pieces = []; st = 0
    for i in range(1, n + 1):
        if i == n or w[i] - w[i - 1] > tol * max(1.0, np.abs(w).max()):
            pieces.append(V[:, st:i]); st = i
    tests = [mats[j] for j in rng.choice(len(mats), size=min(ntest, len(mats)), replace=False)]
    comps = {}
    for Q in pieces:
        chi = tuple(np.round([np.trace(Q.T @ T @ Q) for T in tests], 5) + 0.0)
        comps.setdefault(chi, []).append(Q)
    out = []
    for chi, Qs in comps.items():
        B = np.hstack(Qs); out.append(dict(basis=B, d=Qs[0].shape[1], m=len(Qs), chi=chi))
    return out

def f_of(Y, g):
    return rank(np.hstack([Y, g @ Y])) - Y.shape[1]

def best_cuts(comps, rh, rK, extra=None, top=5):
    """all nonempty proper subsets of components; returns sorted list of (ratio, dimY, cap, demand, subset)"""
    res = []
    k = len(comps)
    for r in range(1, k + 1):
        for sub in itertools.combinations(range(k), r):
            Y = np.hstack([comps[i]['basis'] for i in sub])
            if Y.shape[1] == Y.shape[0]: continue
            c = f_of(Y, rh)
            if c == 0: continue
            dmd = f_of(Y, rK)
            row = (dmd / c, Y.shape[1], c, dmd, sub)
            if extra is not None: row = row + tuple(f_of(Y, e) for e in extra)
            res.append(row)
    res.sort(key=lambda t: -t[0])
    return res[:top], res
