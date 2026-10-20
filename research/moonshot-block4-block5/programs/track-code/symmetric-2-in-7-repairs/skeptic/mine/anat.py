import sys, pickle
import numpy as np
from mycheck import *
dmp = pickle.load(open(sys.argv[1], 'rb')); m = dmp['m']
isfree = make_free_test('exact', m); K = dmp['K']; gens = dmp['gens']
al = np.array(dmp['alpha'], dtype=complex); al[dmp['free_idx'][0]] = 0.37
cv = coefs(dmp, al)
for l in range(dmp['W']):
    rows, cols, ent = closure(dmp, l, isfree)
    _, A = residual(len(rows), len(cols), ent, cv)
    print(f"array {l} closure matrix at a2 = 0.37 (columns: generator index):", [i for (i, k) in cols])
    print(np.round(A.real, 3))
# relative letter between the two kernel-frame terms of g7 (array 0)
T = gens[7]; hs = [h for (l, h) in T if l == 0]
u = np.conj(hs[0]).T @ hs[1]           # h0^-1 h1 (right ratio); conjugate by K h0^-1 for the left ratio at K
uL = K @ u @ np.conj(K).T
for name, M in (("h0^-1 h1", u), ("K h0^-1 h1 K^-1", uL)):
    P = np.eye(16, dtype=complex); o = 0
    while True:
        P = P @ M; o += 1
        if np.allclose(P, np.eye(16)) or o > 64: break
    print(name, "free(E):", isfree(M), "order:", o, "scalar:", np.allclose(M, M[0, 0] * np.eye(16)), "eig:", np.round(np.unique(np.round(np.linalg.eigvals(M), 6)), 3))
# junk letters: array-1 terms of g7 vs g5's terms
print("g7 x1 cells vs g5 cells: ratio letters free:", [[isfree(h7 @ np.conj(h5).T) for (l5, h5) in gens[5]] for (l7, h7) in T if l7 == 1])
