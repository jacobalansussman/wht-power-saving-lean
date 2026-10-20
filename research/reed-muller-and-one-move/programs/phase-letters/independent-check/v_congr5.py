import numpy as np
from vb import *
m = 5; gens = []
for i in range(m):
    for j in range(m):
        if i != j:
            E = np.eye(m, dtype=np.int64); E[i, j] = 1; gens.append(E)
I = np.eye(m, dtype=np.int64); seen = {key(I)}; fr = [I]
while fr:
    nf = []
    for S in fr:
        for E in gens:
            T = (E @ S @ E.T) % 2; k = key(T)
            if k not in seen: seen.add(k); nf.append(T)
    fr = nf
print("m=5: orbit of I under congruence (BFS over elementary generators):", len(seen))
