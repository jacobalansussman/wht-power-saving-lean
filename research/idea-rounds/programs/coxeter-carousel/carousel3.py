"""Recount distinct frame paths of the 192 block-3 Coxeter-type steps with canonical frames (full span)."""
import itertools
import numpy as np
from coxeter import sp_h, mul, inv, lag_dist
from carousel2 import P_block3

def span_key(F):
    vecs = set()
    cols = [F[:, k] for k in range(F.shape[1])]
    for bits in itertools.product([0, 1], repeat=len(cols)):
        v = np.zeros(F.shape[0], dtype=np.uint8)
        for b, col in zip(bits, cols):
            if b:
                v ^= col
        vecs.add(v.tobytes())
    return frozenset(vecs)

m = 3
h1 = sp_h(m, 0)
LZ = np.zeros((6, 3), dtype=np.uint8); LZ[3:, :] = np.eye(3, dtype=np.uint8)
LX = np.zeros((6, 3), dtype=np.uint8); LX[:3, :] = np.eye(3, dtype=np.uint8)
paths, full_loops = {}, 0
for p in P_block3():
    c = mul(p, h1)
    Q = np.eye(6, dtype=np.uint8); fr = []; ok = True
    for j in range(1, m + 1):
        Q = mul(c, Q); F = mul(inv(Q), LZ)
        if lag_dist(F, LZ, m) != j: ok = False; break
        fr.append(span_key(F))
    if ok and lag_dist(F, LX, m) == 0:
        key = tuple(fr)
        paths[key] = paths.get(key, 0) + 1
        # does the walk come back to L_Z after 2m rounds (a closed Petrie polygon)?
        Q2 = Q.copy()
        for j in range(m):
            Q2 = mul(c, Q2)
        if lag_dist(mul(inv(Q2), LZ), LZ, m) == 0:
            full_loops += 1
txt = (f"block 3: Coxeter-type steps 192 -> distinct geodesic frame paths L_Z->L_X: {len(paths)} "
       f"(maximal flags of F_2^3: 21); steps per path: {sorted(set(paths.values()))}; "
       f"steps whose walk returns to L_Z after 2m=6 rounds: {full_loops}")
print(txt); open("out_carousel3.txt", "w").write(txt + "\n")
