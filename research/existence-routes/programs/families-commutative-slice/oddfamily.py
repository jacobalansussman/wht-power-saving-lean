"""oddfamily.py: the second natural family V_m = all odd-weight directions of F_2^m (the vectors with z.z = 1 for the
kernel's form), W_m = the deletion bound (min weight of an odd-trace form on V_m); clique test (Theorem 1) on it.
Lean bit-level enumeration of all 2^(m(m+1)/2) binary quadratic forms (m <= 6 on the Mac).
Also works for the uniform family (all nonzero directions) as a check: argument 'all'."""
import sys, numpy as np
from clique_test import max_clique

m = int(sys.argv[1]); which = sys.argv[2] if len(sys.argv) > 2 else "odd"
pairs = [(j, k) for j in range(m) for k in range(j, m)]
P = len(pairs)
dirs = [z for z in range(1, 1 << m) if (which == "all" or bin(z).count("1") % 2 == 1)]
mon = np.array([sum(1 << b for b, (j, k) in enumerate(pairs) if (z >> j & 1) and (z >> k & 1)) for z in dirs], dtype=np.uint32)
diagmask = sum(1 << pairs.index((j, j)) for j in range(m))
nf = 1 << P
POP = np.array([bin(i).count("1") & 1 for i in range(1 << 16)], dtype=np.uint8)
def par(x):
    return POP[x & 0xFFFF] ^ POP[x >> 16]
best = None
chunk = 1 << 18
wts = np.zeros(nf, dtype=np.int16)
odd = np.zeros(nf, dtype=bool)
for s in range(0, nf, chunk):
    q = np.arange(s, min(nf, s + chunk), dtype=np.uint32)
    odd[s:s + len(q)] = par(q & np.uint32(diagmask)).astype(bool)
    w = np.zeros(len(q), dtype=np.int16)
    for mz in mon:
        w += par(q & mz)
    wts[s:s + len(q)] = w
W = int(wts[odd].min())
tight = np.nonzero(odd & (wts == W))[0]
M = len(dirs)
adj = [[False] * M for _ in range(M)]; live = [False] * M
for qv in tight.astype(np.uint32):
    D = [i for i in range(M) if par(np.uint32(qv) & mon[i])]
    for a in D:
        live[a] = True
        for b in D:
            if a != b: adj[a][b] = True
nodes = [i for i in range(M) if live[i]]
cl = max_clique(adj, nodes)
sav = np.log(m * W / M) / np.log(m)
print(f"m={m} family={which}: M={M} directions, deletion bound W={W} (moves/array {M/W:.3f}, saving if it existed {sav:.3f}); "
      f"tight odd forms {len(tight)}, live {len(nodes)}, max clique {len(cl)} -> {'EXCLUDED' if len(cl) > W else 'passes'}", flush=True)
