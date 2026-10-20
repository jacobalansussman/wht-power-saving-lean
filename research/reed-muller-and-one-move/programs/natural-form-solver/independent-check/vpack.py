import itertools
from vfield import Model
m = 5; M = Model(m); N = 31
DM = {z: M.dmask(z) for z in range(1, N + 1)}
bucket = {}
for S in itertools.combinations(range(1, N + 1), m):
    v = 0
    for z in S: v ^= DM[z]
    bucket.setdefault(v, []).append(frozenset(S))
fr = {k: bucket[M.smask(M.target(k))] for k in range(N)}
def packing(ks):
    W = len(ks); best = 0
    # exhaustive over all 7^W assignments (frame or none) with disjointness, by recursion
    def rec(i, used, cnt):
        nonlocal best
        if cnt + (W - i) <= best: return
        if i == W: best = max(best, cnt); return
        for F in fr[ks[i]]:
            if not (F & used): rec(i + 1, used | F, cnt + 1)
        rec(i + 1, used, cnt)
    rec(0, frozenset(), 0); return best
sets = {"rm-theory 8": (0,2,4,8,18,24,25,27), "additive 8": (0,0,1,9,16,18,21,25), "Frobenius 8": (0,0,0,1,2,4,8,16),
        "N1": (0,8,9,10,14,17,23,28), "N2": (0,4,6,11,14,21,22,23), "rm-theory 9": (0,1,2,5,11,26,27,29,30),
        "C1uC3 10": (1,2,4,8,16,3,6,12,17,24), "slice 8": (0,)*8, "slice 7": (0,)*7}
for name, ks in sets.items():
    print(name, ks, "packing", packing(ks), flush=True)
# pairwise: do any two frames of I intersect?
F0 = fr[0]
print("slice: min pairwise intersection of the 6 frames of I:", min(len(a & b) for a, b in itertools.combinations(F0, 2)))
# additive set: check that alpha^(2k) for k in (1,9,16,18,21,25) and 0 ... form subspace
F = M.F
els = [F.e[(2*k) % 31] for k in (0,1,9,16,18,21,25)]
print("additive elements", els, "closed under xor:", all((a ^ b) in els or a == b for a in els for b in els))
print("T0+T1+T18 elements:", F.e[0] ^ F.e[2] ^ F.e[36 % 31])
