# quot5.py -- block 4 (2295 frames): (a) distance-profile contractions (small k), (b) sampled greedy merging from
# singletons (large k): how many merges before the contraction's routing value rho falls below 41/13?
import sys, random, json, time, collections
import numpy as np
import dpg, quot3

m = 4
t0 = time.time()
n, Dm, OPP, ADJ = quot3.setup(m)
npairs = int(OPP.sum()); nedges = int(ADJ.sum() // 2)
print(f"m=4: {n} frames, {nedges} edges, {npairs} ordered opposite pairs, built {time.time()-t0:.0f}s"); sys.stdout.flush()
best = quot3.profiles(m, n, Dm, OPP, ADJ, int(sys.argv[1]), 3, kmax=300)
print("  distance-profile contractions, best rho per k (k <= 40 shown):")
for k in sorted(best):
    if k <= 40: print(f"   k={k}: rho = {best[k][0]:.5f} (refs {best[k][1]}, largest parts {best[k][2][-5:]})")
sys.stdout.flush()
# (b) sampled greedy merging
rng = random.Random(5)
D = Dm.astype(np.int32); P = OPP.astype(np.int32); C = ADJ.astype(np.int32)
cross = nedges; num = int((P * D).sum())
alive = np.ones(n, dtype=bool)
curve = {n: num / npairs / (cross / nedges)}
steps = int(sys.argv[2]); cand = int(sys.argv[3])
edges = [(i, j) for i in range(n) for j in range(i + 1, n) if ADJ[i, j]]
for s in range(steps):
    bestc = None
    for _ in range(cand):
        while True:
            a, b = edges[rng.randrange(len(edges))]
            # map to current representatives
            a = rep[a] if s else a; b = rep[b] if s else b
            if a != b and C[a, b] > 0: break
        Dn = np.minimum(D, np.minimum(D[:, [a]] + D[[b], :], D[:, [b]] + D[[a], :]))
        numn = int((P * Dn).sum())
        crossn = cross - int(C[a, b])
        val = numn / npairs / (crossn / nedges)
        if bestc is None or val > bestc[0]: bestc = (val, a, b, Dn, numn, crossn)
    val, a, b, Dn, numn, crossn = bestc
    # merge b into a: rows/cols of b become copies, then neutralise b
    P[a, :] += P[b, :]; P[:, a] += P[:, b]; P[b, :] = 0; P[:, b] = 0
    C[a, :] += C[b, :]; C[:, a] += C[:, b]; C[b, :] = 0; C[:, b] = 0; C[a, a] = 0
    D = Dn; D[b, :] = D[a, :]; D[:, b] = D[:, a]
    num = int((P * D).sum()); cross = crossn
    if s == 0: rep = list(range(n))
    rep = [a if r == b else r for r in rep]
    alive[b] = False
    k = int(alive.sum())
    curve[k] = num / npairs / (cross / nedges)
    if s % 10 == 0 or curve[k] < 41 / 13 + 0.02:
        print(f"   merges {s+1}: k = {k}, rho = {curve[k]:.5f} ({time.time()-t0:.0f}s)"); sys.stdout.flush()
    if curve[k] < 3.0: break
json.dump(dict(profiles={k: best[k][0] for k in best}, merge_curve=curve), open("quot5_m4.json", "w"))
