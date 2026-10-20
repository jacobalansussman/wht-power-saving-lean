"""t0_control.py -- control: rebuild the dual polar graph and reproduce v-w3-blueprint's exhaustive counts for single-term
middle shared moves fed by a distance-d term (t-frame: t -> L_Z, inputs -> M at distance d, kernel -> N' complement of M).
Then add the two constraints the skeptic's count did not impose:
  (i)  t must be exactly r = 4 - d moves from the kernel (d(N', L_Z) = 4 - d, forced since d(in, K) = 4), and
  (ii) the shared move is an exact orbit sum: its cell (t-frame) is an H-fixed neighbour N of L_Z with d(N, M) = d + 1
       and d(N, N') = 3 - d (the shared move progresses towards the kernel; d = 2: N adjacent to N').
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp8 import *
from collections import Counter

t0 = time.time()
P = Polar()
print('Lagrangians', len(P.L), 'built in %.1fs' % (time.time() - t0), flush=True)
iZ = P.idx[LZ]
Kc = sigma(G.K)
i1 = P.idx[img_mask(table(inv(Kc)), LZ)]
print('kernel cell distance from inputs:', P.dZ[i1], P.d(iZ, i1))
GR = groups()
D = P.D
for name, (c, cp) in GR.items():
    W = len(c)
    pc = [P.perm(x) for x in c]; pci = [P.perm(inv(x)) for x in c]; pcp = [P.perm(x) for x in cp]
    assert all(p[iZ] == iZ and p[i1] == i1 for p in pc), 'c not in Levi'
    assert all(p[iZ] == iZ and p[i1] == i1 for p in pcp), "c' not in Levi"
    fix = [j for j in P.nb[iZ] if all(p[j] == j for p in pc)]
    print(name, 'H-fixed neighbours of L_Z (exact orbit-sum directions):', fix, flush=True)
    for d in (1, 2, 3):
        Ms = [i for i in range(len(P.L)) if P.dZ[i] == d]
        cnt = Counter(); cex = 0; cexA = 0; cexB = 0; nA = 0; nB = 0
        for Mi in Ms:
            comp = np.nonzero(D[Mi] == 4)[0]
            for Ni in comp:
                CN = [pci[g][Ni] for g in range(W)]
                S = [g for g in range(W) if CN[g] == Ni]
                junk = [g for g in range(W) if g not in S]
                full = [g for g in junk if D[CN[g], Mi] == 4]
                cnt[(len(S), len(full))] += 1
                esc = bool(junk) and not full
                cex += esc
                okA = P.dZ[Ni] == 4 - d
                okB = okA and any(D[N, Mi] == d + 1 and D[N, Ni] == 3 - d for N in fix)
                nA += okA; nB += okB
                cexA += esc and okA; cexB += esc and okB
        print(f'd={d} {name}: pairs {sum(cnt.values())}; (|S|,#full junk) {dict(sorted(cnt.items()))}; escaping {cex};'
              f' with (i): {nA} pairs, escaping {cexA}; with (i)+(ii): {nB} pairs, escaping {cexB}', flush=True)
print('done %.1fs' % (time.time() - t0))
