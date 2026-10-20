"""t2_d3.py -- the skeptic's d = 3 configurations (a shared move F fed by a distance-3 term t, then ONE move of the top
orbit to the kernel; e.g. the second shared move of F O O F O) with the constraints his count left out:
  (i')  t within 2 moves of the kernel: d(N', L_Z) in {1, 2};
  (ii') F's own cell N (t-frame; a neighbour of L_Z) adjacent to the kernel N' (the top orbit's single move reaches it);
        for an exact (or psi-twisted) orbit sum N is H-fixed.
Reports |S|, number of junk cells at full input distance, and the S = {0} escapes (no full-distance junk)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp8 import *
from collections import Counter
P = Polar(); D = P.D; iZ = P.idx[LZ]
GR = groups()
for name, (c, cp) in GR.items():
    W = len(c)
    pc = [P.perm(x) for x in c]; pci = [P.perm(inv(x)) for x in c]
    fix = [j for j in P.nb[iZ] if all(p[j] == j for p in pc)]
    for d in (2, 3):
        Ms = [i for i in range(len(P.L)) if P.dZ[i] == d]
        cntA = Counter(); cntB = Counter(); exA = []; exB = []
        for Mi in Ms:
            for Ni in np.nonzero(D[Mi] == 4)[0]:
                if P.dZ[Ni] not in ((2,) if d == 2 else (1, 2)):
                    continue
                CN = [pci[g][Ni] for g in range(W)]
                S = [g for g in range(W) if CN[g] == Ni]
                junk = [g for g in range(W) if g not in S]
                full = [g for g in junk if D[CN[g], Mi] == 4]
                key = (int(P.dZ[Ni]), len(S), len(full))
                cntA[key] += 1
                if junk and not full: exA.append((Mi, int(Ni)))
                if any(D[N, Ni] == 1 for N in fix):
                    cntB[key] += 1
                    if junk and not full: exB.append((Mi, int(Ni)))
        print(f'{name} d={d}: with (i\'): (d(t,K), |S|, #full junk): {dict(sorted(cntA.items()))}  escapes {len(exA)}')
        print(f'{name} d={d}: with (i\')+(ii\'): {dict(sorted(cntB.items()))}  escapes {len(exB)}', flush=True)
