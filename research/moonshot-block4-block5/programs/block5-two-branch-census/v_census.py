# v_census.py -- skeptic's own code (v-b5-twobranch). Exact over F_2, dual polar graph of Sp(10,2).
# Tests the "three middle frames" census of b5-twobranch PROOF sec. 2 over ALL anchor pairs (P1, Q) of
# several random walks M0 - L0 - N1 (no transitivity assumed), and the gate argument that proves it:
#   M0's gate in the quad Q_I (I = N1 cap A1) is N1; Q ~ M0 => d(gate(Q), N1) <= 1.
# Lagrangian = 1024-bit mask over vectors of F_2^10 (x = bits 0-4, z = bits 5-9).
import random, sys, time
from collections import Counter

def pc(x): return bin(x).count('1')
def omega(u, v): return (pc((u & 31) & (v >> 5)) + pc((u >> 5) & (v & 31))) & 1

def span_mask(vecs):
    S = {0}
    for v in vecs:
        S |= {s ^ v for s in S}
    m = 0
    for s in S: m |= 1 << s
    return m

def elems(m):
    out = []; v = 0
    while m:
        if m & 1: out.append(v)
        m >>= 1; v += 1
    return out

def dist(a, b): return 5 - (pc(a & b).bit_length() - 1)

NBR = {}
def nbrs(L):
    if L in NBR: return NBR[L]
    Le = elems(L); res = set()
    for v in range(1024):
        if (L >> v) & 1: continue
        H = [h for h in Le if omega(v, h) == 0]
        m = 0
        for h in H: m |= (1 << h) | (1 << (h ^ v))
        res.add(m)
    res = list(res); assert len(res) == 62
    NBR[L] = res; return res

def isotropic_check(L):
    Le = elems(L); assert len(Le) == 32
    for a in Le[:8]:
        for b in Le:
            assert omega(a, b) == 0

def run(seed):
    rnd = random.Random(seed)
    L0 = span_mask([1 << (5 + i) for i in range(5)])          # Z-type Lagrangian
    M0 = rnd.choice(nbrs(L0))
    N1 = rnd.choice([N for N in nbrs(L0) if dist(N, M0) == 2])
    isotropic_check(M0); isotropic_check(N1)
    anchors = [P for P in nbrs(M0) if dist(P, L0) == 2 and dist(P, N1) == 3]
    nN1 = nbrs(N1)
    stats = Counter(); gate_viol = 0; ntarg = 0; dPQ = Counter()
    for P1 in anchors:
        branches = [N for N in nN1 if dist(P1, N) == 4]
        targets = set()
        for N in branches:
            for A in nbrs(N):
                if dist(P1, A) == 5: targets.add(A)
        targets = list(targets); ntarg += len(targets)
        commons = {}
        for A in targets:
            assert dist(N1, A) == 2
            C = [N for N in nN1 if dist(N, A) == 1]
            assert len(C) == 3, len(C)
            assert all(dist(P1, N) == 4 for N in C)
            # gate of M0 in the quad of I = N1 cap A1 must be N1: d(M0,X) = 2 + d(N1,X) for X in quad
            assert dist(M0, A) == 4
            commons[A] = C
        for Qa in anchors:
            if Qa == P1: continue
            dPQ[dist(P1, Qa)] += 1
            for A in targets:
                lam = dist(Qa, A)
                pat = tuple(sorted(dist(Qa, N) for N in commons[A]))
                stats[(dist(P1, Qa), lam, pat)] += 1
    return len(anchors), ntarg, dPQ, stats

if __name__ == '__main__':
    seeds = [int(s) for s in sys.argv[1:]] or [1]
    tot = Counter()
    for s in seeds:
        t = time.time()
        na, nt, dPQ, st = run(s)
        print(f"walk seed {s}: {na} anchors, {nt} (P1,target) pairs, d(P1,Q) counts {dict(dPQ)} ({time.time()-t:.0f}s)")
        for k in sorted(st): print("   (d(P1,Q), lam1, d(Q,.) of the 3 common nbrs of N1, A1):", k, st[k])
        for k, v in st.items(): tot[(k[1], k[2])] += v
        sys.stdout.flush()
    print("TOTAL by (lam1, pattern):")
    for k in sorted(tot): print("  ", k, tot[k])
