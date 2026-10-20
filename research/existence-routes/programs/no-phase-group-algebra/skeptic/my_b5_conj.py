#!/usr/bin/env python3
"""my_b5_conj.py (v-e3-rmbuild): block 5, W = 8, M = 31, H = a random GL(5,2)-conjugate of the Singer cycle
(x^5+x^2+1).  With H transitive, |D_Q| = 31 |supp Q| and deletion never fails (|supp Q| >= 8 for Q != 0).  A form is
tight iff |supp Q| = 8 (Q = l_a l_b, two independent linear forms) and sum_j Q(h^-1 e_j) = 1 for all 31 h in H.
Counts the tight forms for many conjugates (the clique test can only bite when tight forms exist)."""
import random, sys
m = 5; N = 32; poly = 0b100101
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 1)


def mulx(v):
    v <<= 1
    if v >> m & 1:
        v ^= poly
    return v


sig = [mulx(v) for v in range(N)]


def lin(cols):
    img = []
    for x in range(N):
        y = 0
        for k in range(m):
            if x >> k & 1:
                y ^= cols[k]
        img.append(y)
    return img


def dot(a, b):
    return bin(a & b).count("1") & 1


ntrials = int(sys.argv[1]) if len(sys.argv) > 1 else 200
hist = {}
for trial in range(ntrials):
    while True:
        g = lin([random.randrange(1, N) for _ in range(m)])
        if len(set(g)) == N:
            break
    gi = [0] * N
    for x in range(N):
        gi[g[x]] = x
    gen = [g[sig[gi[x]]] for x in range(N)] if trial else sig
    H = [list(range(N))]
    for _ in range(30):
        H.append([gen[v] for v in H[-1]])
    Hinv = []
    for h in H:
        inv = [0] * N
        for x in range(N):
            inv[h[x]] = x
        Hinv.append(inv)
    tight = 0
    for a in range(1, N):
        for b in range(a + 1, N):
            ok = True
            for k in range(31):
                if sum(dot(a, Hinv[k][1 << j]) & dot(b, Hinv[k][1 << j]) for j in range(m)) % 2 == 0:
                    ok = False; break
            if ok:
                tight += 1
    hist[tight] = hist.get(tight, 0) + 1
print("block 5, 8 in 31, %d Singer conjugates (first = e3-rmbuild's own): tight-form counts -> multiplicity %s"
      % (ntrials, sorted(hist.items())))
