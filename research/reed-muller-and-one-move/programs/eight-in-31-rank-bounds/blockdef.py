"""blockdef.py (r3-classes): block-diagonal deficiency test.  If moves X split into parts P_1..P_t with every cross
cell E[P_i,P_j] (i != j) a derived zero, then rank E[X,X] = sum rank E[P_i,P_i] <= sum u(P_i), while
L[X,X] = E[X,X] + K[X,X] is unit triangular (rank |X|) and rank K <= W.  So sum_i (|P_i| - u(P_i)) <= W.
Parts used: single moves x with E_xx = 0 (u = 0) and subsets P of an odd form's support D_G (u = s_G).
This file: exploration on one target set (zeros from a linkrec/starrec Z line)."""
import struct, sys, itertools
raw = open(sys.argv[1], 'rb').read(); nf = struct.unpack_from('<i', raw, 0)[0]
mask = struct.unpack_from('<%dI' % nf, raw, 4); wt = struct.unpack_from('<%di' % nf, raw, 4 + 4 * nf)
vid = struct.unpack_from('<%di' % nf, raw, 4 + 8 * nf); o = 4 + 12 * nf
nv = struct.unpack_from('<i', raw, o)[0]; pk = struct.unpack_from('<%dI' % nv, raw, o + 4)
ks = [int(x) for x in sys.argv[2].split(',')]; W = len(ks)
Z = [0] + [int(h, 16) for h in sys.argv[3].split()]
w = [sum(pk[v] >> k & 1 for k in ks) for v in range(nv)]
ALL = (1 << 31) - 1
# symmetric zero adjacency: adj[x] = moves y != x with E_xy = E_yx = 0
adj = [0] * 32
for x in range(1, 32):
    for y in range(1, 32):
        if y != x and (Z[x] >> (y - 1) & 1) and (Z[y] >> (x - 1) & 1): adj[x] |= 1 << (y - 1)
diag0 = sum(1 << (x - 1) for x in range(1, 32) if Z[x] >> (x - 1) & 1)
def popc(m): return bin(m).count('1')
def maxclique(P):
    best = [0]
    def rec(R, P, sz):
        if not P:
            best[0] = max(best[0], sz); return
        if sz + popc(P) <= best[0]: return
        while P:
            if sz + popc(P) <= best[0]: return
            v = (P & -P).bit_length(); P &= P - 1
            rec(R | 1 << (v - 1), P & adj[v], sz + 1)
    rec(0, P, 0); return best[0]
print('zero-diagonal moves', popc(diag0), 'zero pairs', sum(popc(a) for a in adj) // 2, 'clique', maxclique(diag0))
best = 0; cnt = 0
for q in range(1, nf):
    r = w[vid[q]]; s = wt[q] - r
    if r < 2: continue
    D = mask[q]
    # common zero-neighbourhood of all of D (outside D)
    N = ALL & ~D
    for x in range(1, 32):
        if D >> (x - 1) & 1: N &= adj[x]
    c = maxclique(N & diag0)
    d = r + c
    if d > best: best = d; print('form', q, 'wt', wt[q], 'r', r, 's', s, '| common zero-nbhd', popc(N), 'clique there', c, '-> deficiency', d)
    # drop t elements of D to enlarge the neighbourhood: P = D \ T, def = |P| - s
    for t in (1, 2):
        for T in itertools.combinations([x for x in range(1, 32) if D >> (x - 1) & 1], t):
            P = D
            for x in T: P &= ~(1 << (x - 1))
            N2 = ALL & ~P
            for x in range(1, 32):
                if P >> (x - 1) & 1: N2 &= adj[x]
            c2 = maxclique(N2 & diag0 & ~D | (N2 & diag0 & sum(1 << (x - 1) for x in T)))
            d2 = popc(P) - s + c2
            if d2 > best: best = d2; print('form', q, 'minus', T, 'r', r, 's', s, 'clique', c2, '-> deficiency', d2)
print('best deficiency', best, 'need >', W)
