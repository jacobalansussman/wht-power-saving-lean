#!/usr/bin/env python3
"""vclique.py (e3-rmbuild): exact deletion and tight-clique tests for the no-phase game with a cyclic address group.

By Theorem A and Lemma C, a unit "W in M" of the game with address group Z_n = <sigma> (and shifts) gives, at every
binary quadratic form Q, the channel equation  Atil Pi_Q (I + Gamma Pi_Q)^-1 Btil = O_Q (x) I_W  (operators on
C^n x registers).  Read with VIRTUAL moves (t, k) (physical move t seen in channel k, direction sigma^-k z_t) and
VIRTUAL arrays (l, k), this is a commutative-slice system with Mn moves and Wn arrays, where virtual array (l, k) is
odd at Q iff tr(Q o sigma^-k) = 1.  Two exact consequences (ex2-families Lemma 2 and Theorem 1, re-derived for
projector targets in PROOF.md, Lemma V):
  deletion:  |D_Q| >= W * n_odd(Q) for every Q, where |D_Q| = #{(t,k): Q(sigma^-k z_t) = 1};
  clique:    call Q tight if equality holds; virtual moves lying in a common tight support are joined; a set of
             pairwise joined virtual moves, each in some tight support, has size <= W n  (rank of Btil Atil).
Usage: vclique.py m n power W dirs(comma)."""
import sys, itertools

class tf:      # the pure helpers of tf.py (copied, so that no numpy/torch is needed)
    PRIM = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011}
    @staticmethod
    def mulx(v, m):
        v <<= 1
        if (v >> m) & 1:
            v ^= tf.PRIM[m]
        return v
    @staticmethod
    def singer_table(m, power=1):
        t = list(range(1 << m))
        for _ in range(power):
            t = [tf.mulx(v, m) for v in t]
        return t
    @staticmethod
    def order_of(tab):
        N = len(tab); cur = list(range(N)); k = 0
        while True:
            cur = [tab[v] for v in cur]; k += 1
            if cur == list(range(N)):
                return k
    @staticmethod
    def binary_forms(m):
        pairs = [(j, k) for j in range(m) for k in range(j, m)]
        out = []
        for mask in range(1 << len(pairs)):
            row = []
            for z in range(1 << m):
                s_ = 0
                for b, (j, k) in enumerate(pairs):
                    if (mask >> b) & 1 and (z >> j) & 1 and (z >> k) & 1:
                        s_ ^= 1
                row.append(2 * s_)
            out.append(tuple(row))
        return out

m, n, power, W = (int(x) for x in sys.argv[1:5])
dirs = [int(d) for d in sys.argv[5].split(",")]
N = 1 << m
sig = tf.singer_table(m, power) if n > 1 else list(range(N))
if n > 1:
    assert tf.order_of(sig) == n
inv = [0] * N
for v in range(N):
    inv[sig[v]] = v
sk = [list(range(N))]
for k in range(1, n):
    sk.append([inv[v] for v in sk[-1]])
forms = [tuple(x // 2 for x in t) for t in tf.binary_forms(m)]       # 0/1 tables
vm = [(t, k, sk[k][z]) for t, z in enumerate(dirs) for k in range(n)]  # virtual moves with their directions
bad_del = []; tight = []
for qi, Q in enumerate(forms):
    D = [i for i, (t, k, w) in enumerate(vm) if Q[w]]
    nodd = sum(1 for k in range(n) if sum(Q[sk[k][1 << j]] for j in range(m)) % 2)
    if len(D) < W * nodd:
        bad_del.append((qi, len(D), W * nodd))
    elif len(D) == W * nodd and nodd > 0:
        tight.append(set(D))
print("m=%d n=%d W=%d M=%d: virtual moves %d, virtual arrays %d" % (m, n, W, len(dirs), len(vm), W * n))
print("deletion failures:", len(bad_del), bad_del[:5])
print("tight forms:", len(tight))
# joined graph on virtual moves
V = len(vm)
live = set().union(*tight) if tight else set()
adj = [set() for _ in range(V)]
for S in tight:
    for a in S:
        adj[a] |= S
# greedy + exact clique search on the direction level (virtual moves of the same direction are interchangeable)
bydir = {}
for i, (t, k, w) in enumerate(vm):
    if i in live:
        bydir.setdefault(w, []).append(i)
dlist = sorted(bydir)
def joined(w1, w2):
    a = bydir[w1][0]; b = bydir[w2][0]
    return b in adj[a]
selfok = {w: (len(bydir[w]) == 1 or bydir[w][1] in adj[bydir[w][0]]) for w in dlist}
best = []
def extend(cl, cand):
    global best
    wt = sum(len(bydir[w]) if selfok[w] else 1 for w in cl)
    if wt > sum(len(bydir[w]) if selfok[w] else 1 for w in best):
        best = list(cl)
    for i, w in enumerate(cand):
        extend(cl + [w], [u for u in cand[i + 1:] if joined(w, u)])
extend([], dlist)
size = sum(len(bydir[w]) if selfok[w] else 1 for w in best)
print("live directions:", len(dlist), "best clique (directions):", best, "virtual size", size, "bound W n =", W * n)
print("VERDICT:", "EXCLUDED (deletion)" if bad_del else ("EXCLUDED (tight clique %d > %d)" % (size, W * n) if size > W * n else "not excluded by these two tests"))
