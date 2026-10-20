"""paircoin.py -- the multi-term loophole, first look: for a shared move with two content terms t1, t2 (both O1 outputs,
t = D1 a), anchored so that t1's branch 0 reaches the kernel, branch (h, i) sits on the cell [K t1^-1 c_h t_i].
Count, over random pairs: cross cells equal to the kernel [K] (a second route into the kernel from the other term),
and cross cells at full input distance that coincide with t1's own full-distance junk cells (possible cancellation)."""
import sys
from symsolve import *
from collections import Counter
W, H, n, mode = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), sys.argv[4]
R = Run(4, W, H, 1); G = R.G; rng = random.Random(3)
Kc = np.conj(G.K).T
def rs(M):
    return int(round(np.log2(G.rowsupp(M).max())))
def term():
    a = rand_free_b(G, rng, 4)
    return G.mul(G.D1, a) if mode == 'mid' else a
cnt = Counter()
for _ in range(n):
    t1 = term(); t2 = term()
    if mode == 'midrel':
        # second term = the first one with an extra free letter in front (same O1 output, two letters in F's content)
        t2 = G.mul(rand_free_b(G, rng, 4), t1)
    t1i = np.conj(t1).T
    own = {}
    for h in R.Hc[1:]:
        M = G.mul(G.K, t1i, h[0], t1)
        if not G.is_monomial(G.mul(M, Kc)) and rs(M) == 4:
            own[G.frame_key(M)] = 1
    kern = 0; coin = 0
    for h in R.Hc:
        M = G.mul(G.K, t1i, h[0], t2)
        if G.is_monomial(G.mul(M, Kc)):
            kern += 1
        elif rs(M) == 4 and G.frame_key(M) in own:
            coin += 1
    cnt[(kern, coin)] += 1
print(H, mode, 'pairs:', n, ' (cross terms reaching the kernel, cross terms on own full-distance junk):', sorted(cnt.items()), flush=True)
