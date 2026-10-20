# Independent recount of Lemma C's table (units of |B|): P = slots, E = W^2 L, orbit o, gauge, dominance threshold.
from fractions import Fraction as Fr
from math import prod, sqrt
def gb(m, k): return prod(2**(m - i) - 1 for i in range(k)) // prod(2**(i + 1) - 1 for i in range(k))
def n_(m, k): return 2**(k*(k+1)//2) * gb(m, k)
def L(m): return prod(2**i + 1 for i in range(1, m + 1))
def slots(W, n): return sum(W + j for j in range(n)) + W * (W + n)        # step j sees W inputs + j earlier outputs
rows = [(2,1,[1,1],[2]),(2,2,[1]*4,[2,2]),(2,1,[1],[2]),(2,2,[1]*3,[2,2]),(2,3,[1]*5,[2]*3),(3,1,[1]*3,[3]),(3,2,[1]*6,[3,3]),
        (3,2,[1]*5,[3,3]),(3,3,[1]*8,[3]*3),(3,5,[1]*14,[3]*5),(3,5,[1]*10,[3,2,2,2,2]),(4,1,[1]*4,[4]),(4,2,[1]*8,[4,4]),
        (4,2,[1]*7,[4,4]),(4,2,[2,2,2,1],[4,4]),(4,3,[1]*10,[4]*3)]
for m, W, st, tg in rows:
    P = slots(W, len(st)); E = W * W * L(m)
    o = sum((2 - Fr(1, n_(m, tg[t]))) if tg[t] == tg[s] else 2 for t in range(W) for s in range(W))
    gauge = sum(Fr(1, n_(m, r)) for r in st)
    print(f'm={m} W={W} steps={st} targets={tg}: P={P} E={E} P-E={P-E} o={float(o):.3f} P-gauge={float(P-gauge):.3f}')
for m in (2, 3, 4):
    print(f'm={m}: dominance threshold c = sqrt(2+2L)-2 = {sqrt(2 + 2*L(m)) - 2:.3f};  exact smallest c with slots(W,cW) >= W^2 L for W=1..5:',
          [next(n for n in range(1, 10**4) if slots(W, n) >= W*W*L(m)) / W for W in range(1, 6)])
