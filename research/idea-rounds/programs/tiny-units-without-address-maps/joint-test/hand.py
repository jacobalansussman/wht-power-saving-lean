#!/usr/bin/python3
"""joint/cert/hand.py -- agent "cert".  The six static inequalities of the hand proof, recomputed from model_cert.txt,
and an exhaustive check that they alone reject every multiset of K directions (plain sums)."""
import itertools
v = [list(map(int, l.split())) for l in open("model_cert.txt")]
cap = {w: v[15 * (w - 1):15 * w] for w in range(1, 8)}; dem = v[105:120]
B1 = list(range(8, 15))                       # block 7 + x = (x, bit 1)
def ineq(xs):
    U = [7 + x for x in xs]; out = [s for s in B1 if s not in U]
    return [sum(cap[w][s][d] for s in out for d in U) for w in range(1, 8)], sum(dem[s][d] for s in out for d in U)
pairs = [(1, 6), (2, 5), (4, 3)]; triples = [(5, 6, 7), (3, 6, 7), (3, 5, 7)]
six = []
for xs in pairs + triples:
    c, D = ineq(xs); six.append((c, D)); print("U = blocks (x,1) with x in %s: sum_w n_w * %s >= %d * W" % (list(xs), c, D))
print("sum of the 3 pair inequalities:   coefficients %s >= %d W" % ([sum(six[i][0][w] for i in range(3)) for w in range(7)], sum(six[i][1] for i in range(3))))
print("sum of the 3 triple inequalities: coefficients %s >= %d W" % ([sum(six[i][0][w] for i in range(3, 6)) for w in range(7)], sum(six[i][1] for i in range(3, 6))))
print("every capacity of the model is even: %s" % all(x % 2 == 0 for w in cap for r in cap[w] for x in r))
def least(W):
    K = 0
    while True:
        for ms in itertools.combinations_with_replacement(range(1, 8), K):
            n = [ms.count(w) for w in range(1, 8)]
            if all(sum(n[i] * c[i] for i in range(7)) >= W * D for c, D in six): return K, ms
        K += 1
for W in range(1, 10):
    K, ms = least(W)
    S = 0
    while 10 * S < 3 * (7 * W + (7 * W) % 2): S += 1          # pairs, with the parity step
    Kf = S
    while 12 * Kf - 2 * S < 3 * (9 * W + (9 * W) % 2): Kf += 1   # triples, with the parity step, c <= K - S
    print("W = %d: the six inequalities alone need K >= %2d (first multiset passing: %s) ; the formula of the hand proof gives K >= %2d ; 13W/5 = %.1f ; a saving of 1e-1 needs K <= %d" % (W, K, ms, Kf, 2.6 * W, int(2.6883 * W)))
