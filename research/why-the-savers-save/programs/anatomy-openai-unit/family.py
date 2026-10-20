# an-openai: the OpenAI exchange network as a family in its own parameters (exact rational arithmetic).
# Reproduces the paper's integers at h=100 and scans the ground-set size h (triples, 3 stages).
from fractions import Fraction as Fr
from math import comb, log
import sys

def complex_net(h, k=3, t=3):
    """OpenAI complex network (Integer-mult Sec. 3 / explicit-DFT Sec. 2) with k-subsets of [h],
    t stages (t=3 in the paper; the signed exchange needs 3 shears). Neighbours: distinct k-sets with
    even intersection (k odd so that t_T.t_T = 1). Centre wires per invocation: h + 1."""
    v = comb(h, k)
    # neighbours: distinct sets with even intersection size
    d = sum(comb(k, i) * comb(h - k, k - i) for i in range(0, k, 2))  # i even, i<k
    N = v ** t
    m = h ** t
    I = t * v ** (t - 1)
    c = h + 1
    W = 2 * N + I * (v * d + c)
    L = I * c * h
    s = W * m - 2 * N + 2 * L
    Delta = W * m - s
    return dict(h=h, k=k, v=v, d=d, N=N, m=m, I=I, W=W, L=L, s=s, Delta=Delta)

def bit_net(h, t=3):
    v = comb(h, 3)
    d = 3 * comb(h - 3, 2)   # intersection exactly one
    N = v ** t; m = h ** t; I = t * v ** (t - 1); c = h
    W = 2 * N + I * (v * d + c)
    L = I * c * h
    s = W * m - N + 2 * L
    return dict(h=h, v=v, d=d, N=N, m=m, W=W, L=L, s=s, Delta=W * m - s)

def saving(W, m, M):
    return log(m * W / M) / log(m)

if __name__ == "__main__":
    r = complex_net(100)
    assert r["d"] == 147731, r["d"]
    assert r["W"] == 1873807244643542670000, r["W"]
    assert r["s"] == 1873807244636671267308000000, r["s"]
    assert r["Delta"] == 6871402692000000
    eta = Fr(r["Delta"], r["W"] * r["m"])
    assert eta == Fr(73, 19906842167500), eta
    print("complex h=100 reproduces paper: W, s, Delta, eta =", r["W"], r["s"], r["Delta"], eta)
    W_star = 2 ** 71
    lam = Fr(r["m"]) - Fr(r["Delta"], W_star)
    th = log(float(lam)) / log(r["m"])
    print("  padded W*=2^71: 1-theta = %.6e ; unpadded saving ln(mW/M)/ln m = %.6e" %
          (1 - th, -log(1 - float(eta)) / log(r["m"])))
    b = bit_net(100)
    assert b["W"] == 177176569091445000000 and b["s"] == 177176569088785861287000000
    print("bit h=100 reproduces paper: eta_b =", Fr(b["Delta"], b["W"] * b["m"]))
    print()
    print("scan of h (complex, 3-subsets, 3 stages): saving s = ln(mW/M)/ln m")
    best = None
    for h in range(20, 61):
        r = complex_net(h)
        if r["Delta"] <= 0:
            print(" h=%d: Delta=%d (no saving)" % (h, r["Delta"])); continue
        sv = -log(1 - r["Delta"] / (r["W"] * r["m"])) / log(r["m"])
        if best is None or sv > best[1]: best = (h, sv)
        if h <= 32 or h % 10 == 0:
            print(" h=%d m=%d W=%.4e Delta=%.4e Delta/N=%.4f W/N=%.1f saving=%.4e" %
                  (h, r["m"], r["W"], r["Delta"], r["Delta"] / r["N"], r["W"] / r["N"], sv))
    print("best h in scan:", best)
