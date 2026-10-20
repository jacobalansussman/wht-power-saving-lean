#!/usr/bin/env python3
"""v-fp-measure: OpenAI network, own derivation of the raw packing P (and the move count M) from the preprint's
rows (net:schedule), its chronological stage-2 reading, and its frame table (net:labels), with fp-measure's class
model (X_d: own twist; Y, centre, side: twist 0, one class).  Per invocation per stage, every wire's label count on
every edge between consecutive touches; frame dimensions from the table (a = h^(j-1), dim E = a, dim B = a-1,
P a line, Q a tensor of lines).  Compared with fp-measure's closed form for P and with M = Wm - 2N + 2L."""
from math import comb

def run(h):
    v = comb(h, 3); d = comb(h - 3, 3) + 3 * (h - 3); m = h ** 3
    cj = comb(h - 1, 2)                     # triples containing a given coordinate
    P = M = 0
    Pst = []
    for j in (1, 2, 3):
        a = h ** (j - 1)
        fr = {"0": a - 1, "1": (a - 1) * h, "2": (a - 1) * h + 1, "3": a * h, "4": (a - 1) * h,
              "5": (a - 1) * h + h - 1, "6": a * h, "7": a * h,
              "Xin": a, "Yin": a - 1, "Xout": a * h, "Yout": a * h - 1, "src": 0, "sink": m}
        def wire(seq, labels):
            """seq: frame keys in time order; labels[i] = label count on edge seq[i] -> seq[i+1]"""
            p = mm = 0
            for (u, w), n in zip(zip(seq, seq[1:]), labels):
                k = abs(fr[w] - fr[u])
                mm += k; p += k * (n - 1)
            return p, mm
        if j == 1:
            X = wire(["Xin", "2", "3", "6", "7", "Xout"], [1] * 5)
            Y = wire(["Yin", "0", "1", "4", "5", "Yout"], [1, 1, 1, d + 2, 2])
            Cj = wire(["src", "1", "3", "4", "6", "sink"], [1, 1, 1 + cj, 1 + cj, 1])
            Cs = wire(["src", "1", "3", "4", "6", "sink"], [1, 1, 1 + v, 1 + v, 1])
            Sd = wire(["src", "0", "2", "5", "7", "sink"], [1, 1, 2, 2, 1])
        elif j == 2:
            X = wire(["Xin", "2", "3", "6", "7", "Xout"], [1, d + 2, 1, 1, 1])
            Y = wire(["Yin", "0", "1", "4", "5", "Yout"], [2] * 5)
            Cj = wire(["src", "1", "3", "4", "6", "sink"], [1, 1 + cj, 1 + cj, 1, 1])
            Cs = wire(["src", "1", "3", "4", "6", "sink"], [1, 1 + v, 1 + v, 1, 1])
            Sd = wire(["src", "0", "2", "5", "7", "sink"], [1, 2, 2, 1, 1])
        else:
            X = wire(["Xin", "2", "3", "6", "7", "Xout"], [1] * 5)
            Y = wire(["Yin", "0", "1", "4", "5", "Yout"], [2] * 5)
            Cj = wire(["src", "1", "3", "4", "6", "sink"], [1] * 5)
            Cs = Cj
            Sd = wire(["src", "0", "2", "5", "7", "sink"], [1] * 5)
        # X, Y wire edges Xout / Yout of stage j are the next stage's Xin / Yin (zero-length joins); terminal
        # source/sink edges of X, Y: X source frame = U_a (dim 1 = Xin at stage 1), sinks: X full, Y = U_a^perp
        inv = v * v
        p = inv * (v * X[0] + v * Y[0] + h * Cj[0] + Cs[0] + v * d * Sd[0])
        mv = inv * (v * X[1] + v * Y[1] + h * Cj[1] + Cs[1] + v * d * Sd[1])
        P += p; M += mv; Pst.append(p)
    # terminal climbs of the banks: X from Xout(3) = h^3 already full; Y from Yout(3) = h^3 - 1 = U_a^perp: 0
    return P, M, Pst

def closed(h):
    v = comb(h, 3); d = comb(h - 3, 3) + 3 * (h - 3)
    P1 = v ** 3 * ((h - 1) * (2 * d + 1) + 8 * h)
    P2 = v ** 3 * ((h - 1) * (d + 1) + h * (h - 1) + 8 * h + d * h * (h - 1))
    P3 = v ** 3 * h * h * (h - 1)
    N = v ** 3; I = 3 * v * v; L = I * (h + 1) * h; W = 2 * N + I * (v * d + h + 1)
    return [P1, P2, P3], W * h ** 3 - 2 * N + 2 * L

for h in (6, 7, 8, 9, 25, 100):
    P, M, Pst = run(h)
    Pc, Mc = closed(h)
    print("h=%d  P(own)=%d  P(fp-measure)=%d  equal %s  per stage equal %s | M(own)=%d  Wm-2N+2L=%d  equal %s  P/M=%.5f"
          % (h, P, sum(Pc), P == sum(Pc), Pst == Pc, M, Mc, M == Mc, P / M))
