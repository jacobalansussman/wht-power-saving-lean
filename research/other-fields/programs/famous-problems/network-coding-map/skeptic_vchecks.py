#!/usr/bin/env python3
"""v-iso-network skeptic checks (14 Oct 2026). Standard library only; seconds; < 50 MB.

V1. Row 1 for the cyclic DFT over a finite field: F diag(w^{ls}/N) F e_j = e_{-(j+s) mod N} over F_p,
    p = 97 (N | 96 for N = 2,4,8,16,32,3,6,12,24,48), all j, s.
V2. Row 8 specialised to a finite field: H^{(x)k} = ((1+i)/2)^{-k} D C^{(x)k} D over F_p, p = 10009
    (p = 1 mod 4, so i exists), k = 1..6. (Shows the WHT <-> C^{(x)k} identity survives reduction mod p,
    which the finite-field step of row 1 needs for circuits with Gaussian-dyadic constants.)
V3. Row 1 constant: the graph-independent threshold. For a graph of max degree 3, at most 1 + 3(2^r - 1)
    vertices lie within distance r of a vertex, so the average distance from an input to N distinct
    outputs is at least g(N) (computed exactly from that count). Averaging over shifts gives some s with
    sum_j dist >= N g(N). The relay graph of a circuit with |C| gates has <= 8|C| + 2N edges. For the plain
    radix-2 WHT (|C| = N log2 N) the guaranteed ratio is (8|C| + 2N) / (N g(N)); compare with the
    graph-specific ratio ~2 reported in iso-network check C. Also: the log2 N needed for
    (log N)^s to beat a constant factor K at s = 1.68e-3, for K = 2 and K = the guaranteed ratio.
V4. Row 3 cross-check: number of Lagrangian subspaces of F_2^{2m} (= vertices of the dual polar graph
    DSp(2m,2)) is prod_{i=1..m} (2^i + 1): 15 for m = 2, 135 for m = 3, matching the frame-graph sizes
    in proveout NOTE (15 and 135 nodes for label dimension 2 and 3).
"""
import math

def check_V1(p=97):
    out = []
    g = next(x for x in range(2, p) if all(pow(x, (p - 1) // q, p) != 1 for q in (2, 3)))
    for N in (2, 3, 4, 6, 8, 12, 16, 24, 32, 48):
        assert (p - 1) % N == 0
        w = pow(g, (p - 1) // N, p)
        Ninv = pow(N, p - 2, p)
        F = [[pow(w, (a * b) % N, p) for b in range(N)] for a in range(N)]
        ok = True
        for s in range(N):
            d = [(pow(w, (l * s) % N, p) * Ninv) % p for l in range(N)]
            for j in range(N):
                v = [F[l][j] for l in range(N)]               # F e_j
                v = [(v[l] * d[l]) % p for l in range(N)]
                y = [sum(F[a][l] * v[l] for l in range(N)) % p for a in range(N)]
                want = [1 if a == (-(j + s)) % N else 0 for a in range(N)]
                if y != want:
                    ok = False
        out.append((N, ok))
    return out

def check_V2(p=10009, kmax=6):
    assert p % 4 == 1
    i = next(x for x in range(2, p) if pow(x, 2, p) == p - 1)
    inv2 = pow(2, p - 2, p)
    cp, cm = ((1 + i) * inv2) % p, ((1 - i) * inv2) % p
    cpinv = pow(cp, p - 2, p)
    res = []
    for k in range(1, kmax + 1):
        N = 1 << k
        ok = True
        for a in range(N):
            for b in range(N):
                d = bin(a ^ b).count("1")
                Cab = (pow(cp, k - d, p) * pow(cm, d, p)) % p
                rhs = (pow(cpinv, k, p) * pow(i, bin(a).count("1"), p) * Cab * pow(i, bin(b).count("1"), p)) % p
                lhs = 1 if bin(a & b).count("1") % 2 == 0 else p - 1
                if lhs != rhs:
                    ok = False
        res.append((k, ok))
    return res

def g_avg_min(N, deg=3):
    """Minimum possible average distance from a vertex to N distinct other-or-same vertices
    in a graph of max degree deg (ball sizes 1, 1+deg, 1+deg+deg(deg-1), ...)."""
    remaining, r, total, ball_prev = N, 0, 0, 0
    while remaining > 0:
        ball = 1 if r == 0 else 1 + deg * ((deg - 1) ** r - 1) // (deg - 2)
        shell = ball - ball_prev
        take = min(shell, remaining)
        total += take * r
        remaining -= take
        ball_prev = ball
        r += 1
    return total / N

def check_V3():
    rows = []
    for k in (8, 16, 32, 64, 128, 256, 1024):
        N = 2 ** k if k <= 64 else None
        if N is None:
            # asymptotic: g(N) = log2 N - c with c -> log2(3) + ~1 ; use the exact formula in log space
            g = k - (math.log2(3) + 1) + 1e-9  # approximation, documented below
        else:
            g = g_avg_min(N)
        C = k  # gates per coordinate: |C| = N log2 N for the plain WHT
        ratio = (8 * C + 2) / g
        rows.append((k, round(g, 3), round(ratio, 3)))
    return rows

def needed_log2N(K, s):
    # (log2 N)^s >= K  <=>  log2 N >= K^(1/s);  report log2(log2 N)
    return math.log2(K) / s

def lagrangians(m):
    return math.prod(2 ** i + 1 for i in range(1, m + 1))

if __name__ == "__main__":
    print("V1 F diag(w^{ls}/N) F = negated shift over F_97:", check_V1())
    print("V2 H = ((1+i)/2)^-k D C^k D over F_10009:", check_V2())
    print("V3 plain WHT, relay graph: k, guaranteed avg distance g(N), guaranteed ratio (8|C|+2N)/(N g):")
    for row in check_V3():
        print("    ", row)
    s = 1.68e-3
    for K in (2.0, 8.0):
        print("V3 s=%.2e, factor K=%g: need log2 N >= 2^%.0f" % (s, K, needed_log2N(K, s)))
    print("V4 Lagrangians of F_2^{2m}:", [(m, lagrangians(m)) for m in (1, 2, 3, 4, 5)])
