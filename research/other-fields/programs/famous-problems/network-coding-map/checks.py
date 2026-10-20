#!/usr/bin/env python3
"""iso-network small checks (14 Oct 2026). Standard library only; seconds; < 100 MB.

A. Hadamard <-> C^{(x)k} <-> dyadic convolution ("Bluestein for (Z/2)^k"):
   (-1)^{x.y} = i^{|x|} i^{|y|} i^{-|x xor y|}, hence
   H^{(x)k} = ((1+i)/2)^{-k} D C^{(x)k} D,  D = diag(i^{|x|}),  C = (1/2)[[1+i,1-i],[1-i,1+i]].
B. Cyclic Bluestein (DFT_n from one cyclic convolution of length L >= 2n-1 with a FIXED chirp) and
   Rader (DFT_p from one cyclic convolution of length p-1 with a FIXED vector): exact identities.
C. The AFKL / ncc-certificate shift map for the WHT: H diag(H e_s / N) H = shift by s (xor), checked
   over the prime field F_p (p = 10007), and the network cost of the PLAIN radix-2 circuit placed
   in that map: edges of the composite network vs the routing lower bound sum_j dist(in_j, out_{j^s}).
   This measures the constant-factor loss of the L0 -> network-coding map (how far a circuit must
   beat the plain FFT before the map yields a counterexample).
"""
import cmath, itertools, random
from collections import deque

def popcount(x):
    return bin(x).count("1")

# ---------- A ----------
def check_A(kmax=6):
    ok = True
    for x in range(4):
        for y in range(4):
            for k in (2,):
                pass
    for k in range(1, kmax + 1):
        N = 1 << k
        c = ((1 + 1j) / 2) ** (-k)
        # C^{(x)k}[a][b] = prod over bits of C[a_t][b_t]; C[0][0]=C[1][1]=(1+i)/2, C[0][1]=C[1][0]=(1-i)/2
        for a in range(N):
            for b in range(N):
                d = popcount(a ^ b)
                Cab = ((1 + 1j) / 2) ** (k - d) * ((1 - 1j) / 2) ** d
                lhs = (-1) ** popcount(a & b)
                rhs = c * (1j ** popcount(a)) * Cab * (1j ** popcount(b))
                if abs(lhs - rhs) > 1e-9:
                    ok = False
    return ok

# ---------- B ----------
def dft(x):
    n = len(x)
    return [sum(x[j] * cmath.exp(-2j * cmath.pi * j * k / n) for j in range(n)) for k in range(n)]

def cyc_conv(a, b):
    L = len(a)
    return [sum(a[j] * b[(t - j) % L] for j in range(L)) for t in range(L)]

def bluestein(x, L):
    n = len(x)
    assert L >= 2 * n - 1
    w = lambda t: cmath.exp(-1j * cmath.pi * t * t / n)          # exp(-pi i t^2 / n)
    a = [x[j] * w(j) if j < n else 0 for j in range(L)]
    b = [0] * L
    for u in range(-(n - 1), n):
        b[u % L] = 1 / w(u)                                          # FIXED kernel
    c = cyc_conv(a, b)
    return [w(k) * c[k] for k in range(n)]

def prim_root(p):
    for g in range(2, p):
        if len({pow(g, e, p) for e in range(p - 1)}) == p - 1:
            return g

def rader(x):
    p = len(x)
    g = prim_root(p)
    ginv = pow(g, p - 2, p)
    a = [x[pow(g, q, p)] for q in range(p - 1)]                     # input permuted
    b = [cmath.exp(-2j * cmath.pi * pow(ginv, q, p) / p) for q in range(p - 1)]   # FIXED kernel
    c = cyc_conv(a, b)
    X = [0] * p
    X[0] = sum(x)
    for q in range(p - 1):
        X[pow(ginv, q, p)] = x[0] + c[q]
    return X

def check_B():
    random.seed(1)
    worst = 0.0
    for n in range(2, 14):
        x = [complex(random.gauss(0, 1), random.gauss(0, 1)) for _ in range(n)]
        L = 1
        while L < 2 * n - 1:
            L *= 2
        y1, y2 = dft(x), bluestein(x, L)
        worst = max(worst, max(abs(u - v) for u, v in zip(y1, y2)))
    for p in (3, 5, 7, 11, 13):
        x = [complex(random.gauss(0, 1), random.gauss(0, 1)) for _ in range(p)]
        y1, y2 = dft(x), rader(x)
        worst = max(worst, max(abs(u - v) for u, v in zip(y1, y2)))
    return worst

# ---------- C ----------
P = 10007

def wht_modp(v):
    v = list(v)
    N = len(v)
    h = 1
    while h < N:
        for i in range(0, N, 2 * h):
            for j in range(i, i + h):
                a, b = v[j], v[j + h]
                v[j], v[j + h] = (a + b) % P, (a - b) % P
        h *= 2
    return v

def check_C_algebra(k):
    N = 1 << k
    Ninv = pow(N, P - 2, P)
    for s in range(N):
        es = [1 if t == s else 0 for t in range(N)]
        d = [(u * Ninv) % P for u in wht_modp(es)]
        for j in range(N):
            ej = [1 if t == j else 0 for t in range(N)]
            y = wht_modp([(u * dd) % P for u, dd in zip(wht_modp(ej), d)])
            want = [1 if t == (j ^ s) else 0 for t in range(N)]
            if y != want:
                return False
    return True

def composite_graph(k, relays):
    """Undirected graph of WHT circuit -> diagonal -> WHT circuit (radix-2, in place).
    Node ('v', layer, idx): value of coordinate idx after `layer` butterfly layers, layers 0..2k
    (layer k -> diagonal applied: modelled as one extra node per coordinate).
    relays=False: each butterfly output node has the two input nodes as neighbours (gate = node).
    relays=True: degree <= 3 as in the ncc-certificate (each use of a value gets its own relay)."""
    N = 1 << k
    adj = {}
    def add(u, v):
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)
    E = 0
    def node(L, i):
        return ('v', L, i)
    layer = 0
    # first transform: layers 0..k
    for half in range(2):
        for t in range(k):
            h = 1 << t
            L = half * (k + 1) + t
            for i in range(N):
                j = i ^ h
                # output (L+1, i) uses inputs (L, i) and (L, j)
                if relays:
                    for src in (i, j):
                        r = ('r', L, src, i)          # relay for this use of value (L, src)
                        add(node(L, src), r); add(r, node(L + 1, i)); E += 2
                else:
                    add(node(L, i), node(L + 1, i)); add(node(L, j), node(L + 1, i)); E += 2
        if half == 0:
            # diagonal: (k, i) -> (k+1, i)
            for i in range(N):
                add(node(k, i), node(k + 1, i)); E += 1
    return adj, E, node(0, 0)[0:1]

def bfs(adj, src):
    dist = {src: 0}
    dq = deque([src])
    while dq:
        u = dq.popleft()
        for w in adj[u]:
            if w not in dist:
                dist[w] = dist[u] + 1
                dq.append(w)
    return dist

def check_C_cost(k, relays):
    N = 1 << k
    adj, E, _ = composite_graph(k, relays)
    last = 2 * k + 1
    D = [bfs(adj, ('v', 0, j)) for j in range(N)]
    sums = []
    for s in range(N):
        sums.append(sum(D[j][('v', last, j ^ s)] for j in range(N)))
    return E, min(sums), max(sums), sum(sums) / N

if __name__ == "__main__":
    print("A  H^(x)k = ((1+i)/2)^-k D C^(x)k D for k=1..6:", check_A())
    print("B  Bluestein (n=2..13, L=2^ceil) and Rader (p=3..13): max abs error %.2e" % check_B())
    for k in range(1, 6):
        print("C  k=%d: H diag(H e_s/N) H = xor-shift by s over F_%d for all s,j: %s" % (k, P, check_C_algebra(k)))
    print("C  composite network cost vs routing lower bound (rate-1, unit capacity):")
    print("   k  N   relays  edges(code cost)  min_s sum dist  max_s sum dist  ratio edges/max_s")
    for k in range(1, 9):
        for relays in (False, True):
            E, mn, mx, avg = check_C_cost(k, relays)
            print("   %d %4d  %-5s  %8d  %10d  %10d  %.3f" % (k, 1 << k, relays, E, mn, mx, E / mx))
