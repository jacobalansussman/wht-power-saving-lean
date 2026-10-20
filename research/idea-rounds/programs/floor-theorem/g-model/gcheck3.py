# g-model check 3: column permutations free, symmetric witnesses, at the engine's block size m = 3 with f = 3 columns
# (length 512), plus m = 2 with f = 3 and f = 4.  Pure Python.  Ranks modulo one prime p = 998244353 (a rank modulo p
# is a lower bound for the true rank; i -> a square root of -1 mod p).
import itertools, sys
p = int(sys.argv[1]) if len(sys.argv) > 1 else 998244353   # second run: 469762049
s = next(r for r in (pow(g, (p - 1) // 4, p) for g in range(2, 200)) if r * r % p == p - 1)
inv2 = pow(2, p - 2, p); a = (1 + s) * inv2 % p; b = (1 - s) * inv2 % p
def rank(M):
    M = [row[:] for row in M]; r = 0; n = len(M); c = len(M[0])
    for col in range(c):
        k = next((j for j in range(r, n) if M[j][col]), None)
        if k is None: continue
        M[r], M[k] = M[k], M[r]; iv = pow(M[r][col], p - 2, p); pv = [x * iv % p for x in M[r]]; M[r] = pv
        for j in range(n):
            if j != r and M[j][col]:
                fc = M[j][col]; M[j] = [(x - fc * y) % p for x, y in zip(M[j], pv)]
        r += 1
        if r == n: break
    return r
def mm(A, B):
    Bt = list(zip(*B)); return [[sum(x * y for x, y in zip(r, c) if x and y) % p for c in Bt] for r in A]
def kron(A, B): return [[x * y % p for x in ra for y in rb] for ra in A for rb in B]
for m, f in ((2, 3), (2, 4), (3, 3)):
    q = 2 ** m; n = q ** f
    def dirz(z):
        D = [[0] * q for _ in range(q)]
        for x in range(q): D[x][x] = a; D[x ^ z][x] = (D[x ^ z][x] + b) % p
        return D
    dg = list(itertools.product(range(q), repeat=f))
    pairs = [sum(d[i] == d[j] for i in range(f) for j in range(i + 1, f)) for d in dg]
    allq = [int(all(x == d[0] for x in d)) for d in dg]
    p12 = [int(d[0] == d[1]) for d in dg]
    out = []
    for k in range(1, m + 1):
        g = [[int(i == j) for j in range(q)] for i in range(q)]
        for j in range(k): g = mm(dirz(1 << j), g)
        G = [[1]]
        for _ in range(f): G = kron(G, g)
        w = [rank([[G[r][c] * (d[c] - d[r]) % p for c in range(n)] for r in range(n)]) for d in (p12, pairs, allq)]
        out.append((k, w))
        print(" m=%d f=%d n=%d rank-%d block%s: weight under 'first two equal' %d (formula %d), under 'number of equal pairs' %d, under 'all equal' %d"
              % (m, f, n, k, " (= kernel)" if k == m else "", w[0], 2 * (q - 2 ** (m - k)) * q ** (f - 2), w[1], w[2])); sys.stdout.flush()
    for name, idx in (("number of equal pairs", 1), ("all equal", 2)):
        ker = out[-1][1][idx]; best = max(w[idx] for k, w in out[:-1])
        print("   => column permutations free, witness '%s': proper blocks >= W * %d/%d = %.4f W" % (name, ker, best, ker / best)); sys.stdout.flush()
