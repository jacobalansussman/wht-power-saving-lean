#!/usr/bin/env python3
"""span_clo.py m f prime: tests (1) and (2) of span.py for larger sizes, span of the group found by closure; vectors are packed into Python integers (64-bit limbs)."""
import sys, itertools, time, struct
m = int(sys.argv[1]); f = int(sys.argv[2]); p = int(sys.argv[3])
d = 2**m; n = d**f; L = n*n
I = next(x for x in range(2, p) if x*x % p == p-1)
inv2 = pow(2, p-2, p); ipow = [1, I, p-1, (p-I) % p]; MASK = (1 << 64) - 1
t0 = time.time()
def comp(a, b):
    pa, ha = a; pb, hb = b
    return (tuple(pa[pb[x]] for x in range(d)), tuple((hb[x] + ha[pb[x]]) % 4 for x in range(d)))
gens = []
for z in range(1, d): gens.append((tuple(x ^ z for x in range(d)), (0,)*d))
for i in range(m):
    for j in range(m):
        if i != j: gens.append((tuple(x ^ (((x >> j) & 1) << i) for x in range(d)), (0,)*d))
for z in range(d):
    for c in range(2):
        if z == 0 and c == 0: continue
        gens.append((tuple(range(d)), tuple((bin(z & x).count('1') + c) % 2 for x in range(d))))
ident = (tuple(range(d)), (0,)*d)
print(f"m={m} f={f} p={p} dim End(V)={L}  (span by closure under the generators; the group is not listed)", flush=True)
cols = list(itertools.product(range(d), repeat=f)); colidx = {c: k for k, c in enumerate(cols)}
def tpow(b):
    pm, ph = b
    return ([colidx[tuple(pm[x] for x in c)] for c in cols], [ipow[sum(ph[x] for x in c) % 4] for c in cols])
fmt = '<%dQ' % L
pack = lambda v: int.from_bytes(struct.pack(fmt, *v), 'little')
unpack = lambda x: struct.unpack(fmt, x.to_bytes(8*L, 'little'))
class Span:
    def __init__(s): s.basis = []
    def add(s, v):           # v: list of residues
        x = pack(v)
        for k, r in s.basis:
            c = ((x >> (64*k)) & MASK) % p
            if c: x += (p - c) * r
        v = [a % p for a in unpack(x)]
        k = next((j for j, a in enumerate(v) if a), None)
        if k is None: return False
        iv = pow(v[k], p-2, p); s.basis.append((k, pack([a*iv % p for a in v]))); return True
    def dim(s): return len(s.basis)
def kron(A, Bm):
    ra, rb = len(A), len(Bm)
    return [[A[i][j]*Bm[k][l] % p for j in range(ra) for l in range(rb)] for i in range(ra) for k in range(rb)]
a, bb = (1 + I)*inv2 % p, (1 - I)*inv2 % p
C = [[a, bb], [bb, a]]
def kpow(M, e):
    R = [[1]]
    for _ in range(e): R = kron(R, M)
    return R
Kf = kpow(kpow(C, m), f)
A = Span(); T = []
def vec(b):
    rows, vals = tpow(b); v = [0]*L
    for c in range(n): v[rows[c]*n + c] = vals[c]
    return v
A.add(vec(ident)); T.append(ident); q = [ident]
while q:
    x = q.pop()
    for g in gens:
        c = comp(g, x)
        if A.add(vec(c)): T.append(c); q.append(c)
tp = [tpow(b) for b in T]
dA = A.dim()
S = Span(); S.basis = list(A.basis)
inA = not S.add([x for row in Kf for x in row])
S.basis = list(A.basis)
print(f"(1) dim A_f = {dA}   K^f in A_f: {inA}   [{time.time()-t0:.1f}s]", flush=True)
cnt = 0
for rows, vals in tp:
    M = [None]*n
    for c in range(n):
        vc = vals[c]; M[rows[c]] = [vc*x % p for x in Kf[c]]
    if S.add([x for row in M for x in row]): cnt += 1
print(f"(2) dim(A_f + A_f K^f) = {S.dim()} = dim A_f + {cnt};  2 dim A_f = {2*dA};  dim(A K^f ∩ A) = {2*dA - S.dim()};"
      f"  finite-f count: blocks >= W * {cnt}/{dA}   [{time.time()-t0:.1f}s]", flush=True)
