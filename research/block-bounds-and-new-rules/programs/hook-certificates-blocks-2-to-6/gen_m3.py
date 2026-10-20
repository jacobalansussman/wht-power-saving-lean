# pf-logging, 13 Oct 2026.  GENERATOR (untrusted) of the block-3 foundation certificates.
# Module W = V (x) Lambda^4 V, V = C^8 (m = 3, f = 5); S = {(x,T): x in T}, N = span{e_s : s in S}.
# Emits:
#   cert_m3_dem.json : a nonsingular minor (mod p) of Q_W[S^c, S], Q = (-i)^{wt(x^y)} (= K up to a nonzero scalar)
#                      -> dem = rank K_W[S^c,S] >= size of the minor
#   cert_m3_M0.json  : a nonsingular minor of the alpha=0 block M_0 (35 x 35) -> cross-check of Lemma C at m = 3
#   cert_m3_B2.json  : exact Z[i] vectors v in N with B2_W v in N, B2 = (2 dir e1)(2 dir e2) -> f(B2) <= 280 - #vectors
# Pure Python; run with python3 -I gen_m3.py
import itertools, json, sys
from fractions import Fraction as Fr

m, k = 3, 4
n = 1 << m
P = 1000000009                       # prime, = 1 mod 4
def find_I(p):
    for g in range(2, 1000):
        if pow(g, (p - 1) // 2, p) == p - 1:
            I = pow(g, (p - 1) // 4, p); assert I * I % p == p - 1; return I
I = find_I(P)
popc = lambda x: bin(x).count('1')
S  = [(x, T) for x in range(n) for T in itertools.combinations(range(n), k) if x in T]
Sc = [(y, U) for y in range(n) for U in itertools.combinations(range(n), k) if y not in U]
Sidx = {s: i for i, s in enumerate(S)}

# ---------- Gaussian integers as (re, im) ----------
def gm(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def ga(a, b): return (a[0] + b[0], a[1] + b[1])
def gneg(a): return (-a[0], -a[1])
def perm_sign(seq):
    s = 1; seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]: s = -s
    return s
PERMS = [(pm, perm_sign(pm)) for pm in itertools.permutations(range(k))]
def gdet(Mx):                         # Leibniz, Gaussian integers
    tot = (0, 0)
    for pm, sg in PERMS:
        t = (sg, 0)
        for r in range(k): t = gm(t, Mx[r][pm[r]])
        tot = ga(tot, t)
    return tot
def modp(a): return (a[0] + a[1] * I) % P

# Q[y,x] = (-i)^{wt(x^y)}  (K = ((1+i)/2)^3 Q)
POWMI = [(1, 0), (0, -1), (-1, 0), (0, 1)]
def Q(y, x): return POWMI[popc(x ^ y) % 4]
# B2' = (2 dir e1)(2 dir e2): entry (1+i) or (1-i) per bit in {0,1}, zero unless bit 2 of x^y is 0
def B2(y, x):
    d = x ^ y
    if d & 4: return (0, 0)
    t = (1, 0)
    for j in (0, 1): t = gm(t, (1, -1) if (d >> j) & 1 else (1, 1))
    return t
def gW(g, row, col):                  # g_W[(y,U),(x,T)] = g[y,x] det g[U,T]
    (y, U), (x, T) = row, col
    return gm(g(y, x), gdet([[g(u, t) for t in T] for u in U]))

def elim_modp(Mx):
    """Gaussian elimination mod P with column pivoting order; returns list of (row, col) pivots (nonsingular minor)."""
    A = [r[:] for r in Mx]; R = len(A); C = len(A[0]); piv = []; used = [False] * R
    for c in range(C):
        pr = next((r for r in range(R) if not used[r] and A[r][c] % P), None)
        if pr is None: continue
        used[pr] = True; piv.append((pr, c))
        inv = pow(A[pr][c], P - 2, P); rowp = A[pr]
        for r in range(R):
            if r != pr and not used[r] and A[r][c] % P:
                f = A[r][c] * inv % P; rr = A[r]
                A[r] = [(rr[j] - f * rowp[j]) % P for j in range(C)]
    return piv

# ---------- 1. dem: minor of Q_W[S^c, S] ----------
Mdem = [[modp(gW(Q, r, c)) for c in S] for r in Sc]
piv = elim_modp(Mdem)
print(f"rank mod p of Q_W[S^c,S] = {len(piv)}  (dim N = {len(S)})", flush=True)
json.dump({"m": m, "k": k, "p": P, "I": I, "map": "Q = (-i)^wt(x xor y) (K up to a nonzero scalar)",
           "rows": [[Sc[r][0], list(Sc[r][1])] for r, c in piv],
           "cols": [[S[c][0], list(S[c][1])] for r, c in piv]}, open("cert_m3_dem.json", "w"))

# ---------- 2. alpha = 0 block M_0 (Lemma C), minor ----------
def shifted(T, s): return [t ^ s for t in T]          # ordered, not sorted
RS  = [T for T in itertools.combinations(range(n), k) if 0 in T]       # reps (0,T) of S
RSc = [U for U in itertools.combinations(range(n), k) if 0 not in U]   # reps (0,U) of S^c
def Mblock(g, alpha, U, T):
    tot = (0, 0)
    for s in range(n):
        sg = -1 if popc(alpha & s) & 1 else 1
        d = gdet([[g(u, t) for t in shifted(T, s)] for u in U])
        tot = ga(tot, gm((sg, 0), gm(g(0, s), d)))
    return tot
M0 = [[modp(Mblock(Q, 0, U, T)) for T in RS] for U in RSc]
piv0 = elim_modp(M0)
print(f"rank mod p of M_0 = {len(piv0)} ; 8 x that = {8 * len(piv0)}", flush=True)
json.dump({"m": m, "k": k, "p": P, "I": I, "alpha": 0,
           "rows": [list(RSc[r]) for r, c in piv0], "cols": [list(RS[c]) for r, c in piv0]}, open("cert_m3_M0.json", "w"))
for alpha in range(1, n):
    Ma = [[modp(Mblock(Q, alpha, U, T)) for T in RS] for U in RSc]
    print(f"  (information) rank mod p of M_alpha, alpha={alpha}: {len(elim_modp(Ma))}", flush=True)

# ---------- 3. B2 kernel vectors, exact over Q(i), via the blocks ----------
class G:
    __slots__ = ("r", "i")
    def __init__(s, r, i=Fr(0)): s.r, s.i = Fr(r), Fr(i)
    def __add__(s, o): return G(s.r + o.r, s.i + o.i)
    def __sub__(s, o): return G(s.r - o.r, s.i - o.i)
    def __mul__(s, o): return G(s.r * o.r - s.i * o.i, s.r * o.i + s.i * o.r)
    def inv(s):
        d = s.r * s.r + s.i * s.i; return G(s.r / d, -s.i / d)
    def nz(s): return s.r != 0 or s.i != 0
def kernel_exact(A):
    A = [r[:] for r in A]; R = len(A); C = len(A[0]); pivcols = []; row = 0
    for c in range(C):
        pr = next((r for r in range(row, R) if A[r][c].nz()), None)
        if pr is None: continue
        A[row], A[pr] = A[pr], A[row]
        iv = A[row][c].inv(); A[row] = [x * iv for x in A[row]]
        for r in range(R):
            if r != row and A[r][c].nz():
                f = A[r][c]; A[r] = [A[r][j] - f * A[row][j] for j in range(C)]
        pivcols.append(c); row += 1
    free = [c for c in range(C) if c not in pivcols]
    basis = []
    for fc in free:
        v = [G(0) for _ in range(C)]; v[fc] = G(1)
        for r, pc in enumerate(pivcols): v[pc] = G(0) - A[r][fc]
        basis.append(v)
    return basis
from math import lcm
vecs = []
for alpha in range(n):
    Ma = [[G(*Mblock(B2, alpha, U, T)) for T in RS] for U in RSc]
    ker = kernel_exact(Ma)
    for v in ker:
        full = {}
        for j, T in enumerate(RS):
            if not v[j].nz(): continue
            for s in range(n):
                sh = shifted(T, s); srt = tuple(sorted(sh)); eps = perm_sign(sh)
                sg = (-1 if popc(alpha & s) & 1 else 1) * eps
                key = Sidx[(s, srt)]
                full[key] = full.get(key, G(0)) + G(sg) * v[j]
        L = 1
        for c in full.values(): L = lcm(L, c.r.denominator, c.i.denominator)
        vecs.append({str(key): [int(c.r * L), int(c.i * L)] for key, c in full.items() if c.nz()})
    print(f"  B2 block alpha={alpha}: block rank {35 - len(ker)}, kernel {len(ker)}", flush=True)
print(f"B2: {len(vecs)} exact kernel vectors -> f(B2) <= {len(S) - len(vecs)}", flush=True)
json.dump({"m": m, "k": k, "map": "B2' = (2 dir e1)(2 dir e2)", "basis_S": [[x, list(T)] for x, T in S],
           "vectors": vecs}, open("cert_m3_B2.json", "w"))
