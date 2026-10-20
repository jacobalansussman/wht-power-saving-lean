#!/usr/bin/env python3
"""ctr-search / p112.py (own code; numpy/scipy from our pylib; run with python3 -I).
The hub's 112-port design in HAMMING coordinates: the E8 root lines of shape (+-1)^4 0^4 on the 14 weight-4 words
("tetrads") of the extended Hamming [8,4,4] code.  (The 8 frame lines +-2e_i complete E8; they are the removed lines.)
Inner products are in {0, +-2, +-4} (norm 4).  Binary Gram  B_ab = 1 + <r_a, r_b>/2 (mod 2): 1 on the diagonal and on
orthogonal pairs (E-pairs: must be delivered exactly), 0 on 60/120-degree pairs (junk allowed).
Labels u_a in F_2^9 with u_a . u_b = B_ab (factored here; any two factorings differ by an isometry of F_2^9).

Decoder model (pr717 D1, and cx-medium section 2): a decoder is a list of terms (psi_k, phi_k) (scatter, gather
coefficient vectors over the 112 ports) with  M = sum_k psi_k phi_k^T,  M_ab = delta_ab on every E-pair.
Costs of a term with scatter support R and gather support A:
  pair(R, A) = rank_F2 B[R, A]            (Theorem C pairing rank; the hub's D1 'cst')
  star(A)    = rank_F2 U(A) = dim span of the gathered labels   (gcert E2: targets wait at frame 0; cx-medium floor)
"""
import sys, os, itertools
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))

def hamming_tetrads():
    gens = [0b11110000, 0b11001100, 0b10101010, 0b11111111]
    words = set()
    for m in range(16):
        w = 0
        for i in range(4):
            if m >> i & 1: w ^= gens[i]
        words.add(w)
    tet = sorted(tuple(i for i in range(8) if w >> (7 - i) & 1) for w in words if bin(w).count('1') == 4)
    assert len(tet) == 14
    return tet

TET = hamming_tetrads()

def build_lines():
    """112 lines; representative has first nonzero coordinate +1.  Returns lines (tuples), tetrad index per line."""
    L, T = [], []
    for ti, t in enumerate(TET):
        for s in itertools.product((1, -1), repeat=3):
            v = [0] * 8
            v[t[0]] = 1
            for j in range(3): v[t[j + 1]] = s[j]
            L.append(tuple(v)); T.append(ti)
    return L, T

LINES, TETOF = build_lines()
N = len(LINES)
FRAME = [tuple(2 if j == i else 0 for j in range(8)) for i in range(8)]   # the 8 removed lines (orthonormal x2)

def ip(a, b): return sum(x * y for x, y in zip(a, b))

G = np.array([[ip(a, b) for b in LINES] for a in LINES], dtype=np.int64)       # integer Gram (diag 4)
Bm = ((1 + G // 2) % 2).astype(np.uint8)                                       # binary Gram
E = (Bm == 1)                                                                  # exact pairs (incl. diagonal)
JUNK = ~E

# ---------------------------------------------------------------- F2 helpers (ints as bit vectors)
def rank2_rows(rows):
    basis = []
    for v in rows:
        v = int(v)
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)

def dot(x, y): return bin(x & y).count('1') & 1

def labels():
    """u_a as 9-bit ints with u_a.u_b = B_ab; Gaussian factoring B = U U^T over F2 (deterministic)."""
    n = N
    rows = [int(''.join(map(str, Bm[i])), 2) for i in range(n)]
    I, basis = [], []
    for i in range(n):
        w = rows[i]
        for b in basis: w = min(w, w ^ b)
        if w: basis.append(w); I.append(i)
    k = len(I)
    S = [[int(Bm[i][j]) for j in I] for i in I]
    def form(x, y): return sum(x[i] & S[i][j] & y[j] for i in range(k) for j in range(k)) & 1
    vecs = [list(map(int, format(x, f'0{k}b'))) for x in range(1, 1 << k)]
    out = []
    for v in vecs:                                   # greedy orthonormal basis for the form S (deterministic order)
        if form(v, v) == 1 and all(form(v, w) == 0 for w in out):
            out.append(v)
            if len(out) == k: break
    assert len(out) == k
    P = np.array(out, dtype=np.int64).T % 2         # columns orthonormal: P^T S P = I
    # U = B[:,I] P  (then U U^T restricted ... check below); since S = P^-T P^-1, U = B[:,I] P works: U U^T = B[:,I] P P^T B[I,:]
    BI = np.array([[int(Bm[i][j]) for j in I] for i in range(n)], dtype=np.int64)
    U = (BI @ P) % 2
    labs = [int(''.join(map(str, row)), 2) for row in U]
    for i in range(n):
        for j in range(n):
            assert dot(labs[i], labs[j]) == Bm[i][j], 'label factoring failed'
    return labs, k

LABS, K = labels()
assert K == 9 and len(set(LABS)) == N and all(bin(u).count('1') % 2 for u in LABS)

def pair_cost(R, A):
    """rank_F2 B[R, A]"""
    if not R or not A: return 0
    rows = []
    Al = list(A)
    for a in R:
        v = 0
        for b in Al: v = (v << 1) | int(Bm[a][b])
        rows.append(v)
    return rank2_rows(rows)

def star_cost(A): return rank2_rows([LABS[b] for b in A])

def supports(Phi, Psi, tol=0):
    n = len(Phi)
    out = []
    for k in range(n):
        A = [b for b in range(N) if Phi[k][b] != 0]
        R = [a for a in range(N) if Psi[k][a] != 0]
        out.append((R, A))
    return out

def costs(Phi, Psi):
    sp = supports(Phi, Psi)
    pc = [pair_cost(R, A) for R, A in sp]
    sg = [star_cost(A) for R, A in sp]
    ss = [star_cost(R) for R, A in sp]
    return dict(pair=sum(pc), pair_terms=pc, star_gather=sum(sg), star_gather_terms=sg, star_scatter=sum(ss),
                star_min=sum(min(x, y) for x, y in zip(sg, ss)), n=len(Phi))

def check_exact(Phi, Psi):
    """exact rational check: (Psi^T Phi)_ab == delta_ab on every E-pair. Phi, Psi: lists of lists of Fraction/int."""
    n = len(Phi)
    bad = 0
    for a in range(N):
        for b in range(N):
            if not E[a][b]: continue
            s = Fraction(0)
            for k in range(n):
                pa, fb = Psi[k][a], Phi[k][b]
                if pa and fb: s += Fraction(pa) * Fraction(fb)
            if s != (1 if a == b else 0): bad += 1
    return bad == 0, bad

def junk_map(Phi, Psi):
    M = {}
    n = len(Phi)
    for a in range(N):
        for b in range(N):
            if E[a][b]: continue
            s = Fraction(0)
            for k in range(n):
                if Psi[k][a] and Phi[k][b]: s += Fraction(Psi[k][a]) * Fraction(Phi[k][b])
            if s: M[(a, b)] = s
    return M

# ---------------------------------------------------------------- cliques (8 mutually orthogonal lines)
def max_cliques():
    adj = [0] * N
    for i in range(N):
        for j in range(N):
            if i != j and E[i][j]: adj[i] |= 1 << j
    best = [0]; res = []
    def bk(R, P, X, size):
        if P == 0 and X == 0:
            if size > best[0]: best[0] = size; res.clear()
            if size == best[0]: res.append(R)
            return
        if size + bin(P).count('1') < best[0]: return
        u = (P | X).bit_length() - 1
        cand = P & ~adj[u]
        while cand:
            v = cand.bit_length() - 1; cand &= ~(1 << v)
            bk(R | (1 << v), P & adj[v], X & adj[v], size + 1)
            P &= ~(1 << v); X |= 1 << v
    bk(0, (1 << N) - 1, 0, 0)
    return best[0], [[i for i in range(N) if c >> i & 1] for c in res]

# ---------------------------------------------------------------- symmetry: signed coordinate maps preserving the code
LINE_INDEX = {l: i for i, l in enumerate(LINES)}

def canon(v):
    v = tuple(v)
    for x in v:
        if x != 0:
            return (v, 1) if x > 0 else (tuple(-y for y in v), -1)
    raise ValueError

def act(perm, signs):
    """signed permutation g: e_i -> signs[i] e_perm[i].  Returns (pi, eps): g(r_a) = eps[a] r_pi[a], or None."""
    pi, eps = [], []
    for l in LINES:
        w = [0] * 8
        for i in range(8):
            if l[i]: w[perm[i]] = signs[i] * l[i]
        c, s = canon(w)
        if c not in LINE_INDEX: return None
        pi.append(LINE_INDEX[c]); eps.append(s)
    return pi, eps

def code_perms():
    """the 1344 coordinate permutations preserving the tetrads"""
    tset = set(TET)
    out = []
    for p in itertools.permutations(range(8)):
        if all(tuple(sorted(p[i] for i in t)) in tset for t in TET): out.append(p)
    return out

if __name__ == '__main__':
    print('lines', N, 'window k', K, 'labels distinct+odd ok')
    print('E-pairs per row (incl diag):', set(int(x) for x in E.sum(1)), ' junk per row:', set(int(x) for x in JUNK.sum(1)))
    om, cl = max_cliques()
    print('omega', om, 'max cliques', len(cl))
    P = code_perms()
    print('code permutations', len(P))
