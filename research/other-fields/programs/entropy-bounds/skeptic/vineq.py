# vineq.py -- v-nc-gain: test nc-gain's Ingleton and DFZ (65) row encodings (nclp.ingleton_row, nclp.dfz65_row).
#  (i) Ingleton must hold (row <= 0) on random linear rank functions (GF(2), GF(3)) and FAIL on the Vamos matroid with
#      A,B,C,D = the four pairs (value +1);  (ii) DFZ (65) must hold on random GF(3) linear rank functions (char != 2).
import sys, random, itertools
import numpy as np
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'lp-and-layouts'))
import nclp
def rank_mod(rows, p):
    M = [list(r) for r in rows]; r = 0
    if not M: return 0
    ncol = len(M[0])
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; inv = pow(M[r][c], p - 2, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
    return r
def linear_h(n, p, dim, rng, maxgen=2):
    gens = [[[rng.randrange(p) for _ in range(dim)] for _ in range(rng.randint(1, maxgen))] for _ in range(n)]
    h = np.zeros(1 << n)
    for mask in range(1, 1 << n):
        rows = [v for i in range(n) if mask >> i & 1 for v in gens[i]]
        h[mask] = rank_mod(rows, p)
    return h
ident = np.arange(1 << 8)
rng = random.Random(5); worst_ing = -1e9; worst_dfz = -1e9
for t in range(150):
    p = (2, 3)[t % 2]; h = linear_h(4, p, 4, rng)
    for A, B, C, D in itertools.permutations([1, 2, 4, 8]):
        worst_ing = max(worst_ing, nclp.row_value(nclp.ingleton_row(ident, A, B, C, D), h))
for t in range(120):
    h = linear_h(7, 3, rng.choice([3, 4, 5]), rng)
    L = [1 << i for i in range(7)]
    for _ in range(30):
        rng.shuffle(L); worst_dfz = max(worst_dfz, nclp.row_value(nclp.dfz65_row(ident, L), h))
print(f"Ingleton on 150 random GF(2)/GF(3) rank functions, all orderings: max row value {worst_ing} (must be <= 0)")
print(f"DFZ (65) on 120 random GF(3) rank functions x 30 letter orders: max row value {worst_dfz} (must be <= 0)")
# Vamos: 8 elements a a' b b' c c' d d' -> bits 0..7; pairs A={0,1}, B={2,3}, C={4,5}, D={6,7}
pairs = [0b11, 0b1100, 0b110000, 0b11000000]
dep = {pairs[0] | pairs[1], pairs[0] | pairs[2], pairs[0] | pairs[3], pairs[1] | pairs[2], pairs[1] | pairs[3]}
h = np.zeros(256)
for mask in range(1, 256):
    k = bin(mask).count('1'); h[mask] = 3 if mask in dep else min(k, 4)
print("Vamos, Ingleton(A,B;C,D) row value:", nclp.row_value(nclp.ingleton_row(ident, *pairs), h), "(must be +1: a violation)")
# Fano matroid (GF(2)^3, 7 points): DFZ (65) is a characteristic != 2 inequality, so some labelling should violate it
pts = [[(v >> i) & 1 for i in range(3)] for v in range(1, 8)]
hF = np.zeros(128)
for mask in range(1, 128): hF[mask] = rank_mod([pts[i] for i in range(7) if mask >> i & 1], 2)
best = max(nclp.row_value(nclp.dfz65_row(ident, [1 << i for i in perm]), hF) for perm in itertools.permutations(range(7)))
print("Fano matroid, DFZ (65) max row value over all 5040 labellings:", best, "(> 0 expected for a char-not-2 inequality)")
hF3 = np.zeros(128)
pts3 = [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]]
for mask in range(1, 128): hF3[mask] = rank_mod([pts3[i] for i in range(7) if mask >> i & 1], 3)
best3 = max(nclp.row_value(nclp.dfz65_row(ident, [1 << i for i in perm]), hF3) for perm in itertools.permutations(range(7)))
print("same 7 vectors over GF(3) (non-Fano), DFZ (65) max row value:", best3, "(must be <= 0)")
