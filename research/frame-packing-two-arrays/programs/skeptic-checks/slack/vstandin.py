#!/usr/bin/env python3
"""v-slack own check of the block-2 stand-in 2 in 3 (ex2-certnum sx/out/scheme_2in3.json, read-only).

(1) Exact evaluation in V (x) V with rho(g) = g (x) g, own arithmetic (Fractions): outputs equal rho(D1 D0) in their
    block and 0 in the other.
(2) Membership of every mixing matrix in A = End(U6) + End(U4) + End(L2): commutes with the swap and maps the
    diagonal span U4 into itself and the symmetric off-diagonal span U6 into itself (block structure).
(3) Depth: every move acts on input sources only, so each output's coefficient on its input lies in
    span(B) + span(B R B). In C[G] every group element of B and of B R B has frame distance <= 1 from L0
    (computed here from matrices), while K = D1 D0 has frame distance 2. Since group elements are a basis of C[G],
    NO lift of this scheme to the engine game is a formal identity (independent of Theorem 2W).
(4) Parity: M = 3, W m = 4, deficit 1 (odd).
"""
from fractions import Fraction as F
import json, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from vframes import matmul, ident, dirz, ONE, ZERO, cadd, cmul, sp_action, lagrangians, dist

import os
PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'skeptic-checks', 'slack', 'scheme_2in3.json')

def parse(M): return [[(F(e[0]), F(e[1])) for e in r] for r in M]

def kron(A, B):
    n, p = len(A), len(B)
    return [[cmul(A[i // p][j // p], B[i % p][j % p]) for j in range(n * p)] for i in range(n * p)]

def zeros(r, c): return [[ZERO] * c for _ in range(r)]

def madd(A, B): return [[cadd(a, b) for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]

def main():
    d = json.load(open(PATH))
    W, s = d['W'], d['s']
    # X0: x -> x ^ 1 on C^4 with x = 2*bit1 + bit0 ; D0 = (1+i)/2 I + (1-i)/2 X0 ; D1 with bit1
    D0 = dirz(2, 1); D1 = dirz(2, 2)
    R = kron(D0, D0); Kr = kron(matmul(D1, D0), matmul(D1, D0))
    srcs = []
    for l in range(W):
        S = zeros(16, 16 * W)
        for i in range(16): S[i][16 * l + i] = ONE
        srcs.append(S)
    for j, row in enumerate(d['moves']):
        acc = zeros(16, 16 * W)
        for i, M in enumerate(row):
            acc = madd(acc, matmul(parse(M), srcs[i]))
        srcs.append(matmul(R, acc))
    ok_out = []
    for t, row in enumerate(d['outputs']):
        acc = zeros(16, 16 * W)
        for i, M in enumerate(row):
            acc = madd(acc, matmul(parse(M), srcs[i]))
        good = True
        for l in range(W):
            blk = [r[16 * l:16 * (l + 1)] for r in acc]
            want = Kr if l == t else zeros(16, 16)
            if blk != want: good = False
        ok_out.append(good)
    # (2) membership in A: index e_a e_b <-> 4a + b
    def in_A(M):
        n = 16
        # commutes with swap
        sw = lambda k: 4 * (k % 4) + (k // 4)
        for i in range(n):
            for j in range(n):
                if M[i][j] != M[sw(i)][sw(j)]: return False
        diag = {4 * a + a for a in range(4)}
        # symmetric part: maps U4 (diag) into U4: columns in diag have support in diag rows
        for j in diag:
            for i in range(n):
                if i not in diag and M[i][j] != ZERO: return False
        # maps symmetric off-diagonal vectors e_ab + e_ba into span with no diag component
        for a in range(4):
            for b in range(a + 1, 4):
                j1, j2 = 4 * a + b, 4 * b + a
                for i in diag:
                    if cadd(M[i][j1], M[i][j2]) != ZERO: return False
        return True
    mats = [parse(M) for row in d['moves'] for M in row] + [parse(M) for row in d['outputs'] for M in row]
    all_in_A = all(in_A(M) for M in mats)
    # (3) depth and frames
    moves_on_inputs_only = all(all(i < W or all(e == ['0', '0'] for r in M for e in r) for i, M in enumerate(row))
                               for row in d['moves'])
    L0, Ls = lagrangians(2)
    aD0 = sp_action(2, D0); aK = sp_action(2, matmul(D1, D0))
    img = lambda a, L: frozenset(a(v) for v in L)
    res = {'outputs_exact_in_VxV': ok_out, 'all_mixing_matrices_in_A_block_algebra': all_in_A,
           'moves_act_on_inputs_only (depth 1)': moves_on_inputs_only,
           'frame distance of R = D0': dist(2, img(aD0, L0), L0),
           'frame distance of K = D1 D0': dist(2, img(aK, L0), L0),
           'M': s, 'Wm': 2 * W, 'deficit': 2 * W - s, 'deficit_parity': 'odd' if (2 * W - s) % 2 else 'even'}
    return res

if __name__ == '__main__':
    r = main()
    for k, v in r.items(): print(k, ':', v)
    json.dump(r, open(__file__.replace('.py', '.json'), 'w'), indent=1)
