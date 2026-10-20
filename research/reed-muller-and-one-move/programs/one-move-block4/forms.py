#!/usr/bin/env python3
"""forms.py (w8-rmexact): block-m binary quadratic forms, GL(m,2) action, the virtual-target distribution of
Lemma C (targets h h^T, h in H), and per-orbit support size s(F) and odd fraction a(F).
Conventions: direction x in 1..2^m-1 (bit i = x_i); a form Q is 10 bits q_ij (i<=j); Q(x) = sum q_ii x_i + sum q_ij x_i x_j;
a symmetric matrix S is 10 bits; <Q,S> = sum q_ij S_ij (i<=j); <Q, x x^T> = Q(x)."""
import itertools, sys
from fractions import Fraction as Fr

def setup(m):
    idx = [(i, j) for i in range(m) for j in range(i, m)]
    pos = {p: b for b, p in enumerate(idx)}
    nf = 1 << len(idx)
    dirs = list(range(1, 1 << m))
    def bit(x, i): return (x >> i) & 1
    def Qval(Q, x):
        v = 0
        for b, (i, j) in enumerate(idx):
            if (Q >> b) & 1 and bit(x, i) and bit(x, j): v ^= 1
        return v
    supp = [frozenset(x for x in dirs if Qval(Q, x)) for Q in range(nf)]
    def matmul(A, B):   # matrices as lists of row bitmasks
        n = len(A); out = []
        for i in range(n):
            r = 0
            for k in range(n):
                if bit(A[i], k): r ^= B[k]
            out.append(r)
        return out
    def transpose(A):
        n = len(A); return [sum(bit(A[r], c) << r for r in range(n)) for c in range(n)]
    def enc_sym(S):
        v = 0
        for b, (i, j) in enumerate(idx):
            if bit(S[i], j): v |= 1 << b
        return v
    def rank2(A):
        rows = list(A); r = 0; n = len(rows)
        for c in range(m):
            p = next((i for i in range(r, n) if bit(rows[i], c)), None)
            if p is None: continue
            rows[r], rows[p] = rows[p], rows[r]
            for i in range(n):
                if i != r and bit(rows[i], c): rows[i] ^= rows[r]
            r += 1
        return r
    # GL(m,2): all invertible matrices (rows as bitmasks); h acts on directions by x -> h x
    GL = [list(rows) for rows in itertools.product(range(1, 1 << m), repeat=m) if rank2(list(rows)) == m]
    def apply(h, x):   # (h x)_i = row_i . x
        return sum((bin(h[i] & x).count('1') & 1) << i for i in range(m))
    return dict(m=m, idx=idx, nf=nf, dirs=dirs, Qval=Qval, supp=supp, GL=GL, apply=apply, matmul=matmul,
                transpose=transpose, enc_sym=enc_sym, rank2=rank2)

def pair(Q, S):
    return bin(Q & S).count('1') & 1

def main():
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    D = setup(m); GL = D['GL']; supp = D['supp']; nf = D['nf']
    print("m=%d |GL|=%d forms=%d" % (m, len(GL), nf))
    # virtual targets: h h^T over h in GL
    from collections import Counter
    T = Counter(D['enc_sym'](D['matmul'](h, D['transpose'](h))) for h in GL)
    print("distinct targets h h^T:", len(T), "multiplicities", set(T.values()))
    a = [Fr(sum(c for S, c in T.items() if pair(Q, S)), len(GL)) for Q in range(nf)]
    # orbits of forms under GL: Q -> Q o h^-1 ; use support sets: supp(Q o h^-1) = h(supp Q)
    key = {supp[Q]: Q for Q in range(nf)}
    assert len(key) == nf
    orb = [-1] * nf; orbits = []
    for Q in range(nf):
        if orb[Q] >= 0: continue
        k = len(orbits); mem = set()
        for h in GL:
            Q2 = key[frozenset(D['apply'](h, x) for x in supp[Q])]
            mem.add(Q2)
        for Q2 in mem: orb[Q2] = k
        orbits.append(sorted(mem))
    for k, mem in enumerate(orbits):
        Q = mem[0]
        print("orbit %d: size %3d  |supp| = %2d  a(F) = %s (odd fraction over targets)  tr-odd for target I: %d/%d" % (
            k, len(mem), len(supp[Q]), a[Q], sum(pair(Q2, D['enc_sym']([1 << i for i in range(m)])) for Q2 in mem), len(mem)))
if __name__ == '__main__':
    main()
