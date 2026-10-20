#!/usr/bin/env python3
"""Streaming cells in the commutative slice, character picture (exact, pure Python).

Slice at block m (Theorem A of e3-rmbuild): group algebra of Sym_m(F_2); a move in direction z is the element
[z z^T]; the target is [I]. Characters are binary quadratic forms Q: the move acts as s_z(Q) = (-1)^{Q(z)}, the
target as eps(Q) = (-1)^{sum_i Q(e_i)}.

A streaming cell: p input registers, r state registers, scratch; a straight-line word of ops run once per period:
  ("mv", reg, z)          paid move: reg *= s_z(Q)
  ("ad", dst, src, c)     free: dst += c * src
  ("sc", reg, c)          free: reg *= c
Per form Q the cell is a numeric matrix K(Q) on (inputs, states) -> (outputs, next states).
Transfer function H(z) = M + z C (I - z D)^{-1} B. Correct iff its Markov parameters are eps(Q) at lag delta
and 0 elsewhere (checking lags 0..delta+r+2 suffices by Cayley-Hamilton).
"""
from fractions import Fraction as F
from itertools import combinations, product
import random, sys

def forms(m):
    pairs = [(i, j) for i in range(m) for j in range(i, m)]
    for bits in range(1 << len(pairs)):
        yield {pr: (bits >> k) & 1 for k, pr in enumerate(pairs)}

def qval(Q, z, m):
    v = 0
    for (i, j), b in Q.items():
        if b and z[i] and z[j]:
            v ^= 1
    return v

def eps(Q, m):
    return -1 if sum(Q[(i, i)] for i in range(m)) % 2 else 1

def run_cell(word, nreg, Q, m, start):
    """Run the word on register vector 'start' (list of Fractions); return final registers."""
    R = list(start)
    for op in word:
        if op[0] == "mv":
            if qval(Q, op[2], m):
                R[op[1]] = -R[op[1]]
        elif op[0] == "ad":
            R[op[1]] += op[3] * R[op[2]]
        elif op[0] == "sc":
            R[op[1]] *= op[2]
    return R

def cell_matrix(cell, Q, m):
    """Columns: inputs then states; rows: outputs then next states."""
    word, nreg, ins, sts, outs, nsts = cell
    cols = ins + sts
    K = []
    for c in cols:
        start = [F(0)] * nreg
        start[c] = F(1)
        R = run_cell(word, nreg, Q, m, start)
        K.append([R[o] for o in outs + nsts])
    # transpose to rows
    return [[K[c][r] for c in range(len(cols))] for r in range(len(outs + nsts))]

def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def markov(cell, Q, m, nlags):
    word, nreg, ins, sts, outs, nsts = cell
    p, r = len(ins), len(sts)
    K = cell_matrix(cell, Q, m)
    M = [row[:p] for row in K[:p]]
    C = [row[p:] for row in K[:p]]
    B = [row[:p] for row in K[p:]]
    D = [row[p:] for row in K[p:]]
    H = [M]
    if r == 0:
        return H + [[[F(0)] * p for _ in range(p)] for _ in range(nlags - 1)]
    X = B
    for _ in range(nlags - 1):
        H.append(matmul(C, X))
        X = matmul(D, X)
    return H

def check(cell, m, delta, label, verbose=True):
    word, nreg, ins, sts, outs, nsts = cell
    p, r = len(ins), len(sts)
    nl = delta + r + 3
    bad = 0
    nq = 0
    for Q in forms(m):
        nq += 1
        e = eps(Q, m)
        H = markov(cell, Q, m, nl)
        for j, Hj in enumerate(H):
            want = e if j == delta else 0
            for a in range(p):
                for b in range(p):
                    if Hj[a][b] != (want if a == b else 0):
                        bad += 1
                        break
                else:
                    continue
                break
    paid = sum(1 for op in word if op[0] == "mv")
    if verbose:
        print(f"  {label}: block {m}, p={p}, r={r}, paid/period={paid}, delay={delta}: "
              f"{'CORRECT at all' if bad == 0 else 'WRONG at ' + str(bad) + ' of'} {nq} forms")
    return bad == 0

def e(i, m):
    return tuple(1 if k == i else 0 for k in range(m))

def plain_pipeline(m):
    # regs: 0 = input, 1..m-1 = states, m = output scratch; stage j applies move e_j
    # next state s_1 = mv_1(x); s_j = mv_j(s_{j-1}); output y = mv_m(s_{m-1})
    word = []
    nreg = 2 * m + 1
    x, S, N, y = 0, list(range(1, m)), list(range(m, 2 * m - 1)), 2 * m - 1
    # output: copy s_{m-1} into y, move e_{m-1}
    if m == 1:
        word += [("ad", y, x, F(1)), ("mv", y, e(0, m))]
        return (word, nreg, [x], [], [y], [])
    word += [("ad", y, S[-1], F(1)), ("mv", y, e(m - 1, m))]
    for j in range(m - 2, 0, -1):          # N[j] = mv_j(S[j-1])
        word += [("ad", N[j], S[j - 1], F(1)), ("mv", N[j], e(j, m))]
    word += [("ad", N[0], x, F(1)), ("mv", N[0], e(0, m))]
    return (word, nreg, [x], S, [y], N)

def plain_block(m):
    word = [("mv", 0, e(i, m)) for i in range(m)]
    return (word, 1, [0], [], [0], [])

def with_scalar_iir(m, a):
    """Plain pipeline, with a FREE scalar IIR pre-filter 1/(1 - a z) and its FIR inverse (1 - a z) after.
    Class (ii): no paid move inside any loop. Extra states: u (pre-filter), v_prev (post-filter memory)."""
    word, nreg, ins, sts, outs, nsts = plain_pipeline(m)
    u, unew, vprev, vnew, xin, yout = nreg, nreg + 1, nreg + 2, nreg + 3, nreg + 4, nreg + 5
    pre = [("ad", unew, xin, F(1)), ("ad", unew, u, F(a)),       # unew = x + a*u   (u = state)
           ("ad", ins[0], unew, F(1))]                            # pipeline input := unew
    post = [("ad", vnew, outs[0], F(1)),                          # v = pipeline output
            ("ad", yout, vnew, F(1)), ("ad", yout, vprev, F(-a))]  # y = v - a*v_prev
    return (pre + word + post, nreg + 6, [xin], sts + [u, vprev], [yout], nsts + [unew, vnew])

def with_paid_loop_hidden(m):
    """Plain pipeline plus a loop register L with a PAID move inside the loop (class iii):
    L_new = x + mv_1(L). L feeds nothing observable: a hidden (unobservable) IIR mode."""
    word, nreg, ins, sts, outs, nsts = plain_pipeline(m)
    L, Lnew = nreg, nreg + 1
    loop = [("ad", Lnew, L, F(1)), ("mv", Lnew, e(0, m)), ("ad", Lnew, ins[0], F(1))]
    return (word + loop, nreg + 2, ins, sts + [L], outs, nsts + [Lnew])

def span_contains_I(dirs, m):
    """Is I in the F_2-span of {z z^T : z in dirs}? (vectors over the m(m+1)/2 coordinates i<=j)"""
    pairs = [(i, j) for i in range(m) for j in range(i, m)]
    def vec(z):
        return sum(((z[i] & z[j]) << k) for k, (i, j) in enumerate(pairs))
    basis = []
    for z in dirs:
        v = vec(z)
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    t = sum((1 << k) for k, (i, j) in enumerate(pairs) if i == j)
    for b in sorted(basis, reverse=True):
        t = min(t, t ^ b)
    return t == 0

if __name__ == "__main__":
    print("E1. Calibration: plain cells are exact streaming units; a cell one move short is not.")
    for m in range(1, 6):
        check(plain_block(m), m, 0, "plain block cell (r=0)")
        check(plain_pipeline(m), m, m - 1, "plain FIR pipeline")
    for m in (2, 3, 4):
        cell = plain_pipeline(m)
        word = [op for op in cell[0] if not (op[0] == "mv" and op[2] == e(m - 1, m))]
        check((word,) + cell[1:], m, m - 1, "negative control (last stage's move deleted)")
    print("E2. Slice consistency: c < m moves per period can never carry the target [I] (rank argument).")
    for m in (3, 4, 5):
        nz = [z for z in product((0, 1), repeat=m) if any(z)]
        cnt = tot = 0
        for c in range(1, m):
            for dirs in combinations(nz, c):
                tot += 1
                cnt += span_contains_I(dirs, m)
        full = sum(span_contains_I(d, m) for d in combinations(nz, m))
        print(f"  block {m}: direction sets with c < m whose span holds I: {cnt} of {tot};"
              f" with c = m: {full} of {sum(1 for _ in combinations(nz, m))}")
    print("E3. Class (ii): free scalar IIR inside the cell (a = 3/7) changes nothing; class (iii) hidden paid loop.")
    for m in (2, 3, 4):
        check(with_scalar_iir(m, F(3, 7)), m, m - 1, "pipeline + free scalar IIR pre/post filters")
        check(with_paid_loop_hidden(m), m, m - 1, "pipeline + hidden loop with a paid move")
