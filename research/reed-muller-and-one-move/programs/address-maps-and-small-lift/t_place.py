"""t_place.py: every placement of the Singer group relative to the target, at block m (rm-theory, claims F6, P1-P3).

Fix the Singer cycle sigma: x -> alpha x in self-dual coordinates. A placement is a target T (symmetric, invertible,
non-alternating m x m over F_2); conjugating H instead of moving T gives the same games (v-e3-rmbuild's sweep).
Lemma V with H = <sigma> (order n = 2^m - 1), M physical moves (any directions), W arrays:
   |D_Q| = M wt(Q)   (each physical move is n virtual moves covering every nonzero direction once),
   n_odd(Q) = #{k : <Q, sigma^-k T sigma^-kT> = 1},
   deletion: M wt(Q) >= W n_odd(Q) for every Q;  tight: equality with |D_Q| > 0.
Outputs, per placement class (orbits of Singer and Frobenius on targets): rho(T) = max_Q n_odd(Q)/wt(Q) (deletion
passes iff M/W >= rho), the number of forms violating / tight for the listed (W, M), and which class holds the trace
placement T = I and e3-rmbuild's polynomial-basis placement.
"""
import sys, itertools
from fractions import Fraction
import numpy as np
from gf import GF
from t_field import mons, rank2


def bits_of_sym(S, ML):
    v = 0
    for idx, (i, j) in enumerate(ML):
        if S[i][j]: v |= 1 << idx
    return v


def matmul2(A, B):
    m = len(A); p = len(B[0]); k = len(B)
    return [[sum(A[i][t] & B[t][j] for t in range(k)) & 1 for j in range(p)] for i in range(m)]


def transpose(A):
    return [list(r) for r in zip(*A)]


def main(m, shapes):
    F = GF(m); q, n = F.q, F.n
    fb = F.self_dual_basis(1)
    ML = mons(m); d = len(ML)
    co = {x: F.coords(x, fb) for x in range(q)}
    # Singer matrix (columns = coords of alpha * f_j) and its inverse power list
    def mult_matrix(a):
        cols = [co[F.mul(a, f)] for f in fb]
        return [[cols[j][i] for j in range(m)] for i in range(m)]
    Sinv = [mult_matrix(F.alpha(-k)) for k in range(n)]
    Frob = (lambda: [[co[F.mul(f, f)][i] for f in fb] for i in range(m)])()
    # all forms as ints; weights
    Z = np.array([[co[x][i] & co[x][j] for (i, j) in ML] for x in range(1, q)], dtype=np.int64)   # n x d
    zmask = np.array([int(sum(int(Z[r, c]) << c for c in range(d))) for r in range(n)], dtype=np.int64)
    Qs = np.arange(1 << d, dtype=np.int64)
    pop = np.zeros(1 << d, dtype=np.int8)
    for b in range(d): pop ^= ((np.arange(1 << d) >> b) & 1).astype(np.int8)    # parity table
    wt = np.zeros(1 << d, dtype=np.int64)
    for r in range(n): wt += pop[Qs & zmask[r]]
    # targets
    targets = []
    for bits in range(1 << d):
        S = [[0] * m for _ in range(m)]
        for idx, (i, j) in enumerate(ML):
            if (bits >> idx) & 1: S[i][j] = S[j][i] = 1
        if not any(S[i][i] for i in range(m)): continue           # alternating
        if rank2([int("".join(map(str, r)), 2) for r in S], m) != m: continue
        targets.append(bits)
    tset = set(targets)
    print(f"=== m = {m}: {len(targets)} targets (invertible non-alternating symmetric); forms 2^{d}; n = {n}")
    def sym_from_bits(bits):
        S = [[0] * m for _ in range(m)]
        for idx, (i, j) in enumerate(ML):
            if (bits >> idx) & 1: S[i][j] = S[j][i] = 1
        return S
    def act(Ag, bits):   # A T A^T
        S = sym_from_bits(bits)
        return bits_of_sym(matmul2(matmul2(Ag, S), transpose(Ag)), ML)
    # classes under Singer x Frobenius
    seen = {}; classes = []
    for t in targets:
        if t in seen: continue
        orb = set(); frontier = [t]
        while frontier:
            u = frontier.pop()
            if u in orb: continue
            orb.add(u)
            frontier.append(act(Sinv[1], u)); frontier.append(act(Frob, u))
        for u in orb: seen[u] = len(classes)
        classes.append(sorted(orb))
    I = [[int(i == j) for j in range(m)] for i in range(m)]
    Ibits = bits_of_sym(I, ML)
    # polynomial-basis placement: target = sum_j c(alpha^j) c(alpha^j)^T
    Pcols = [co[F.alpha(j)] for j in range(m)]
    Tp = [[0] * m for _ in range(m)]
    for v in Pcols:
        for i in range(m):
            for j in range(m): Tp[i][j] ^= v[i] & v[j]
    Tpbits = bits_of_sym(Tp, ML)
    print(f"   {len(classes)} classes under Singer x Frobenius; sizes {sorted(len(c) for c in classes)}")
    print(f"   trace placement T = I in class {seen[Ibits]}; e3-rmbuild polynomial placement in class {seen[Tpbits]}")
    tot_fail = {s: 0 for s in shapes}; tot_tight = {s: 0 for s in shapes}
    for ci, cl in enumerate(classes):
        T = cl[0]
        S = sym_from_bits(T)
        masks = [bits_of_sym(matmul2(matmul2(Sinv[k], S), transpose(Sinv[k])), ML) for k in range(n)]
        nodd = np.zeros(1 << d, dtype=np.int64)
        for mk in masks: nodd += pop[Qs & mk]
        nz = wt > 0
        ratios = nodd[nz] / wt[nz]
        rho_idx = np.argmax(ratios)
        qi = np.nonzero(nz)[0][rho_idx]
        rho = Fraction(int(nodd[qi]), int(wt[qi]))
        vals = sorted(set(nodd[1:].tolist()))
        line = f"   class {ci:3d} size {len(cl):4d}: rho = {rho} = {float(rho):.4f} (at wt {wt[qi]}, n_odd {nodd[qi]}); n_odd values {vals[:12]}{'...' if len(vals) > 12 else ''}"
        for (W, M) in shapes:
            lhs = M * wt; rhs = W * nodd
            fail = int(np.sum(lhs[1:] < rhs[1:])); tight = int(np.sum((lhs[1:] == rhs[1:]) & (nodd[1:] > 0)))
            line += f" | {W} in {M}: viol {fail} tight {tight}"
            if fail: tot_fail[(W, M)] += len(cl)
            if tight: tot_tight[(W, M)] += len(cl)
        print(line)
    for s in shapes:
        print(f"   {s[0]} in {s[1]}: targets excluded by deletion {tot_fail[s]} of {len(targets)} "
              f"({tot_fail[s]}/{len(targets)} = {tot_fail[s] / len(targets):.4f}); targets with a tight form {tot_tight[s]}")


if __name__ == "__main__":
    m = int(sys.argv[1])
    shapes = [tuple(map(int, a.split(","))) for a in sys.argv[2:]]
    main(m, shapes)
