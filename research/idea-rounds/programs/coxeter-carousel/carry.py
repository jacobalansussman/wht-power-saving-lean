"""Carry identities behind the Hamming-weight-phasing verdict (track wildcard, wht44/ideas, 16 Oct 2026).

In the transform basis (conjugate by J = normalised Walsh matrix on all m*f bits), OpenAI's kernel C becomes the
phase S = diag(1, i).  So, with y in F_2^(m x f) (row j = the f bits of direction j):
  paid move along j  ->  i^{w_j(y)},   w_j = weight of row j
  target C^(x)mf     ->  i^{|y|}
Hamming-weight phasing writes i^w = i^(w mod 2) * (-1)^(floor(w/2) mod 2) and (floor(w/2) mod 2) = e_2 mod 2
(the carry).  We check exactly, on every y, for m = 2, 3 and f = 1..4:
  (a) J C_j J = diag(i^{w_j})                                   (numerically, on real arrays)
  (b) i^{w_j} = [(1+i)/2 + (1-i)/2 (-1)^{p_j}] * (-1)^{e2(row j)}   (p_j = parity of row j; first factor = 2 Paulis)
  (c) i^{|y|} = [(1+i)/2 + (1-i)/2 (-1)^{p}] * (-1)^{e2(all bits)}
  (d) (-1)^{e2(all)} = prod_j (-1)^{e2(row j)} * prod_{j<j'} (-1)^{p_j p_j'}, and each (-1)^{p_j p_j'} is the
      combination (1 + Z_j + Z_j' - Z_j Z_j')/2 of four free Paulis.
So every paid move is (free) x (row carry) and the target is (free) x (product of the m row carries) x (free).
"""
import itertools
import numpy as np
from coxeter import apply_dir

def e2(bits):
    w = int(sum(bits))
    return (w * (w - 1) // 2) % 2

def check(m, f):
    n = m * f
    ok = {"b": True, "c": True, "d": True}
    for y in itertools.product([0, 1], repeat=n):
        rows = [y[j * f:(j + 1) * f] for j in range(m)]
        p = [sum(r) % 2 for r in rows]
        for j in range(m):
            w = sum(rows[j])
            lhs = 1j ** w
            rhs = ((1 + 1j) / 2 + (1 - 1j) / 2 * (-1) ** p[j]) * (-1) ** e2(rows[j])
            ok["b"] &= abs(lhs - rhs) < 1e-12
        W = sum(y)
        lhs = 1j ** W
        rhs = ((1 + 1j) / 2 + (1 - 1j) / 2 * (-1) ** (W % 2)) * (-1) ** e2(y)
        ok["c"] &= abs(lhs - rhs) < 1e-12
        prod = 1
        for j in range(m):
            prod *= (-1) ** e2(rows[j])
        for j in range(m):
            for k in range(j + 1, m):
                Zj, Zk = (-1) ** p[j], (-1) ** p[k]
                pauli4 = (1 + Zj + Zk - Zj * Zk) / 2
                ok["d"] &= abs(pauli4 - (-1) ** (p[j] * p[k])) < 1e-12
                prod *= pauli4
        ok["d"] &= abs(prod - (-1) ** e2(y)) < 1e-12
    return ok

def check_a(m, f, rng):
    H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    C = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])
    shape = (2 ** f,) * m
    x = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    worst = 0.0
    for j in range(m):
        t = x.copy()
        for k in range(m):
            t = apply_dir(t, H, k, m, f)
        t = apply_dir(t, C, j, m, f)
        for k in range(m):
            t = apply_dir(t, H, k, m, f)
        # expected: multiply x[y] by i^{w_j(y)}; index of direction j is the axis j value (f bits)
        wj = np.array([bin(v).count("1") for v in range(2 ** f)])
        phase = (1j ** wj).reshape([2 ** f if a == j else 1 for a in range(m)])
        worst = max(worst, np.max(np.abs(t - phase * x)))
    return worst

if __name__ == "__main__":
    rng = np.random.default_rng(5)
    out = []
    for m in (2, 3):
        for f in (1, 2, 3, 4):
            ok = check(m, f)
            out.append(f"m={m} f={f}: (a) max error {check_a(m, f, rng):.1e}; (b) {ok['b']}; (c) {ok['c']}; (d) {ok['d']}")
    txt = "\n".join(out)
    print(txt)
    open("out_carry.txt", "w").write(txt + "\n")
