"""zx4_xw_sumfree.py -- a2-zx (meta2, 14 Oct 2026).  P1/P2 in the everything-free game, sum-free model.

Free wrappers now use XW letters on the tensor-encoded register (W arrays x m rows x f columns):
  CS in every column of one array (rows 0,1): i^(x0 x1) = T0 T1 CNOT T1^dag CNOT;
  Toffoli / CCZ in every column (block >= 3): CCZ on rows 0,1,2;
  chain map chi on one row of one array: column c ^= column c-1 (CNOTs from wire (r, c-1) to (r, c)), and E
  (row 0 col c ^= row 1 col c-1);
  plus the engine's CNOT, CZ, S, X, Z.
Paid: XPhase(1/2) on every wire (plain schedule).  Circuit U = A (I (x) K) B with A, B random XW wrappers.
The Lagrangian invariant does not apply (non-Clifford), but the support count does: every gate except HAD and
XPhase(odd/2) is monomial, so (#HAD + #sqrtX) >= log2(column support) = mWf for ANY circuit equal to U.
Check: PyZX full_reduce -> extract -> basic_optimization (XPhase expanded) never goes below mWf; tensor check n <= 10.
"""
import sys, random, statistics
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pyzx as zx
from pyzx.circuit import Circuit
from fractions import Fraction
from zx1_sumfree import noX, counts, plain


def wire(m, f, w, r, c):
    return (w * f + c) * m + r


def xw_layer(m, W, f, k, rng):
    n = m * W * f
    c = Circuit(n)
    for _ in range(k):
        t = rng.randrange(7)
        w = rng.randrange(W)
        if t == 0:      # CS in every column of array w
            for col in range(f):
                a, b = wire(m, f, w, 0, col), wire(m, f, w, 1, col)
                c.add_gate("ZPhase", a, Fraction(1, 4)); c.add_gate("ZPhase", b, Fraction(1, 4))
                c.add_gate("CNOT", a, b); c.add_gate("ZPhase", b, Fraction(7, 4)); c.add_gate("CNOT", a, b)
        elif t == 1 and m >= 3:   # CCZ in every column
            for col in range(f):
                c.add_gate("CCZ", wire(m, f, w, 0, col), wire(m, f, w, 1, col), wire(m, f, w, 2, col))
        elif t == 2 and f >= 2:   # chi on row r (in place: later columns first)
            r = rng.randrange(m)
            for col in range(f - 1, 0, -1):
                c.add_gate("CNOT", wire(m, f, w, r, col - 1), wire(m, f, w, r, col))
        elif t == 3 and f >= 2:   # E: row 0 col c ^= row 1 col c-1
            for col in range(f - 1, 0, -1):
                c.add_gate("CNOT", wire(m, f, w, 1, col - 1), wire(m, f, w, 0, col))
        elif t == 4:              # column-uniform CNOT between rows (possibly across arrays)
            w2 = rng.randrange(W); r1 = rng.randrange(m); r2 = rng.randrange(m)
            if (w2, r2) != (w, r1):
                for col in range(f): c.add_gate("CNOT", wire(m, f, w, r1, col), wire(m, f, w2, r2, col))
        elif t == 5:              # S on a row
            r = rng.randrange(m)
            for col in range(f): c.add_gate("ZPhase", wire(m, f, w, r, col), Fraction(1, 2))
        else:                     # X on a row
            r = rng.randrange(m)
            for col in range(f): c.add_gate("NOT", wire(m, f, w, r, col))
    return c


def main():
    rng = random.Random(1410)
    rows = []
    print("m W f n | case | T-count in -> after | HAD+sqrtX after full_reduce/extract/basic | mWf | tensor")
    for m in (2, 3, 4):
        for W in (1, 2):
            for f in (1, 2, 3):
                n = m * W * f
                if n > 24: continue
                for s in range(3):
                    A = xw_layer(m, W, f, 6, rng); B = xw_layer(m, W, f, 6, rng)
                    c = Circuit(n); c.add_circuit(B); c.add_circuit(plain(m, W, f)); c.add_circuit(A)
                    t0 = zx.tcount(c)
                    try:
                        g = c.to_graph(); zx.simplify.full_reduce(g); t1 = zx.tcount(g)
                        c2 = zx.optimize.basic_optimization(noX(zx.extract_circuit(g.copy())))
                        h, sx, ng = counts(c2)
                        ok = ""
                        if n <= 10: ok = "T" if zx.compare_tensors(c, c2, preserve_scalar=False) else "F"
                        rows.append((m, W, f, h + sx, m * W * f, ok))
                        print("%d %d %d %2d | xw%d | T %3d -> %3d | %3d | %3d | %s" % (m, W, f, n, s, t0, t1, h + sx, m * W * f, ok), flush=True)
                    except Exception as e:
                        print("%d %d %d %2d | xw%d | ERR %s" % (m, W, f, n, s, str(e)[:60]), flush=True)
    below = [r for r in rows if r[3] < r[4]]
    bad = [r for r in rows if r[5] == "F"]
    ratios = [r[3] / r[4] for r in rows]
    print("\nSUMMARY: runs %d; below mWf: %d; tensor mismatches: %d; ratio min %.3f median %.3f max %.3f"
          % (len(rows), len(below), len(bad), min(ratios), statistics.median(ratios), max(ratios)))


if __name__ == "__main__":
    main()
