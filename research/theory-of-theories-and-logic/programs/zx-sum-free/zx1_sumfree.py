"""zx1_sumfree.py -- a2-zx (meta2, 14 Oct 2026).  Prediction P2: calibration of PyZX on SUM-FREE diagrams.

Tensor encoding: W arrays, block m, f columns -> n = m*W*f qubits.  The paid move D_z (one direction, every column,
one array) is the X-spider with phase pi/2 (sqrt X) on f wires; in PyZX: XPhase(q, 1/2) on each of the f wires.
Plain schedule = one XPhase(1/2) on every wire: mW paid moves, mWf sqrt-X gates.

Invariant (Bruhat / Lagrangian distance; the frame distance of the dual polar graph):
   d(U) = rank_F2 of the X-part of {U Z_q U^dag : q}  =  n - dim(L_Z cap U L_Z).
Free (Z-frame-preserving) gates leave d unchanged on either side; one sqrt-X or H gate changes d by at most 1.
So paid moves >= d(U)/f for any sum-free circuit, and d(I_W (x) K) = mWf.

For each (m, W, f): plain circuit, and circuits A (I (x) K) B with random free wrappers A, B (CNOT, CZ, S, X, Z on
random wires).  Then PyZX passes: (i) basic_optimization; (ii) full_reduce -> extract_circuit -> basic_optimization;
(iii) teleport_reduce -> to circuit -> basic_optimization.  Report HAD counts, sqrt-X counts, d before/after, and a
tensor check (compare_tensors) when n <= 10.
"""
import sys, random, time
import numpy as np
import pyzx as zx
from pyzx.circuit import Circuit
from pyzx.circuit.gates import HAD, ZPhase, XPhase, CNOT, CZ, NOT, Z, S, SWAP
from fractions import Fraction


def tableau_d(circ):
    """rank over F2 of the X-part of the images of Z_q under the circuit (signs ignored)."""
    n = circ.qubits
    # rows: images of Z_q; columns x[0..n), z[0..n)
    x = np.zeros((n, n), dtype=np.uint8)
    z = np.eye(n, dtype=np.uint8)
    for g in circ.to_basic_gates().gates:
        nm = g.name
        if nm == "HAD":
            q = g.target; x[:, q], z[:, q] = z[:, q].copy(), x[:, q].copy()
        elif nm in ("ZPhase", "S", "Z", "T"):
            ph = Fraction(g.phase) % 2 if hasattr(g, "phase") else None
            if ph is None:
                ph = {"S": Fraction(1, 2), "Z": Fraction(1)}[nm]
            if ph.denominator > 2: raise ValueError("non-Clifford phase")
            if ph.denominator == 2: z[:, g.target] ^= x[:, g.target]
        elif nm in ("XPhase", "NOT"):
            ph = Fraction(g.phase) % 2 if nm == "XPhase" else Fraction(1)
            if ph.denominator > 2: raise ValueError("non-Clifford phase")
            if ph.denominator == 2: x[:, g.target] ^= z[:, g.target]
        elif nm == "CNOT":
            c, t = g.control, g.target
            x[:, t] ^= x[:, c]; z[:, c] ^= z[:, t]
        elif nm == "CZ":
            a, b = g.control, g.target
            z[:, a] ^= x[:, b]; z[:, b] ^= x[:, a]
        elif nm == "SWAP":
            a, b = g.control, g.target
            x[:, [a, b]] = x[:, [b, a]]; z[:, [a, b]] = z[:, [b, a]]
        else:
            raise ValueError("gate " + nm)
    return rank2(x)


def rank2(M):
    M = M.copy() % 2; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]: piv = i; break
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
        if r == rows: break
    return r


def free_layer(n, k, rng):
    c = Circuit(n)
    for _ in range(k):
        t = rng.randrange(5)
        a = rng.randrange(n); b = rng.randrange(n)
        while b == a and n > 1: b = rng.randrange(n)
        if t == 0 and n > 1: c.add_gate("CNOT", a, b)
        elif t == 1 and n > 1: c.add_gate("CZ", a, b)
        elif t == 2: c.add_gate("ZPhase", a, Fraction(1, 2))
        elif t == 3: c.add_gate("NOT", a)
        else: c.add_gate("ZPhase", a, Fraction(1))
    return c


def plain(m, W, f):
    n = m * W * f
    c = Circuit(n)
    for q in range(n): c.add_gate("XPhase", q, Fraction(1, 2))
    return c


def counts(c):
    b = c.to_basic_gates()
    h = sum(1 for g in b.gates if g.name == "HAD")
    sx = sum(1 for g in b.gates if g.name == "XPhase" and Fraction(g.phase) % 1 == Fraction(1, 2))
    return h, sx, len(b.gates)


def noX(c):
    """expand XPhase(a) as HAD ZPhase(a) HAD: PyZX 0.10.7 basic_optimization returned an UNEQUAL circuit on one
    input containing XPhase gates (out_zx1_v1.txt, row m=3 W=1 f=1 scr2, teleport pass; reproduced in zx1_dbg2.py);
    with XPhase expanded first it is correct."""
    o = Circuit(c.qubits)
    for g in c.to_basic_gates().gates:
        if g.name == "XPhase":
            o.add_gate("HAD", g.target); o.add_gate("ZPhase", g.target, g.phase); o.add_gate("HAD", g.target)
        else: o.add_gate(g)
    return o


def run_pass(name, c):
    if name == "basic":
        return zx.optimize.basic_optimization(noX(c))
    if name == "full":
        g = c.to_graph(); zx.simplify.full_reduce(g)
        c2 = zx.extract_circuit(g.copy())
        return zx.optimize.basic_optimization(noX(c2))
    if name == "teleport":
        g = c.to_graph(); g2 = zx.simplify.teleport_reduce(g)
        c2 = Circuit.from_graph(g2)
        return zx.optimize.basic_optimization(noX(c2))
    raise ValueError(name)


def main():
    rng = random.Random(20261014)
    print("m W f n | case | d_before | pass: HAD sqrtX d_after tensor_ok")
    rows = []
    for m in (2, 3, 4):
        for W in (1, 2, 3):
            for f in (1, 2, 3):
                n = m * W * f
                if n > 36: continue
                cases = [("plain", plain(m, W, f))]
                for s in range(3):
                    A = free_layer(n, 4 * n, rng); B = free_layer(n, 4 * n, rng)
                    c = Circuit(n); c.add_circuit(B); c.add_circuit(plain(m, W, f)); c.add_circuit(A)
                    cases.append(("scr%d" % s, c))
                for cname, c in cases:
                    d0 = tableau_d(c)
                    line = ["%d %d %d %2d | %-5s | d=%3d (mWf=%3d) |" % (m, W, f, n, cname, d0, m * W * f)]
                    for pname in ("basic", "full", "teleport"):
                        t0 = time.time()
                        try:
                            c2 = run_pass(pname, c)
                            h, sx, ng = counts(c2)
                            d1 = tableau_d(c2)
                            ok = ""
                            if n <= 10:
                                ok = "T" if zx.compare_tensors(c, c2, preserve_scalar=False) else "F"
                            line.append(" %s: H=%3d sX=%3d d=%3d %s" % (pname, h, sx, d1, ok))
                            rows.append((m, W, f, cname, pname, d0, h, sx, d1, ok))
                        except Exception as e:
                            line.append(" %s: ERR %s" % (pname, str(e)[:40]))
                    print("".join(line), flush=True)
    # summary
    bad = [r for r in rows if r[8] != r[5] or r[5] != r[0] * r[1] * r[2]]
    below = [r for r in rows if r[6] + r[7] < r[0] * r[1] * r[2]]
    tens = [r for r in rows if r[9] == "F"]
    print("\nSUMMARY: rows %d; invariant d != mWf or changed: %d; (HAD + sqrtX) < mWf: %d; tensor mismatches: %d"
          % (len(rows), len(bad), len(below), len(tens)))
    import statistics
    for pname in ("basic", "full", "teleport"):
        ratios = [(r[6] + r[7]) / (r[0] * r[1] * r[2]) for r in rows if r[4] == pname and r[3] != "plain"]
        rp = [(r[6] + r[7]) / (r[0] * r[1] * r[2]) for r in rows if r[4] == pname and r[3] == "plain"]
        print("  %-8s scrambled: (HAD+sqrtX)/mWf min %.3f median %.3f max %.3f | plain: min %.3f max %.3f"
              % (pname, min(ratios), statistics.median(ratios), max(ratios), min(rp), max(rp)))


if __name__ == "__main__":
    main()
