#!/usr/bin/env python3
"""phi_check.py (e3-rmbuild): exact check of Theorem A's lift phi on real engine words (block 3), using the
theory guide's exact formal-content engine hb.py (contents in C[G], Gaussian-rational coefficients, group elements
as exact matrices).

A word of C[Sym_m x| H] (moves [zz^T], address maps) is lifted to an engine word: each move becomes dir(z) followed
by the free element c_z = ((1-i)[1] + (1+i)[X^z])/2 (copy to scratch, shift, mix), and at the end each delivering
register is multiplied by prod_j c_{e_j}^{-1} = prod_j c_{e_j} X^{e_j}.  The engine word must then deliver exactly
[K] y_l (one term, coefficient 1) in every delivering register.  Controls that must FAIL: the same words without the
correction (raw moves), and a direction set whose Gram matrix is not I.
"""
import sys
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import hb
from hb import GQ
from fractions import Fraction as Fr

m = 3
E = hb.Engine(m)
kK = E.key(E.K)
HALF_1mi = GQ(Fr(1, 2), Fr(-1, 2))
HALF_1pi = GQ(Fr(1, 2), Fr(1, 2))


def c_steps(r, s, z):
    return [("mix", s, [(r, GQ(1))]), ("free", s, "X%d" % z),
            ("mix", r, [(r, HALF_1mi), (s, HALF_1pi)]), ("mix", s, [])]


def lift(sym_steps, deliver, scratch, corrected=True):
    out = []
    for st in sym_steps:
        if st[0] == "smove":
            _, r, z = st
            out.append(("move", r, z))
            if corrected:
                out += c_steps(r, scratch, z)
        else:
            out.append(st)
    if corrected:
        for r in deliver:
            for j in range(m):
                out += c_steps(r, scratch, 1 << j)
                out.append(("free", r, "X%d" % (1 << j)))
    return out


def is_unit(steps, W, R, deliver):
    w = hb.Word(E, W, R)
    _, C, _ = w.run(steps)
    ok = True
    for l, r in enumerate(deliver):
        c = C[r]
        good = (len(c) == 1 and (l, kK) in c and c[(l, kK)].a == 1 and c[(l, kK)].b == 0)
        ok = ok and good
    return ok, [len(C[r]) for r in deliver]


def gname(cols):
    """address map name for hb: images of all x given images of basis vectors"""
    img = []
    for x in range(1 << m):
        y = 0
        for k in range(m):
            if (x >> k) & 1:
                y ^= cols[k]
        img.append(y)
    return "G" + "".join("%x" % t for t in img)


def inv_cols(cols):
    img = {}
    for x in range(1 << m):
        y = 0
        for k in range(m):
            if (x >> k) & 1:
                y ^= cols[k]
        img[y] = x
    return [img[1 << k] for k in range(m)]


tests = []
# T1 plain, one array
tests.append(("T1 plain e1 e2 e3 (W=1)", [("smove", 0, 1), ("smove", 0, 2), ("smove", 0, 4)], 1, 2, [0], True))
# T2 the four non-coordinate directions (Gram I), one array
tests.append(("T2 dirs 3,5,6,7 (W=1)", [("smove", 0, z) for z in (3, 5, 6, 7)], 1, 2, [0], True))
# T3 U_3 = 2 in 7 as a partition, two arrays
tests.append(("T3 U_3 2 in 7 (W=2)", [("smove", 0, z) for z in (1, 2, 4)] + [("smove", 1, z) for z in (3, 5, 6, 7)],
              2, 3, [0, 1], True))
# T4 address map: y -> P_h y, moves h e_j... wait: moves in directions z with Gram(h^-1 z) = I, then P_h^-1
cols = [3, 6, 7]           # h e_1 = 011, h e_2 = 110, h e_3 = 111 (invertible)
hi = inv_cols(cols)
T4 = [("free", 0, gname(cols))] + [("smove", 0, z) for z in cols] + [("free", 0, gname(hi))]
tests.append(("T4 address map h, moves h e_j, then h^-1 (W=1)", T4, 1, 2, [0], True))
# T5 address map with a non-orthonormal seen set: moves e1,e2,e3 seen through h^-1 (Gram h^-1 h^-T != I): must fail
T5 = [("free", 0, gname(cols))] + [("smove", 0, z) for z in (1, 2, 4)] + [("free", 0, gname(hi))]
tests.append(("T5 control: Gram != I (must fail)", T5, 1, 2, [0], False))
# T6 control: dirs 3,5,6 only (Gram != I): must fail
tests.append(("T6 control: dirs 3,5,6 (must fail)", [("smove", 0, z) for z in (3, 5, 6)], 1, 2, [0], False))


# T7 channels + mixing + address maps: u = y1 + P_h y2, w = y1 - P_h y2 (h = cyclic coordinate permutation, orthogonal),
# both moved by e1, e2, e3; readout (u + w)/2 = K' y1 and P_h^-1 (u - w)/2 = K' y2.  Registers: y1=0, y2=1, u=2, w=3, s=4.
ccols = [2, 4, 1]          # h e1 = e2, h e2 = e3, h e3 = e1
chi = inv_cols(ccols)
T7 = [("free", 1, gname(ccols)),
      ("mix", 2, [(0, GQ(1)), (1, GQ(1))]), ("mix", 3, [(0, GQ(1)), (1, GQ(-1))])]
T7 += [("smove", r, z) for r in (2, 3) for z in (1, 2, 4)]
T7 += [("mix", 0, [(2, GQ(Fr(1, 2))), (3, GQ(Fr(1, 2)))]), ("mix", 1, [(2, GQ(Fr(1, 2))), (3, GQ(Fr(-1, 2)))]),
       ("free", 1, gname(chi))]
tests.append(("T7 two channels, mixing, address map (W=2)", T7, 2, 5, [0, 1], True))
# T8 control: the same with a non-orthogonal h (Gram(h^-1 e_j) != I): array 2 must fail
T8 = [st if not (st[0] == "free" and st[2] == gname(ccols)) else ("free", 1, gname([3, 6, 7])) for st in T7]
T8 = [st if not (st[0] == "free" and st[2] == gname(chi)) else ("free", 1, gname(inv_cols([3, 6, 7]))) for st in T8]
tests.append(("T8 control: non-orthogonal h (must fail)", T8, 2, 5, [0, 1], False))

bad = 0
for name, st, W, R, deliver, expect in tests:
    ok, sizes = is_unit(lift(st, deliver, R - 1), W, R, deliver)
    raw_ok, _ = is_unit(lift(st, deliver, R - 1, corrected=False), W, R, deliver)
    flag = "OK" if ok == expect else "UNEXPECTED"
    if ok != expect:
        bad += 1
    print("%-48s lifted word unit: %-5s (expected %-5s) terms %s | raw word (no c_z) unit: %s  %s"
          % (name, ok, expect, sizes, raw_ok, flag))
print("unexpected:", bad)
