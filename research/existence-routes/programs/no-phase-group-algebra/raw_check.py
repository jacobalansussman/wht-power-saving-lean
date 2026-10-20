#!/usr/bin/env python3
"""raw_check.py (e3-rmbuild): run the lifted calibration words of phi_check.py on RAW arrays at f = 1..6 digit
columns (block 3: arrays of 8^f complex entries, random inputs), and compare every delivering register with the kernel
applied column by column.  A matrix letter acts as the f-fold product M (x) ... (x) M, applied factor by factor."""
import sys
import numpy as np
sys.argv = [sys.argv[0]]
src = open("phi_check.py").read()
exec(src.split("bad = 0\nfor name, st")[0])          # builds E, tests, lift(), gname() ...


def apply_f(M, v, f):
    """apply M on every one of the f tensor factors of v (shape (8,)*f)."""
    for ax in range(f):
        v = np.moveaxis(np.tensordot(M, v, axes=([1], [ax])), 0, ax)
    return v


def run_raw(steps, W, R, f, rng):
    shape = (1 << m,) * f
    regs = [np.zeros(shape, dtype=complex) for _ in range(R)]
    ys = [rng.standard_normal(shape) + 1j * rng.standard_normal(shape) for _ in range(W)]
    for l in range(W):
        regs[l] = ys[l].copy()
    for st in steps:
        if st[0] == "move":
            regs[st[1]] = apply_f(E.dir[st[2]], regs[st[1]], f)
        elif st[0] == "free":
            regs[st[1]] = apply_f(E.letter(st[2]), regs[st[1]], f)
        elif st[0] == "mix":
            new = np.zeros(shape, dtype=complex)
            for s_, mu in st[2]:
                new = new + mu.cx() * regs[s_]
            regs[st[1]] = new
    return regs, ys


rng = np.random.default_rng(7)
worst = 0.0; bad = 0
for name, st, W, R, deliver, expect in tests:
    if not expect:
        continue
    word = lift(st, deliver, R - 1)
    for f in range(1, 7):
        regs, ys = run_raw(word, W, R, f, rng)
        err = max(np.max(np.abs(regs[r] - apply_f(E.K, ys[l], f))) / np.max(np.abs(ys[l]))
                  for l, r in enumerate(deliver))
        worst = max(worst, err)
        if err > 1e-9:
            bad += 1
            print("MISMATCH", name, f, err)
    print("%-48s f = 1..6: max relative error %.2e" % (name, worst))
print("mismatches:", bad)
