#!/usr/bin/python3
"""regmachine.py -- the game the Lean engine certifies, transcribed by hand from the Lean definitions.

Agent label: w-engine, 10 October 2026.  NOT Lean and not checked by Lean: a Python transcription, written to
test my reading of the definitions.  Each definition names the Lean file and line it copies.
REL = the third revision of this repository (Lean file and line references below are to it)

Two simulators of the same word language:
  * numeric  : `gwalk` on raw arrays for a given number f of digit columns (what Lean's GLiveKernel quantifies over);
  * formal   : every array holds a rational combination of terms "(a matrix of the group) applied to input l"
               (the register machine).  By the reduction checked in OUT/w-engine.md section 3 a word passes the
               numeric test for EVERY f exactly when it passes the formal test.

All matrix entries that occur are Gaussian numbers with a small power of 2 in the denominator, and the test data
are Gaussian integers, so double-precision arithmetic is exact here; equality is tested exactly.
"""
import sys
from fractions import Fraction

import numpy as np

I = 1j
CA = (1 + 1j) / 2          # REL/OAI/Computability/FourierTransform/BinaryFrames.lean:67   ca = (1+I)/2
CB = (1 - 1j) / 2          # BinaryFrames.lean:68                                           cb = (1-I)/2

# A vector of Space alpha = alpha -> ZMod 2 (BinaryBasis.lean:10-11) is a bit mask; coordinate i is bit i.


def dot(x, y):             # BinaryBasis.lean:26
    return bin(x & y).count("1") & 1


def shift(m, z):           # BinaryFrames.lean:137-138   shift z i j = 1 iff i + z = j
    n = 1 << m
    a = np.zeros((n, n), complex)
    for i in range(n):
        a[i, i ^ z] = 1
    return a


def dirm(m, z):            # BinaryFrames.lean:152       dir z = ca . 1 + cb . shift z
    return CA * np.eye(1 << m) + CB * shift(m, z)


def tint(b):               # BinaryFrames.lean:12        tint i = 1 if i = 0 else I
    return 1 if b == 0 else I


def phase_mat(m, z, c):    # REL/Work/GFrame/Engine/Letters.lean:42   diagonal (tint (dot z x + c))
    return np.diag([tint((dot(z, x) + c) & 1) for x in range(1 << m)]).astype(complex)


def linear_map(m, cols):   # an additive bijection of the address space (Letters.lean:31-33): x -> sum of cols[i]
    pi = []
    for x in range(1 << m):
        y = 0
        for i in range(m):
            if (x >> i) & 1:
                y ^= cols[i]
        pi.append(y)
    assert sorted(pi) == list(range(1 << m)), "not a bijection"
    return pi


def perm_mat(m, pi):       # Letters.lean:39              permMat G x y = 1 iff G.pi x = y
    n = 1 << m
    a = np.zeros((n, n), complex)
    for x in range(n):
        a[x, pi[x]] = 1
    return a


def kernel(m):             # BinaryFrames.lean:92 (kernel = digitProd fun _ => littleC), littleC at :69,
    n = 1 << m             # digitProd at ProductMatrix.lean:14-16 (entry = product over the coordinates)
    a = np.ones((n, n), complex)
    for x in range(n):
        for y in range(n):
            for i in range(m):
                a[x, y] *= CA if ((x >> i) & 1) == ((y >> i) & 1) else CB
    return a


def block_mat(m, zs):      # REL/Work/Block/Split.lean:125  blockMat z = (List.ofFn fun i => dir (z i)).prod
    a = np.eye(1 << m, dtype=complex)
    for z in zs:
        a = a @ dirm(m, z)
    return a


def f2_rank(vs):
    basis = []
    for v in vs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def mat_act(f, a, x):      # DirectionalWords.lean:17-18   matAct f M x = (digitProd (fun _ : Fin f => M)) *v x
    for ax in range(f):    # x has shape (2^m,)*f; the f-fold product of M acts on every axis
        x = np.moveaxis(np.tensordot(a, x, axes=([1], [ax])), 0, ax)
    return x


# ---- letters -------------------------------------------------------------------------------------------------
# ("mix", g)             BMove.localM g, g a rational R x R matrix        Word.lean:28    free, no condition
# ("shift", r, z)        BMove.shift r z                                  Word.lean:29    free, no condition
# ("block", r, [z...])   BMove.block r rk z ind                           Word.lean:30    paid: rank rk
# ("perm", r, cols)      GMove.perm r G, G additive bijection             Letters.lean:49 free
# ("phase", r, z, c)     GMove.phase r z c                                Letters.lean:50 free

def letter_matrix(m, letter):
    k = letter[0]
    if k == "shift":
        return shift(m, letter[2])
    if k == "block":
        return block_mat(m, letter[2])
    if k == "perm":
        return perm_mat(m, linear_map(m, letter[2]))
    if k == "phase":
        return phase_mat(m, letter[2], letter[3])
    raise ValueError(k)


def proper(m, word):       # Word.lean:55-57 and Letters.lean:66: a block needs 1 <= rk < u; plus `ind` (Word.lean:30)
    for letter in word:
        if letter[0] == "block":
            zs = letter[2]
            if not (1 <= len(zs) < m):
                return False
            if f2_rank(zs) != len(zs):
                return False
    return True


def ranks(word):           # the paid ranks, in order (Word.lean:44-46, Letters.lean:64)
    return [len(letter[2]) for letter in word if letter[0] == "block"]


def gwalk(m, f, nroles, word, data):   # Letters.lean:151-158; one array per role, shape (2^m,)*f
    data = list(data)
    for letter in word:
        if letter[0] == "mix":         # Word.lean:198: point f (actPoint g); DirectionalWords.lean:39-40, 64
            g = letter[1]
            data = [sum(complex(g[r][j]) * data[j] for j in range(nroles)) for r in range(nroles)]
        else:                          # roleAct f r M  (DirectionalWords.lean:31-32)
            r = letter[1]
            data[r] = mat_act(f, letter_matrix(m, letter), data[r])
    return data


def live_kernel_numeric(m, nroles, live, word, fs=(0, 1, 2, 3), trials=2, seed=1):
    """GLiveKernel (Letters.lean:176-178) tested on random Gaussian-integer inputs: scratch roles start at zero,
    every live role must end with the kernel applied to its own input.  Returns {f: True/False}."""
    rng = np.random.default_rng(seed)
    kk = kernel(m)
    out = {}
    for f in fs:
        ok = True
        shape = (1 << m,) * f
        for _ in range(trials):
            data = []
            for r in range(nroles):
                if r in live:
                    data.append(rng.integers(-3, 4, shape).astype(complex) + 1j * rng.integers(-3, 4, shape))
                else:
                    data.append(np.zeros(shape, complex))
            res = gwalk(m, f, nroles, word, data)
            for r in live:
                if not np.array_equal(res[r], mat_act(f, kk, data[r])):
                    ok = False
        out[f] = ok
    return out


# ---- the formal register machine -------------------------------------------------------------------------------
MATS = {}


def key(a):
    b = a * 8
    br, bi = np.rint(b.real), np.rint(b.imag)
    assert np.array_equal(br, b.real) and np.array_equal(bi, b.imag), "entry outside Z[i]/8"
    k = br.astype(np.int16).tobytes() + bi.astype(np.int16).tobytes()
    MATS.setdefault(k, a)
    return k


def formal_start(m, nroles, live):
    rows = [dict() for _ in range(nroles)]
    one = key(np.eye(1 << m, dtype=complex))
    for l, r in enumerate(live):
        rows[r][(l, one)] = Fraction(1)
    return rows


def formal_step(m, nroles, rows, letter):
    if letter[0] == "mix":
        g = letter[1]
        new = []
        for r in range(nroles):
            acc = {}
            for j in range(nroles):
                c = Fraction(g[r][j])
                if c:
                    for t, v in rows[j].items():
                        acc[t] = acc.get(t, 0) + c * v
            new.append({t: v for t, v in acc.items() if v != 0})
        return new
    r = letter[1]
    a = letter_matrix(m, letter)
    new = list(rows)
    new[r] = {(l, key(a @ MATS[k])): v for (l, k), v in rows[r].items()}
    return new


def formal_run(m, nroles, live, word, trace=False):
    rows = formal_start(m, nroles, live)
    states = [rows]
    for letter in word:
        rows = formal_step(m, nroles, rows, letter)
        states.append(rows)
    kk = key(kernel(m))
    ok = all(rows[r] == {(l, kk): Fraction(1)} for l, r in enumerate(live))
    return (ok, states) if trace else ok


def coherent(rows):
    """Is this state a state of the equal-frame game for SOME choice of start matrices S_l?
    That is: are there matrices M_r (one per non-empty row) and S_l with  row r = sum_l c_(r,l) [M_r S_l^-1] ?
    Returns (True, None) or (False, reason)."""
    rowmaps = []
    for r, row in enumerate(rows):
        d = {}
        for (l, k) in row:
            if l in d:
                return False, "array %d holds two different matrices applied to input %d" % (r, l)
            d[l] = MATS[k]
        rowmaps.append(d)
    s, mr = {}, {}
    todo = [r for r, d in enumerate(rowmaps) if d]
    while todo:
        progress = False
        for r in list(todo):
            d = rowmaps[r]
            known = [l for l in d if l in s]
            if not known:
                continue
            mr[r] = d[known[0]] @ s[known[0]]
            for l in d:
                cand = np.linalg.inv(d[l]) @ mr[r]
                if l in s:
                    if not np.allclose(s[l], cand):
                        return False, "array %d: its terms on inputs %d and %d need different start matrices" % (
                            r, known[0], l)
                else:
                    s[l] = cand
            todo.remove(r)
            progress = True
        if not progress:                      # a new connected piece: fix one start matrix freely
            r = todo[0]
            l0 = next(iter(rowmaps[r]))
            s[l0] = np.eye(rowmaps[r][l0].shape[0], dtype=complex)
    return True, None


def show_row(row, names):
    if not row:
        return "0"
    parts = []
    for (l, k), v in sorted(row.items(), key=lambda t: (t[0][0], names.get(t[0][1], "?"))):
        c = "" if v == 1 else ("-" if v == -1 else "(%s)" % v)
        parts.append("%s[%s]v%d" % (c, names.get(k, "g%04x" % (hash(k) & 0xffff)), l))
    return " + ".join(parts)
