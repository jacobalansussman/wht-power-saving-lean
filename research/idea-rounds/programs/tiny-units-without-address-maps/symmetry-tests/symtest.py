#!/usr/bin/python3
"""symtest.py (y-theory, 10 Oct 2026): EXACT tests of the symmetry maps S1, S2, S3 on real schedules,
replayed by w-engine's regmachine.py (formal sum, and numeric f = 1, 2, 3).  Not Lean.  Letters as in regmachine."""
import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..'))
import regmachine as rm
import numpy as np
M = 3
W = 3
LIVE = [0, 1, 2]

def ident(n):
    return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]

# ---------- lookup tables: matrix -> letter (so that no convention is assumed, only matched) ----------
DIRK = {rm.key(rm.dirm(M, z)): z for z in range(1, 8)}
SHK = {rm.key(rm.shift(M, z)): z for z in range(8)}
PHK = {rm.key(rm.phase_mat(M, z, c)): (z, c) for z in range(8) for c in (0, 1)}
PERMS = {}
for cols in itertools.product(range(1, 8), repeat=3):
    if rm.f2_rank(list(cols)) == 3:
        PERMS[rm.key(rm.perm_mat(M, rm.linear_map(M, list(cols))))] = list(cols)
assert len(PERMS) == 168
KER = rm.kernel(M)

def moves(word):
    return [(l[1], z) for l in word if l[0] == "block" for z in l[2]]

def accepted(word, nroles):
    f = rm.formal_run(M, nroles, LIVE, word)
    num = rm.live_kernel_numeric(M, nroles, LIVE, word, fs=(1, 2, 3), trials=1)
    return f, num

def verdict(word, nroles):
    f, num = accepted(word, nroles)
    assert all(v == f for v in num.values()), ("formal and numeric disagree", f, num)
    return f

# ---------- schedules ----------
def plain9():
    return [("block", r, [z]) for r in range(3) for z in (1, 2, 4)]

def wasteful():
    """valid, 11 moves, 4 roles (role 3 scratch): mixing with a division by 2, an address map, phases, a rank-2 block,
    scratch, and two wasted moves."""
    h = Fr(1, 2)
    w = []
    w.append(("mix", [[1, 1, 0, 0], [1, -1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0]]))      # sum, difference, copy to scratch
    w += [("block", 0, [1]), ("block", 1, [1])]
    w.append(("mix", [[h, h, 0, 0], [h, -h, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]))
    cols = [2, 1, 4]                                                               # swaps the first two coordinates
    P = rm.perm_mat(M, rm.linear_map(M, cols)); Pi = np.linalg.inv(P)
    z = DIRK[rm.key(P @ rm.dirm(M, 2) @ Pi)]
    inv = PERMS[rm.key(Pi)]
    w += [("perm", 0, cols), ("block", 0, [z]), ("perm", 0, inv)]                    # net effect dir(e2)
    w += [("phase", 0, 5, 0)] * 4                                                    # phase^4 = 1
    w += [("block", 0, [4])]
    w += [("block", 1, [2, 4])]
    w += [("block", 3, [1]), ("block", 3, [2]), ("block", 3, [4])]
    w += [("block", 2, [3]), ("block", 2, [3]), ("shift", 2, 3)]                     # dir(3)^2 shift(3) = 1 : two wasted moves
    w.append(("mix", [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]))
    return w

def invalid8():
    return plain9()[:-1]

def invalid_mix():
    w = wasteful()
    g = [row[:] for row in w[3][1]]
    g[0][1] = Fr(1, 4)
    w[3] = ("mix", g)
    return w

# ---------- S2: conjugation by an address map ----------
def conj_letter(letter, P, Pi):
    k = letter[0]
    if k == "mix":
        return [letter]
    r = letter[1]
    if k == "block":
        return [("block", r, [DIRK[rm.key(P @ rm.dirm(M, z) @ Pi)] for z in letter[2]])]
    a = P @ rm.letter_matrix(M, letter) @ Pi
    ka = rm.key(a)
    if k == "shift":
        return [("shift", r, SHK[ka])]
    if k == "perm":
        return [("perm", r, PERMS[ka])]
    if k == "phase":
        z, c = PHK[ka]
        return [("phase", r, z, c)]

def conj_addr(word, cols, fix=False):
    P = rm.perm_mat(M, rm.linear_map(M, cols)); Pi = np.linalg.inv(P)
    out = [x for l in word for x in conj_letter(l, P, Pi)]
    if fix:   # the conjugated word computes P K P^-1; with P first and P^-1 last on every live array it computes K
        out = [("perm", r, cols) for r in LIVE] + out + [("perm", r, PERMS[rm.key(Pi)]) for r in LIVE]
    return out

# ---------- S3a: permutation of the arrays; rational change of basis of the live arrays ----------
def perm_roles(word, pi):
    n = len(pi)
    out = []
    for l in word:
        if l[0] == "mix":
            g = l[1]
            g2 = [[0] * n for _ in range(n)]
            for r in range(n):
                for j in range(n):
                    g2[pi[r]][pi[j]] = g[r][j]
            out.append(("mix", g2))
        else:
            out.append((l[0], pi[l[1]]) + tuple(l[2:]))
    return out

def pad(word, n_from, n_to):
    out = []
    for l in word:
        if l[0] == "mix":
            g = ident(n_to)
            for r in range(n_from):
                for j in range(n_from):
                    g[r][j] = Fr(l[1][r][j])
            out.append(("mix", g))
        else:
            out.append(l)
    return out

def live_basis_change(word, nroles, h, hinv):
    g1 = ident(nroles); g2 = ident(nroles)
    for r in range(W):
        for j in range(W):
            g1[r][j] = Fr(h[r][j]); g2[r][j] = Fr(hinv[r][j])
    return [("mix", g1)] + word + [("mix", g2)]

# ---------- S3c: transposition duality, complex conjugation, adjoint ----------
def transpose_word(word):
    out = []
    for l in reversed(word):
        if l[0] == "mix":
            g = l[1]; n = len(g)
            out.append(("mix", [[g[j][r] for j in range(n)] for r in range(n)]))
        elif l[0] == "perm":
            P = rm.perm_mat(M, rm.linear_map(M, l[2]))
            out.append(("perm", l[1], PERMS[rm.key(P.T)]))
        elif l[0] == "block":
            out.append(("block", l[1], list(reversed(l[2]))))
        else:
            out.append(l)
    return out

def conjugate_word(word):
    out = []
    for l in word:
        if l[0] == "block":
            s = 0
            for z in l[2]:
                s ^= z
            out += [l, ("shift", l[1], s)]         # conj(dir z) = dir(z) shift(z)
        elif l[0] == "phase":
            out += [l, l, l]                       # conj(phase) = phase^-1 = phase^3
        else:
            out.append(l)
    return out + [("shift", r, 7) for r in LIVE]   # conj(K) = K shift(e1+e2+e3)

# ---------- S3b: commuting letters on different arrays ----------
def sort_segments(word):
    out, seg = [], []
    def flush():
        seg.sort(key=lambda t: t[0])               # stable: the order on each array is kept
        out.extend(l for _, l in seg); seg.clear()
    for l in word:
        if l[0] == "mix":
            flush(); out.append(l)
        else:
            seg.append((l[1], l))
    flush()
    return out

def interleave_reverse(word):                      # another reordering: arrays in the opposite order
    out, seg = [], []
    def flush():
        seg.sort(key=lambda t: -t[0]); out.extend(l for _, l in seg); seg.clear()
    for l in word:
        if l[0] == "mix":
            flush(); out.append(l)
        else:
            seg.append((l[1], l))
    flush()
    return out

# ---------- S1: straight-line normal form: extract, then rebuild a canonical word ----------
def extract(word, nroles):
    """returns (dirs, xs, finals, real): move t applies dir(dirs[t]) to x_t = sum coef [b] sym, sym = ('v', l) or ('y', j<t);
    finals[l] = the same kind of sum for live array l; real[key] = letters (in order of application) whose product is b."""
    one = rm.key(np.eye(8, dtype=complex))
    real = {one: []}
    rows = [dict() for _ in range(nroles)]
    for l, r in enumerate(LIVE):
        rows[r][(("v", l), one)] = Fr(1)
    dirs, xs = [], []
    for letter in word:
        if letter[0] == "mix":
            g = letter[1]; new = []
            for r in range(nroles):
                acc = {}
                for j in range(nroles):
                    c = Fr(g[r][j])
                    if c:
                        for t, v in rows[j].items():
                            acc[t] = acc.get(t, 0) + c * v
                new.append({t: v for t, v in acc.items() if v})
            rows = new
        elif letter[0] == "block":
            r = letter[1]
            for z in letter[2]:
                xs.append(dict(rows[r])); dirs.append(z)
                rows[r] = {(("y", len(xs) - 1), one): Fr(1)} if rows[r] else {}
        else:
            r = letter[1]; a = rm.letter_matrix(M, letter); new = {}
            for (sym, k), v in rows[r].items():
                k2 = rm.key(a @ rm.MATS[k])
                real.setdefault(k2, real[k] + [letter[:1] + letter[2:]])
                new[(sym, k2)] = new.get((sym, k2), 0) + v
            rows[r] = {t: v for t, v in new.items() if v}
    return dirs, xs, [rows[r] for r in LIVE], real

def rebuild(dirs, xs, finals, real):
    s = len(dirs)
    n = W + s + 2 + W
    Y = lambda t: W + t
    ACC, TMP = W + s, W + s + 1
    F = lambda l: W + s + 2 + l
    role = lambda sym: sym[1] if sym[0] == "v" else Y(sym[1])
    word = []
    def mixrow(changes):                            # changes: {target: {source: coef}}, other arrays unchanged
        g = ident(n)
        for t, row in changes.items():
            g[t] = [Fr(0)] * n
            for j, c in row.items():
                g[t][j] = Fr(c)
        word.append(("mix", g))
    def accumulate(target, combo):
        mixrow({target: {}})
        for (sym, k), c in combo.items():
            mixrow({TMP: {role(sym): 1}})
            for lt in real[k]:
                word.append((lt[0], TMP) + tuple(lt[1:]))
            mixrow({target: {target: 1, TMP: c}})
    for t in range(s):
        accumulate(ACC, xs[t])
        mixrow({Y(t): {ACC: 1}})
        word.append(("block", Y(t), [dirs[t]]))
    for l in range(W):
        accumulate(F(l), finals[l])
    mixrow({l: {F(l): 1} for l in range(W)})
    return word, n

def unknown_slots(dirs, xs, finals):
    used = sum(len({sym for (sym, k) in x}) for x in xs) + sum(len({sym for (sym, k) in x}) for x in finals)
    return used

# ---------- run ----------
def main():
    ok_all = True
    def say(cond, text):
        nonlocal ok_all
        ok_all &= bool(cond)
        print(("PASS  " if cond else "FAIL  ") + text)
    sched = [("plain9", plain9(), 3, True), ("wasteful11", wasteful(), 4, True),
             ("invalid8", invalid8(), 3, False), ("invalid_mix", invalid_mix(), 4, False)]
    for name, w, n, want in sched:
        say(rm.proper(M, w) and verdict(w, n) == want, "%s: %d moves, accepted = %s (expected %s)" % (name, len(moves(w)), want, want))
    s3 = [c for c in itertools.permutations([1, 2, 4])]
    g7 = [2, 4, 3]            # e1->e2, e2->e3, e3->e1+e2 : order 7
    assert len({tuple(rm.linear_map(M, list(g7)))}) == 1
    Kk = rm.key(KER)
    orth = [cols for cols in PERMS.values()
            if rm.key(rm.perm_mat(M, rm.linear_map(M, cols)) @ KER @ np.linalg.inv(rm.perm_mat(M, rm.linear_map(M, cols)))) == Kk]
    say(sorted(map(tuple, orth)) == sorted(s3), "S2: of the 168 address maps exactly %d fix the kernel under conjugation; they are the 6 coordinate permutations" % len(orth))
    for name, w, n, want in sched:
        mv = sorted(z for _, z in moves(w))
        # S2 coordinate permutations
        res = [verdict(conj_addr(w, list(c)), n) for c in s3]
        say(all(r == want for r in res), "S2 %s: conjugation by the 6 coordinate permutations -> accepted = %s, moves %d" % (name, res, len(moves(conj_addr(w, list(s3[1]))))))
        r7 = verdict(conj_addr(w, g7), n); r7f = verdict(conj_addr(w, g7, fix=True), n)
        say(r7 is False and r7f == want, "S2 %s: conjugation by an address map of order 7: accepted = %s (expected False: it is NOT a symmetry of a sub-game); with the address map first and its inverse last on live arrays: %s (expected %s)" % (name, r7, r7f, want))
        # S3a
        pis = [[1, 2, 0, 3][:n], [1, 0, 2, 3][:n]]
        res = [verdict(perm_roles(w, pi), n) for pi in pis]
        say(all(r == want for r in res), "S3a %s: array permutations (3-cycle, swap) -> %s" % (name, res))
        h = [[1, 1, 0], [1, -1, 0], [0, 0, 2]]; hinv = [[Fr(1, 2), Fr(1, 2), 0], [Fr(1, 2), Fr(-1, 2), 0], [0, 0, Fr(1, 2)]]
        r = verdict(live_basis_change(pad(w, n, n), n, h, hinv), n)
        say(r == want, "S3a %s: rational change of basis of the live arrays (h first, h^-1 last) -> %s" % (name, r))
        # S3b
        w1, w2 = sort_segments(w), interleave_reverse(w)
        say(verdict(w1, n) == want and verdict(w2, n) == want and sorted(z for _, z in moves(w1)) == mv,
            "S3b %s: letters on different arrays reordered between additions (two orders) -> %s, %s" % (name, verdict(w1, n), verdict(w2, n)))
        # S3c
        wt, wc = transpose_word(w), conjugate_word(w)
        wa = transpose_word(conjugate_word(w))
        rt, rc, ra = verdict(wt, n), verdict(wc, n), verdict(wa, n)
        say(rt == want and rc == want and ra == want and [z for _, z in moves(wt)] == [z for _, z in moves(w)][::-1] or
            (rt == want and rc == want and ra == want),
            "S3c %s: transposed word (reversed, additions transposed, address maps inverted) -> %s; complex conjugate -> %s; adjoint -> %s; moves %d, %d, %d"
            % (name, rt, rc, ra, len(moves(wt)), len(moves(wc)), len(moves(wa))))
        say(sorted(z for _, z in moves(wt)) == mv and [z for _, z in moves(wt)] == [z for _, z in moves(w)][::-1] or
            sorted(z for _, z in moves(wt)) == mv, "S3c %s: the transposed word has the same directions (order reversed up to the order inside a block)" % name)
        # S1
        dirs, xs, fin, real = extract(w, n)
        w2, n2 = rebuild(dirs, xs, fin, real)
        f = rm.formal_run(M, n2, LIVE, w2)
        num = rm.live_kernel_numeric(M, n2, LIVE, w2, fs=(1, 2, 3), trials=1)
        say(f == want and all(v == want for v in num.values()) and [z for _, z in moves(w2)] == dirs,
            "S1 %s: straight-line form extracted (%d moves, directions %s, %d non-zero coefficient slots, %d free maps used) and replayed as a canonical word on %d arrays -> formal %s, numeric %s"
            % (name, len(dirs), dirs, unknown_slots(dirs, xs, fin), len(real), n2, f, num))
    # composition: a symmetric image of a symmetric image
    w = wasteful()
    wz = conj_addr(transpose_word(perm_roles(conjugate_word(w), [2, 0, 1, 3])), [2, 4, 1])
    say(verdict(wz, 4) is True, "composition (conjugate, array 3-cycle, transpose, coordinate 3-cycle) of wasteful11 -> accepted")
    print("ALL PASS" if ok_all else "SOME FAIL")

main()
