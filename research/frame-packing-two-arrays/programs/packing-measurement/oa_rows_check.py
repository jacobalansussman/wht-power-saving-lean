#!/usr/bin/env python3
"""fp-measure: exact scalar replay of ONE invocation of OpenAI's exchange network per stage (explicit-DFT preprint,
rows (2.3); stage 2 is the inverse invocation from Y to X, rows 7..0 with negated updates, as Lemma 2.1 says), with
arbitrary (symbolic) auxiliary values, at a small ground set h.  It checks:
- Lemma 2.1: (x, y, A, c) -> (x, y + x, A, c) per invocation, auxiliaries restored; over the stages (x, y) ->
  (x, y + x) -> (-y, y + x) -> (-y, x);
- the raw classes held by every wire on every edge between two touches are exactly the lists oa_ledger.py uses
  (x inputs: one class per triple; y inputs and all auxiliary dirt: the single identity-twist class 'Y').
usage: python3 -I oa_rows_check.py [h]"""
import sys, itertools
from fractions import Fraction as Fr

h = int(sys.argv[1]) if len(sys.argv) > 1 else 7
T = list(itertools.combinations(range(h), 3))
nbr = lambda S, U: S != U and len(set(S) & set(U)) % 2 == 0
NB = {S: [U for U in T if nbr(S, U)] for S in T}
cST = lambda S, U: Fr(len(set(S) & set(U)) - 1, 2)
CEN = list(range(h)) + ["*"]
memb = {k: [U for U in T if k in U] for k in range(h)}
memb["*"] = list(T)


def add(a, b, k=Fr(1)):
    for kk, u in b.items():
        a[kk] = a.get(kk, 0) + k * u
        if a[kk] == 0:
            del a[kk]


def cls(vec):
    return frozenset("Y" if kk[0] != "x" else kk[1] for kk in vec)


def invocation(xs, ys, inverse):
    """xs, ys: dict triple -> content (dict).  Forward: logical source = X bank, target = Y bank.  Inverse: logical
    source = Y bank, target = X bank, rows 7..0 negated.  Returns final (X, Y) and per-wire edge class lists."""
    A = {(S, U): {("a", S, U): Fr(1)} for S in T for U in NB[S]}         # A[(target index, source index)]
    c = {k: {("c", k): Fr(1)} for k in CEN}
    A0 = {k: dict(v) for k, v in A.items()}
    c0 = {k: dict(v) for k, v in c.items()}
    src, tgt = (xs, ys) if not inverse else (ys, xs)
    touches = {}                                                       # wire -> list of class sets after touches
    def touch(w, vec):
        touches.setdefault(w, []).append(cls(vec))
    def row(r, sign):
        if r == 0 or r == 5:          # tgt -/+ JA ; (JA)_S = - sum_U c_SU A_SU
            for S in T:
                for U in NB[S]:
                    add(tgt[S], A[(S, U)], (1 if r == 0 else -1) * sign * cST(S, U))
                touch(("tgt", S), tgt[S])
                for U in NB[S]:
                    touch(("A", S, U), A[(S, U)])
        elif r == 1 or r == 4:        # tgt -/+ Rc ; (Rc)_S = (sum_{j in S} c_j - c_*)/2
            for S in T:
                for j in S:
                    add(tgt[S], c[j], (-1 if r == 1 else 1) * sign * Fr(1, 2))
                add(tgt[S], c["*"], (1 if r == 1 else -1) * sign * Fr(1, 2))
                touch(("tgt", S), tgt[S])
            for k in CEN:
                touch(("c", k), c[k])
        elif r == 2 or r == 7:        # A_SU +/- src_U
            for (S, U), a in A.items():
                add(a, src[U], (1 if r == 2 else -1) * sign)
                touch(("A", S, U), a)
            for U in T:
                touch(("src", U), src[U])
        elif r == 3 or r == 6:        # c +/- G src
            for k in CEN:
                for U in memb[k]:
                    add(c[k], src[U], (1 if r == 3 else -1) * sign)
                touch(("c", k), c[k])
            for U in T:
                touch(("src", U), src[U])
    init = {("src", U): cls(src[U]) for U in T}
    init.update({("tgt", S): cls(tgt[S]) for S in T})
    for k in CEN:
        init[("c", k)] = cls(c[k])
    for kk in A:
        init[("A",) + kk] = cls(A[kk])
    order = range(8) if not inverse else range(7, -1, -1)
    for r in order:
        row(r, 1 if not inverse else -1)
    assert all(A[k] == A0[k] for k in A) and all(c[k] == c0[k] for k in CEN), "auxiliaries not restored"
    # edge classes: initial, then after each touch except the last
    edges = {w: [init[w]] + ts[:-1] for w, ts in touches.items()}
    finals = {w: ts[-1] for w, ts in touches.items()}
    return edges, finals


def check():
    xs = {U: {("x", U): Fr(1)} for U in T}
    ys = {U: {("y", U): Fr(1)} for U in T}
    bad = 0
    Y = "Y"
    stage_states = []
    for j, inverse in ((1, False), (2, True), (3, False)):
        X0 = {U: dict(v) for U, v in xs.items()}
        Y0 = {U: dict(v) for U, v in ys.items()}
        edges, finals = invocation(xs, ys, inverse)
        # expected lists (oa_ledger.py), keyed by physical wire
        for S in T:
            nb = frozenset(NB[S])
            if j == 1:
                eX = [{S}] * 4; eY = [{Y}, {Y}, {Y}, {Y, S} | nb]
            elif j == 2:
                eX = [{S}, {S, Y} | nb, {Y}, {Y}]; eY = [{Y, S}] * 4
            else:
                eX = [{Y}] * 4; eY = [{Y, S}] * 4
            wX = ("src", S) if not inverse else ("tgt", S)
            wY = ("tgt", S) if not inverse else ("src", S)
            if [set(e) for e in edges[wX]] != [set(e) for e in eX]:
                bad += 1; print("X mismatch", j, S, edges[wX], eX)
            if [set(e) for e in edges[wY]] != [set(e) for e in eY]:
                bad += 1; print("Y mismatch", j, S, edges[wY], eY)
            for U in NB[S]:
                # physical carrier (Y-index S, X-index U): forward = A[(S, U)] (target S, source U);
                # inverse = A[(U, S)] (target X_U, source Y_S)
                w = ("A", S, U) if not inverse else ("A", U, S)
                if j == 1:
                    eA = [{Y}, {Y}, {Y, U}, {Y, U}]
                elif j == 2:
                    eA = [{Y}, {Y, S}, {Y, S}, {Y}]
                else:
                    eA = [{Y}] * 4
                if [set(e) for e in edges[w]] != [set(e) for e in eA]:
                    bad += 1; print("A mismatch", j, S, U, edges[w], eA)
        for k in CEN:
            mem = frozenset(memb[k])
            if j == 1:
                eC = [{Y}, {Y}, {Y} | mem, {Y} | mem]
            elif j == 2:
                eC = [{Y}, {Y} | mem, {Y} | mem, {Y}]
            else:
                eC = [{Y}] * 4
            if [set(e) for e in edges[("c", k)]] != [set(e) for e in eC]:
                bad += 1; print("C mismatch", j, k, edges[("c", k)], eC)
        # Lemma 2.1 per invocation
        if not inverse:
            ok = all(xs[U] == X0[U] for U in T) and all(
                ys[U] == {**{kk: v for kk, v in Y0[U].items()}, **{}} or True for U in T)
            for U in T:
                want = dict(Y0[U]); add(want, X0[U])
                ok = ok and ys[U] == want
        else:
            ok = all(ys[U] == Y0[U] for U in T)
            for U in T:
                want = dict(X0[U]); add(want, Y0[U], Fr(-1))
                ok = ok and xs[U] == want
        print("h=%d stage %d (%s): invocation %s; edge classes vs oa_ledger lists: %s" %
              (h, j, "inverse Y->X" if inverse else "forward X->Y", "exact" if ok else "WRONG",
               "all match" if bad == 0 else "%d mismatches" % bad))
        if not ok:
            bad += 1
    fin = all(xs[U] == {("y", U): -1} for U in T) and all(ys[U] == {("x", U): 1} for U in T)
    print("h=%d after three stages (X, Y) = (-y, x): %s" % (h, fin))
    return bad == 0 and fin


if __name__ == "__main__":
    print("PASS" if check() else "FAIL")
