#!/usr/bin/env python3
"""checks.py (p6-prove, 15 Oct 2026): exact checks with sl-combinatorics' ledger engine (imported, not modified).

For each word: is it an exact unit (formal identity)?  M?  And which hypotheses hold at EVERY moment:
  single-term  : every entry (register, input) has at most one term;
  single-label : every entry's terms share one image in Sp(2m,2) (Pauli twins allowed);
  frame-coherent in some gauge (ledger.find_coherent_gauge), no-frame-packing gauge (ledger.find_nopack_gauge).
Cases:
  plain W=2 (calibration); C(2m)O3 order-3 control in 2 registers (calibration: outside every class here);
  SPLIT1: W = 1 wasteful unit with a moved same-input split at two frames (refutes the structural claim);
  HOLO3:  W = 3, NO scratch, single-term contents, holonomy, M = 3m (the episode lemma fails at W = 3);
"""
import sys, os
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "ledger"))
from ledger import Engine, run, find_coherent_gauge, find_nopack_gauge, _step, ONE, GQ, mul2, key2, gq
import numpy as np

def scan(E, W, R, word):
    C = [dict() for _ in range(R)]
    for r in range(W): C[r] = {(r, E.ID): ONE}
    st_ok = lab_ok = True; maxterms = 1; maxlabels = 1
    def check(C):
        nonlocal st_ok, lab_ok, maxterms, maxlabels
        for c in C:
            per = {}
            for (l, gid) in c: per.setdefault(l, []).append(gid)
            for l, gs in per.items():
                maxterms = max(maxterms, len(gs))
                if len(gs) > 1: st_ok = False
                labs = {key2(E.syms[g]) for g in gs}
                maxlabels = max(maxlabels, len(labs))
                if len(labs) > 1: lab_ok = False
    check(C)
    for stp in word:
        C = _step(E, C, stp); check(C)
    return st_ok, lab_ok, maxterms, maxlabels

def thmB(E, W, R, word):
    """Theorem B hypothesis: every MOVED register has <= 1 term per input, and if it holds two inputs their
    group elements differ by a FREE right factor (level of g1^-1 g2 == 0)."""
    if W != 2: return None
    C = [dict() for _ in range(R)]
    for r in range(W): C[r] = {(r, E.ID): ONE}
    ok = True; pairs = 0; twists = set()
    for stp in word:
        if stp[0] == "move":
            c = C[stp[1]]; per = {}
            for (l, gid) in c: per.setdefault(l, []).append(gid)
            if any(len(v) > 1 for v in per.values()): ok = False
            elif len(per) == 2:
                g1, g2 = per[0][0], per[1][0]
                t, _ = E.mul(E.inv(g1)[0], g2)
                pairs += 1; twists.add(key2(E.syms[t]))
                if E.level(E.syms[t]) != 0: ok = False
        C = _step(E, C, stp)
    return ok, pairs, len(twists)

def thmC(E, W, R, word):
    """Theorem C hypothesis with h_l = level for every input: after every paid move the moved register's terms
    (all inputs) sit at one level.  Also the co-level version h_l = m - d(., L_T)."""
    C = [dict() for _ in range(R)]
    for r in range(W): C[r] = {(r, E.ID): ONE}
    Kinv = E.inv(E.Kid)[0]
    ok_lv = ok_co = True
    for stp in word:
        C = _step(E, C, stp)
        if stp[0] == "move":
            c = C[stp[1]]
            lv = {E.level(E.syms[g]) for (l, g) in c}
            # distance of position lambda(g) = g^-1 L0 from L_T = K^-1 L0 equals level of g K^-1 ... use g * K^-1
            co = {E.level(E.syms[E.mul(g, Kinv)[0]]) for (l, g) in c}
            if len(lv) > 1: ok_lv = False
            if len(co) > 1: ok_co = False
    return ok_lv, ok_co

def true_unit(E, W, R, word):
    """unit check in the TRUE every-f model: terms keyed (input, gid, s) = i^s * mats[gid] (ledger folds i^s into
    the coefficient, which can fake cancellations); delivering register r must end as one term (r, Kid, any s)
    with a unit coefficient (the scalar matrix i^s is a free letter)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from t_chain import tstep
    C = [dict() for _ in range(R)]
    for r in range(W): C[r] = {(r, E.ID, 0): ONE}
    for stp in word:
        if stp[0] == "gate":
            regs, A = stp[1], stp[2]; old = [C[r] for r in regs]
            for i, r in enumerate(regs):
                new = {}
                for j in range(len(regs)):
                    mu = gq(A[i][j])
                    for kk, co in old[j].items(): new[kk] = new.get(kk, GQ(0)) + mu * co
                C[r] = {k: v for k, v in new.items() if not v.iszero()}
        else:
            C = tstep(E, C, stp)
    for r in range(W):
        c = C[r]
        if len(c) != 1: return False
        (l, g, s), co = next(iter(c.items()))
        if l != r or g != E.Kid or not (co.a * co.a + co.b * co.b == 1): return False
    return True

def report(name, E, W, R, word):
    res = run(E, W, R, word)
    st, lab, mt, ml = scan(E, W, R, word)
    tb = thmB(E, W, R, word)
    tc = thmC(E, W, R, word)
    gc, why = find_coherent_gauge(E, W, R, word)
    gn, whyn = find_nopack_gauge(E, W, R, word)
    tu = true_unit(E, W, R, word)
    print("== %-8s m=%d W=%d registers=%d  unit=%s (true every-f model: %s)  M=%d (W*m=%d)"
          % (name, E.m, W, R, res["unit"], tu, res["M"], W * E.m))
    print("   single-term always: %s (max terms per entry %d); single-label always: %s (max labels per entry %d)"
          % (st, mt, lab, ml))
    print("   frame-coherent gauge: %s [%s]; no-frame-packing gauge: %s [%s]"
          % (gc is not None, why, gn is not None, whyn))
    print("   Theorem C hypothesis (every move output at one level): %s; (at one distance from L_T): %s" % tc)
    if tb is not None:
        print("   Theorem B hypothesis (moved registers single-term, two-input ones with free twist): %s;"
              " packed moves %d; distinct twists %d" % tb)
    return res

def plain(m, W):
    return [("move", r, 1 << j) for j in range(m) for r in range(W)]

def main():
    out = {}
    for m in (2, 3):
        E = Engine(m)
        report("plain", E, 2, 2, plain(m, 2))
        # gauge exchange (sl-combinatorics t_known): 2 registers, single-term, inside Cor 2' and Theorem A
        w = [("move", 0, 1), ("mix", 0, [(0, 1, None), (1, 1, None)])]
        for j in range(1, m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        w += [("mix", 0, [(0, 1, None), (1, -1, None)]), ("move", 1, 1)]
        report("gaugeX", E, 2, 2, w)
        if m >= 3:
            # order-3 control (o3-build C(2m)O3; word copied from sl-combinatorics t_known.py)
            from fractions import Fraction as Fr
            s3 = E.addr([2, 4, 1] + [1 << j for j in range(3, m)]); s32 = E.addr([4, 1, 2] + [1 << j for j in range(3, m)])
            w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s3)])]
            for j in range(m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
            h = Fr(1, 2)
            w += [("mix", 1, [(0, h, None), (1, -h, None), (0, -h, s3), (1, h, s3), (0, h, s32), (1, -h, s32)]),
                  ("mix", 0, [(0, 1, None), (1, -1, None)])]
            report("C%dO3" % (2 * m), E, 2, 2, w)
        # SPLIT1: W = 1, registers c=0, t=1 (copy of x), s=2
        X1 = "X1"
        w = [("mix", 1, [(0, 1, None)]), ("mix", 2, [(0, 1, None)]), ("move", 2, 1),
             ("mix", 0, [(0, 1, None), (2, 1, None)]),          # c = (1 + D1) x : two frames, one input
             ("move", 0, 1),                                     # packed move: c = (D1 + X1) x
             ("mix", 0, [(0, 1, None), (1, -1, X1)])]            # c = D1 x
        w += [("move", 0, 1 << j) for j in range(1, m)]
        report("SPLIT1", E, 1, 3, w)
        # HOLO3: W = 3, no scratch; s = swap of the first two coordinates (commutes with K)
        s = E.addr([2, 1] + [1 << j for j in range(2, m)])
        w = [("mix", 0, [(0, 1, None), (1, 1, None)]),           # R0 = x0 + x1
             ("free", 1, s),                                     # R1 = s x1
             ("mix", 1, [(1, 1, None), (2, 1, None)]),           # R1 = s x1 + x2
             ("mix", 2, [(2, 1, None), (0, 1, None)])]           # R2 = x0 + x1 + x2   (holonomy: s vs 1)
        w += plain(m, 3)
        w += [("mix", 2, [(2, 1, None), (0, -1, None)]),         # R2 = K x2
              ("mix", 1, [(1, 1, None), (2, -1, None)]),         # R1 = s K x1
              ("free", 1, s),                                    # R1 = K x1   (s^2 = 1)
              ("mix", 0, [(0, 1, None), (1, -1, None)])]         # R0 = K x0
        report("HOLO3", E, 3, 3, w)

if __name__ == "__main__":
    main()
