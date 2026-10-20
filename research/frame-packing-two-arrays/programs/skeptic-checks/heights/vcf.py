#!/usr/bin/env python3
"""vcf.py (v-p6-prove): the literal reading of NOTE.md's "for a common-frame design, suitable heights exist exactly
when every start frame is opposite every end frame" is false without PROOF.md's restriction to heights H o gauge.
Plain W = 2 schedule with sl-alg Definition 2 twists k_1 = e, k_2 = D1: common-frame (every row times diag(k) is one
group element), source a_s = lam(k_s), sink a*_s = k_s^-1 L_T; d(a_2, a*_1) = m - 1 < m (a crossing overlap), yet the
level height system makes every move output homogeneous (vunits.out: plain, system exists)."""
from veng import Eng, run
for m in (2, 3, 4):
    E = Eng(m)
    k = [E.ID, E.letters["D1"]]
    word = [("move", r, 1 << j) for j in range(m) for r in range(2)]
    ok, M, hist = run(E, 2, 2, word)
    coherent = all(len({E.mulk(g, k[l]) for (l, g) in c}) <= 1 for C in hist for c in C)
    a = [E.pos(kk) for kk in k]
    LT = E.pos(E.K)
    astar = [E.span([E.apply_symp(E.symp(kk), v) for v in LT]) for kk in k]   # k^-1 L_T as a position: lam(K k) = k^-1 lam(K)
    astar = [E.pos(E.mulk(E.K, kk)) for kk in k]
    print("m=%d unit=%s M=%d common-frame with twists (e, D1): %s; d(a_s, a*_t) = %s"
          % (m, ok, M, coherent, {(s, t): E.d(a[s], astar[t]) for s in range(2) for t in range(2)}))
