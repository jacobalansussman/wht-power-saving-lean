#!/usr/bin/env python3
"""t_lemma6.py (pk2-orbits): computed check of PROOF.md Lemma 6.  For random trunk skeletons (2 in 2m - 1), for
every prefix state of each input and every slot of the branch that carries that input as junk, test whether some free chart mu makes the injected prefix walker land EXACTLY on the junk element at the end of
the branch (mu = w^-1 J nu0^-1 must lie in P).  Lemma 6 says never.  Also counts frame-level landings (Lemma 3.2),
which exist."""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "heights")); sys.path.insert(0, os.path.join(HERE, "..", "frame-game")); sys.path.insert(0, HERE)
from fpa import *
import design as D, orb
for m in (3, 4, 5):
    S = Sp(m); rng = random.Random(600 + m)
    exact = frame = tests = designs = 0
    for tr in range(40):
        P = D.Prog(S, m, rng, dup=0)
        if not P.ok: continue
        st = P.run(); designs += 1
        for (branch, junk_input, prefix) in (("A", 1, "R2"), ("B", 0, "R1")):
            n = len(P.regs[branch])
            if n == 0: continue
            fin = st[P.out[(branch, n - 1)]]
            junk = [g for (l, g) in fin if l == junk_input]
            # prefix states of the junk input: the input itself and its prefix outputs
            pre = [S.I] + [next(iter(st[P.out[(prefix, k)]]))[1] for k in range(len(P.regs[prefix]))]
            for k in range(n):
                w = S.I
                for u in P.regs[branch][k:]: w = S.mul(S.transv(u), w)
                for g0 in pre:
                    for J in junk:
                        tests += 1
                        mu = S.mul(S.inv(w), S.mul(J, S.inv(g0)))
                        if S.in_P(mu): exact += 1
                        if S.rank_rows(list(S.pos(g0)) + list(S.pos(J))) - m == S.level(w): frame += 1
    print(f"m={m}: {designs} trunk skeletons, {tests} (slot, prefix walker, junk) tests: exact landings {exact} "
          f"(Lemma 6: 0), frame-level landings {frame}", flush=True)
