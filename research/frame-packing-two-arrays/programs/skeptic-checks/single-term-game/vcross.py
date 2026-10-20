#!/usr/bin/env python3
"""vcross.py (v-p6-refute): cross-check p6-refute's stg.py functions against my vstg.py on the SAME game, built once
from my vgroups (so both use one position numbering and one twist list).
  (1) closure: stg.closure == vstg.close on every successor of every state of my BFS up to a depth;
  (2) moves: {closure(s') : s' in stg.successors(s)} union {s}  ==  {close(s') : s' in my unpruned moves(s)} union {s};
  (3) goal_in_one: stg.goal_in_one(s, level) == brute force over my unpruned moves, for every state up to a depth,
      with the real target and with an artificial LOWER target (a level-(m-1) free orbit) so that goals are frequent;
  (4) bfs3's level filter: whenever brute force finds a goal from a closed state c, c passes the filter.
Usage: vcross.py GAME DEPTH"""
import sys
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "single-term-game"))
import stg
import vstg

def build(name):
    s = vstg.S(vstg.make(name)); g = s.g
    free = [i for i, p in enumerate(s.tp) if p[g.raw] == g.raw]
    G = stg.Game(g.N, g.nb, g.lv, g.target, s.tp, free, name=name)
    G.tmask = sum(1 << x for x in g.target)
    return s, G

def brute_goal(s, st, tmask):
    for s2 in s.moves(st):
        c = s.close(*s2)
        if (c[0] & tmask) and (c[1] & tmask): return True
    return False

def main(name, D):
    s, G = build(name); g = s.g
    assert G.inv == s.inv
    # artificial lower target: one free orbit at level m-1
    low = next(x for x in range(g.N) if g.lv[x] == g.m - 1)
    lowset = {s.tp[a][low] for a in range(len(s.tp)) if s.tp[a][g.raw] == g.raw}
    lowmask = sum(1 << x for x in lowset)
    st0 = s.close(1 << g.raw, 1 << g.raw, ())
    seen = {st0}; layer = [st0]
    nst = ncl = nmv = ng = 0; bad = []; goals = {"real": 0, "low": 0}; filt_bad = 0
    for M in range(1, D + 1):
        new = []
        for st in layer:
            nst += 1
            mine = set()
            for s2 in s.moves(st):
                a = s.close(*s2); b = stg.closure(G, *s2); ncl += 1
                if a != b: bad.append(("closure", st, s2, a, b))
                mine.add(a)
            theirs = {stg.closure(G, *s2) for (k, s2) in stg.successors(G, st)}
            nmv += 1
            if (mine | {st}) != (theirs | {st}):
                bad.append(("moves", st, len(mine - theirs - {st}), len(theirs - mine - {st})))
            for lab, tm, lvl in (("real", G.tmask, g.m), ("low", lowmask, g.m - 1)):
                G.tmask = tm
                bg = brute_goal(s, st, tm)
                sg = stg.goal_in_one(G, st, lvl)
                ng += 1
                if bg != sg: bad.append(("goal_in_one", lab, st, bg, sg))
                if bg: goals[lab] += 1
                if bg and lab == "real":
                    # bfs3 filter (computed as in stg.bfs3 on the closed state)
                    Q = [-1, -1]
                    for x in range(g.N):
                        if (st[0] >> x) & 1: Q[0] = max(Q[0], g.lv[x])
                        if (st[1] >> x) & 1: Q[1] = max(Q[1], g.lv[x])
                    for (t, xs) in st[2]:
                        for x in range(g.N):
                            if (xs >> x) & 1: Q[0] = max(Q[0], g.lv[x]); Q[1] = max(Q[1], g.lv[s.tp[t][x]])
                    if Q[0] < g.m - 1 or Q[1] < g.m - 1: filt_bad += 1
            G.tmask = sum(1 << x for x in g.target)
            for a in mine:
                if a not in seen: seen.add(a); new.append(a)
        layer = new
    print("%s depth %d: states %d, closure comparisons %d, move-set comparisons %d, goal tests %d; goals found real %d low %d;"
          " filter misses %d; DISAGREEMENTS %d" % (name, D, nst, ncl, nmv, ng, goals["real"], goals["low"], filt_bad, len(bad)))
    for b in bad[:5]: print("  ", b)

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
