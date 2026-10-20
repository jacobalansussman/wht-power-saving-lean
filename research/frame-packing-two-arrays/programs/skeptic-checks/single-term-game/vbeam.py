#!/usr/bin/env python3
"""vbeam.py (v-p6-refute): adversarial beam search for a violation of (I) (L1 + L2 > M) in the single-term game, in
games too big to exhaust (H(4,3) at depths 6-7: the block-4 Hamming sub-game, where p6-refute stopped at 5).
Uses my solver's moves and closure (vstg.S).  Score = L1 + L2 + w * (best locked potential) + small random tie-break,
where the locked potential of a component is (highest input-1 level) + (highest input-2 level) in it.
A found violation is printed with its whole play.  Usage: vbeam.py GAME DEPTH WIDTH W SEED"""
import sys, random
import vstg

def main(name, D, width, w, seed):
    rng = random.Random(seed)
    s = vstg.S(vstg.make(name)); g = s.g; N = g.N; tp = s.tp
    lvmax = lambda A: max(g.lv[x] for x in range(N) if (A >> x) & 1)
    def pot(L):
        best = 0
        par = {}
        def f(a):
            while par.get(a, a) != a: a = par[a]
            return a
        for (t, xs) in L:
            for x in range(N):
                if (xs >> x) & 1:
                    ru, rw = f(x), f(N + tp[t][x])
                    if ru != rw: par[ru] = rw
        hi = {}
        for (t, xs) in L:
            for x in range(N):
                if (xs >> x) & 1:
                    r = f(x); h = hi.setdefault(r, [0, 0])
                    h[0] = max(h[0], g.lv[x]); h[1] = max(h[1], g.lv[tp[t][x]])
        for h in hi.values(): best = max(best, h[0] + h[1])
        return best
    st0 = s.close(1 << g.raw, 1 << g.raw, ())
    beam = [(st0, [])]; seen = {st0}
    for M in range(1, D + 1):
        cand = []
        for (st, hist) in beam:
            for s2 in s.moves(st):
                c = s.close(*s2)
                if c in seen: continue
                seen.add(c)
                v = lvmax(c[0]) + lvmax(c[1])
                if v > M:
                    print("VIOLATION at M=%d: L1+L2=%d state %s" % (M, v, c), flush=True)
                    for h in hist + [c]: print("   ", h)
                    return
                cand.append((v + w * pot(c[2]) + 1e-3 * rng.random(), c, hist + [c]))
        cand.sort(key=lambda z: -z[0])
        beam = [(c, h) for (_, c, h) in cand[:width]]
        best = cand[0] if cand else None
        print("%s M=%d candidates %d, best score %.2f, best L1+L2 %d, best potential %d" % (
            name, M, len(cand), best[0], lvmax(best[1][0]) + lvmax(best[1][1]), pot(best[1][2])), flush=True)
    print("no violation found (beam %d, w %.2f, seed %d)" % (width, w, seed))

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]))
