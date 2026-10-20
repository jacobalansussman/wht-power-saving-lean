#!/usr/bin/env python3
"""vstg.py (v-p6-refute): my own exhaustive solver of the single-term two-input game (Lemma M + Lemma G, permissive
closure), written independently of stg.py.  No move pruning, no symmetry reduction; global dedup (each state is
recorded at its least depth, which is the strongest test of (I)).  Moves are generated from group elements:
  single:   e in an available coset of input l, a in H  ->  position pos(a e)
  copy:     lineage (t, xs): x in xs, a in H            ->  add pos(a rep(x)) to xs
  created:  e1 in rep of an available x0, e2 ANY element of an available coset y0 of input 2, a in H
            -> new lineage (twist = right action of e1^-1 e2, {pos(a e1)})
  composite: path of >= 3 edges in a locked tree component from (1,u) to (2,w), twist = product along the path
            -> new lineage at pos(a rep(u))
State: (A1 bitmask, A2 bitmask, sorted tuple of lineages (twist permutation id, xs bitmask)).
Usage: vstg.py GAME DEPTH [maxstates]"""
import sys, time
import vgroups as VG

def make(name):
    if ":" in name:
        nm, ci = name.split(":")
        gens, Pg = VG.TOYS[nm]
        return VG.all_classes(gens, Pg, nm)[int(ci)]
    return {"D5": lambda: VG.dihedral(5), "D7": lambda: VG.dihedral(7), "D9": lambda: VG.dihedral(9),
            "B3": lambda: VG.hyperoct(3), "B4": lambda: VG.hyperoct(4), "S4Z3": VG.s4z3, "J42": VG.johnson42,
            "H23": lambda: VG.hamming(2), "H33": lambda: VG.hamming(3), "H43": lambda: VG.hamming(4)}[name]()

class S:
    def __init__(self, g):
        self.g = g; N = g.N
        # twist permutations of positions (right action), deduplicated
        self.tp = []; self.tid = {}; self.elem_t = []
        for e in range(g.n):
            p = g.ract[e]
            if p not in self.tid: self.tid[p] = len(self.tp); self.tp.append(p)
            self.elem_t.append(self.tid[p])
        # left action of paid letters on positions: from position x (rep) by a in H
        self.step = [sorted({g.pos[g.M[a][g.reps[x]]] for a in g.H}) for x in range(N)]
        # sanity: same as g.nb
        assert self.step == [sorted(set(v)) for v in g.nb]
        # composition of twist permutations: comp[i][j] = perm i after perm j (x -> i(j(x)))
        T = len(self.tp)
        self.comp = [[self.tid[tuple(self.tp[i][self.tp[j][x]] for x in range(N))] for j in range(T)] for i in range(T)]
        self.inv = [self.tid[tuple(sorted(range(N), key=lambda x: self.tp[i][x]))] for i in range(T)]
        self.tmask = sum(1 << x for x in g.target)
        # created pairs: for x0, y0 the set of twists t with t(x0) = y0, generated from group elements
        self.cre = {}
        for x0 in range(N):
            e1 = g.reps[x0]
            for e2 in range(g.n):
                t = self.elem_t[g.M[g.I[e1]][e2]]
                self.cre.setdefault((x0, g.pos[e2]), set()).add(t)
        for k in self.cre: self.cre[k] = sorted(self.cre[k])

    def close(self, A1, A2, lins):
        N = self.g.N; tp = self.tp
        while True:
            par = list(range(2 * N))
            def f(a):
                while par[a] != a:
                    par[a] = par[par[a]]; a = par[a]
                return a
            ed = []
            for (t, xs) in lins:
                x = 0; m = xs
                while m:
                    if m & 1:
                        u, w = x, N + tp[t][x]; ed.append((u, w))
                        ru, rw = f(u), f(w)
                        if ru != rw: par[ru] = rw
                    m >>= 1; x += 1
            if not ed: return (A1, A2, lins)
            nv = {}; ne = {}; av = {}
            for (u, w) in ed:
                r = f(u); ne[r] = ne.get(r, 0) + 1
                s = nv.setdefault(r, set()); s.add(u); s.add(w)
            un = set()
            for r, vs in nv.items():
                a = any((v < N and (A1 >> v) & 1) or (v >= N and (A2 >> (v - N)) & 1) for v in vs)
                if a or ne[r] >= len(vs): un.add(r)
            if not un: return (A1, A2, lins)
            for r in un:
                for v in nv[r]:
                    if v < N: A1 |= 1 << v
                    else: A2 |= 1 << (v - N)
            nl = []
            for (t, xs) in lins:
                keep = 0; x = 0; m = xs
                while m:
                    if m & 1 and f(x) not in un: keep |= 1 << x
                    m >>= 1; x += 1
                if keep: nl.append((t, keep))
            lins = tuple(sorted(nl))

    def moves(self, st):
        A1, A2, lins = st; N = self.g.N; tp = self.tp; out = []
        bitsA1 = [x for x in range(N) if (A1 >> x) & 1]; bitsA2 = [x for x in range(N) if (A2 >> x) & 1]
        for x in bitsA1:
            for x2 in self.step[x]: out.append((A1 | (1 << x2), A2, lins))
        for y in bitsA2:
            for y2 in self.step[y]: out.append((A1, A2 | (1 << y2), lins))
        for i, (t, xs) in enumerate(lins):
            for x in range(N):
                if (xs >> x) & 1:
                    for x2 in self.step[x]:
                        nl = list(lins); nl[i] = (t, xs | (1 << x2)); out.append((A1, A2, tuple(sorted(nl))))
        for x0 in bitsA1:
            for y0 in bitsA2:
                for t in self.cre.get((x0, y0), ()):
                    for x2 in self.step[x0]:
                        out.append((A1, A2, tuple(sorted(lins + ((t, 1 << x2),)))))
        # composites
        adj = {}
        for (t, xs) in lins:
            for x in range(N):
                if (xs >> x) & 1:
                    u, w = x, N + tp[t][x]
                    adj.setdefault(u, []).append((w, t)); adj.setdefault(w, []).append((u, t))
        for u0 in [v for v in adj if v < N]:
            seen = {u0: (0, 0)}; stck = [u0]       # twist id 0 must be the identity
            while stck:
                v = stck.pop(); tv, d = seen[v]
                for (nb, t) in adj[v]:
                    if nb in seen: continue
                    nt = self.comp[t][tv] if v < N else self.comp[self.inv[t]][tv]
                    seen[nb] = (nt, d + 1); stck.append(nb)
                    if nb >= N and d + 1 >= 3:
                        for x2 in self.step[u0]:
                            out.append((A1, A2, tuple(sorted(lins + ((nt, 1 << x2),)))))
        return out

def lvmax(g, A):
    return max(g.lv[x] for x in range(g.N) if (A >> x) & 1)

def bfs(name, D, maxstates=3000000, report=print):
    g = make(name); s = S(g)
    assert s.tp[0] == tuple(range(g.N))
    st0 = s.close(1 << g.raw, 1 << g.raw, ())
    seen = {st0}; layer = [st0]; first_goal = None; rows = []
    for M in range(1, D + 1):
        t0 = time.time(); new = []; mx = -1; bad = 0; goal = 0
        for st in layer:
            for s2 in s.moves(st):
                c = s.close(*s2)
                if c in seen: continue
                seen.add(c); new.append(c)
                v = lvmax(g, c[0]) + lvmax(g, c[1]); mx = max(mx, v)
                if v > M: bad += 1
                if (c[0] & s.tmask) and (c[1] & s.tmask): goal += 1
        layer = new
        if goal and first_goal is None: first_goal = M
        rows.append((M, len(new), mx, bad, goal))
        report("%s M=%d new=%d max(L1+L2 among new)=%d (I) %s goal %d [%.1fs]" % (name, M, len(new), mx,
               "holds" if not bad else "VIOLATED %d" % bad, goal, time.time() - t0))
        sys.stdout.flush()
        if len(seen) > maxstates: report("  cap"); break
    return rows

if __name__ == "__main__":
    bfs(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 3000000)
