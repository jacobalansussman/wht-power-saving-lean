#!/usr/bin/env python3
"""vgroups.py (v-p6-refute, 15 Oct 2026): my own small permutation groups and frame games, independent of toylab.
Convention: permutations are tuples, (a*b)[i] = a[b[i]] (a after b).  Positions = right cosets P g; a paid move is
left multiplication by an element of the paid double coset P h P; twists act on positions by right multiplication."""
import itertools
from collections import deque

def mul(a, b): return tuple(a[i] for i in b)
def inv(a):
    r = [0] * len(a)
    for i, x in enumerate(a): r[x] = i
    return tuple(r)

def closure_group(gens):
    n = len(gens[0]); e = tuple(range(n))
    seen = {e: 0}; elems = [e]; q = deque([e])
    while q:
        g = q.popleft()
        for s in gens:
            h = mul(s, g)
            if h not in seen: seen[h] = len(elems); elems.append(h); q.append(h)
    return elems, seen

def cyc(n, *cs):
    p = list(range(n))
    for c in cs:
        for i in range(len(c)): p[c[i]] = c[(i + 1) % len(c)]
    return tuple(p)

class VGame:
    def __init__(self, gens, Pgens, h, name):
        self.name = name
        self.elems, self.idx = closure_group(gens)
        self.n = len(self.elems)
        Pel, _ = closure_group(Pgens) if Pgens else ([tuple(range(len(gens[0])))], None)
        self.P = [self.idx[p] for p in Pel]
        E = self.elems
        self.M = [[self.idx[mul(a, b)] for b in E] for a in E]       # M[a][b] = a*b
        self.I = [self.idx[inv(a)] for a in E]
        # right cosets P g
        self.pos = [-1] * self.n; self.reps = []
        for g in range(self.n):
            if self.pos[g] < 0:
                c = len(self.reps); self.reps.append(g)
                for p in self.P: self.pos[self.M[p][g]] = c
        self.N = len(self.reps)
        hs = h if isinstance(h, list) else [h]
        H = set()
        for hh in hs:
            hi = self.idx[hh]
            H |= {self.M[self.M[p][hi]][q] for p in self.P for q in self.P}
        H |= {self.I[x] for x in H}
        assert not (H & set(self.P))
        self.H = sorted(H)
        self.nb = [sorted({self.pos[self.M[x][self.reps[c]]] for x in self.H}) for c in range(self.N)]
        lv = [-1] * self.N; lv[self.pos[0]] = 0; fr = [self.pos[0]]
        while fr:
            nf = []
            for c in fr:
                for c2 in self.nb[c]:
                    if lv[c2] < 0: lv[c2] = lv[c] + 1; nf.append(c2)
            fr = nf
        self.lv = lv; self.m = max(lv)
        assert min(lv) >= 0
        self.raw = self.pos[0]
        # right action of an element t on positions
        self.ract = [tuple(self.pos[self.M[self.reps[c]][t]] for c in range(self.N)) for t in range(self.n)]
        top = [c for c in range(self.N) if lv[c] == self.m]
        # target orbit: images of one top position under right multiplication by P
        self.target = {self.ract[p][top[0]] for p in self.P}

def dihedral(n):
    r = tuple((i + 1) % n for i in range(n)); s = tuple((-i) % n for i in range(n))
    return VGame([r, s], [s], r, "D%d" % n)

def hyperoct(m):
    n = 2 * m
    gens = [cyc(n, [0, 1])] + [cyc(n, [2 * i, 2 * i + 2], [2 * i + 1, 2 * i + 3]) for i in range(m - 1)]
    return VGame(gens, gens[1:], cyc(n, [0, 1]), "B%d" % m)

def s4z3():
    gens = [cyc(4, [0, 1]), cyc(4, [0, 1, 2, 3])]
    # P = <(0 1 2)>; paid class representative chosen so that the frame graph is connected
    for h in [cyc(4, [0, 3]), cyc(4, [0, 1, 3]), cyc(4, [2, 3])]:
        try:
            g = VGame(gens, [cyc(4, [0, 1, 2])], h, "S4Z3")
            return g
        except AssertionError:
            continue

def johnson42():
    gens = [cyc(4, [0, 1]), cyc(4, [0, 1, 2, 3])]
    return VGame(gens, [cyc(4, [0, 1]), cyc(4, [2, 3])], cyc(4, [1, 2]), "J42")

def hamming(k):
    n = 3 * k; gens = []
    for i in range(k): gens += [cyc(n, [3 * i, 3 * i + 1]), cyc(n, [3 * i, 3 * i + 1, 3 * i + 2])]
    Pg = [cyc(n, [3 * i + 1, 3 * i + 2]) for i in range(k)]
    return VGame(gens, Pg, [cyc(n, [3 * i, 3 * i + 1]) for i in range(k)], "H%d3" % k)

if __name__ == "__main__":
    for g in [dihedral(7), dihedral(9), hyperoct(3), hyperoct(4), s4z3(), johnson42(), hamming(2), hamming(3)]:
        print(g.name, "|G|", g.n, "|P|", len(g.P), "frames", g.N, "deg", sorted({len(x) for x in g.nb}), "m", g.m,
              "target orbit", len(g.target))

# ------------------------------------------------------------------ every paid class of a (G, P) pair (my own enumeration)
def all_classes(gens, Pgens, name):
    base = VGame.__new__(VGame)
    elems, idx = closure_group(gens); n = len(elems)
    Pel, _ = closure_group(Pgens)
    P = {idx[p] for p in Pel}
    Mt = lambda a, b: idx[mul(elems[a], elems[b])]
    seen = {}; reps = []
    for x in range(n):
        if x in seen: continue
        k = len(reps); reps.append(x)
        for p in P:
            for q in P: seen[Mt(Mt(p, x), q)] = k
    out = []; used = set()
    for k, h in enumerate(reps):
        if h in P or k in used: continue
        ki = seen[idx[inv(elems[h])]]; used |= {k, ki}
        try:
            g = VGame(gens, Pgens, [elems[h]], "%s:c%d" % (name, len(out)))
        except AssertionError:
            out.append(None); continue
        out.append(g)
    return out

TOYS = {
    "A5Z2": ([cyc(5, [0, 1, 2]), cyc(5, [0, 1, 2, 3, 4])], [cyc(5, [0, 1], [2, 3])]),
    "A5Z3": ([cyc(5, [0, 1, 2]), cyc(5, [0, 1, 2, 3, 4])], [cyc(5, [0, 1, 2])]),
    "A5Z5": ([cyc(5, [0, 1, 2]), cyc(5, [0, 1, 2, 3, 4])], [cyc(5, [0, 1, 2, 3, 4])]),
    "A5V4": ([cyc(5, [0, 1, 2]), cyc(5, [0, 1, 2, 3, 4])], [cyc(5, [0, 1], [2, 3]), cyc(5, [0, 2], [1, 3])]),
    "S5D5": ([cyc(5, [0, 1]), cyc(5, [0, 1, 2, 3, 4])], [cyc(5, [0, 1, 2, 3, 4]), cyc(5, [1, 4], [2, 3])]),
    "S4Z2": ([cyc(4, [0, 1]), cyc(4, [0, 1, 2, 3])], [cyc(4, [0, 1])]),
    "S4Z3all": ([cyc(4, [0, 1]), cyc(4, [0, 1, 2, 3])], [cyc(4, [0, 1, 2])]),
    "D7all": ([tuple((i + 1) % 7 for i in range(7)), tuple((-i) % 7 for i in range(7))], [tuple((-i) % 7 for i in range(7))]),
    "J63": ([cyc(6, [0, 1]), cyc(6, [0, 1, 2, 3, 4, 5])], [cyc(6, [0, 1]), cyc(6, [0, 1, 2]), cyc(6, [3, 4]), cyc(6, [3, 4, 5])]),
}
TOYS.update({
    "S5S3": ([cyc(5, [0, 1]), cyc(5, [0, 1, 2, 3, 4])], [cyc(5, [0, 1, 2]), cyc(5, [0, 1])]),
    "A633": ([cyc(6, [0, 1, 2]), cyc(6, [1, 2, 3, 4, 5])], [cyc(6, [0, 1, 2]), cyc(6, [3, 4, 5]), cyc(6, [0, 1], [3, 4])]),
})
