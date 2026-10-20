#!/usr/bin/env python3
"""vchain.py (v-p6-prove): random exact test of Theorem C's chain, with my own engine (true model) and with
HETEROGENEOUS height systems (a different 1-Lipschitz height per input), including heights with negative values and
h(L_T) < m (quantitative form).  p6-prove tested only h = level.

A paid move is kept only if its output register is homogeneous for the chosen heights.  E_k = span of eps(w) over kept
outputs of degree k (E_0 = C^W).  Checked after every kept move:
  (a) eps(w) in E_{k'-1} + E_{k'} + E_{k'+1}  (sum before adding w), and for k' >= 1 also in E_{k'-1};
  (b) E_k subset of E_{k-1} for every k >= 1;  (c) sum_{k>=1} dim E_k <= moves.
Also (d): for every register content after every step, its homogeneous components' eps lie in the matching E_k
(the gradedness of V that the proof uses), i.e. eps(c^(k)) in E_k for every k.
"""
import sys, random
from fractions import Fraction as Fr
from veng import Eng, step, Q

def rank(vecs):
    rows = [list(v) for v in vecs if any(not x.z() for x in v)]
    r = 0; nc = len(rows[0]) if rows else 0
    for c in range(nc):
        p = next((i for i in range(r, len(rows)) if not rows[i][c].z()), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        a = rows[r][c]; n = a.abs2(); inv = Q(a.a / n, -a.b / n)
        rows[r] = [x * inv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and not rows[i][c].z():
                f = rows[i][c]; rows[i] = [x + Q(-1) * f * y for x, y in zip(rows[i], rows[r])]
        r += 1
    return r

def inspan(vs, w): return rank(vs + [w]) == rank(vs)

def main(seed):
    rng = random.Random(seed)
    stats_all = []
    for m in (2, 3):
        E = Eng(m); L0 = E.L0(); LT = E.pos(E.K)
        lev = lambda x: E.d(L0, x); col = lambda x: m - E.d(x, LT)
        H = {"level": lev, "colevel": col, "mid": lambda x: (lev(x) + col(x)) // 2,
             "vee": lambda x: abs(lev(x) - 1) - 1, "cap1": lambda x: min(lev(x), 1)}
        systems = [("level,colevel", ["level", "colevel"]), ("colevel,mid,level", ["colevel", "mid", "level"]),
                   ("vee,level", ["vee", "level"]), ("mid,cap1,vee", ["mid", "cap1", "vee"])]
        free = [n for n in E.letters if not n.startswith("D")]
        poscache = {}
        def hgt(l, g, hs):
            if g not in poscache: poscache[g] = E.pos(g)
            return H[hs[l]](poscache[g])
        for (sname, hs) in systems:
            W = len(hs); R = W + 2
            st = dict(moves=0, rejected=0, multi_input_kept=0, fail_a=0, fail_b=0, fail_c=0, fail_d=0, maxdeg=0, mindeg=0)
            ntr = 30 if m == 2 else 10
            for _ in range(ntr):
                C = [dict() for _ in range(R)]
                for r in range(W): C[r] = {(r, E.ID): Q(1)}
                Ek = {0: [[Q(1) if i == j else Q(0) for i in range(W)] for j in range(W)]}
                mv = 0
                for _ in range(30):
                    t = rng.random()
                    if t < 0.45:
                        r = rng.randrange(R); z = rng.randrange(1, 1 << m)
                        new = step(E, C, ("move", r, z))
                        if not new[r]: continue
                        degs = {hgt(l, g, hs) for (l, g) in new[r]}
                        if len(degs) != 1: st["rejected"] += 1; continue
                        k = degs.pop(); st["maxdeg"] = max(st["maxdeg"], k); st["mindeg"] = min(st["mindeg"], k)
                        if len({l for (l, g) in new[r]}) > 1: st["multi_input_kept"] += 1
                        w = [Q(0)] * W
                        for (l, g), co in new[r].items(): w[l] = w[l] + co
                        nb = Ek.get(k - 1, []) + Ek.get(k, []) + Ek.get(k + 1, [])
                        if not inspan(nb, w): st["fail_a"] += 1
                        if k >= 1 and not inspan(Ek.get(k - 1, []), w): st["fail_a"] += 1
                        Ek.setdefault(k, []).append(w); C = new; mv += 1; st["moves"] += 1
                        for kk in Ek:
                            if kk >= 1 and rank(Ek.get(kk - 1, []) + Ek[kk]) != rank(Ek.get(kk - 1, [])): st["fail_b"] += 1
                        if sum(rank(Ek[kk]) for kk in Ek if kk >= 1) > mv: st["fail_c"] += 1
                    elif t < 0.8:
                        r = rng.randrange(R)
                        terms = [(s, rng.choice([1, -1, 2, Fr(1, 2), Q(0, 1)]), rng.choice([None, None, rng.choice(free)]))
                                 for s in rng.sample(range(R), rng.randint(1, min(3, R)))]
                        C = step(E, C, ("mix", r, terms))
                    else:
                        r = rng.randrange(R); C = step(E, C, ("free", r, rng.choice(free)))
                    # (d) gradedness: eps of each homogeneous component of each register lies in E_k
                    for c in C:
                        comp = {}
                        for (l, g), co in c.items():
                            k = hgt(l, g, hs); v = comp.setdefault(k, [Q(0)] * W); v[l] = v[l] + co
                        for k, v in comp.items():
                            if any(not x.z() for x in v) and not inspan(Ek.get(k, []), v): st["fail_d"] += 1
            print("m=%d heights=%-18s W=%d %s" % (m, sname, W, st), flush=True)

if __name__ == "__main__":
    for sd in (sys.argv[1:] or ["1"]):
        main(int(sd))
