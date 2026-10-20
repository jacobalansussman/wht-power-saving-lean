#!/usr/bin/env python3
"""t_chain.py (p6-prove): random exact test of the two mechanical claims inside Theorem C's proof.

Random words (ledger engine, exact) at blocks 2 and 3 with W = 2 or 3 inputs and scratch, in which a paid move is
kept only if its output register is level-homogeneous (Theorem C's hypothesis with h_l = level).  Generators of V
are the move outputs; E_k = span of eps(w) over outputs w of level k, E_0 = C^W.  Checked after every move:
  (a) eps(w) lies in E_{k'-1} (before adding w), k' = level of w >= 1;     (b) E_k is contained in E_{k-1}, k >= 1;
  (c) sum_k>=1 dim E_k <= number of moves so far.
Also counted: moves whose input register was NOT homogeneous (allowed by the theorem; output-only hypothesis).

TRUE formal model: ledger.py normalises group elements modulo {1, i, -1, -i} and folds the phase into the
coefficient.  In the every-f model a scalar MATRIX produced by letters is a group element (it acts as s^f), so here
every term carries its own power of i: key (input, gid, s) = the group element i^s * mats[gid].  The augmentation
eps (each group element -> 1) is computed on these true terms.  (Run 1 of this script used ledger coefficients and
showed spurious chain failures at block 3; see the run record, not published.)
"""
import sys, random
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "ledger"))
from ledger import Engine, ONE, GQ, UNITS, gq
from fractions import Fraction as Fr

def rank(vecs):
    """rank over Q(i) of a list of vectors with GQ entries (Gaussian elimination, exact)"""
    rows = [list(v) for v in vecs if any(not x.iszero() for x in v)]
    r = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(rows)):
            if not rows[i][c].iszero(): piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = rows[r][c].inv()
        rows[r] = [x * inv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and not rows[i][c].iszero():
                f = rows[i][c]; rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r

LOG = {}
def logi(ph):
    for k, u in enumerate(UNITS):
        if u == ph: return k
    raise ValueError(ph)

def tstep(E, C, st):
    """one step in the true model; contents {(l, gid, s): GQ}"""
    C = [dict(c) for c in C]
    def act(name, c):
        out = {}
        for (l, gid, s), co in c.items():
            nid, ph = E.lmul(name, gid); kk = (l, nid, (s + logi(ph)) % 4)
            out[kk] = out.get(kk, GQ(0)) + co
        return {k: v for k, v in out.items() if not v.iszero()}
    kind = st[0]
    if kind in ("move", "free"):
        r = st[1]; name = ("D%d" % st[2]) if kind == "move" else st[2]
        C[r] = act(name, C[r])
    elif kind == "mix":
        r, terms = st[1], st[2]; new = {}
        for (sreg, mu, name) in terms:
            mu = gq(mu); src = C[sreg] if name is None else act(name, C[sreg])
            for kk, co in src.items(): new[kk] = new.get(kk, GQ(0)) + mu * co
        C[r] = {k: v for k, v in new.items() if not v.iszero()}
    return C

def eps(E, c, W):
    v = [GQ(0) for _ in range(W)]
    for (l, gid, s), co in c.items(): v[l] = v[l] + co
    return v

def level_set(E, c): return {E.level(E.syms[g]) for (l, g, s) in c}

def trial(E, W, R, steps, rng, stats):
    C = [dict() for _ in range(R)]
    for r in range(W): C[r] = {(r, E.ID, 0): ONE}
    Ek = {0: [[GQ(1) if i == j else GQ(0) for i in range(W)] for j in range(W)]}
    moves = 0
    free = [n for n in E.letters if not n.startswith("D")]
    for _ in range(steps):
        t = rng.random()
        if t < 0.45:
            r = rng.randrange(R); z = rng.randrange(1, 1 << E.m)
            new = tstep(E, C, ("move", r, z))
            if not new[r]: continue
            ls = level_set(E, new[r])
            if len(ls) != 1: stats["rejected"] += 1; continue
            if len(level_set(E, C[r])) > 1: stats["inhom_input"] += 1
            k = ls.pop(); w = eps(E, new[r], W)
            if k >= 1:
                below = Ek.get(k - 1, [])
                if rank(below + [w]) != rank(below): stats["fail_a"] += 1
            Ek.setdefault(k, []).append(w)
            C = new; moves += 1; stats["moves"] += 1
            for kk in sorted(Ek):
                if kk >= 1 and rank(Ek.get(kk - 1, []) + Ek[kk]) != rank(Ek.get(kk - 1, [])): stats["fail_b"] += 1
            tot = sum(rank(Ek[kk]) for kk in Ek if kk >= 1)
            if tot > moves: stats["fail_c"] += 1
        elif t < 0.75:
            r = rng.randrange(R); terms = [(s, rng.choice([1, -1, 2, Fr(1, 2)]), rng.choice([None, None, rng.choice(free)]))
                                           for s in rng.sample(range(R), rng.randint(1, min(3, R)))]
            C = tstep(E, C, ("mix", r, terms))
        else:
            r = rng.randrange(R); C = tstep(E, C, ("free", r, rng.choice(free)))
    stats["maxsum"] = max(stats["maxsum"], sum(rank(Ek[kk]) for kk in Ek if kk >= 1))

if __name__ == "__main__":
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    for m in (2, 3):
        E = Engine(m)
        for W in (2, 3):
            stats = dict(moves=0, rejected=0, inhom_input=0, fail_a=0, fail_b=0, fail_c=0, maxsum=0)
            for _ in range(60 if m == 2 else 25):
                trial(E, W, W + 2, 40, rng, stats)
            print("m=%d W=%d %s" % (m, W, stats), flush=True)
