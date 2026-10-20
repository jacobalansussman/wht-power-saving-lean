#!/usr/bin/env python3
"""t_cutcheck.py (sl-combinatorics): direct check of Lemma 4.2 (cut condition) on actual words.
For each exact unit with a no-frame-packing gauge: map every live move to an edge of the frame graph Gamma_m (the
moved register's single frame class before -> after), check each move is an edge (distance exactly 1), and check
c(delta Q) >= #{a : exactly one of S_a, T_a in Q} for EVERY vertex set Q at block 2 (16,383 sets) and for 3,000
random sets plus all balls and vector sets at block 3.  Also checks M >= W m whenever W <= 2 (Theorem 4.3).
usage: python3 t_cutcheck.py [seed]"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ledger import *
from cutconst import setup
import t_random

def lag_of(E, lam):
    """bitset over vector indices of the Lagrangian lam^-1 L0"""
    m = E.m; Li = inv2(lam)
    cols = [Li[:, m + j] for j in range(m)]
    out = 0
    for mask in range(1 << m):
        v = np.zeros(2 * m, dtype=np.int64)
        for j in range(m):
            if (mask >> j) & 1: v ^= cols[j]
        idx = 0
        for k in range(2 * m):
            if v[k] & 1: idx |= 1 << k
        out |= 1 << idx
    return out

def check(E, W, R, word, labs_index, D, Qs, name):
    h, why = find_nopack_gauge(E, W, R, word)
    res = run(E, W, R, word)
    if not res["unit"]: return "not a unit"
    if h is None: return "skipped (%s)" % why
    # replay, recording the moved register's frame class before and after each move
    C = [dict() for _ in range(R)]
    for r in range(W): C[r] = {(r, E.ID): ONE}
    V = len(labs_index); cap = {}
    for st in word:
        if st[0] == "move":
            r = st[1]
            if C[r]:
                (l, gid) = next(iter(C[r]))
                lam0 = mul2(E.syms[gid], h[l])
                nid, _ = E.lmul("D%d" % st[2], gid)
                lam1 = mul2(E.syms[nid], h[l])
                u, v = labs_index[lag_of(E, lam0)], labs_index[lag_of(E, lam1)]
                assert D[u, v] == 1, "a move is not an edge"
                e = (min(u, v), max(u, v)); cap[e] = cap.get(e, 0) + 1
        C = ledger_step(E, C, st)
    S = [labs_index[lag_of(E, h[a])] for a in range(W)]
    T = [labs_index[lag_of(E, mul2(E.Ksym, h[a]))] for a in range(W)]
    assert all(D[S[a], T[a]] == E.m for a in range(W))
    bad = 0
    for q in Qs:
        crossing = sum(c for (u, v), c in cap.items() if q[u] != q[v])
        dem = sum(1 for a in range(W) if q[S[a]] != q[T[a]])
        if crossing < dem: bad += 1
    ok43 = (W > 2) or res["M"] >= W * E.m
    return "cut condition violations %d of %d sets; M=%d, Wm=%d, Thm 4.3 %s" % (bad, len(Qs), res["M"], W * E.m, "ok" if ok43 else "VIOLATED")

import ledger
ledger_step = ledger._step

def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    rng = random.Random(seed)
    for m in (2, 3):
        E = Engine(m)
        labs, L0, D = setup(m)
        labs_index = {x: i for i, x in enumerate(labs)}
        V = len(labs)
        if m == 2:
            Qs = [np.array([(mk >> v) & 1 for v in range(V)], dtype=bool) for mk in range(1, 1 << (V - 1))]
        else:
            Qs = [np.array([rng.random() < p for _ in range(V)], dtype=bool) for p in (0.1, 0.3, 0.5) for _ in range(1000)]
            for k in range(m): Qs.append(D[0] <= k)
            for w in range(1, 1 << (2 * m)): Qs.append(np.array([(x >> w) & 1 for x in labs], dtype=bool))
        cases = []
        cases.append(("plain W=2", 2, 2, plain_word(m, 2, rng)))
        cases.append(("pair trick", 2, 2, [("gate", [0, 1], [[1, 1], [1, -1]])] + plain_word(m, 2, rng) +
                      [("gate", [0, 1], [[Fr(1, 2), Fr(1, 2)], [Fr(1, 2), Fr(-1, 2)]])]))
        w = [("move", 0, 1), ("mix", 0, [(0, 1, None), (1, 1, None)])]
        for j in range(1, m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        w += [("mix", 0, [(0, 1, None), (1, -1, None)]), ("move", 1, 1)]
        cases.append(("gauge-exchange", 2, 2, w))
        w = [("mix", 2, [(0, 1, None)]), ("mix", 3, [(1, 1, None)])]
        w += [("move", 2, 1 << j) for j in range(m)] + [("move", 3, 1 << j) for j in range(m)]
        w += [("mix", 0, [(2, 1, None)]), ("mix", 1, [(3, 1, None)])]
        cases.append(("copy-move-copy", 2, 4, w))
        for k in range(6):
            W = rng.choice([2, 3])
            while True:
                A = [[GQ(rng.randint(-2, 2), rng.randint(-1, 1)) for _ in range(W)] for _ in range(W)]
                Ai = t_random.inv_gq_matrix(A)
                if Ai is not None: break
            cases.append(("random gate W=%d #%d" % (W, k), W, W, [("gate", list(range(W)), A)] + plain_word(m, W, rng) +
                          [("gate", list(range(W)), Ai)]))
        if m == 3:
            s = E.addr([2, 4, 1]); s2 = E.addr([4, 1, 2]); hh = Fr(1, 2)
            w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s)])]
            for j in range(m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
            w += [("mix", 1, [(0, hh, None), (1, -hh, None), (0, -hh, s), (1, hh, s), (0, hh, s2), (1, -hh, s2)]),
                  ("mix", 0, [(0, 1, None), (1, -1, None)])]
            cases.append(("order-3 control C6O3", 2, 2, w))
        for (name, W, R, w) in cases:
            print("m=%d %-24s W=%d R=%d: %s" % (m, name, W, R, check(E, W, R, w, labs_index, D, Qs, name)), flush=True)

if __name__ == "__main__":
    main()
