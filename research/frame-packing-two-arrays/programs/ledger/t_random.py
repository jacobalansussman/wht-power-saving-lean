#!/usr/bin/env python3
"""t_random.py (sl-combinatorics): (1) the raw ledger identity on random words (non-units) at blocks 2-4 in random
gauges; (2) random exact units of the packing type (joint lineages with an odd-order coordinate permutation s and the
readout (1+s)^-1 = (1 - s + s^2)/2), W = 2 and 3, random move orders and interleavings; (3) random coherent units of
the pair-trick type (random invertible Gaussian-integer gate, plain moves, inverse gate).
usage: python3 t_random.py [seed] [nwords]"""
import sys, os, json, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ledger import *

def random_word(E, R, L, rng, frees):
    w = []
    for _ in range(L):
        t = rng.random()
        r = rng.randrange(R)
        if t < 0.4: w.append(("move", r, rng.randrange(1, 1 << E.m)))
        elif t < 0.65: w.append(("free", r, rng.choice(frees)))
        else:
            k = rng.randint(1, 3); terms = []
            for _ in range(k):
                terms.append((rng.randrange(R), rng.choice([1, -1, Fr(1, 2), Fr(-1, 2), GQ(0, 1), GQ(1, 1), 2]),
                              rng.choice([None, None, rng.choice(frees)])))
            w.append(("mix", r, terms))
    return w

def inv_gq_matrix(A):
    n = len(A); M = [[gq(A[i][j]) for j in range(n)] + [GQ(1) if i == j else GQ(0) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next((i for i in range(c, n) if not M[i][c].iszero()), None)
        if p is None: return None
        M[c], M[p] = M[p], M[c]
        iv = M[c][c].inv(); M[c] = [x * iv for x in M[c]]
        for i in range(n):
            if i != c and not M[i][c].iszero():
                f = M[i][c]; M[i] = [M[i][j] - f * M[c][j] for j in range(2 * n)]
    return [row[n:] for row in M]

def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    nwords = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    rng = random.Random(seed)
    stats = dict(raw_words=0, raw_bad=0, pack_units=0, pack_bad=0, pack_nonunit=0, gate_units=0, gate_bad=0,
                 gate_nonunit=0)
    examples = []
    for m in (2, 3, 4):
        E = Engine(m)
        frees = [n for n in E.letters if not n.startswith("D")]
        frees += [E.addr([1 << ((j + 1) % m) for j in range(m)])]          # a cyclic coordinate shift
        for W in (2, 3):
            R = W + 1
            for _ in range(nwords if m < 4 else nwords // 3):
                w = random_word(E, R, rng.randint(4, 14), rng, frees)
                g = [random_sp(E, rng) for _ in range(W)] if rng.random() < 0.7 else gauge_trivial(E, W)
                res = run(E, W, R, w, gauge=g)
                stats["raw_words"] += 1
                if not res["raw_ok"]: stats["raw_bad"] += 1; examples.append(("raw_bad", m, W, w))
        if m >= 3:
            cyc = []
            for trip in itertools.combinations(range(m), 3):
                a, b, c = trip
                cols = [1 << j for j in range(m)]; cols[a], cols[b], cols[c] = 1 << b, 1 << c, 1 << a
                s = E.addr(cols)
                cols2 = [1 << j for j in range(m)]; cols2[a], cols2[b], cols2[c] = 1 << c, 1 << a, 1 << b
                cyc.append((s, E.addr(cols2)))
            for _ in range(nwords // 4):
                s, s2 = rng.choice(cyc)
                W = rng.choice([2, 3])
                # registers 0..W-1: r0 = x0 + x1 (+ x2), r1 = x0 - s x1, r2 = x0 - s x2 (W=3)
                w = []
                if W == 2:
                    w += [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s)])]
                else:
                    # r0 = x0 + x1 + x2 ; r1 = x0 + x1 + x2 - (1+s) x1 ; r2 = x0 + x1 + x2 - (1+s) x2
                    w += [("mix", 0, [(0, 1, None), (1, 1, None), (2, 1, None)]),
                          ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s)]),
                          ("mix", 2, [(0, 1, None), (2, -1, None), (2, -1, s)])]
                moves = []
                for r in range(W):
                    order = list(range(m)); rng.shuffle(order)
                    moves += [("move", r, 1 << j) for j in order]
                rng.shuffle(moves)
                w += moves
                h = Fr(1, 2)
                inv1 = lambda r: [(0, h, None), (r, -h, None), (0, -h, s), (r, h, s), (0, h, s2), (r, -h, s2)]
                if W == 2:
                    w += [("mix", 1, inv1(1)), ("mix", 0, [(0, 1, None), (1, -1, None)])]
                else:
                    w += [("mix", 1, inv1(1)), ("mix", 2, inv1(2)), ("mix", 0, [(0, 1, None), (1, -1, None), (2, -1, None)])]
                for gi in range(3):
                    g = gauge_trivial(E, W) if gi == 0 else [random_sp(E, rng) for _ in range(W)]
                    res = run(E, W, W, w, gauge=g)
                    if not res["unit"]: stats["pack_nonunit"] += 1; break
                    stats["pack_units"] += 1
                    if not (res["raw_ok"] and res["unit_ok"]): stats["pack_bad"] += 1
                    if gi == 0 and len(examples) < 12:
                        examples.append(("pack", m, W, dict(M=res["M"], P=res["P"], B=res["B"], D=res["D"], Z=res["Z"])))
        for _ in range(nwords // 4):
            W = rng.choice([2, 3])
            while True:
                A = [[GQ(rng.randint(-2, 2), rng.randint(-1, 1)) for _ in range(W)] for _ in range(W)]
                Ai = inv_gq_matrix(A)
                if Ai is not None: break
            w = [("gate", list(range(W)), A)] + plain_word(m, W, rng) + [("gate", list(range(W)), Ai)]
            for gi in range(2):
                g = gauge_trivial(E, W) if gi == 0 else [random_sp(E, rng) for _ in range(W)]
                res = run(E, W, W, w, gauge=g)
                if not res["unit"]: stats["gate_nonunit"] += 1; break
                stats["gate_units"] += 1
                if not (res["raw_ok"] and res["unit_ok"]): stats["gate_bad"] += 1
        print("m=%d done: %s" % (m, stats), flush=True)
    print("FINAL", json.dumps(stats))
    for e in examples[:12]: print("example", e[0], e[1], e[2], e[3] if e[0] == "pack" else "")
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "t_random_s%d.json" % seed), "w") as f:
        json.dump(dict(stats=stats, examples=[str(e) for e in examples]), f, indent=1)

if __name__ == "__main__":
    main()
