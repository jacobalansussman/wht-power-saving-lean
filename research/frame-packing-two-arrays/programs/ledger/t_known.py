#!/usr/bin/env python3
"""t_known.py (sl-combinatorics): the ledger on known / constructed exact units at blocks 1-4.
Every case: exact unit check, raw identity, unit-form identity, in the trivial gauge, in the coherent gauge (if one
exists) and in random symplectic gauges.  usage: python3 -I t_known.py [seed]"""
import sys, os, json, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ledger import *

def report(name, E, W, R, word, rng, ngauges=4, start=None, deliver=None, out=None):
    rows = []
    g_coh, why = find_coherent_gauge(E, W, R, word, start=start)
    gauges = [("trivial", gauge_trivial(E, W))]
    if g_coh is not None: gauges.append(("coherent", g_coh))
    for k in range(ngauges): gauges.append(("random%d" % k, [random_sp(E, rng) for _ in range(W)]))
    for gname, g in gauges:
        res = run(E, W, R, word, start=start, deliver=deliver, gauge=g)
        row = dict(case=name, m=E.m, W=W, R=R, gauge=gname, M=res["M"], unit=res["unit"], raw_ok=res["raw_ok"],
                   unit_ok=res.get("unit_ok"), N=res.get("N"), X=res.get("X"), D=res["D"], Z=res["Z"],
                   B=res.get("B"), P=res["P"], I0=res["I0"], maxlab=res["maxlab"])
        rows.append(row)
    print("== %s (m=%d, W=%d, R=%d, M=%d): %s" % (name, E.m, W, R, rows[0]["M"], why))
    for r in rows:
        print("   %-9s unit=%s raw_ok=%s unit_ok=%s | N=%s X=%s D=%s Z=%s B=%s P=%s I0=%s maxlab=%s"
              % (r["gauge"], r["unit"], r["raw_ok"], r["unit_ok"], r["N"], r["X"], r["D"], r["Z"], r["B"], r["P"],
                 r["I0"], r["maxlab"]))
    if out is not None: out.append(dict(case=name, coherence=why, rows=rows))
    return rows

def free_fix(E, gid, target, maxdepth=4, names=None):
    """BFS over free letters: a word b with b . gid = scalar * target"""
    if names is None: names = [n for n in E.letters if not n.startswith("D")]
    from collections import deque
    seen = {gid: []}; dq = deque([gid])
    while dq:
        g = dq.popleft()
        if g == target: return seen[g]
        if len(seen[g]) >= maxdepth: continue
        for n in names:
            h, _ = E.lmul(n, g)
            if h not in seen: seen[h] = seen[g] + [n]; dq.append(h)
    return None

def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    rng = random.Random(seed)
    out = []
    for m in (1, 2, 3, 4):
        E = Engine(m)
        for W in (1, 2, 3):
            report("plain", E, W, W, plain_word(m, W, rng), rng, out=out)
        # detour: dir z twice = X^z, fixed by the free shift
        w = plain_word(m, 2, rng); z = rng.randrange(1, 1 << m)
        w = w[:2] + [("move", 0, z), ("move", 0, z), ("free", 0, "X%d" % z)] + w[2:]
        report("detour(z=%d)" % z, E, 2, 2, w, rng, out=out)
        # pair trick: gate at the raw frame, plain, inverse gate at the target frame
        A = [[1, 1], [1, -1]]; Ai = [[Fr(1, 2), Fr(1, 2)], [Fr(1, 2), Fr(-1, 2)]]
        report("pair trick", E, 2, 2, [("gate", [0, 1], A)] + plain_word(m, 2, rng) + [("gate", [0, 1], Ai)], rng, out=out)
        # copy-move-copy back (overwrites, births and deaths)
        w = [("mix", 2, [(0, 1, None)]), ("mix", 3, [(1, 1, None)])]
        w += [("move", 2, 1 << j) for j in range(m)] + [("move", 3, 1 << j) for j in range(m)]
        w += [("mix", 0, [(2, 1, None)]), ("mix", 1, [(3, 1, None)])]
        report("copy-move-copy", E, 2, 4, w, rng, out=out)
        if m >= 2:
            # gauge exchange: line source for input 1, paid back by one descent (coherent gauge) / packing (trivial)
            w = [("move", 0, 1), ("mix", 0, [(0, 1, None), (1, 1, None)])]
            for j in range(1, m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
            w += [("mix", 0, [(0, 1, None), (1, -1, None)]), ("move", 1, 1)]
            report("gauge-exchange", E, 2, 2, w, rng, out=out)
            # a lateral move: after K, dir(e0+e1) keeps the level; undo with free letters
            gK = E.Kid
            g2, _ = E.lmul("D3", gK)
            fix = free_fix(E, g2, gK)
            if fix is not None:
                w = plain_word(m, 1, rng) + [("move", 0, 3)] + [("free", 0, n) for n in fix]
                report("lateral(D3 after K; fix %s)" % fix, E, 1, 1, w, rng, out=out)
            else:
                print("== lateral: no free fix within depth 4 at m=%d" % m)
        if m >= 3:
            # the order-3 control (o3-build C8O3 type): joint lineages K(x1 + x2), K(x1 - s x2), order-3 halving readout
            s = E.addr([2, 4, 1] + [1 << j for j in range(3, m)]); s2 = E.addr([4, 1, 2] + [1 << j for j in range(3, m)])
            w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s)])]
            for j in range(m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
            h = Fr(1, 2)
            w += [("mix", 1, [(0, h, None), (1, -h, None), (0, -h, s), (1, h, s), (0, h, s2), (1, -h, s2)]),
                  ("mix", 0, [(0, 1, None), (1, -1, None)])]
            report("order-3 control C%dO3" % (2 * m), E, 2, 2, w, rng, out=out)
            # same with the +-1 halving (pair trick written as mixes): coherent
            w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -2, None)])]
            for j in range(m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
            w += [("mix", 1, [(0, h, None), (1, -h, None)]), ("mix", 0, [(0, 1, None), (1, -1, None)])]
            report("+-1 halving control C%dH" % (2 * m), E, 2, 2, w, rng, out=out)
        if m >= 2:
            # wasteful packing gadget: one move advances two labels, later paid back
            w = [("mix", 2, [(1, 1, None)]), ("mix", 0, [(0, 1, None), (1, 1, None)]), ("move", 0, 1),
                 ("mix", 0, [(0, 1, None), (2, -1, None)]), ("move", 0, 2), ("move", 2, 2),
                 ("mix", 0, [(0, 1, None), (2, 1, None)]), ("move", 1, 1), ("move", 1, 2),
                 ("mix", 0, [(0, 1, None), (1, -1, None)])]
            for j in range(2, m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
            report("packing gadget", E, 2, 3, w, rng, out=out)
    bad = [(c["case"], r["gauge"]) for c in out for r in c["rows"] if not (r["raw_ok"] and (r["unit_ok"] in (True, None)))]
    nonunit = [c["case"] for c in out if not c["rows"][0]["unit"]]
    print("\nSUMMARY: %d cases, %d (case, gauge) rows; identity violations: %s; non-units: %s"
          % (len(out), sum(len(c["rows"]) for c in out), bad, nonunit))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "t_known.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)

if __name__ == "__main__":
    main()
