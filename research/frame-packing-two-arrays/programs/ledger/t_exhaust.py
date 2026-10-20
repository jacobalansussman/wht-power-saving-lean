#!/usr/bin/env python3
"""t_exhaust.py (sl-combinatorics): EVERY word of length <= L over a fixed alphabet at block m (default m = 2,
W = 2 inputs, R registers), run through the ledger in the trivial gauge and in one fixed random gauge.  Checks the raw
identity on every word (every prefix of a longer word is itself a shorter word), lists every exact unit found, its
ledger and its coherence class.  usage: python3 t_exhaust.py m R L [seed]"""
import sys, os, json, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ledger import *

def alphabet(E, R):
    m = E.m
    al = []
    for r in range(R):
        for z in range(1, 1 << m): al.append(("move", r, z))
    frees = ["X%d" % ((1 << m) - 1), "P%d_0" % ((1 << m) - 1)]
    if m >= 2: frees.append(E.addr([2, 1] + [1 << j for j in range(2, m)]))
    for r in range(R):
        for n in frees: al.append(("free", r, n))
    coefs = [(1, 1), (1, -1), (1, GQ(0, 1)), (-2, 1), (Fr(-1, 2), Fr(1, 2)), (0, 1)]
    for r in range(R):
        for s in range(R):
            if s == r: continue
            for (a, b) in coefs:
                terms = ([(r, a, None)] if a != 0 else []) + [(s, b, None)]
                al.append(("mix", r, terms))
    return al

def main():
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    R = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    L = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 7
    W = 2
    rng = random.Random(seed)
    E = Engine(m)
    al = alphabet(E, R)
    g_rand = [random_sp(E, rng) for _ in range(W)]
    t0 = time.time()
    n = 0; bad = []; units = []
    for length in range(0, L + 1):
        for word in itertools.product(al, repeat=length):
            word = list(word)
            for gname, g in (("trivial", None), ("random", g_rand)):
                res = run(E, W, R, word, gauge=g)
                n += 1
                if not res["raw_ok"] or (res["unit"] and not res["unit_ok"]): bad.append((gname, word))
                if gname == "trivial" and res["unit"]:
                    gc, why = find_coherent_gauge(E, W, R, word)
                    units.append(dict(word=str(word), M=res["M"], P=res["P"], B=res["B"], D=res["D"], Z=res["Z"],
                                      coherence=why))
        print("length %d done: %d runs, %d bad, %d units so far, %.0f s" % (length, n, len(bad), len(units),
              time.time() - t0), flush=True)
    byM = {}
    for u in units: byM.setdefault((u["M"], u["coherence"].split(" (")[0]), []).append(u)
    print("alphabet size", len(al), "| runs", n, "| identity violations", len(bad))
    for k in sorted(byM): print("  units with M=%d, %s: %d" % (k[0], k[1], len(byM[k])))
    minM = min((u["M"] for u in units), default=None)
    print("  smallest M among units:", minM, "(W m = %d)" % (W * m))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "t_exhaust_m%d_R%d_L%d.json" % (m, R, L)), "w") as f:
        json.dump(dict(m=m, R=R, L=L, alphabet=len(al), runs=n, bad=[str(b) for b in bad[:20]], nbad=len(bad),
                       units=units[:400], nunits=len(units)), f, indent=1)

if __name__ == "__main__":
    main()
