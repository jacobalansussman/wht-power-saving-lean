"""verify_ctl.py -- o3-build: full verification pipeline on the order-3 halving-end control word (C8O3, 6-move prefix
+ 2 end moves = a plain 2-in-8 in the basis y = (x1 + x2, x1 - s x2)): exact end solve, formal re-check, raw complex
arrays at f = 1..5 (f >= 6 needs >= 16^6 entries per array; the formal identity covers every f).  All letters are
column-uniform (address maps / progress letters), so the automaton test for column-mixing letters does not apply."""
import sys, os, random, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o3lib import *
import search as S
from end2 import end2, patterns, check_solution, raw_check_end2
w = W4(3); rng = random.Random(3)
B, sn = S.tmpl("C8O3", w, rng)
G = w.G
Es = S.E_pool(w, B.gens, 4)
s = w.sig[sn]
for i6 in range(len(Es)):
    al = patterns(G, s)["O3"]
    st, sol, nc = end2(G, B.gens, [Es[i6], Es[i6]], al, 3000)
    if st == "yes":
        print(f"C8O3 control (s = {sn}): O3 end with E6 = E7 = E[{i6}]: exact solve yes ({nc} columns); "
              f"formal re-check {check_solution(G, B.gens, [Es[i6], Es[i6]], al, sol)}", flush=True)
        dens = sorted({Fraction(x).denominator for (_, x) in sol})
        print("   denominators in the solution:", dens, flush=True)
        r = raw_check_end2(w, B.prog, B.gens, [Es[i6], Es[i6]], al, sol, fs=(1, 2, 3, 4, 5))
        print("   raw complex arrays, max relative error by f:", {k: float('%.2e' % v) for k, v in r.items()}, flush=True)
        priv = w.end_private(B.gens)[1]
        print("   same prefix, private end in the basis x:", priv, flush=True)
        json.dump({"s": sn, "prog": B.prog, "raw": r, "priv_x": priv}, open("out_verify_ctl.json", "w"))
        break
