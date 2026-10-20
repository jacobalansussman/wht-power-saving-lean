# e3-boost calib.py: calibrate the boosting screen before any claim about an open shape.
# (a) relations known to exist must PASS every certificate; (b) relations proved impossible must be EXCLUDED by a
# rigorous certificate; (c) cases the price screen cannot see are listed as blind; (d) the booster test on made-up
# rules: the window test flags Rs and R2b, certificate A removes them; (e) the job-model exponent of Lean-pipeline
# units equals the pipeline's own figure.
import json, gzip
from fractions import Fraction as F
from math import log
import certlib as C

def show(name, tasks, calls, m=None, expect=None):
    v = C.mk(tasks, calls)
    ver, bad, res = C.verdict(v, m)
    w, at = C.window(v)
    rt = C.roots(v)
    worst = min(res, key=lambda r: float(r[3]) / max(1e-300, max(abs(float(x[3])) for x in res)))
    tag = "" if expect is None else ("  OK" if ver.startswith(expect) else "  <-- CALIBRATION FAILURE")
    print(f"{name}: D(0)={C.D(v,0):.4g} D(1)={C.D(v,1):.4g} roots in (0,1]={rt[:6]} min D on window={w:.4g} at {at:.4f}"
          f" | {ver}" + (f" by {[(n, b, s) for (n, b, s, _) in bad][:6]}" if bad else "") + tag)

print("(a) known relations: must PASS")
for k in range(2, 9): show(f"  plain [{k}] <= {k}[1]", {k: 1}, {1: k}, expect="PASSES")
for a, b in [(1, 1), (1, 2), (2, 2), (2, 3), (3, 5)]:
    show(f"  Kronecker [{a+b}] <= [{a}]+[{b}]", {a + b: 1}, ({a: 2} if a == b else {a: 1, b: 1}), expect="PASSES")
for m in range(2, 9): show(f"  A5 [1] <= 2[{m}] (register {m})", {1: 1}, {m: 2}, expect="PASSES")

def hub(name, W, m, calls, ref):
    v = C.mk({m: W}, calls)
    lo, hi = 0.99, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if C.D(v, mid) < 0: hi = mid
        else: lo = mid
    ver, bad, res = C.verdict(v, m)
    print(f"  {name}: W={W} at block {m}; J-balance root tau = {lo:.12f}, saving 1 - tau = {1 - lo:.11e} "
          f"(pipeline / earlier value {ref}); certificates: {ver}")
p = json.load(open(__import__("os").environ["T74E_JSON"]))
hub("T74E (hub #1839 word, our chain; lg-a1 price file)", p["W"], p["m"], {int(r): n for r, n in p["profile"].items()},
    f"root_lower {p['exact_rate'].get('root_lower', 'n/a') if isinstance(p.get('exact_rate'), dict) else 'n/a'}; "
    f"dft_numerator {p['dft_numerator']} -> DFT saving {p['dft_numerator']/1e10:.7e}")
L1042 = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}; I1042 = {1: 2249, 2: 1, 3: 400, 7: 3, 16: 47, 48: 4850}
calls = {r: 144 * n for r, n in L1042.items()}
for r, n in I1042.items(): calls[r] = calls.get(r, 0) + n
hub("#1042 (hub; ex-laser reprice_hub.py control)", 14336, 73, calls, "stages.py 1.674989e-3")

print("(b) proved impossible: must be EXCLUDED rigorously")
show("  engine 2 in 5 at block 3 (proved impossible)", {3: 2}, {1: 5}, expect="EXCLUDED")
show("  W2: 2 arrays, 3 rank-1 + 1 rank-2 at block 3", {3: 2}, {1: 3, 2: 1}, expect="EXCLUDED")
show("  W2: 2 arrays, 1 rank-1 + 2 rank-2 at block 3", {3: 2}, {1: 1, 2: 2}, expect="EXCLUDED")
show("  Rs 2[2] <= [1]+[3] (spreading; certificate A)", {2: 2}, {1: 1, 3: 1}, expect="EXCLUDED")
show("  R2 2[2] <= 1.05[1]+[3]", {2: 2}, {1: 1.05, 3: 1}, expect="EXCLUDED")
show("  R2b 2[4] <= [1]+[7]", {4: 2}, {1: 1, 7: 1}, expect="EXCLUDED")
for k, j in [(3, 1), (3, 2), (4, 1), (4, 3), (5, 2), (6, 1)]:
    show(f"  spreading 2[{k}] <= [{k-j}]+[{k+j}]", {k: 2}, {k - j: 1, k + j: 1}, expect="EXCLUDED")
show("  hook-free check: 1 array in 2.75 moves at block 3 (hook 2.8)", {3: 4}, {1: 11}, expect="EXCLUDED")
show("  b56: 3 arrays with four rank-2 blocks at block 6 (lift-excluded, sampled)", {6: 3}, {2: 4})
show("  block 2: 2 in 3 (engine; plain is optimal at block 2: certificate B is perfect there)", {2: 2}, {1: 3},
     expect="EXCLUDED")

print("(c) blind by construction (price screen sees no program):")
show("  1+4 in 10 at block 3 (exactly hook-tight; excluded only by ms-mixed Theorem T, a program-level argument)",
     {3: 1, 2: 4}, {1: 10})
show("  slice 3 in 11 at block 4 (Theorem 17 holds only in the commutative slice; engine certificates allow it)",
     {4: 3}, {1: 11})
print("  absorb4 limit (residual 2/3) and px-t514 (wrong at f = 6): not relations; a price screen cannot see them.")

print("(d) made-up rules of ex2-duality: the window test flags them, certificate A removes them")
for name, t, c in [("R1 2[3] <= 5[1]", {3: 2}, {1: 5}), ("Rs", {2: 2}, {1: 1, 3: 1}),
                   ("R2b", {4: 2}, {1: 1, 7: 1})]:
    v = C.mk(t, c); w, at = C.window(v)
    print(f"  {name}: D(1)={C.D(v,1):.3f}, min D on [0.9, 0.9983] = {w:.4f} -> "
          f"{'booster-shaped (D(1) >= 0, negative in window)' if C.D(v,1) >= 0 and w < 0 else ('unit' if C.D(v,1) < 0 else 'no')}")
def theta_set(rels, n=100000):
    return [i / n for i in range(n + 1) if all(C.D(C.mk(t, c), i / n) >= -1e-12 for t, c in rels)]
for name, rels in [("R1 alone", [({3: 2}, {1: 5})]), ("Rs alone", [({2: 2}, {1: 1, 3: 1})]),
                   ("R1 + Rs", [({3: 2}, {1: 5}), ({2: 2}, {1: 1, 3: 1})])]:
    S = theta_set(rels)
    print(f"  {name}: max Theta = {max(S):.4f}")
