# ex2-duality calib.py: calibration of the repaired duality on known cases.
# Part A (exact-ish arithmetic, instant): for a set R of generating relations of the job semiring J
#   (each relation: tasks {size: weight} served by resources {size: weight}), the spectral points of <R> are
#   phi_theta([n]) = n^theta for theta in Theta(R) = {theta in [0,1] : sum_res w r^theta >= sum_task w k^theta for all R}.
#   The duality says the best exponent reachable from R is max Theta(R).  We compare max Theta with the exponents the
#   project computed independently (ex-laser Theorem 1 tau_c, the engine's moment root).
# Part B (LP, seconds): ex-laser's fluid LP (grouped_lp.py, imported read-only) on
#   (i) its grouped test rule and the rank-1 rule it reduces to by recycling (prediction: identical LP values);
#   (ii) a mixed rule whose grouped paid rank (2) is NOT a target rank (prediction from the duality: exponent -> tau_c).
import sys
from math import log
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "tau-theorem-and-laser"))

def D(rel, t):
    tasks, res = rel
    return sum(w * r ** t for r, w in res.items()) - sum(w * k ** t for k, w in tasks.items())

def theta_set(rels, grid=20000):
    ok = [all(D(r, i / grid) >= -1e-12 for r in rels) for i in range(grid + 1)]
    # report max theta and whether the set is an interval [0, max]
    mx = max(i for i in range(grid + 1) if ok[i]) / grid
    interval = all(ok[: int(mx * grid) + 1])
    return mx, interval

KRON = [({a + b: 1}, ({a: 2} if a == b else {a: 1, b: 1})) for a in range(1, 6) for b in range(1, 6)]   # [a+b] <= [a] + [b]
PLAIN = [({m: 1}, {1: m}) for m in range(2, 9)]                                    # [m] <= m [1]
A5 = [({1: 1}, {m: 2}) for m in range(2, 9)]                                       # [1] <= 2 [m] (lower boundedness)
BASE = KRON + PLAIN + A5

cases = [
    ("base only (plain, Kronecker, lower boundedness)", BASE, None),
    ("+ 1+4 in 10 at block 3 (rank-1 moves; ex-laser tau_c 0.8811)", BASE + [({3: 1, 2: 4}, {1: 10})], 0.8811),
    ("+ grouped test rule {3:1,2:4} <= {1:6,2:2} (ex-laser tau_c 0.8214)", BASE + [({3: 1, 2: 4}, {1: 6, 2: 2})], 0.8214),
    ("+ its recycled rank-1 form {3:1,2:2} <= {1:6}", BASE + [({3: 1, 2: 2}, {1: 6})], 0.8214),
    ("+ grouped mixed {4:1,3:2} <= {1:5,2:2} (paid rank 2 not a target rank)", BASE + [({4: 1, 3: 2}, {1: 5, 2: 2})], None),
    ("+ pure grouped 2[3] <= 3[1] + [2] (engine moment root)", BASE + [({3: 2}, {1: 3, 2: 1})], None),
    ("+ hypothetical engine 2 in 5 at block 3 (proved impossible as a finite unit)", BASE + [({3: 2}, {1: 5})], None),
]
print("Part A: max Theta(R) = best exponent the duality allows from R; saving = 1 - max Theta")
for name, rels, ref in cases:
    mx, itv = theta_set(rels)
    print(f"  {name}: max Theta = {mx:.4f} (saving {1 - mx:.4f}); Theta is [0, max]: {itv}"
          + (f"; independent value {ref}" if ref else ""))
print("  2 in 5 at block 3 is in the ASYMPTOTIC closure of J iff theta* <= log_3(5/2) =", round(log(2.5) / log(3), 4),
      "; known: theta* <= 1 - 1.68e-3 and theta* >= 0 only, so the duality cannot refuse it.")

if len(sys.argv) > 1 and sys.argv[1] == "lp":
    from grouped_lp import solve, tau_c
    Ns = [int(x) for x in sys.argv[2:]] or [60, 120, 240, 480]
    lp_cases = [("grouped test rule {3:1,2:4} <= {1:6,2:2}", {3: 1, 2: 4}, {1: 6, 2: 2}),
                ("recycled rank-1 form {3:1,2:2} <= {1:6}", {3: 1, 2: 2}, {1: 6}),
                ("grouped mixed {4:1,3:2} <= {1:5,2:2}", {4: 1, 3: 2}, {1: 5, 2: 2})]
    print("Part B: fluid LP (ex-laser grouped_lp.solve), n =", Ns)
    for name, Wt, P in lp_cases:
        tc = tau_c(Wt, P)
        costs = [solve(N, Wt, P) for N in Ns]
        print(f"  {name}: tau_c {tc:.4f}; LP cost {[round(c, 6) for c in costs]}; cost/n^tau_c "
              f"{[round(c / N ** tc, 4) for N, c in zip(Ns, costs)]}; local slopes "
              f"{[round(log(costs[i + 1] / costs[i]) / log(Ns[i + 1] / Ns[i]), 4) for i in range(len(Ns) - 1)]}", flush=True)
