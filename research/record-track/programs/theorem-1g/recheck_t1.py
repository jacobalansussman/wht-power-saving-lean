# recheck_t1.py (rl-prove): re-run ex-laser's literal lockstep simulation (Lemma 2) on its four small cases and the
# exact balance/Phi check of Theorem 1 Step 1-2 at L = 30, 200 (imports ex-laser's scripts read-only from
# research/existence-routes/programs/tau-theorem-and-laser/).
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "existence-routes", "programs", "tau-theorem-and-laser"))
from lockstep_sim import run
from frontier_flow import flow, tau_c
for (W3, W2, M, L) in ((1, 4, 10, 12), (1, 4, 10, 30), (1, 1, 4, 30), (1, 4, 10, 54), (7, 38, 88, 12)):
    st, S, lv = run(L, W3, W2, M)
    print(f"lockstep {W3}+{W2} in {M}, L={L}: classes {st['classes']}, sources used == Source: {st['src_used'] == S}, leaves == M u(0): {st['leaves'] == lv}", flush=True)
for (W3, W2, M) in ((1, 4, 10), (7, 38, 88)):
    tc = tau_c(W3, W2, M)
    for L in (30, 200):
        U, src, leaves, bad = flow(L, W3, W2, M, exact=True)
        print(f"balance {W3}+{W2} in {M}, L={L}: unbalanced {len(bad)}, leaves {leaves}, Phi(src) {sum(float(w)*s**tc for s, w in src.items()):.10f}, "
              f"sources in ({min(src)}, {max(src)}] vs (L/3, L] = ({L/3:.1f}, {L}]", flush=True)
