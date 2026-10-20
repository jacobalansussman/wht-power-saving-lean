# e3-boost check_survivor.py: the smallest block-6 survivor against every certificate in certlib (all statuses,
# including the projections of every symmetric family from projall.py and the sampled lifts), the walk condition,
# its balance and window, and its margins.
import certlib as C
base = [(n, M, P, s) for (n, M, P, s) in list(C.PROFILES) if 3 <= M <= 7 and s not in ("sampled", "conj")
        and not n.startswith("P")]
for (n, M, P, s) in base:
    cur = list(P)
    for M2 in range(M + 1, 13):
        cur = C.project(cur)
        if max(cur) > 0: C.add(f"Q{n}to{M2}", M2, cur, "argued2")
ALL = C.RIGOROUS + ("argued", "argued2", "sampled", "conj")
for tasks, calls, reg0 in [({2: 1, 5: 4}, {3: 2, 4: 1, 6: 2}, 6), ({2: 1, 5: 1, 6: 1}, {3: 2, 7: 1}, 7),
                           ({2: 1, 5: 1, 6: 2}, {3: 1, 4: 2, 8: 1}, 8)]:
    v = C.mk(tasks, calls)
    for reg in range(reg0, 9):
        res = C.certify(v, reg, ALL)
        bad = [r for r in res if r[3] < 0]
        tight = sorted(set(f"{n}x{b}" for (n, b, s, mg) in res if mg == 0))
        print(f"tasks {tasks} calls {calls} on a register of {reg} rows: {len(res)} certificate tests, failed: "
              f"{[(n, b, s, float(mg)) for (n, b, s, mg) in bad][:5]}; tight: {tight[:12]}")
    print("   walk condition:", C.walk_ok(v, reg0), "; D(0) =", C.D(v, 0), "; D(1) =", C.D(v, 1), "; roots:", C.roots(v),
          "; D(0.9) =", round(C.D(v, 0.9), 6), "; D(0.9983) =", round(C.D(v, 0.9983), 8))
