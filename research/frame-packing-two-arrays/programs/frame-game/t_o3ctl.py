"""t_o3ctl.py: calibration in the frame game: the order-3 control (r1 = y1 + y2, r2 = y1 - s y2, K on both; s a
3-cycle of coordinates, which commutes with K) must pass as 2 in 2m; the raw-gauge ledger P, B is printed;
and the same word with r2's K cut by one move must fail."""
from fpa import *
for m in (3, 4):
    S = Sp(m)
    cols = [1 << ((k + 1) % 3) if k < 3 else 1 << k for k in range(m)]    # e0->e1->e2->e0
    s3 = S.addr(cols)
    assert S.mul(s3, S.K) == S.mul(S.K, s3)
    y1, y2 = inp(S, 0), inp(S, 1)
    r1 = cadd((1, y1), (1, y2)); r2 = cadd((1, y1), (-1, lmul(S, s3, y2)))
    states = [y1, y2]; P = 0
    for r in (r1, r2):
        c = r
        for j in range(m):
            labs = {g for (l, g) in c}
            P += len(labs) - 1
            c = lmul(S, S.T(1 << j), c); states.append(c)
    res = [end_solve(S, states, {(l, S.K): Fr(1)}) for l in (0, 1)]
    print(f"m={m}: order-3 control, end solves {[r[0] for r in res]}; raw-gauge packing P = {P} (expected m = {m})")
    for l, r in enumerate(res):
        ok, diff = verify_combo(S, states, r[1]["solution"], {(l, S.K): Fr(1)})
        print(f"   output {l}: solution of {len(r[1]['solution'])} terms re-verified exactly: {ok}")
    st2 = states[:-1]
    print(f"   r2 cut by one move: {[end_solve(S, st2, {(l, S.K): Fr(1)})[0] for l in (0, 1)]}")
