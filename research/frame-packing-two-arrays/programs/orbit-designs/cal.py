#!/usr/bin/env python3
"""cal.py (pk2-orbits): calibration of the orbit pipeline (psi, linearised chart screen, heightfilter shim, exact).
1. plain 2 in 2m (m = 3, 4, 5): must pass psi + exact; heightfilter must find a height system (cannot save).
2. order-3 control as an orbit design: r1 = y1 + y2, r2 = y1 - sigma y2, K on both, sigma in the two-frame
   stabiliser P cap K^-1 P K.  Prediction: exact YES iff 1 + sigma' (sigma' = K sigma K^-1) is invertible in
   C[<sigma'>], i.e. iff ord(sigma) is odd; even order must FAIL (1 + sigma' is a zero divisor).  Also the 3-cycle
   address map (o3-build's word).
3. cut: r2's K cut by one move must fail; heightfilter on the control finds a height system (it does not save).
"""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "heights")); sys.path.insert(0, os.path.join(HERE, "..", "frame-game")); sys.path.insert(0, HERE)
from fpa import *
import orb

def run_regs(S, regs):
    states = [inp(S, 0), inp(S, 1)]
    for c in regs:
        for j in range(S.m):
            c = lmul(S, S.T(1 << j), c); states.append(c)
    return states

def verdicts(S, F, states, W=2):
    ps = "".join("1" if psi_screen(S, states, l) else "0" for l in range(W))
    ex = []
    for l in range(W):
        r = end_solve(S, states, {(l, S.K): Fr(1)})
        if r[0] == "yes":
            ok, _ = verify_combo(S, states, r[1]["solution"], {(l, S.K): Fr(1)})
            ex.append("Y%d" % len(r[1]["solution"]) if ok else "?")
        else: ex.append(r[0])
    hv, bb = orb.height_verdict(S, F, states, W)
    return ps, ex, hv, bb

rng = random.Random(7)
for m in (3, 4, 5):
    S = Sp(m); F = orb.FShim(S)
    st = run_regs(S, [inp(S, 0), inp(S, 1)])
    print(f"m={m} plain 2 in {2*m}: psi={verdicts(S, F, st)[0]} exact={verdicts(S, F, st)[1]} "
          f"height system={verdicts(S, F, st)[2]} Theorem-C bound={verdicts(S, F, st)[3]}", flush=True)
    st2 = st[:-1]
    v = verdicts(S, F, st2)
    print(f"   plain cut by one move: psi={v[0]} exact={v[1]}", flush=True)
    # order-3 address map (o3-build)
    if m >= 3:
        cols = [1 << ((k + 1) % 3) if k < 3 else 1 << k for k in range(m)]
        s3 = S.addr(cols)
        st = run_regs(S, [cadd((1, inp(S, 0)), (1, inp(S, 1))), cadd((1, inp(S, 0)), (-1, lmul(S, s3, inp(S, 1))))])
        v = verdicts(S, F, st)
        print(f"   o3 control (3-cycle address map): psi={v[0]} exact={v[1]} height system={v[2]} bound={v[3]}", flush=True)
        v = verdicts(S, F, st[:-1])
        print(f"   o3 control cut by one move: psi={v[0]} exact={v[1]}", flush=True)
    # two-frame stabiliser sigma of various orders
    seen = {}
    for trial in range(400):
        sg = orb.twoframe(S, rng)
        n = orb.order(S, sg, cap=40)
        if n is None or n < 2 or seen.get(n, 0) >= 2: continue
        seen[n] = seen.get(n, 0) + 1
        st = run_regs(S, [cadd((1, inp(S, 0)), (1, inp(S, 1))), cadd((1, inp(S, 0)), (-1, lmul(S, sg, inp(S, 1))))])
        v = verdicts(S, F, st)
        print(f"   two-frame sigma of order {n}: psi={v[0]} exact={v[1]} (prediction: {'YES' if n % 2 else 'NO'}) "
              f"height system={v[2]}", flush=True)
    sys.stdout.flush()
