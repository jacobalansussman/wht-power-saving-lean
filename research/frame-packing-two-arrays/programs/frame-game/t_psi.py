"""t_psi.py: the psi screen on calibration words and on single-trunk designs."""
import random
from fpa import *
import design as D
for m in (3, 4):
    S = Sp(m)
    # plain
    states = [inp(S, 0), inp(S, 1)]
    for l in (0, 1):
        c = inp(S, l)
        for j in range(m): c = lmul(S, S.T(1 << j), c); states.append(c)
    print(f"m={m} plain 2 in 2m: psi screen", [psi_screen(S, states, l) for l in (0, 1)],
          "| cut by one move:", [psi_screen(S, states[:-1], l) for l in (0, 1)])
    rng = random.Random(m)
    tally = {}
    for tr in range(300):
        P = D.Prog(S, m, rng)
        if not P.ok: continue
        P.run()
        k = tuple(psi_screen(S, P.states, l) for l in (0, 1))
        tally[k] = tally.get(k, 0) + 1
    print(f"m={m} single trunk, no injections (300 designs): psi screen outcomes {tally}")

# the order-3 control under psi
for m in (3, 4):
    S = Sp(m)
    cols = [1 << ((k + 1) % 3) if k < 3 else 1 << k for k in range(m)]
    s3 = S.addr(cols)
    y1, y2 = inp(S, 0), inp(S, 1)
    states = [y1, y2]
    for r in (cadd((1, y1), (1, y2)), cadd((1, y1), (-1, lmul(S, s3, y2)))):
        c = r
        for j in range(m): c = lmul(S, S.T(1 << j), c); states.append(c)
    print(f"m={m} order-3 control: psi screen", [psi_screen(S, states, l) for l in (0, 1)],
          "| cut by one move:", [psi_screen(S, states[:-1], l) for l in (0, 1)])
