"""vhit.py (v-fp-build): rebuild fp-build's closest design (remote-out/hit_m4_s102_4.pkl) with my own frame-game code
(vfg.py), recompute its psi residual, and test whether the psi 'pass' is attained or only approached as the trunk's
y1 weight 1 + c0 -> 0 (a degenerate limit in which the trunk packs nothing).
Semantics replayed from fp-build/design.py Prog.run (read, not imported)."""
import sys
import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)); sys.path.insert(0, _HERE)
import numpy as np
from vfg import FG, cadd
from loadhit import load

d = load(sys.argv[1] if len(sys.argv) > 1 else
         _os.path.join(_HERE, "..", "..", "..", "data", "frame-game", "design-runs", "hit_m4_s102_4.pkl"))
m = d["m"]; S = FG(m)
p2 = S.from_rows(d["p2"])
inj = [(rg, sl, src, S.from_rows(mu)) for (rg, sl, src, mu) in d["inj"]]
ident = S.I
print("injections whose chart is the identity:", sum(1 for x in inj if x[3] == ident), "of", len(inj))

def build(c):
    states = [{(0, S.I): 1.0 + 0j}, {(1, S.I): 1.0 + 0j}]
    out = {}; selfinj = []
    for reg in d["order"]:
        if reg == "R1": cur = dict(states[0])
        elif reg == "R2": cur = dict(states[1])
        elif reg == "Rt":
            c1 = states[out[("R1", d["a1"] - 1)]] if d["a1"] else states[0]
            c2 = states[out[("R2", d["a2"] - 1)]] if d["a2"] else states[1]
            cur = cadd((1, c1), (1, S.lmul(p2, c2)))
        else:
            cur = dict(states[out[("Rt", len(d["regs"]["Rt"]) - 1)]])
        for k, u in enumerate(d["regs"][reg]):
            for i, (rg, sl, src, mu) in enumerate(inj):
                if rg == reg and sl == k:
                    add = S.lmul(mu, states[src])
                    # is the injected content proportional to the register's current content? (a pure rescaling)
                    keys = set(kk for kk, v in add.items() if abs(v) > 0)
                    ck = set(kk for kk, v in cur.items() if abs(v) > 1e-300)
                    selfinj.append((i, keys <= ck))
                    cur = cadd((1, cur), (c[i], add))
            cur = S.lmul(S.tv(u), cur)
            states.append(cur); out[(reg, k)] = len(states) - 1
    return states, selfinj

def psi_res(states):
    keys = sorted({(l, S.pos(g)) for st in states for (l, g) in st} | {(l, S.pos(S.K)) for l in (0, 1)}, key=str)
    ki = {k: i for i, k in enumerate(keys)}
    A = np.zeros((len(keys), len(states)), complex)
    for r, st in enumerate(states):
        for (l, g), v in st.items(): A[ki[(l, S.pos(g))], r] += v
    out = []
    for l in (0, 1):
        b = np.zeros(len(keys), complex); b[ki[(l, S.pos(S.K))]] = 1
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        out.append((float(np.linalg.norm(A @ x - b)), float(np.max(np.abs(x)))))
    return out

c = np.array(d["c"], complex)
states, selfinj = build(c)
print("injections that only add the register's own current walkers (rescalings):",
      [i for i, s in selfinj if s], "->", sum(1 for i, s in selfinj if s), "of", len(selfinj))
print("effective scalars: trunk y1 weight 1+c0 =", np.round(1 + c[0], 9),
      "| A slot 0: 1+c5+c6+c11 =", np.round(1 + c[5] + c[6] + c[11], 4), " landing c1 =", np.round(c[1], 4))
print("psi residual at fp-build's scalars (per output: residual, max |readout coeff|):", psi_res(states))
print("trunk walkers (state 5):", {(l, S.pos(g)): np.round(v, 9) for (l, g), v in states[5].items()})
for t in (1.0, 1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 0.0):
    cc = c.copy(); cc[0] = -1 + t * np.exp(0.7j)
    r = psi_res(build(cc)[0])
    print(f"  |1+c0| = {t:.0e}: psi residuals (out0, out1) = ({r[0][0]:.3e}, {r[1][0]:.3e}); max readout coeff {max(r[0][1], r[1][1]):.2e}")
# other scalars random, c0 near -1: is the small residual driven by c0 alone?
rng = np.random.default_rng(5)
for trial in range(3):
    cc = rng.standard_normal(len(c)) + 1j * rng.standard_normal(len(c)); cc[0] = -1 + 1e-6
    r = psi_res(build(cc)[0])
    print(f"  random other scalars, |1+c0| = 1e-6: psi residuals ({r[0][0]:.3e}, {r[1][0]:.3e})")
cc = c.copy(); cc[0] = -1
st0 = build(cc)[0]
print("with c0 = -1 exactly, input-0 walkers outside the raw input state:",
      sum(1 for st in st0[2:] for (l, g), v in st.items() if l == 0 and abs(v) > 1e-12))
# exact support check: which stored states touch the psi coordinates (0, X1) [trunk's y1 position] and (1, L1)?
from fractions import Fraction as Fr
import random as _r
rr = _r.Random(3)
cq = [Fr(rr.randint(1, 97), rr.randint(1, 97)) * rr.choice([1, -1]) for _ in range(len(c))]
stq, _ = build(cq)
X1 = [S.pos(g) for (l, g) in stq[5] if l == 0][0]; L1 = S.pos(S.K)
def psi_of(st):
    out = {}
    for (l, g), v in st.items():
        k = (l, S.pos(g)); out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}
touch = [(r, psi_of(st).get((0, X1), 0), psi_of(st).get((1, L1), 0)) for r, st in enumerate(stq)]
print("states touching (0, X1) or (1, L1) [state, coeff at (0,X1), coeff at (1,L1)], exact, random rational scalars:",
      [(r, str(a), str(b)) for r, a, b in touch if a != 0 or b != 0], " (1+c0 =", 1 + cq[0], ")")
print("level of X1:", S.dist(S.L0, X1), "; input-0 walkers by state:", [sum(1 for (l, g) in st if l == 0) for st in stq])
