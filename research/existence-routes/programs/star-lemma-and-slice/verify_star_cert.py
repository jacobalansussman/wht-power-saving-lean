#!/usr/bin/env python3
"""verify_star_cert.py (e3-slice6): independent check of the block-5 star certificate (M >= 93/20 W).

Rebuilds the binary quadratic forms on F_2^5 itself (no code from ex2-families), checks the MEANING of every row of the
certificate, and adds the rows up in exact fractions.
  STAR-L row :  2 c_x + 2 c_y - w_P1 - w_P2 - w_P3 <= 0, where {x, y} is a line and P1, P2, P3 are the three planes
               through it inside an even 3-flat F that is the support of an EVEN form (E[F,F] nonsingular, Theorem 10b);
               the planes partition F minus the line (star lemma, k = 0).
  bound:w row:  w_P - sum_{z in S} c_z <= -1 with S the support of an ODD form containing P (Theorem 10c), or
               w_P - sum_{z in P} c_z <= 0.
Variables: c_z = moves in direction z+1 per array; w_P = rank E[P,P] / W.
Conclusion checked: sum_i lambda_i row_i reads  sum_z beta_z c_z <= -93/20 with every beta_z >= -1 and every w cancelled,
so sum_z c_z >= 93/20.
Usage: verify_star_cert.py cert.json vars.json
"""
import sys, json, itertools
from fractions import Fraction as Fr

m = 5
pts = range(1, 1 << m)
bit = lambda z, j: (z >> j) & 1
pairs = [(j, k) for j in range(m) for k in range(j, m)]
supports = {}                      # frozenset(direction indices z-1) -> set of traces of forms with that support
for coef in itertools.product((0, 1), repeat=len(pairs)):
    S = frozenset(z - 1 for z in pts if sum(c * bit(z, j) * bit(z, k) for c, (j, k) in zip(coef, pairs)) % 2)
    tr = sum(c for c, (j, k) in zip(coef, pairs) if j == k) % 2
    supports.setdefault(S, set()).add(tr)
odd_supports = {S for S, t in supports.items() if 1 in t}
even_supports = {S for S, t in supports.items() if 0 in t}

def is_flat(idx, dim):
    P = [z + 1 for z in idx]
    if len(P) != 1 << dim: return False
    a = P[0]; D = {p ^ a for p in P}
    return all((u ^ v) in D for u in D for v in D)

cert = json.load(open(sys.argv[1])); var = json.load(open(sys.argv[2]))
key = {int(k): v for k, v in var.items()}
total = {}; rhs_total = Fr(0); bad = []
for r in cert:
    lam = Fr(r["lambda"]); co = {int(j): Fr(v).limit_denominator(10) for j, v in r["coeffs"].items()}; rhs = Fr(r["rhs"]).limit_denominator(10)
    if lam < 0: bad.append(("negative lambda", r["row"]))
    cs = {key[j][1]: v for j, v in co.items() if key[j][0] == "c"}
    ws = {frozenset(key[j][1]): v for j, v in co.items() if key[j][0] == "w"}
    others = [key[j] for j in co if key[j][0] not in ("c", "w")]
    if others: bad.append(("unexpected variable", r["row"], others))
    if r["kind"] == "STAR-L":
        ok = (rhs == 0 and len(cs) == 2 and all(v == 2 for v in cs.values()) and len(ws) == 3 and all(v == -1 for v in ws.values()))
        line = frozenset(cs)
        planes = list(ws)
        F = frozenset().union(*planes) if planes else frozenset()
        ok = ok and all(is_flat(sorted(P), 2) and line <= P for P in planes) and is_flat(sorted(F), 3)
        ok = ok and F in even_supports
        # the three planes through the line inside F partition F minus the line
        ok = ok and sorted(len(P - line) for P in planes) == [2, 2, 2] and len(F - line) == 6 and \
            len(frozenset().union(*[P - line for P in planes])) == 6
        if not ok: bad.append(("STAR-L meaning", r["row"]))
    elif r["kind"] == "bound:w":
        ok = len(ws) == 1 and list(ws.values())[0] == 1 and all(v == -1 for v in cs.values())
        P = list(ws)[0]; S = frozenset(cs)
        ok = ok and is_flat(sorted(P), 2)
        if rhs == -1: ok = ok and S in odd_supports and P <= S
        elif rhs == 0: ok = ok and S == P
        else: ok = False
        if not ok: bad.append(("bound:w meaning", r["row"]))
    else:
        bad.append(("unknown kind", r["kind"], r["row"]))
    for j, v in co.items(): total[j] = total.get(j, Fr(0)) + lam * v
    rhs_total += lam * rhs
w_left = {j: v for j, v in total.items() if key[j][0] == "w" and v != 0}
beta = {key[j][1]: v for j, v in total.items() if key[j][0] == "c"}
print("rows", len(cert), "| meaning failures", len(bad), bad[:3])
print("w variables left over:", len(w_left))
print("combined: sum beta_z c_z <=", rhs_total, "; min beta_z =", min(beta.get(z, Fr(0)) for z in range(31)))
if not bad and not w_left and all(beta.get(z, Fr(0)) >= -1 for z in range(31)):
    print(f"VERIFIED: every commutative-slice unit at block 5 has M >= {-rhs_total} W = {float(-rhs_total):.6f} W")
