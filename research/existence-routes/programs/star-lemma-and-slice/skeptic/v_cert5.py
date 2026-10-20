#!/usr/bin/env python3
"""v_cert5.py (v-e3-slice6): my own check of e3-slice6's block-5 certificate (Theorem E: M >= 93/20 W).
Reads only the certificate rows (star_cert_b5.json) and the variable names (star_vars_b5.json); rebuilds all 2^15
binary quadratic forms on F_2^5 myself; checks that every row is a VALID inequality for the true values
  c_x = C_x / W (moves in direction x per array, >= 0),  w_P = rank E[P,P] / W (>= 0)
namely  bound:w   w_P <= sum_{x in G} c_x - 1   with G = supp of an ODD form, P a plane inside G  (Thm 10(c));
        STAR-L    2 (c_x + c_y) <= w_P1 + w_P2 + w_P3  with {x,y} a line, P_j the three planes through it, and
                  P1 u P2 u P3 = supp of an EVEN form (Thm 10(b) + star lemma with k = 0, a >= 0 dropped);
then adds lambda_i * row_i in exact fractions and checks the sum gives  -sum_x c_x <= -93/20  using only c, w >= 0.
"""
import json, sys
from fractions import Fraction as Fr
d = sys.argv[1]
cert = json.load(open(d + '/star_cert_b5.json')); var = json.load(open(d + '/star_vars_b5.json'))
m = 5; pts = list(range(1, 32))
mons = [(j, k) for j in range(m) for k in range(j, m)]
monmask = [sum(1 << z for z in pts if (z >> j) & 1 and (z >> k) & 1) for (j, k) in mons]
odd, even = set(), set()
mask = 0; tr = 0
for g in range(1 << len(mons)):
    if g:
        b = (g & -g).bit_length() - 1; mask ^= monmask[b]; tr ^= (mons[b][0] == mons[b][1])
    (odd if tr else even).add(mask)
print("forms:", 1 << len(mons), " distinct odd supports", len(odd), " even supports", len(even))
def dirs_of_c(i): return var[str(i)][1] + 1     # c index -> direction z (indices are z-1)
def plane_of_w(i):
    kind, L = var[str(i)]; assert kind == 'w'; return frozenset(x + 1 for x in L)
def is_affine_plane(S):
    S = sorted(S)
    if len(S) != 4: return False
    a = S[0]; V = {x ^ a for x in S}
    return 0 in V and all((u ^ v) in V for u in V for v in V) and 0 not in S
def is_line(S): return len(S) == 2 and 0 not in S
nbad = 0; y = {}; total_rhs = Fr(0)
for r in cert:
    lam = Fr(r['lambda']); assert lam >= 0
    co = {int(k): Fr(v).limit_denominator(1) for k, v in r['coeffs'].items()}
    rhs = Fr(r['rhs']).limit_denominator(1)
    cvars = {i: a for i, a in co.items() if var[str(i)][0] == 'c'}
    wvars = {i: a for i, a in co.items() if var[str(i)][0] == 'w'}
    others = [i for i in co if var[str(i)][0] not in ('c', 'w')]
    ok = not others
    if r['kind'] == 'bound:w':
        ok = ok and len(wvars) == 1 and list(wvars.values()) == [1] and rhs == -1 and all(a == -1 for a in cvars.values())
        G = sum(1 << dirs_of_c(i) for i in cvars); P = plane_of_w(list(wvars)[0])
        ok = ok and G in odd and is_affine_plane(P) and all((G >> z) & 1 for z in P)
    elif r['kind'] == 'STAR-L':
        ok = ok and rhs == 0 and len(cvars) == 2 and all(a == 2 for a in cvars.values()) and len(wvars) == 3 and all(a == -1 for a in wvars.values())
        L = {dirs_of_c(i) for i in cvars}; Ps = [plane_of_w(i) for i in wvars]
        ok = ok and is_line(L) and all(is_affine_plane(P) and L <= P for P in Ps) and len(set(Ps)) == 3
        U = set().union(*Ps); Fm = sum(1 << z for z in U)
        ok = ok and len(U) == 8 and Fm in even
        # the three planes through L inside the 3-flat U are exactly these (U is a 3-flat: even support of weight 8)
    else:
        ok = False
    if not ok: nbad += 1; print("BAD ROW", r)
    for i, a in co.items(): y[i] = y.get(i, 0) + lam * a
    total_rhs += lam * rhs
cneg = {i: v for i, v in y.items() if var[str(i)][0] == 'c'}
wneg = {i: v for i, v in y.items() if var[str(i)][0] == 'w' and v < 0}
print("rows", len(cert), "invalid rows", nbad)
print("combined c coefficients: min", min(cneg.values()), "max", max(cneg.values()), "(need >= -1; c vars present", len(cneg), "of 31)")
print("combined w coefficients < 0:", len(wneg), "; max w coefficient", max((v for i, v in y.items() if var[str(i)][0] == 'w'), default=0))
print("sum lambda*rhs =", total_rhs)
if nbad == 0 and not wneg and all(v >= -1 for v in cneg.values()) and len(cneg) == 31:
    print("VERIFIED: sum_x c_x >=", -total_rhs, "=", float(-total_rhs))
else: print("NOT VERIFIED")
