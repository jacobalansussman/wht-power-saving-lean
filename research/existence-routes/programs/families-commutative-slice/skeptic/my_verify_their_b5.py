# My own re-check of ex2-families' rc_cert_b5.json using only the variable meanings in rc_vars_b5.json and my own
# geometry (geo.py): every row must be a valid inequality for every slice unit (by Theorem 10 + rank subadditivity),
# and the lambda-combination must give sum_z c_z >= 155/36 exactly.
import json, sys
from fractions import Fraction as Fr
import geo
src = "../"
V = json.load(open(src + "rc_vars_b5.json")); C = json.load(open(src + "rc_cert_b5.json"))
m = 5
def dset(idxs): return frozenset(i + 1 for i in idxs)          # c index i <-> direction i + 1
def piece_cells(var):
    t = V[var]
    if t[0] == "w": H = dset(t[1]); return {(x, y) for x in H for y in H}, H, H
    if t[0] == "v": R, Cc = dset(t[1]), dset(t[2]); return {(x, y) for x in R for y in Cc}, R, Cc
    if t[0] == "y": x = t[1] + 1; return {(x, x)}, frozenset({x}), frozenset({x})
    if t[0] == "z": x, y = t[1] + 1, t[2] + 1; return {(x, y)}, frozenset({x}), frozenset({y})
    raise ValueError(t)
errs = 0; kinds = {}
tot = {}; const = Fr(0)
for row in C:
    lam = Fr(row["lambda"]); assert lam > 0
    co = {k: Fr(v).limit_denominator(4) for k, v in row["coeffs"].items()}
    assert all(abs(float(co[k]) - row["coeffs"][k]) < 1e-12 for k in co)
    rhs = Fr(row["rhs"]).limit_denominator(4)
    cs = {k: v for k, v in co.items() if V[k][0] == "c"}
    ps = {k: v for k, v in co.items() if V[k][0] != "c"}
    kind = row["kind"]; kinds[kind] = kinds.get(kind, 0) + 1
    ok = True
    if kind.startswith("bound"):
        # p - sum_G c <= -1 : G odd form support containing R u C ; or size bound p - sum_R c <= 0
        (pk, pv), = ps.items(); ok &= pv == 1
        cells, R, Cc = piece_cells(pk)
        G = dset(int(V[k][1]) for k in cs)
        ok &= all(v == -1 for v in cs.values())
        if rhs == -1:
            ok &= geo.is_form_support(G, m) and geo.trace(G, m) == 1 and (R | Cc) <= G
        elif rhs == 0:
            ok &= (G == R or G == Cc)
        else:
            ok = False
    elif kind == "R10":
        F = dset(int(V[k][1]) for k in cs)
        ok &= all(v == 1 for v in cs.values()) and rhs == 0
        ok &= geo.is_form_support(F, m) and geo.trace(F, m) == 0
        cover = {}
        for k, v in ps.items():
            cells, R, Cc = piece_cells(k)
            mult = -v
            if V[k][0] == "v":   # shared variable for (R, C) and (C, R): coefficient -2 = both orientations
                ok &= mult == 2
                allc = list(cells) + [(y, x) for (x, y) in cells]
            else:
                ok &= mult == 1
                allc = list(cells)
            for cell in allc: cover[cell] = cover.get(cell, 0) + 1
        ok &= set(cover) == {(x, y) for x in F for y in F} and all(v == 1 for v in cover.values())
    else:
        ok = False
    if not ok:
        errs += 1; print("BAD ROW", row["row"], kind)
    for k, v in co.items(): tot[k] = tot.get(k, Fr(0)) + lam * v
    const += lam * rhs
pc = [v for k, v in tot.items() if V[k][0] != "c"]
cc = [tot.get(k, Fr(0)) for k in V if V[k][0] == "c"]
print("rows", len(C), kinds, "invalid rows", errs)
print("piece coefficients in the sum: min", min(pc), "max", max(pc), "(need >= 0)")
print("c coefficients: min", min(cc), "(need >= -1)")
print("constant (sum lambda * rhs):", const, "=> sum c >= ", -const, "=", float(-const))
