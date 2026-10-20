# verify_cert_b5.py: independent exact check of rc_cert_b5.json (Theorem 18: block-5 commutative slice, M >= 155W/36).
# Rebuilds the binary quadratic forms on F_2^5 itself (no refl.py), checks the meaning of every certificate row from its
# coefficients, and sums lambda * row exactly.
#   bound rows:  piece - c(T) <= 0 (T = a size set inside the piece)  or  piece - c(S) <= -1 (S = an ODD support
#                containing the piece's directions);
#   R10 rows:    c(F) - w_H - w_H' - 2 (v_KK' + v_KL' + v_LK' + v_LL') <= 0 with F an EVEN 3-flat, H, H' complementary
#                planes of F, K, L (K', L') parallel lines of H (H') with K u K', K u L', L u K', L u L' planes.
import json, itertools
from fractions import Fraction as Fr
m = 5
pairs = [(j, k) for j in range(m) for k in range(j, m)]
odd_supps, even_supps = set(), set()
for bits in range(1 << len(pairs)):
    q = [(bits >> i) & 1 for i in range(len(pairs))]
    tr = sum(q[i] for i, (j, k) in enumerate(pairs) if j == k) % 2
    S = frozenset(z for z in range(1, 32) if sum(q[i] * ((z >> j) & 1) * ((z >> k) & 1) for i, (j, k) in enumerate(pairs)) % 2)
    (odd_supps if tr else even_supps).add(S)
names = json.load(open("rc_vars_b5.json"))           # index -> repr of key
key = {int(i): k for i, k in names.items()}
def dirs_of(k):                                       # directions (1..31) covered by a piece variable
    t = k[0]
    if t == "y": return {k[1] + 1}
    if t == "z": return {k[1] + 1, k[2] + 1}
    if t == "w": return {z + 1 for z in k[1]}
    if t == "v": return {z + 1 for z in k[1]} | {z + 1 for z in k[2]}
def isflat(P, dim):
    P = sorted(P); a = P[0]; D = {p ^ a for p in P}
    return len(P) == 2 ** dim and all((u ^ v) in D for u in D for v in D)
cert = json.load(open("rc_cert_b5.json"))
lin = {}; const = Fr(0); errors = 0
for item in cert:
    lam = Fr(item["lambda"]); co = {int(j): Fr(v).limit_denominator(1) for j, v in item["coeffs"].items()}
    rhs = Fr(item["rhs"]).limit_denominator(1)
    if lam < 0: errors += 1
    cs = {key[j][1] + 1: v for j, v in co.items() if key[j][0] == "c"}
    ps = {j: v for j, v in co.items() if key[j][0] != "c"}
    if item["kind"].startswith("bound"):
        ok = len(ps) == 1 and list(ps.values())[0] == 1 and all(v == -1 for v in cs.values())
        (j,) = ps.keys(); U = dirs_of(key[j]); T = set(cs)
        if rhs == 0: ok = ok and T <= U and len(T) >= 1                       # size bound
        elif rhs == -1: ok = ok and frozenset(T) in odd_supps and U <= T      # slack bound
        else: ok = False
        if key[j][0] == "w": ok = ok and isflat(U, 2)
        if key[j][0] == "v":
            K = {z + 1 for z in key[j][1]}; Kp = {z + 1 for z in key[j][2]}
            ok = ok and isflat(K, 1) and isflat(Kp, 1) and isflat(K | Kp, 2)
        if not ok: errors += 1; print("bad bound row", item["row"])
    elif item["kind"] == "R10":
        F = {z for z, v in cs.items() if v == 1}
        ok = all(v == 1 for v in cs.values()) and frozenset(F) in even_supps and isflat(F, 3)
        ws = [j for j, v in ps.items() if key[j][0] == "w"]; vs = [j for j, v in ps.items() if key[j][0] == "v"]
        ok = ok and len(ws) == 2 and all(ps[j] == -1 for j in ws) and len(vs) == 4 and all(ps[j] == -2 for j in vs)
        H, Hp = dirs_of(key[ws[0]]), dirs_of(key[ws[1]])
        ok = ok and H | Hp == F and not (H & Hp)
        cells = set()
        for j in vs:
            X = {z + 1 for z in key[j][1]}; Y = {z + 1 for z in key[j][2]}
            ok = ok and ((X <= H and Y <= Hp) or (X <= Hp and Y <= H))
            cells |= {(a, b) for a in X for b in Y} | {(b, a) for a in X for b in Y}
        # pieces must tile F x F exactly: H x H, H' x H', and the cross cells (both orders) once each
        allcells = {(a, b) for a in H for b in H} | {(a, b) for a in Hp for b in Hp} | cells
        ok = ok and len(cells) == 2 * len(H) * len(Hp) and allcells == {(a, b) for a in F for b in F}
        if not ok: errors += 1; print("bad R10 row", item["row"])
    else:
        errors += 1; print("unexpected kind", item["kind"])
    for j, v in co.items(): lin[j] = lin.get(j, Fr(0)) + lam * v
    const += lam * rhs
left = {j: v for j, v in lin.items() if key[j][0] != "c" and v != 0}
beta = {j: v for j, v in lin.items() if key[j][0] == "c"}
ok_c = all(v >= -1 for v in beta.values())
print(f"rows {len(cert)}; row errors {errors}; piece variables left uncancelled {len(left)}; c-coefficients >= -1: {ok_c}; "
      f"sum lambda * rhs = {const}")
if errors == 0 and not left and ok_c:
    print(f"VERIFIED: sum_z c_z >= {-const} = {float(-const):.6f}, i.e. M >= {-const} W for every block-5 slice unit")
