# verify_cert_b4.py: independent exact check of rc_cert_b4.json (Theorem 17: block-4 commutative slice, M >= 4W).
# Does NOT use refl.py or the LP code: rebuilds binary quadratic forms on F_2^4 itself, re-derives every inequality from
# its label, checks each label's claims (odd/even, supports, flats, membership), and sums lambda * row in exact rationals.
# Rows (W = 1; c_z = moves in direction z per array; y_{F,x}, z_{F,a,b} rank bounds):
#   Sslack(F, x, S): y_{F,x} - c(S) <= -1      requires: S = support of an ODD form, x in S, x in supp F
#   Cslack(F, a, b, S): z_{F,a,b} - c(S) <= -1 requires: S = support of an ODD form, a, b in S and in supp F
#   R(F): c(supp F) - sum_x y_{F,x} - 2 sum_{a<b} z_{F,a,b} <= 0   requires: F EVEN, supp F a 2-flat (4 points, XOR 0)
import json, itertools
from fractions import Fraction as Fr
m = 4
pairs = [(j, k) for j in range(m) for k in range(j, m)]
forms = {}
for bits in range(1 << len(pairs)):
    q = [(bits >> i) & 1 for i in range(len(pairs))]
    trace = sum(q[i] for i, (j, k) in enumerate(pairs) if j == k) % 2
    supp = frozenset(z - 1 for z in range(1, 16)
                     if sum(q[i] * ((z >> j) & 1) * ((z >> k) & 1) for i, (j, k) in enumerate(pairs)) % 2)
    forms[bits] = (trace, supp)
odd_supps = {s for (t, s) in forms.values() if t == 1}
# the LP code indexes forms in itertools.product order over the pair coefficients (first pair = most significant bit)
prod_index = list(itertools.product([0, 1], repeat=len(pairs)))
def form_of(f):
    q = prod_index[f]
    bits = sum(q[i] << i for i in range(len(pairs)))
    return forms[bits]
cert = json.load(open("rc_cert_b4.json"))
lin = {}; const = Fr(0); errors = 0
def addterm(k, v):
    lin[k] = lin.get(k, Fr(0)) + v
for item in cert:
    row, lam = item["row"], Fr(item["lambda"])
    if lam < 0: errors += 1; print("negative lambda", row)
    kind = row[0]
    if kind == "Sslack":
        _, f, x, S = row; S = frozenset(S); tF, sF = form_of(f)
        ok = (S in odd_supps) and (x in S) and (x in sF) and tF == 0
        if not ok: errors += 1; print("bad row", row)
        addterm(("y", f, x), lam)
        for z in S: addterm(("c", z), -lam)
        const += -lam
    elif kind == "Cslack":
        _, f, a, b, S = row; S = frozenset(S); tF, sF = form_of(f)
        ok = (S in odd_supps) and a in S and b in S and a in sF and b in sF and tF == 0
        if not ok: errors += 1; print("bad row", row)
        addterm(("z", f, min(a, b), max(a, b)), lam)
        for z in S: addterm(("c", z), -lam)
        const += -lam
    elif kind == "R":
        _, f = row; tF, sF = form_of(f)
        pts = [z + 1 for z in sF]
        flat = len(pts) == 4 and (pts[0] ^ pts[1] ^ pts[2] ^ pts[3]) == 0
        if not (tF == 0 and flat): errors += 1; print("bad R row", row)
        for x in sF: addterm(("c", x), lam); addterm(("y", f, x), -lam)
        for a, b in itertools.combinations(sorted(sF), 2): addterm(("z", f, a, b), -2 * lam)
    else:
        errors += 1; print("unknown row", row)
leftover = {k: v for k, v in lin.items() if k[0] != "c" and v != 0}
beta = {k[1]: v for k, v in lin.items() if k[0] == "c"}
ok_c = all(v >= -1 for v in beta.values())
print(f"rows {len(cert)}; label errors {errors}; rank-bound variables left uncancelled: {len(leftover)}; "
      f"c-coefficients all >= -1: {ok_c}; constant: sum lambda * rhs = {const}")
if errors == 0 and not leftover and ok_c and const == -4:
    print("VERIFIED: sum_z beta_z c_z <= -4 with beta_z >= -1, hence sum_z c_z >= 4, i.e. M >= 4W for every unit")
