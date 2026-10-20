# sym_cert_b4.py: average the block-4 certificate over the 24 coordinate permutations (the LP is invariant under them),
# re-verify the averaged certificate exactly (rows keyed by supports, forms rebuilt independently), and print its
# structure: which even planes, how many rows, how many distinct weights.
import json, itertools
from fractions import Fraction as Fr
import numpy as np
from refl import quad_forms
m = 4
vals, tr = quad_forms(m); V = vals[:, 1:]
supp = lambda f: frozenset(int(z) + 1 for z in np.nonzero(V[f])[0])        # directions 1..15
odd_supps = {supp(q) for q in range(len(V)) if tr[q] == 1}
even_supps = {supp(q) for q in range(len(V)) if tr[q] == 0}
cert = json.load(open("rc_cert_b4.json"))
def perm_dir(p, z): return sum(((z >> j) & 1) << p[j] for j in range(m))
acc = {}
for p in itertools.permutations(range(m)):
    P = lambda X: frozenset(perm_dir(p, z) for z in X)
    for it in cert:
        r = it["row"]; lam = Fr(it["lambda"]) / 24
        F = P(supp(r[1]))
        if r[0] == "Sslack": key = ("S", F, perm_dir(p, r[2] + 1), P({z + 1 for z in r[3]}))
        elif r[0] == "Cslack":
            a, b = sorted((perm_dir(p, r[2] + 1), perm_dir(p, r[3] + 1))); key = ("C", F, a, b, P({z + 1 for z in r[4]}))
        else: key = ("R", F)
        acc[key] = acc.get(key, Fr(0)) + lam
# exact verification
lin = {}; const = Fr(0); err = 0
for key, lam in acc.items():
    if key[0] == "S":
        _, F, x, S = key
        if not (S in odd_supps and x in S and x in F and F in even_supps): err += 1
        lin[("y", F, x)] = lin.get(("y", F, x), 0) + lam
        for z in S: lin[("c", z)] = lin.get(("c", z), 0) - lam
        const -= lam
    elif key[0] == "C":
        _, F, a, b, S = key
        if not (S in odd_supps and a in S and b in S and a in F and b in F and F in even_supps): err += 1
        lin[("z", F, a, b)] = lin.get(("z", F, a, b), 0) + lam
        for z in S: lin[("c", z)] = lin.get(("c", z), 0) - lam
        const -= lam
    else:
        _, F = key
        pts = sorted(F)
        if not (F in even_supps and len(pts) == 4 and pts[0] ^ pts[1] ^ pts[2] ^ pts[3] == 0): err += 1
        for x in pts: lin[("c", x)] = lin.get(("c", x), 0) + lam; lin[("y", F, x)] = lin.get(("y", F, x), 0) - lam
        for a, b in itertools.combinations(pts, 2): lin[("z", F, a, b)] = lin.get(("z", F, a, b), 0) - 2 * lam
left = [k for k, v in lin.items() if k[0] != "c" and v != 0]
beta = {k[1]: v for k, v in lin.items() if k[0] == "c"}
print(f"averaged certificate: {len(acc)} rows; label errors {err}; uncancelled {len(left)}; min beta {min(beta.values())}; "
      f"constant {const}")
print("VERIFIED (symmetric certificate): sum c >= 4" if (err == 0 and not left and min(beta.values()) >= -1 and const == -4)
      else "NOT verified")
Rs = sorted({k[1] for k in acc if k[0] == "R"}, key=sorted)
print("even planes used:", [sorted(F) for F in Rs])
print("distinct weights on R rows:", sorted({str(acc[("R", F)]) for F in Rs}))
print("beta by direction (coefficient of c_z is -beta):", {z: str(-v) for z, v in sorted(beta.items())})
json.dump([{"key": [k[0], sorted(k[1])] + [x if isinstance(x, int) else sorted(x) for x in k[2:]], "lambda": str(v)}
           for k, v in acc.items()], open("rc_cert_b4_sym.json", "w"))
