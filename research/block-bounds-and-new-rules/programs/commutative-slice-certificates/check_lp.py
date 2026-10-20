# pf-logging, 14 Oct 2026.  INDEPENDENT CHECKER (pure Python, exact rationals; no code shared with gen_lp.py) of a
# mixed-strategy certificate for g-play's commutative-slice test at block size m, valid for every W.
# Classes are computed from the functions themselves (as in check_matroid.py):
#   lift(z)(y) = z.y mod 2 (as 0/1);  kappa(y) = popcount(y) mod 4;
#   Phi(q) = ( q(e_j) - q(0) mod 2 )_j  then  ( (q(0) - q(e_j) - q(e_l) + q(e_j+e_l)) mod 4 ) / 2  for j < l  (asserted even).
# Checks: mu is a probability distribution (exact), every y in its support has y.Phi(kappa) = 1, and for every direction
# z != 0:  E_mu[ y.Phi(lift z) ] <= v.
# Consequence (PROOF.md section 5b): for ANY multiset z_1..z_n with n v < W there is y in the support with
# #{i : y.Phi(z_i) = 1} <= W - 1, so the multiset fails the test.     Usage: python3 -I check_lp.py CERT.json
import sys, json
from fractions import Fraction
C = json.load(open(sys.argv[1])); m = int(C["m"]); v = Fraction(C["v"]); M = 1 << m
mu = {int(y): Fraction(q) for y, q in C["mu"].items()}
pairs = [(j, l) for j in range(m) for l in range(j + 1, m)]
def Phi(q):
    bits = [(q(1 << j) - q(0)) % 2 for j in range(m)]
    for (j, l) in pairs:
        d = (q(0) - q(1 << j) - q(1 << l) + q((1 << j) | (1 << l))) % 4
        assert d % 2 == 0
        bits.append(d // 2)
    return sum(b << i for i, b in enumerate(bits))
dot = lambda a, b: bin(a & b).count('1') & 1
PHk = Phi(lambda y: bin(y).count('1') % 4)
PHz = {z: Phi(lambda y, z=z: bin(z & y).count('1') % 2) for z in range(1, M)}
ok_dist = all(q >= 0 for q in mu.values()) and sum(mu.values()) == 1
ok_kappa = all(dot(y, PHk) == 1 for y in mu)
worst = max(sum((q for y, q in mu.items() if dot(y, PHz[z])), Fraction(0)) for z in range(1, M))
ok_v = worst <= v
print(f"m={m}: support {len(mu)}, distribution exact: {ok_dist}, all y.Phi(kappa)=1: {ok_kappa}, "
      f"max_z E[y.Phi(z)] = {worst} <= v = {v}: {ok_v}")
allok = ok_dist and ok_kappa and ok_v
print(f"RESULT: {'PASS' if allok else 'FAIL'} -> at block {m}, every multiset of n directions with n < {1/v} W fails the test")
sys.exit(0 if allok else 1)
