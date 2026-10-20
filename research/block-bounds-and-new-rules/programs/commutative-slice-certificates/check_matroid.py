# pf-logging, 13 Oct 2026.  INDEPENDENT CHECKER of the commutative-slice exclusion certificates (no shared code with
# gen_matroid.py).  Classes are computed from the functions themselves:
#   lift(z)(y) = (z.y mod 2) as the integer 0 or 1;   kappa(y) = popcount(y) mod 4   (y in F_2^m)
#   Phi(q) = ( q(e_j) - q(0) mod 2 )_{j}  followed by  ( (q(0) - q(e_j) - q(e_l) + q(e_j+e_l)) mod 4 ) / 2  for j < l
# (asserting that the second differences are even).  For each multiset (canonical order: all nondecreasing
# sequences of n directions from 1..2^m-1, enumerated here by our own recursion) the certificate gives y; we check
#   y.Phi(kappa) = 1   and   #{ positions i : y.Phi(z_i) = 1 } <= W - 1.
# Then deleting those positions leaves directions whose classes all lie in the hyperplane y = 0, which does not
# contain Phi(kappa): kappa is not in their span, and the multiset fails the test of g-play/matroid.py.
# Multisets marked 0xFFFF (no certificate) are checked to PASS by brute force (all deletions, all subset sums).
# Usage: python3 -I check_matroid.py m W n CERTFILE
import sys
from array import array
m, W, n = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
M = 1 << m
pairs = [(j, l) for j in range(m) for l in range(j + 1, m)]
def Phi(q):
    v = []
    for j in range(m): v.append((q(1 << j) - q(0)) % 2)
    for (j, l) in pairs:
        d = (q(0) - q(1 << j) - q(1 << l) + q((1 << j) | (1 << l))) % 4
        assert d % 2 == 0
        v.append(d // 2)
    return sum(b << i for i, b in enumerate(v))
lift = lambda z: (lambda y: bin(z & y).count('1') % 2)
PHz = {z: Phi(lift(z)) for z in range(1, M)}
PHk = Phi(lambda y: bin(y).count('1') % 4)
dot = lambda a, b: bin(a & b).count('1') & 1
cert = array('H'); cert.frombytes(open(sys.argv[4], 'rb').read())
def multisets(n, lo):
    if n == 0: yield (); return
    for z in range(lo, M):
        for rest in multisets(n - 1, z): yield (z,) + rest
def passes_bruteforce(ms):
    from itertools import combinations
    for J in combinations(range(len(ms)), W - 1):
        rest = [PHz[ms[i]] for i in range(len(ms)) if i not in J]
        span = {0}
        for v in rest: span |= {s ^ v for s in span}
        if PHk not in span: return False
    return True
cnt = ok = unc = unc_ok = 0
for idx, ms in enumerate(multisets(n, 1)):
    cnt += 1
    if idx >= len(cert): break
    y = cert[idx]
    if y == 0xFFFF:
        unc += 1; unc_ok += passes_bruteforce(ms); continue
    if dot(y, PHk) == 1 and sum(dot(y, PHz[z]) for z in ms) <= W - 1: ok += 1
complete = (cnt == len(cert))
print(f"m={m} W={W} n={n}: multisets enumerated {cnt}, certificate entries {len(cert)}, complete: {complete}")
print(f"  exclusion certificates valid: {ok}; entries without certificate: {unc} (of which pass the test by brute force: {unc_ok})")
allok = complete and ok + unc == cnt and unc_ok == unc
print(f"RESULT: {'PASS' if allok else 'FAIL'} -> {ok} of {cnt} multisets EXCLUDED by certificate; {unc} pass the test")
sys.exit(0 if allok else 1)
