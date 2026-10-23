# v-pf-induction-compute: independent recomputation of PROOF.md Theorems 3-5 numbers.
# Uses N_r(s) closed form (verified above by brute force for m <= 7, and proved in Lemma 2).
import sys
from fractions import Fraction as Fr
def closed(m, r):
    d = {}
    for s, c in ((r, (2**(m-r) - 1) + (4**m - 2**(m+r))),
                 (r-1, 2*(2**m - 2**(m-r)) + (2**r - 1) * 2**(r-1) * 2**(m-r)),
                 (r-2, (2**r - 1) * (2**(r-1) - 1) * 2**(m-r))):
        if c: d[s] = d.get(s, 0) + c
    assert sum(d.values()) == 4**m - 1
    return d
base = [0, 27456, 48576, 68320, 86592]
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
a = list(base); tops = {4: Fr(base[4], base[1])}; profs = {}
for m in range(5, MMAX + 1):
    b = [0] * (m + 1)
    for r in range(1, m + 1):
        b[r] = sum(c * a[s] for s, c in closed(m, r).items() if 1 <= s <= m - 1)
    a = b; tops[m] = Fr(a[m], a[1])
    if m in (5, 6, 7, 30, 60, MMAX): profs[m] = [Fr(a[r], a[1]) for r in range(m + 1)]
print("T_5 =", tops[5], float(tops[5]), " claimed 834241/278850:", tops[5] == Fr(834241, 278850))
print("T_6 =", tops[6], float(tops[6]), " claimed 1188222/399685:", tops[6] == Fr(1188222, 399685))
print("T_7 =", tops[7], float(tops[7]), " claimed 2361819/799370:", tops[7] == Fr(2361819, 799370))
mn = min(tops[m] for m in range(5, MMAX + 1)); arg = min(range(5, MMAX + 1), key=lambda m: tops[m])
print(f"exact tops m=5..{MMAX}: min {float(mn):.15f} at m={arg}; >= 2.9447: {mn >= Fr(29447, 10000)}")
print("monotone decreasing in m:", all(tops[m+1] <= tops[m] for m in range(5, MMAX)))
print("tops at m = 10, 20, 40, 60, 100, MMAX:", [f"{float(tops[m]):.12f}" for m in (10, 20, 40, 60, 100, MMAX) if m <= MMAX])
for m, R in profs.items():
    print(f"block {m}: profile nondecreasing: {all(R[r] <= R[r+1] for r in range(1, m))}; R_r <= r: {all(R[r] <= r for r in range(1, m+1))}; "
          f"subadditive: {all(R[r+s] <= R[r] + R[s] for r in range(1, m) for s in range(1, m - r + 1))}" if m <= 60 else
          f"block {m}: profile nondecreasing: {all(R[r] <= R[r+1] for r in range(1, m))}")
# domination inequalities (my own derivation of probabilities from the deficit chain)
viol = 0
for k in range(1, 201):
    P = 4**k - 1
    for j in range(0, k):
        r = k - j; N = closed(k, r)
        up = N.get(r - 2, 0); down = N.get(r, 0)
        if Fr(up, P) > Fr(1, 2**(j + 1)): viol += 1
        if j >= 1 and Fr(down, P) < 1 - Fr(1, 2**j): viol += 1
print("domination violations k<=200:", viol)
# corrupted dominating chain (up 2^-(j+2)) must show violations
cv = sum(1 for k in range(2, 41) for j in range(0, k) if Fr(closed(k, k - j).get(k - j - 2, 0), 4**k - 1) > Fr(1, 2**(j + 2)))
print("corrupted dominating chain (up = 2^-(j+2)): violations", cv)
# stationary law and c*_M for M = 30, 60 (normaliser bounded above)
J = 100; w = [Fr(1)]
for j in range(1, J + 1): w.append(w[-1] / (2**j - 1))
Zup = sum(w) + 2 * w[J]
for M in (30, 60):
    R = profs[M]
    c = sum(w[j] * R[M - j] for j in range(0, M)) / Zup
    print(f"c*_{M} = {float(c):.15f} (>= 2.9447: {c >= Fr(29447, 10000)}); T_{M} = {float(tops[M]):.15f}")
# S_inf closed forms
for m0 in (2, 3, 4, 5):
    S = Fr(1)
    for k in range(m0 + 1, 300): S *= Fr(4**k - 3 * 2**(k - 1) - 1, 4**k - 1)
    print(f"S_inf(m0={m0}): product to 300 = {float(S):.15f}, closed form {float(Fr(2*(2**m0-1), 2**(m0+1)+1)):.15f}, "
          f"ceiling m0/S_inf = {float(m0 / Fr(2*(2**m0-1), 2**(m0+1)+1)):.4f}")
# excluded units at blocks 5, 6
for m in (5, 6):
    T = tops[m]
    print(f"block {m}: 3 in 8 excl {8 < 3*T}; 4 in 11 excl {11 < 4*T}; 5 in 14 excl {14 < 5*T}; "
          f"max W with 3W-1 excluded = {max(W for W in range(1, 2000) if 3*W - 1 < W*T)}")
