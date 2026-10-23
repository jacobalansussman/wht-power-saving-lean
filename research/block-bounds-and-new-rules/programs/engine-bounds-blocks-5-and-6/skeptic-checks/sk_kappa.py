# v-pf-induction-compute: kappa_m from the standard dual-polar eigenvalues theta_j = 2(2^(m-j)-1) - (2^j-1)
# (BCN 9.4, q = 2, e = 1), cosines by the three-term recursion; checks theta_j are eigenvalues (end condition),
# and that max_j (1-u_m)/(1-u_1) = 2 - 2^(1-m) at j = 1, for m = 2..40 (extends the proof's m <= 12).
from fractions import Fraction as Fr
for m in range(2, 41):
    b = [2**(i+1) * (2**(m-i) - 1) for i in range(m + 1)]; c = [2**i - 1 for i in range(m + 1)]
    k = b[0]; a = [k - b[i] - c[i] for i in range(m + 1)]
    best = None
    for j in range(1, m + 1):
        th = 2 * (2**(m-j) - 1) - (2**j - 1)
        u = [Fr(1), Fr(th, k)]
        for i in range(1, m): u.append((th * u[i] - c[i] * u[i-1] - a[i] * u[i]) / b[i])
        assert c[m] * u[m-1] + a[m] * u[m] == th * u[m]
        v = (1 - u[m]) / (1 - u[1])
        if best is None or v > best[0]: best = (v, j)
    assert best[0] == 2 - Fr(2, 2**m) and best[1] == 1, (m, best)
print("kappa_m = 2 - 2^(1-m), attained at theta_1 = 2^m - 3, for every m = 2..40: True")
