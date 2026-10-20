# Star-lemma recursion (e3-slice6, Theorem C/D). E nonsingular on AG(n,2) with rank <= s on every affine hyperplane.
# For a d-flat Psi, g = 2^(n-d) - 1 (d+1)-flats through Psi partition the rest; a nonsingular k x k block inside Psi
# forces (g-2) k <= g s_{d+1} - 2^(d+1).  So rank on d-flats <= s_d = floor((g s_{d+1} - 2^(d+1))/(g-2)), capped at 2^d.
# Contradiction when s_d < 0, or when s_d <= d - 1 (Lemma 11: every s_d + 2 <= d + 1 points lie in a d-flat).
import math
def excluded(n, s):
    sd = s
    for d in range(n - 2, -1, -1):
        g = 2 ** (n - d) - 1
        sd = min(2 ** d, (g * sd - 2 ** (d + 1)) // (g - 2))
        if sd < 0 or sd <= d - 1:
            return True, d
    return False, None
print(" m  n  |F|  max s excluded  => U_m excluded for W >=   (Cor. 15 had W > 0.9|F|)  fraction")
for n in range(2, 13):
    m = n + 2; F = 2 ** n
    smax = max([s for s in range(0, F) if excluded(n, s)[0]], default=-1)
    Wmin = F - smax
    bar = (2 ** m - 1) / m ** 0.9
    sav = lambda W: math.log(m * W / (2 ** m - 1)) / math.log(m)
    openab = [W for W in range(1, F + 1) if W >= bar and W < Wmin]
    print(f"{m:2d} {n:2d} {F:5d} {smax:6d} {Wmin:10d}  (old {math.floor(0.9*F)+1:5d})  {Wmin/F:.4f}  "
          f"bar needs W >= {math.ceil(bar)}; open above bar: {('%d..%d' % (openab[0], openab[-1])) if openab else 'none'}"
          + (f" (savings {sav(openab[0]):.3f}..{sav(openab[-1]):.3f})" if openab else "")
          + f"; saving at W = {Wmin}: {sav(Wmin):.3f}")
