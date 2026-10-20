# b56-note: arithmetic for NOTE.md (exact fractions; nothing heavy). Recomputes every number the table quotes.
from fractions import Fraction as F
from math import log, comb, ceil
def sav(m, L): return log(m / float(L)) / log(m)
def unit(m, W, M): return log(m * W / M) / log(m)
rows = {
 2: [("both", F(2))],
 3: [("both", F(14, 5))],
 4: [("both", F(86592, 27456)), ("engine-lift", F(137700, 50595))],
 5: [("both", F(970920, 438480)), ("engine-lift", F(11738925, 4404165)), ("engine-proj", F(9486, 3575))],
 6: [("both", F(122976, 60512)), ("engine-proj", F(405993, 153725))],
}
for m, lst in rows.items():
    for g, L in lst:
        print(f"m={m} {g}: L={L} = {float(L):.4f}  cap={sav(m, L):.4f}")
        ex, op = [], []
        for W in range(1, 7):
            Mmin = ceil(L * W)          # M >= L W required
            op.append(f"{W} in {Mmin} ({unit(m, W, Mmin):.4f})" if Mmin < m * W else f"{W}: none below plain")
        print("   smallest NOT excluded per W:", "; ".join(op))
print()
for m in range(2, 11):
    n = 2 ** m
    hook = F(4 * (n - 1), n + 2); coord = F(4 * (n - 1), n)
    kmax = max(k for k in range(1, n // 2 + 1) if (k - 1) ** 2 < n)
    arg = F(2 * (n - 1), n - kmax + 1)
    # check cap identity and psi bound at kmax
    cap = F(n, 2) * comb(n - 2, kmax - 1); dem = kmax * comb(n, kmax) - comb(n, kmax - 1)
    assert F(dem) / cap == arg and cap == F(kmax * (n - kmax) * comb(n, kmax), 2 * (n - 1))
    s_h = sav(m, hook) if hook < m else 0.0; s_c = sav(m, coord) if coord < m else 0.0
    print(f"m={m}: hook ceiling {hook} = {float(hook):.4f} (s* {s_h:.4f}); coord ceiling {coord} = {float(coord):.4f} "
          f"(s* {s_c:.4f}); argued hook k={kmax}: {arg} = {float(arg):.4f} (cap {sav(m, arg):.4f}); bar m^0.9 {m**0.9:.3f}")
print()
for m, W, M in [(3,5,14),(3,3,8),(4,2,7),(4,3,10),(4,3,11),(4,4,14),(5,3,11),(5,2,8),(5,3,12),(5,4,16),(6,2,8),(6,3,12),(6,4,16),(5,3,8),(6,3,8),(5,2,5),(6,2,5),(5,3,7),(6,3,7)]:
    print(f"m={m}: {W} in {M}: value {unit(m, W, M):.4f}  M/W={M/W:.4f}")
