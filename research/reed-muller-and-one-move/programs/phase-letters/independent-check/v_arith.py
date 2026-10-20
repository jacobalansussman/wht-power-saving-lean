from fractions import Fraction as F
import math
for m in range(2, 41):
    W = 2 ** (m - 2); M = 2 ** m - 1
    assert F(M, W) == 2 * (2 - F(2) ** (1 - m))
print("2 T_m(m) = (2^m-1)/2^(m-2) for m = 2..40: ok")
print("4 in 15: lambda <", F(4 * 4, 15) - 1)
for W in (8, 9, 10):
    print(f"{W} in 31: lambda <", F(5 * W, 31) - 1, float(F(5 * W, 31) - 1))
bar = 4 / math.exp(0.1 * math.log(4))
print("bar at block 4: M/W <", bar)
for W in range(1, 6):
    print(W, [M for M in range(1, 40) if F(41 * W, 13) <= M and M / W < bar])
for (W, M) in ((2, 7), (3, 10), (3, 11), (4, 13), (4, 14), (4, 15), (5, 16), (5, 17)):
    print(f"{W} in {M}: saving {math.log(4 * W / M) / math.log(4):.4f}")
for (W, M) in ((8, 31), (9, 31), (10, 31)):
    print(f"{W} in {M}: saving {math.log(5 * W / M) / math.log(5):.4f}, per array {M / W:.4f}")
print("2.561:", F(105, 41), float(F(105, 41)), "5*41/35 =", F(5 * 41, 35), "< 6")
