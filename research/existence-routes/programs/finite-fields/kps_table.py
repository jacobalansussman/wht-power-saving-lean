"""kps_table.py -- ex2-modp: how many primes the 'solutions mod many primes' route needs (PROOF.md Prop. 5).
Krick-Pardo-Sombra 2001, Theorem 1 (read on arXiv math/9911094, p. 2): f_1..f_s in Z[x_1..x_n] with no common complex
zero, d = max deg, h = max log|coef|  =>  a = sum g_i f_i with a in Z \ {0}, h(a) <= 4n(n+1)d^n (h + log s + (n+7) log(n+1) d).
Every prime p with X(F_p-bar) nonempty divides a.  So F_p-bar points at more than h(a)/log(P) distinct primes >= P prove
X(C) nonempty.  (Masser-Philippon: h(a) >= d^(n-1) h in the worst case, so the d^n growth is real.)"""
import math
def need(n, d, h, s, P=2**30):
    H = 4*n*(n+1)*d**n*(h + math.log(s) + (n+7)*math.log(n+1)*d)
    return H, H/math.log(P)
print(f"{'case':42s} {'n':>3s} {'d':>2s} {'s':>4s} {'h(a) bound':>11s} {'primes >= 2^30':>15s} {'primes >= 2^62':>15s}")
for name, n, d, s in [('toy: one alpha, one beta (3-row core)', 2, 4, 3),
                      ('plain chart, joint 2 arrays (cal_plain4F)', 7, 4, 4),
                      ('border chart 1_0_1_0, one array', 5, 4, 8),
                      ('border chart 0_0_1_1 (61 x 13 + 2 alphas)', 15, 4, 62),
                      ('relaxation hit (12 alphas, 9 betas)', 21, 5, 40),
                      ('core 133 x 37 + 3 alphas', 40, 4, 134)]:
    H, k = need(n, d, math.log(2), s)
    H2, k2 = need(n, d, math.log(2), s, 2**62)
    print(f"{name:42s} {n:3d} {d:2d} {s:4d} {H:11.3g} {k:15.3g} {k2:15.3g}")
