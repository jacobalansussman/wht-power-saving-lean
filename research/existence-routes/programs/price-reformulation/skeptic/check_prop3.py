#!/usr/bin/env python3
"""v-ex2-reform: independent recheck of ex2-reform Proposition 3 arithmetic, Part A, and a test of the
consistency of hypothesis H_theta (OpenAI's Theorem 4.1 laws with the tensor law replaced by homogeneity).
Plain python + numpy, exact where it matters (fractions / integer matrices)."""
import math, itertools
from fractions import Fraction

print("== 1. Proposition 3 arithmetic (own code, scan over integer d) ==")
def lower(d, th): return (2*d)**th/2 - 4
def upper(d, th, beta): return beta*(d/2 + 1) + 2
def beta_open(th): return 3*2**(-th)
def best_gap(th, beta=None, dmax=200000):
    if beta is None: beta = beta_open(th)
    bd, bg = None, -1e300
    for d in range(1, dmax+1):
        g = lower(d, th) - upper(d, th, beta)
        if g > bg: bd, bg = d, g
    return bd, bg
# theta = 1: first d with contradiction
d1 = next(d for d in range(1, 1000) if lower(d,1) > upper(d,1,beta_open(1)))
print(f"theta=1: first contradiction at d={d1}; paper: d>30")
# theta_c by bisection on max_d gap (integer d up to 2e5 is enough below 0.97)
lo, hi = 0.90, 0.97
for _ in range(50):
    mid = (lo+hi)/2
    if best_gap(mid, dmax=20000)[1] > 0: hi = mid
    else: lo = mid
print(f"theta_c (integer d <= 20000) = {hi:.5f}; best d there = {best_gap(hi, dmax=20000)[0]}")
for th in [0.9, 0.95, 0.96]:
    d, g = best_gap(th, dmax=20000)
    print(f"theta={th}: OpenAI beta={beta_open(th):.4f}; max_d(lower-upper) = {g:+.4f} at d = {d}; value at d=120: {lower(120,th)-upper(120,th,beta_open(th)):+.3f}")
# needed beta at 0.9
lo, hi = 0.5, 2.0
for _ in range(60):
    mid = (lo+hi)/2
    if best_gap(0.9, beta=mid, dmax=20000)[1] > 0: lo = mid
    else: hi = mid
bn = lo; dn = best_gap(0.9, beta=bn*0.9999, dmax=20000)[0]
print(f"theta=0.9: needed dense-pair price < {bn:.5f} (OpenAI-under-homogeneity {beta_open(0.9):.5f}); optimum d = {dn}")
print(f"   candidate 'gap factors' at d={dn}: 2/beta_needed = {2/bn:.4f}; d/(needed cheap side) = {dn/(bn*(dn/2+1)+2):.4f};"
      f" (2d)^0.1 = {(2*dn)**0.1:.4f}; beta_open/beta_needed = {beta_open(0.9)/bn:.4f}; OpenAI at theta=1: 4/3 = {4/3:.4f}")

print()
print("== 2. Is H_theta consistent?  Test matrix K = U (+) I_2 (width 4) ==")
def kron(A, B):
    return [[A[i//len(B)][j//len(B[0])]*B[i%len(B)][j%len(B[0])] for j in range(len(A[0])*len(B[0]))] for i in range(len(A)*len(B))]
def eye(n): return [[int(i==j) for j in range(n)] for i in range(n)]
U = [[1,0],[1,1]]; I2 = eye(2)
K = [[1,0,0,0],[1,1,0,0],[0,0,1,0],[0,0,0,1]]
K2 = kron(K, K)
order = sorted(range(16), key=lambda t: ((t//4)>=2, (t%4)>=2, t))
B = [[K2[order[i]][order[j]] for j in range(16)] for i in range(16)]
blocks = [kron(U,U), kron(U,I2), kron(I2,U), eye(4)]
BD = [[0]*16 for _ in range(16)]; o = 0
for blk in blocks:
    for i in range(4):
        for j in range(4): BD[o+i][o+j] = blk[i][j]
    o += 4
print("K(x)K is permutation-similar to (U(x)U) (+) (U(x)I2) (+) (I2(x)U) (+) I4 :", B == BD)
for th in [1.0, 0.95, 0.9]:
    pU = 1.0
    pUU = 4 * 2**th * pU / 2          # homogeneity: p(X^(x)2)/|X^(x)2| = 2^th p(X)/|X|, |U|=2
    pK = pU + 0.0                      # exact direct-sum law, p(I2)=0
    via_sums = pUU + 2*pU + 2*pU + 0.0 # sum law on the block form; p(U(x)I2)=p(U(+)U)=2 p(U)
    via_hom = 16 * 2**th * pK / 4      # homogeneity applied to K
    print(f"theta={th}: p(K(x)K) by sum law = {via_sums:.5f}; by homogeneity = {via_hom:.5f}; difference = {via_sums-via_hom:+.5f}")
print("=> for theta<1 the two laws of H_theta give different prices for one matrix: H_theta is inconsistent.")
print("   (in general: p(A(+)I_q)^(x)2 forces 2^(1+theta) = 2^theta + 2, i.e. theta = 1)")

print()
print("== 3. What a theta-price obeys instead (OpenAI's own amplification sector formula) ==")
# Natural per-coordinate theta-price of A (+) I_{Q-q}: E[(Bin(j,p))^theta] ~ (p j)^theta  => factor p^theta, not p
for th in [0.9]:
    for (q,Q) in [(2,4),(2,16),(2,2**20)]:
        p = q/Q
        print(f"theta={th}: a call A(+)I on q={q} of Q={Q} coords: sum law charges p = {p:.3e} per coordinate,"
              f" dilution (OpenAI 2.2, g rho^alpha) charges p^theta = {p**th:.3e}; ratio {p**th/p:.3g}")

print()
print("== 4. Part A recheck (root of W m^z = sum r^z), with status of each shape ==")
cases = [
 ("2 in 3 at block 2", 2,2,[1]*3, "stand-in (ex-dimension representation game), not an engine unit"),
 ("3 in 5 at block 2", 3,2,[1]*5, "stand-in, not an engine unit"),
 ("W=2 block 3 ranks 1,1,1,1,1", 2,3,[1]*5, "engine 2 in 5 at block 3: PROVED IMPOSSIBLE (BRIEF item 3)"),
 ("W=2 block 3 ranks 2,1,1,1", 2,3,[2,1,1,1], "engine 2 in 5 at block 3: PROVED IMPOSSIBLE"),
 ("W=2 block 3 ranks 2,2,1", 2,3,[2,2,1], "engine 2 in 5 at block 3: PROVED IMPOSSIBLE"),
 ("W=3 block 3 eight rank 1", 3,3,[1]*8, "w-engine: 'a formula for an object nobody has'"),
 ("W=3 block 3 ranks 2,2,2,2", 3,3,[2]*4, "w-engine: 'a formula for an object nobody has'"),
 ("5 in 14 at block 3", 5,3,[1]*14, "hypothetical shape"),
 ("2 in 7 at block 4", 2,4,[1]*7, "hypothetical shape"),
 ("3 in 10 at block 4", 3,4,[1]*10, "hypothetical shape"),
 ("W=2 block 4 ranks 2,2,2,1", 2,4,[2,2,2,1], "hypothetical shape"),
 ("plain 4 in 20 at block 5", 4,5,[1]*20, "the plain schedule: a real unit"),
]
def root(W, m, rs):
    lo, hi = 1e-12, 1.0
    f = lambda z: W*m**z - sum(r**z for r in rs)
    for _ in range(200):
        mid = (lo+hi)/2
        if f(mid) > 0: hi = mid
        else: lo = mid
    return hi
for name, W, m, rs, st in cases:
    z = root(W, m, rs)
    print(f"  {name:32s} saving 1-z = {1-z:.4f}   [{st}]")
# mixed targets
def root_mixed(T, rs):
    lo, hi = 1e-12, 1.0
    f = lambda z: sum(W*m**z for W,m in T) - sum(r**z for r in rs)
    for _ in range(200):
        mid = (lo+hi)/2
        if f(mid) > 0: hi = mid
        else: lo = mid
    return hi
for name, T, rs in [("1+4 in 10 (mixed)", [(1,3),(4,2)], [1]*10), ("7+38 in 88 (mixed)", [(7,3),(38,2)], [1]*88)]:
    print(f"  {name:32s} saving 1-z = {1-root_mixed(T, rs):.4f}   [ex-laser example shape]")

print()
print("== 5. A finite win (A, b, g) priced two ways: H_theta's sum law vs OpenAI's own amplification (Lemma 2.2) ==")
# word for A^(x)b on Q = q^b coordinates with g gate calls A (+) I_{Q-q}
# sum law + homogeneity: q^(b-1) b^theta <= g  -> refutes theta-prices for theta > ln(g/q^(b-1))/ln b
# OpenAI amplification (sector formula, g rho^alpha < 1, rho = 1/(b q^(b-1))): exponent alpha = ln g / ln(b q^(b-1))
for (q,b,g) in [(2,2,3),(2,3,11),(2,4,31),(3,2,5)]:
    th_sum = math.log(g/q**(b-1))/math.log(b)
    alpha = math.log(g)/math.log(b*q**(b-1))
    print(f"  q={q} b={b} g={g}: sum-law threshold {th_sum:.4f} (saving {1-th_sum:.4f}) vs amplification exponent {alpha:.4f} (saving {1-alpha:.4f})")
print("  => with the sum law a sparse gate call is priced too cheaply; H_theta prices disagree with the exponent the")
print("     win actually yields, so they are not theta-prices of Proposition 1.")
