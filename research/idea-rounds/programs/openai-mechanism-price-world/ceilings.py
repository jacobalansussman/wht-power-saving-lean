"""Ceilings for any quantitative version of OpenAI's reflection contradiction.

Pure arithmetic (no search). Run: python3 -I ceilings.py

Inputs read from 08-reflection.tex and GraphBlocks.lean:
  lower bound  p(R) >= (d - 4) w,  upper bound p(R) <= (3d/4 + 7/2) w,
  so the cheap side computes M (tensor price d w) at ratio r(d) = 3/4 + 7.5/d,
  w = 4^(d+1) coordinates, M = Gamma (x) Gamma lives on 2d+2 binary directions.

Three ways of turning a ratio r at size w into a saving:
  fixed  : OpenAI's amplification for a fixed matrix, s = ln(1/r) / ln(b q^(b-1)),
           here taken as U^(x)2d (q = 2, b = 2d)  [same order as ln(1/r)/ln w]
  width  : s = ln(1/r) / ln w  (dilution scale, the skeptic's charge in v-ex2-reform)
  uniform: a unit right for every f at block m = 2d+2, s = ln(1/r) / ln m
"""
import math

def r_of(d, c=7.5, pair=1.5):
    # cheap side per coordinate: pair*(d/2+1) + 2 + 4 = (pair/2) d + (pair + 6); OpenAI: 3d/4 + 7.5
    return (pair / 2 * d + (pair + 6 if c is None else c)) / d

def scan(label, fn, ds=range(31, 2001)):
    best = max(((fn(d), d) for d in ds if r_of(d) < 1), default=None)
    return best

print("OpenAI's family, dense pair 3/2, overhead 7.5 w")
for name, den in [("fixed U^(x)2d", lambda d: math.log(2 * d) + (2 * d - 1) * math.log(2)),
                  ("width ln w", lambda d: (d + 1) * math.log(4)),
                  ("uniform block 2d+2", lambda d: math.log(2 * d + 2))]:
    vals = [(math.log(1 / r_of(d)) / den(d), d) for d in range(31, 3001) if r_of(d) < 1]
    s, d = max(vals)
    print(f"  {name:20s} max s = {s:.4e} at d = {d}  (r = {r_of(d):.4f})")
for d in (31, 32, 45, 60, 80, 250):
    r = r_of(d)
    print(f"  d={d:4d} r={r:.5f} fixed={math.log(1/r)/(math.log(2*d)+(2*d-1)*math.log(2)):.3e}"
          f" width={math.log(1/r)/((d+1)*math.log(4)):.3e} uniform={math.log(1/r)/math.log(2*d+2):.3e}")

print("\nIf the overhead 7.5 w could be cut to c w (same 3/2 dense pair), width-scale ceiling:")
for c in (7.5, 4, 2, 1, 0.5):
    vals = [(math.log(1 / ((0.75 * d + c) / d)) / ((d + 1) * math.log(4)), d) for d in range(2, 3001) if 0.75 * d + c < d]
    s, d = max(vals)
    print(f"  c={c:4}: max s = {s:.4e} at d = {d}")

print("\nIf a dense pair cost 1 instead of 3/2 (cheap side d/2 + 7 per coordinate):")
for name, den in [("width ln w", lambda d: (d + 1) * math.log(4)),
                  ("uniform block 2d+2", lambda d: math.log(2 * d + 2))]:
    vals = [(math.log(d / (0.5 * d + 7)) / den(d), d) for d in range(15, 5001) if 0.5 * d + 7 < d]
    s, d = max(vals)
    print(f"  {name:20s} max s = {s:.4e} at d = {d}")

print("\nStreaming in the block: uniform ceiling if the block is k*(2d+2) directions")
for k in (1, 10**3, 10**4, 10**5, 10**6):
    vals = [(math.log(1 / r_of(d)) / math.log(k * (2 * d + 2)), d) for d in range(31, 3001) if r_of(d) < 1]
    s, d = max(vals)
    print(f"  k={k:>8}: max s = {s:.4e} at d = {d};  fixed-matrix with k: {max(math.log(1/r_of(dd))/(k*(2*dd+2)*math.log(2)) for dd in range(31,3001) if r_of(dd)<1):.3e}")

print("\nLean packing (InventoryValidity.scalable_square_win): T=(2W-1)u calls vs 2Wu, so r = 1-1/(2W),")
print("q = W u, b = 2, s = ln(1/r)/ln(2 W u) ~ 1/(2W ln(2Wu)); W >= K*B + I >= 5*B with B = sum of type widths.")
for d in (32, 60):
    w = 4 ** (d + 1)
    W = 5 * w  # most favourable: R is the widest type, Z = H = 1, u ~ W
    s = -math.log1p(-1 / (2 * W)) / math.log(2 * W * W)
    print(f"  d={d}: B>=w={w:.3e}, W>={W:.3e}: s <= {s:.3e}")

print("\nFixed-matrix versus uniform amplification of the same ratio (record-like unit, block 73, s=1.7173e-3):")
m, s = 73, 1.7173e-3
r = m ** (-s)
print(f"  r = {r:.6f}; uniform gives {math.log(1/r)/math.log(m):.4e}; fixed (q=2,b=m) gives "
      f"{math.log(1/r)/(math.log(m)+(m-1)*math.log(2)):.4e}; factor {(math.log(m)+(m-1)*math.log(2))/math.log(m):.1f}")
