# an-openai: price OpenAI's complex network two ways.
#  (A) OpenAI's own pricing (explicit-DFT Thm 2.6): every residual dimension is one directional call of size f:
#      lambda = M/W, saving = ln(mW/M)/ln m  (unpadded; the paper pads W to 2^71).
#  (B) Edge-as-one-child pricing (the hub's moment equation, as in hm-gate/gateprice.py's header):
#      every edge with residual dimension r is ONE call of size r*f (Phi_E is address-conjugate to C^{⊗r});
#      tau solves sum_edges (r_e/m)^tau = W, saving b = 1 - tau.
# Residual lists from the paper's table (verified by labels_check.py at h = 7, 8, 9).
from math import comb, log
def pieces(h):
    v = comb(h, 3); d = comb(h - 3, 3) + 3 * (h - 3); m = h ** 3
    out = {}  # size -> multiplicity
    def add(r, n):
        if r > 0 and n > 0: out[r] = out.get(r, 0) + n
    for j in (1, 2, 3):
        a = h ** (j - 1); f = h ** (3 - j); inv = v * v
        for _wire in ("X", "Y"):
            add((a - 1) * (h - 1), inv * v); add(h - 1, inv * v)
        c = h + 1
        add((a - 1) * h, inv * c); add(h, 3 * inv * c); add(a * h * (f - 1), inv * c)
        sd = v * d
        add(a - 1, inv * sd); add((a - 1) * (h - 1) + 1, inv * sd); add(h - 2, inv * sd); add(1, inv * sd)
        add(a * h * (f - 1), inv * sd)
    N = v ** 3; I = 3 * v * v
    W = 2 * N + I * (v * d + h + 1)
    return out, W, m
def tau_root(out, W, m):
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        s = sum(n * (r / m) ** mid for r, n in out.items())
        if s > W: lo = mid
        else: hi = mid
    return (lo + hi) / 2
if __name__ == "__main__":
    print(" h        m          W          saving(A: per dimension)   saving(B: edge = one child)")
    for h in (22, 23, 24, 25, 26, 28, 30, 40, 50, 100):
        out, W, m = pieces(h)
        M = sum(r * n for r, n in out.items())
        sA = log(m * W / M) / log(m)
        t = tau_root(out, W, m)
        print("%3d %9d %11.4e   %+.4e                 %+.4e" % (h, m, W, sA, 1 - t))
