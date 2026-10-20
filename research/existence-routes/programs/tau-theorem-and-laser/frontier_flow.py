# frontier_flow.py (meta/ex-laser, 14 Oct 2026).
# Explicit schedule for a MIXED unit at block 3 (W3 arrays get the 3c-bit kernel, W2 arrays a 2c-bit rank-2 target,
# M paid moves = c-bit jobs), in the per-entry job model of meta/literature/tau_schedule_lp.py.
# Construction (PROOF.md section 2): lattice columns c = 3^a 2^b with 3c <= L; u(a,b) = binom(a+b,a) p^a q^b
# (p = W3/M, q = W2/M) units at column c; then every size c with 3c <= L is produced exactly as often as consumed
# (Pascal), size 1 leaves M*u(0,0) = M leaf jobs, and the sizes in (L/3, L] must be supplied ("frontier sources").
# A job of size n is handled by: give each frontier size s a share w(s) of its batch, split that share as s + (n-s),
# send the s-part into the flow, and handle the (n-s)-part the same way (recursively).
# Checks: (1) exact rational balance of the flow at every lattice size; (2) Phi-conservation
# sum_s Source(s) s^tau_c = M (leaves); (3) T(n)/n^tau_c stays bounded; (4) comparison with the LP optimum.
import sys
from fractions import Fraction as Fr
from math import comb, log

def tau_c(W3, W2, M, W1=0):
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if W3 * 3 ** mid + W2 * 2 ** mid + W1 > M: hi = mid
        else: lo = mid
    return lo

def lattice(L):
    pts = []
    a = 0
    while 3 ** a <= L:
        b = 0
        while 3 ** a * 2 ** b <= L:
            pts.append((a, b)); b += 1
        a += 1
    return pts

def flow(L, W3, W2, M, exact=False):
    """Return (units dict, sources dict size->weight, leaves) for the frontier construction with cutoff L."""
    p = Fr(W3, M) if exact else W3 / M
    q = Fr(W2, M) if exact else W2 / M
    size = lambda a, b: 3 ** a * 2 ** b
    U = {}
    for (a, b) in lattice(L):
        if 3 * size(a, b) <= L:
            U[(a, b)] = comb(a + b, a) * p ** a * q ** b
    prod, cons = {}, {}
    for (a, b), u in U.items():
        c = size(a, b)
        prod[c] = prod.get(c, 0) + M * u
        cons[3 * c] = cons.get(3 * c, 0) + W3 * u
        cons[2 * c] = cons.get(2 * c, 0) + W2 * u
    src, leaves, bad = {}, 0, []
    for s in set(prod) | set(cons):
        d = prod.get(s, 0) - cons.get(s, 0)
        if s == 1:
            leaves = d
        elif s in prod:
            if (d != 0) if exact else abs(d) > 1e-12: bad.append((s, d))
        else:
            src[s] = -d
    return U, src, leaves, bad

def T_table(nmax, W3, W2, M, tau):
    """T[n]: per-entry cost (leaf jobs) of the recursive frontier schedule; plain cost n is always allowed."""
    T = [0.0, 1.0]
    best_L = [0, 0]
    cache = {}
    for n in range(2, nmax + 1):
        best, bL = float(n), 0
        for L in range(3, n + 1):
            if L not in cache:
                U, src, leaves, bad = flow(L, W3, W2, M)
                tot = sum(src.values())
                cache[L] = [(s, w / tot) for s, w in src.items()], leaves / tot
            dist, lv = cache[L]
            val = lv + sum(w * T[n - s] for s, w in dist)
            if val < best: best, bL = val, L
        T.append(best); best_L.append(bL)
    return T, best_L

if __name__ == "__main__":
    print("check (1)+(2): exact balance and Phi conservation", flush=True)
    for (W3, W2, M) in ((1, 4, 10), (1, 1, 4), (7, 38, 88), (1, 6, 14)):
        tc = tau_c(W3, W2, M)
        for L in (30, 200, 3000):
            U, src, leaves, bad = flow(L, W3, W2, M, exact=True)
            phi = sum(float(w) * s ** tc for s, w in src.items())
            print(f"  {W3}+{W2} in {M}, L={L}: units at {len(U)} columns, sources at {len(src)} sizes "
                  f"(range {min(src)}..{max(src)}), unbalanced sizes {len(bad)}, leaves {leaves}, "
                  f"sum Source*s^tau_c = {phi:.10f}", flush=True)
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 480
    print(f"check (3): T(n)/n^tau_c for the recursive frontier schedule, n <= {nmax}", flush=True)
    for (W3, W2, M) in ((1, 4, 10), (1, 1, 4), (7, 38, 88)):
        tc = tau_c(W3, W2, M)
        T, bL = T_table(nmax, W3, W2, M, tc)
        pts = [n for n in (60, 120, 240, 480, 960, 1920, 3840) if n <= nmax]
        print(f"  {W3}+{W2} in {M}: tau_c = {tc:.4f}; " +
              "; ".join(f"n={n}: T={T[n]:.2f}, T/n^tc={T[n]/n**tc:.4f}, L={bL[n]}" for n in pts), flush=True)
        mx = max(T[n] / n ** tc for n in range(2, nmax + 1))
        print(f"     max over 2<=n<={nmax} of T/n^tc = {mx:.4f}", flush=True)
