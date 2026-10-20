# v-e3-boost vbig.py: (1) my own high-precision roots of the J-balance for the hub words T74E and #1042 (decimal);
# (2) distance-set certificates at large M via a tridiagonal recursion (column vectors), applied to the hub words, to
# test the validity of distance-set certificates with S = parity sets on REAL units.
import json, sys
from decimal import Decimal as Dc, getcontext
getcontext().prec = 50
import vcert
def root(W, m, calls):
    def D(t): return sum(Dc(n) * (Dc(r).ln() * t).exp() for r, n in calls.items()) - Dc(W) * (Dc(m).ln() * t).exp()
    lo, hi = Dc("0.99"), Dc(1)
    for _ in range(140):
        mid = (lo + hi) / 2
        if D(mid) < 0: hi = mid
        else: lo = mid
    return lo
p = json.load(open(__import__("os").environ["T74E_JSON"]))
T74 = (p["W"], p["m"], {int(r): n for r, n in p["profile"].items()})
L1042 = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}; I1042 = {1: 2249, 2: 1, 3: 400, 7: 3, 16: 47, 48: 4850}
c1042 = {r: 144 * n for r, n in L1042.items()}
for r, n in I1042.items(): c1042[r] = c1042.get(r, 0) + n
H1042 = (14336, 73, c1042)
for name, (W, m, calls) in (("T74E", T74), ("#1042", H1042)):
    t = root(W, m, calls)
    print(f"{name}: W={W} m={m} sum calls rank = {sum(r*n for r,n in calls.items())} vs W*m = {W*m}; "
          f"root tau = {str(t)[:20]}  saving 1 - tau = {str(1 - t)[:18]}", flush=True)
print("T74E root_lower =", str(p["exact_rate"])[:120])

def colvec(m, j):
    a, b, c, k = vcert.iarray(m); n = m + 1
    def Lv(x):  # (L x)_i = sum_t L[i][t] x_t,  L[t-1][t] = b[t-1], L[t][t] = a[t], L[t+1][t] = c[t+1]
        y = [0] * n
        for t in range(n):
            if x[t] == 0: continue
            if t - 1 >= 0: y[t - 1] += b[t - 1] * x[t]
            y[t] += a[t] * x[t]
            if t + 1 < n: y[t + 1] += c[t + 1] * x[t]
        return y
    e = [0] * n; e[j] = 1
    V = [e, Lv(e)]
    for i in range(1, m):
        X = Lv(V[i]); Y = [X[r] - a[i] * V[i][r] - b[i - 1] * V[i - 1][r] for r in range(n)]
        assert all(y % c[i + 1] == 0 for y in Y)
        V.append([y // c[i + 1] for y in Y])
    return V   # V[i][d] = p^d_{ij}
def dprof(M, S, ds):
    cols = {j: colvec(M, j) for j in S}
    ks = {j: cols[j][j][0] for j in S}          # p^0_{jj} = k_j
    tot = sum(ks.values())
    return {d: tot - sum(cols[j][i][d] for i in S for j in S) for d in ds}
if __name__ == "__main__":
    # sanity vs vcert at M = 6
    import itertools
    for S in [(1,), (1, 3, 5), (0, 2, 4, 6), (2, 3)]:
        assert [dprof(6, S, range(1, 7))[d] for d in range(1, 7)] == vcert.dist_profile(6, S)
    print("column recursion = matrix recursion at M = 6: True", flush=True)
    for name, (W, m, calls) in (("T74E", T74), ("#1042", H1042)):
        ds = sorted(set(calls) | {m})
        for S in [tuple(range(1, m + 1, 2)), tuple(range(0, m + 1, 2)), (1,), (0, 1), tuple(range(1, m + 1))]:
            P = dprof(m, S, ds)
            margin = sum(n * P[r] for r, n in calls.items()) - W * P[m]
            print(f"  {name} Dist{m} S={'odd' if S == tuple(range(1, m+1, 2)) else 'even' if S == tuple(range(0, m+1, 2)) else S if len(S) < 5 else '1..m'}: "
                  f"margin/(W*P[m]) = {margin / (W * P[m]):.6f}  -> {'PASS' if margin >= 0 else 'FAIL'}", flush=True)
