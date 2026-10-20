# v-e3-boost vsurv.py: the survivors against MY certificates (A, C_k, distance sets with all S, the new line family with
# all T), on registers M = m .. m+4 (spare rows, b = 1) and tensor scaling b = 2 (M = 2m); balances and walk condition.
import sys, re, math
from itertools import combinations
from functools import lru_cache
import vcert, vlines

def G(n, k):
    if k < 0 or k > n: return 0
    r = 1
    for i in range(k): r = r * (2 ** (n - i) - 1) // 1
    den = 1
    for i in range(k): den *= 2 ** (i + 1) - 1
    return r // den
def Aprof(M): return [2 ** M - 2 ** (M - d) for d in range(1, M + 1)]
def Cprof(M, k): return [G(M, k) - G(M - d, k) for d in range(1, M + 1)]
@lru_cache(None)
def pm(M): return vcert.pmat(M)
def Dprof(M, S):
    p, ks = pm(M); tot = sum(ks[i] for i in S)
    return [tot - sum(p[d][i][j] for i in S for j in S) for d in range(1, M + 1)]
@lru_cache(None)
def linecounts(M): return [vlines.counts(M, d) for d in range(1, M + 1)]
def Lprof(M, T):
    out = []
    for J in linecounts(M):
        inV = sum(c for (tu, tv), c in J.items() if tv in T)
        both = sum(c for (tu, tv), c in J.items() if tv in T and tu in T)
        out.append(inV - both)
    return out
def D(v, t): return sum(c * d ** t for d, c in v.items())
def walk_ok(v, M):   # my own: reachable distances by the exact-step polygon recursion, per sub-multiset of calls
    calls = [d for d, c in v.items() if c > 0 for _ in range(c)]
    tasks = [d for d, c in v.items() if c < 0]
    reach = set()
    for r in range(1, len(calls) + 1):
        for R in set(combinations(sorted(calls), r)):
            # exact recursion x_i in [|r_i - x_{i-1}|, min(M, r_i + x_{i-1})], any order: take union over orders' interval
            cur = {0}
            for ri in R:
                cur = {y for x in cur for y in range(abs(ri - x), min(M, ri + x) + 1)}
            reach |= cur
    return all(k in reach for k in tasks)
def tests(v, m, Mextra=4):
    res = []
    for M in range(m, m + Mextra + 1):
        res.append((f"A{M}", sum(c * Aprof(M)[d - 1] for d, c in v.items())))
        for k in range(1, M + 1):
            P = Cprof(M, k); res.append((f"C{M}k{k}", sum(c * P[d - 1] for d, c in v.items())))
        for r in range(1, M + 2):
            for S in combinations(range(M + 1), r):
                P = Dprof(M, S); res.append((f"Dist{M}S{S}", sum(c * P[d - 1] for d, c in v.items())))
        TY = vlines.TYPES
        for r in range(1, 6):
            for T in combinations(TY, r):
                P = Lprof(M, set(T)); res.append((f"Lines{M}T{T}", sum(c * P[d - 1] for d, c in v.items())))
    M = 2 * m
    res.append((f"A{M}x2", sum(c * Aprof(M)[2 * d - 1] for d, c in v.items())))
    for k in range(1, M + 1):
        P = Cprof(M, k); res.append((f"C{M}k{k}x2", sum(c * P[2 * d - 1] for d, c in v.items())))
    for i in range(M + 1):
        for j in range(i, M + 1):
            P = Dprof(M, tuple(range(i, j + 1))); res.append((f"Dist{M}S{i}-{j}x2", sum(c * P[2 * d - 1] for d, c in v.items())))
    return res
def parse(line):
    t = re.search(r"tasks (\{[^}]*\}) calls (\{[^}]*\})", line)
    tasks, calls = eval(t.group(1)), eval(t.group(2))
    v = {}
    for k, w in tasks.items(): v[k] = v.get(k, 0) - w
    for r, p in calls.items(): v[r] = v.get(r, 0) + p
    return {d: c for d, c in v.items() if c}
if __name__ == "__main__":
    named = {"S6": {2: -1, 5: -4, 3: 2, 4: 1, 6: 2}, "S7": {2: -1, 5: -1, 6: -1, 3: 2, 7: 1},
             "S8": {2: -1, 5: -1, 6: -2, 3: 1, 4: 2, 8: 1}}
    for name, v in named.items():
        m = max(v)
        F = sum(v.values()); P1 = sum(c * d for d, c in v.items())
        # window: sign of D on a fine grid and the root by bisection
        g = [i / 100000 for i in range(1, 100000)]
        neg = [t for t in g if D(v, t) < 0]
        lo, hi = min(neg) - 1e-5, min(neg)
        for _ in range(60):
            mid = (lo + hi) / 2
            if D(v, mid) < 0: hi = mid
            else: lo = mid
        contiguous = len(neg) == 100000 - 1 - round(min(neg) * 100000) + 1
        res = tests(v, m)
        bad = [r for r in res if r[1] < 0]
        tight = [r[0] for r in res if r[1] == 0]
        print(f"{name} v={v}: F={F} plain={P1} D'(1)={sum(c*d*math.log(d) for d,c in v.items()):.5f} root={hi:.5f} "
              f"negative on ({hi:.4f},1) contiguous={contiguous} D(0.9)={D(v,0.9):.6f} D(0.9983)={D(v,0.9983):.3e}; "
              f"walk={walk_ok(v, m)}; {len(res)} of my tests, failed {bad[:5]}; tight {len(tight)} e.g. {tight[:6]}", flush=True)
    # every survivor of out_enum_678_proj.txt
    lines = [l for l in open("../out_enum_678_proj.txt") if "tasks {" in l]
    nb = 0
    for l in lines:
        v = parse(l); m = max(v)
        res = tests(v, m, Mextra=1); bad = [r for r in res if r[1] < 0]
        if bad: nb += 1; print("KILLED:", v, bad[:4], flush=True)
    print(f"{len(lines)} listed survivors re-tested with my certificates (M = m, m+1, b = 2): killed {nb}")
