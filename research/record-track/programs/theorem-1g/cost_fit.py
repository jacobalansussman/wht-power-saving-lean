# cost_fit.py (bound-check/rl-vprove-compute, 14 Oct 2026). Skeptic's own cost recursion for the schedule of
# rl-prove Theorem 1G (dilated lockstep tree + top split + exit calls); all ranks 1 = ex-laser Theorem 1.
# Merged (by (depth, rho) with path multiplicity) computation written independently; cross-checked against
# explicit_sim.py counts at small L (printed). Cost per entry, weights normalised to source total 1 (lambda = 1/#src):
#   flow(n)  = cu * lambda * #units + sum_src lambda T(n - s) + sum_exit lambda T(e)          (L = n)
#   T(n)     = min(n, flow(n))   ["best"]      or   flow(n) whenever the tree exists  ["forced"]
# cu = per-unit overhead in array-passes, given on the command line as a multiple of W (0 = paid structure only;
# 3 = one regroup pass + F = 2W free register-passes per unit).
# Reports T/n^tau, T/(n^tau ln^kappa n), local slopes ln(T(2n)/T(n))/ln 2, and least-squares slopes.
import sys
from math import log, factorial, prod, ceil
from collections import defaultdict

def tau_c(Wt, rounds):
    Pr = rounds if isinstance(rounds, dict) else {r: rounds.count(r) for r in set(rounds)}
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if sum(n * r ** mid for r, n in Pr.items()) > sum(w * k ** mid for k, w in Wt.items()): lo = mid
        else: hi = mid
    return (lo + hi) / 2

class Rule:
    def __init__(self, Wt, rounds):
        self.Wt = Wt; self.K = sorted(Wt); self.m = max(self.K); self.rounds = rounds
        self.Pr = defaultdict(int)
        if isinstance(rounds, dict): self.Pr.update(rounds)
        else:
            for r in rounds: self.Pr[r] += 1
        self.W = sum(Wt.values()); self.tau = tau_c(Wt, rounds)
        self._u = {}
    def u(self, a):
        if a not in self._u:
            self._u[a] = factorial(sum(a)) // prod(factorial(x) for x in a) * prod(self.Wt[k] ** x for k, x in zip(self.K, a))
        return self._u[a]
    def c(self, a): return prod(k ** x for k, x in zip(self.K, a))
    def pts(self, rho, size, L):
        K, m = self.K, self.m; out = []
        def rec(i, a, cc, left):
            if i == len(K) - 1:
                c2 = cc * K[i] ** left
                if m * rho * c2 <= L: out.append((tuple(a + [left]), c2))
                return
            for e in range(left + 1):
                c2 = cc * K[i] ** e
                if m * rho * c2 > L: break
                rec(i + 1, a + [e], c2, left - e)
        rec(0, [], 1, size); return out
    def flow(self, L):
        K, m = self.K, self.m
        if m > L: return None
        G = 0
        while self.pts(1, G + 1, L): G += 1
        src = defaultdict(int); ex = defaultdict(int); units = 0; unbalanced = 0
        layer = {1: 1}                       # rho -> multiplicity at depth d
        for d in range(G + 1):
            for rho, N in layer.items():
                for a, ca in self.pts(rho, G - d, L): units += N * self.u(a)
            if d == 0:
                for a, ca in self.pts(1, G, L):
                    for i, k in enumerate(K): src[k * ca] += self.Wt[k] * self.u(a)
            if d == G:
                for rho, N in layer.items():
                    if m * rho <= L:
                        for r, pr in self.Pr.items(): ex[rho * r] += N * pr
                break
            nxt = defaultdict(int)
            for rho, N in layer.items():
                par = {a: self.u(a) for a, ca in self.pts(rho, G - d, L)}
                for r, pr in self.Pr.items():
                    rc = rho * r; nxt[rc] += N * pr
                    dem = defaultdict(int)
                    for a, ca in self.pts(rc, G - d - 1, L):
                        for i, k in enumerate(K):
                            b = list(a); b[i] += 1; dem[tuple(b)] += self.Wt[k] * self.u(a)
                    for b, x in dem.items():
                        if b in par:
                            if par[b] != x: unbalanced += 1
                        else:
                            s = rc * self.c(b); src[s] += N * pr * x
                    for b in par:
                        if b not in dem: unbalanced += 1
            layer = nxt
        return dict(G=G, src=src, ex=ex, units=units, unbalanced=unbalanced)

def run(name, Wt, rounds, Nmax, cu_list=(0,)):
    R = Rule(Wt, rounds); tau = R.tau; m = R.m; kmin = min(R.K); rmax = max(R.Pr)
    theta = log(rmax) / log(kmin) if rmax > 1 else 0.0
    eta = min((1 - (1 - t) ** tau) / t ** tau for t in [1 / m + (1 - 1 / m) * i / 20000 for i in range(20001)])
    kappa = 0 if theta == 0 else max(0, ceil(log(2 / eta) / log(1 / theta))) if theta < 1 else float('inf')
    print(f"== {name}: tau_c={tau:.6f} theta={theta:.4f} eta_m={eta:.4f} kappa={kappa} W={R.W} P={sum(R.Pr.values())}", flush=True)
    flows = {}
    bad = 0; maxexit_ratio = 0.0; src_band_bad = 0; phi_dev = 0.0
    for n in range(1, Nmax + 1):
        f = R.flow(n)
        if f is None: continue
        S = sum(f['src'].values()); bad += f['unbalanced']
        for s in f['src']:
            if not (n / m < s <= n): src_band_bad += 1
        ps = sum(x * s ** tau for s, x in f['src'].items()); pe = sum(x * e ** tau for e, x in f['ex'].items())
        phi_dev = max(phi_dev, abs(pe / ps - 1))
        if theta > 0: maxexit_ratio = max(maxexit_ratio, max(f['ex']) / (rmax * (n / m) ** theta))
        flows[n] = ([(s, x / S) for s, x in f['src'].items()], [(e, x / S) for e, x in f['ex'].items()], f['units'] / S, ps / S)
    print(f"   structure over n<= {Nmax}: unbalanced points {bad}; sources outside (n/m, n]: {src_band_bad}; "
          f"max |Phi_exit/Phi_src - 1| = {phi_dev:.2e}; max exit / (r_max (n/m)^theta) = {maxexit_ratio:.4f}", flush=True)
    res = {}
    for cuw in cu_list:
        cu = cuw * R.W
        for mode in ("best", "forced"):
            T = [0.0] * (Nmax + 1)
            for n in range(1, Nmax + 1):
                if n in flows:
                    srcs, exs, uu, phat = flows[n]
                    v = cu * uu + sum(w * T[n - s] for s, w in srcs) + sum(w * T[e] if e < n else float('inf') for e, w in exs)
                    T[n] = min(n, v) if mode == "best" else v
                else:
                    T[n] = float(n)
            res[(cu, mode)] = T
            line = [f"   cu={cu} {mode}:"]
            for n in [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]:
                if n <= Nmax:
                    lk = max(1.0, log(n)) ** (kappa if kappa != float('inf') else 0)
                    line.append(f"n={n} T/n^t={T[n]/n**tau:.4f} T/(n^t ln^k)={T[n]/n**tau/lk:.3g}")
            print("  ".join(line), flush=True)
            sl = []
            n = 100
            while 2 * n <= Nmax:
                # average over a window to smooth lattice effects
                lo = sum(T[i] for i in range(n, n + n // 10)) / (n // 10); hi = sum(T[i] for i in range(2 * n, 2 * n + n // 5) if i <= Nmax) / len([i for i in range(2 * n, 2 * n + n // 5) if i <= Nmax])
                sl.append(f"{n}->{2*n}:{log(hi/lo)/log(2):.4f}"); n *= 2
            print(f"     local slopes (window-averaged) vs tau={tau:.4f}: " + " ".join(sl), flush=True)
            # least squares over [Nmax/16, Nmax]
            xs = [log(i) for i in range(Nmax // 16, Nmax + 1)]; ys = [log(T[i]) for i in range(Nmax // 16, Nmax + 1)]
            mx = sum(xs) / len(xs); my = sum(ys) / len(ys)
            b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
            print(f"     LSQ slope of ln T on ln n over [{Nmax//16},{Nmax}]: {b:.4f}  (tau={tau:.4f}; diff {b - tau:+.4f}); "
                  f"max_n T/n^tau = {max(T[i]/i**tau for i in range(1, Nmax+1)):.4f}", flush=True)
    return R, flows, res

RULES = {
    "D": ("D rank-1 control 1+4 in 10 (Thm 1)", {3: 1, 2: 4}, [1] * 10),
    "E": ("E rank-1 7+38 in 88 (Thm 1)", {3: 7, 2: 38}, [1] * 88),
    "A": ("A {4:1,3:2} r=[1,1,1,1,2,2]", {4: 1, 3: 2}, [1, 1, 1, 1, 2, 2]),
    "B": ("B {7:1,4:1} r=[1,1,1,1,2,2,2]", {7: 1, 4: 1}, [1, 1, 1, 1, 2, 2, 2]),
    "C": ("C {5:1,4:3} r=[1]*6+[2]*2+[3]*2", {5: 1, 4: 3}, [1] * 6 + [2] * 2 + [3] * 2),
    "H": ("H {3:2} r=[2,1,2]", {3: 2}, [2, 1, 2]),
    "J": ("J {5:2} r=[4,1,4] theta=.861", {5: 2}, [4, 1, 4]),
    "T": ("T74 mixed what-if rule {73:43556,25:22748}, paid {1:1560432,2:336336,3:278724,4:84912,5:54612,8:5328}",
          {73: 43556, 25: 22748}, {1: 1560432, 2: 336336, 3: 278724, 4: 84912, 5: 54612, 8: 5328}),
    "X": ("X NON-separated {3:1,2:4} r=[1]*6+[2]*2", {3: 1, 2: 4}, [1, 2, 1, 1, 2, 1, 1, 1]),
}

if __name__ == "__main__":
    which = sys.argv[1]; Nmax = int(sys.argv[2]); cus = tuple(float(x) for x in sys.argv[3].split(",")) if len(sys.argv) > 3 else (0.0,)
    name, Wt, rl = RULES[which]
    run(name, Wt, rl, Nmax, cus)
