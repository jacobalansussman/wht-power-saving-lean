# v-pf-induction-logic, 13 Oct 2026.  Test of PROOF.md Proposition 7(c), "general tensor products:
# T(F (x) G) <= 2 max(T(F), T(G))", on families whose dimension VARIES with the frame (allowed by the
# Averaging Lemma: "no equal dimensions" needed).  Cut families: F_U = span{e} if U in C, else 0.
# Then delta_F(U,V) = [U in C, V not in C], T(F) = P_m(cut C)/P_m... (directed cut ratio), and
# F_C (x) F_D = F_(C cap D).  Block m = 2 (15 frames) and m = 3 (135 frames, random search).
import itertools, random, sys
from fractions import Fraction as Fr

def lagrangians(m):
    n = 2 * m
    def form(u, v):
        a, b = u >> m, u & ((1 << m) - 1); c, d = v >> m, v & ((1 << m) - 1)
        return (bin(a & d).count('1') + bin(b & c).count('1')) & 1
    vecs = range(1, 1 << n)
    out = set()
    def span(gens):
        s = {0}
        for g in gens: s |= {x ^ g for x in s}
        return frozenset(s)
    def rec(gens, sp):
        if len(gens) == m: out.add(sp); return
        for v in vecs:
            if v in sp: continue
            if all(form(v, g) == 0 for g in gens):
                if gens and v < max(gens): continue
                rec(gens + [v], span(gens + [v]))
    rec([], frozenset({0}))
    return sorted(out, key=lambda s: sorted(s))

def setup(m):
    L = lagrangians(m); N = len(L)
    dist = [[m - (len(L[i] & L[j]).bit_length() - 1) for j in range(N)] for i in range(N)]
    return L, N, dist

def ratio(C, N, dist, m):
    num = {1: 0, m: 0}; cnt = {1: 0, m: 0}
    for i in range(N):
        for j in range(N):
            d = dist[i][j]
            if d in cnt:
                cnt[d] += 1
                if i in C and j not in C: num[d] += 1
    if num[1] == 0: return None
    return Fr(num[m], cnt[m]) / Fr(num[1], cnt[1])

if __name__ == '__main__':
    m = 2
    L, N, dist = setup(m)
    print(f"m={m}: {N} frames; distance counts", {d: sum(row.count(d) for row in dist) for d in range(m + 1)})
    # all cut ratios at m = 2
    R = {}
    for mask in range(1, (1 << N) - 1):
        C = frozenset(i for i in range(N) if mask >> i & 1)
        R[C] = ratio(C, N, dist, m)
    vals = sorted(set(R.values()))
    print("distinct cut ratios:", [str(v) for v in vals[:6]], "...", [str(v) for v in vals[-4:]])
    best = None
    sets = list(R.items())
    lowsets = [(C, r) for C, r in sets if r is not None]
    # for each target intersection with high ratio, find two supersets with low ratio meeting exactly in it
    byratio = sorted(lowsets, key=lambda t: t[1])
    for T, rT in sorted(lowsets, key=lambda t: -t[1])[:40]:
        sups = [(C, r) for C, r in lowsets if T <= C and C != T]
        sups.sort(key=lambda t: t[1])
        for (C1, r1), (C2, r2) in itertools.combinations(sups[:300], 2):
            if C1 & C2 != T: continue
            gain = rT / max(r1, r2)
            if best is None or gain > best[0]:
                best = (gain, T, rT, C1, r1, C2, r2)
    gain, T, rT, C1, r1, C2, r2 = best
    print(f"best: T(F)={r1} T(G)={r2} T(F(x)G)={rT}  ratio T(F(x)G)/max = {gain} = {float(gain):.4f}  (claim: <= 2)")
    print("   |C1|, |C2|, |C1 cap C2| =", len(C1), len(C2), len(T))
