#!/usr/bin/env python3
"""v_mixed.py (v-e3-slice6): my own mixed-bound star recursion for U_6 + Z, to check the claim that (63,63) at W=14
is excluded, and to re-test the 95 open W=13 cases. Per even weight-16 support F avoiding Z:
 sigma_A (hyperplanes) = min(8, min slack of odd 16-supports G through A); planes k <= sum_{3 hyp through P} sigma - 8,
 also k <= min sigma_A and <= 4; lines 5k <= sum_{7 planes} - 4; points 13k <= sum_{15 lines} - 2. Any negative bound,
 or no perfect matching in the allowed entry pattern (E_xx needs point bound >= 1, E_xy needs line bound >= 1),
 makes E[F,F] singular, contradicting Theorem 10(b).
"""
import itertools
m = 6; pts = list(range(1, 64))
mons = [(j, k) for j in range(m) for k in range(j, m)]
monmask = [sum(1 << z for z in pts if (z >> j) & 1 and (z >> k) & 1) for (j, k) in mons]
odd, even = set(), set(); mask = 0; tr = 0
for g in range(1 << len(mons)):
    if g:
        b = (g & -g).bit_length() - 1; mask ^= monmask[b]; tr ^= (mons[b][0] == mons[b][1])
    if bin(mask).count('1') == 16: (odd if tr else even).add(mask)
odd = sorted(odd); even = sorted(even)
def subflats(Fp, d):
    a = Fp[0]; D = [p ^ a for p in Fp]; out = set()
    for vs in itertools.combinations([x for x in D if x], d):
        S = {0}
        for v in vs: S |= {s ^ v for s in S}
        if len(S) != 1 << d: continue
        for t in D: out.add(frozenset(p for p, dd in zip(Fp, D) if (dd ^ t) in S))
    return sorted(out, key=sorted)
def tomask(S): return sum(1 << z for z in S)
geo = {}
for F in even:
    Fp = [z for z in pts if (F >> z) & 1]
    geo[F] = (Fp, subflats(Fp, 3), subflats(Fp, 2), subflats(Fp, 1))
def pm(nodes, allowed):
    match = {}
    def aug(u, seen):
        for v in nodes:
            if allowed(u, v) and v not in seen:
                seen.add(v)
                if v not in match or aug(match[v], seen): match[v] = u; return True
        return False
    return all(aug(u, set()) for u in nodes)
def excluded(Z, W):
    for F in even:
        if any((F >> z) & 1 for z in Z): continue
        Fp, H3, P2, L1 = geo[F]
        sig = {}
        for A in H3:
            Am = tomask(A)
            c = [16 + sum((G >> z) & 1 for z in Z) - W for G in odd if G & Am == Am]
            sig[A] = min(c + [8])
        sp = {P: min([4, sum(sig[A] for A in H3 if P <= A) - 8] + [sig[A] for A in H3 if P <= A]) for P in P2}
        if min(sp.values()) < 0: return 'neg plane'
        sl = {L: min([2, (sum(sp[P] for P in P2 if L <= P) - 4) // 5] + [sp[P] for P in P2 if L <= P]) for L in L1}
        if min(sl.values()) < 0: return 'neg line'
        spt = {x: min([1, (sum(sl[L] for L in L1 if x in L) - 2) // 13] + [sl[L] for L in L1 if x in L]) for x in Fp}
        if min(spt.values()) < 0: return 'neg point'
        if not pm(Fp, lambda u, v: spt[u] >= 1 if u == v else sl[frozenset((u, v))] >= 1): return 'no matching'
    return None
print("(63,63) W=14:", excluded((63, 63), 14))
print("calibration: U_6 W=13 (Z empty):", excluded((), 13), "; U_6 W=12:", excluded((), 12))
opens = [((63,), 13)] + [(Z, 13) for Z in itertools.combinations_with_replacement(pts, 2) if 63 in Z or Z[0] ^ Z[1] == 63]
res = [(Z, W, excluded(Z, W)) for Z, W in opens]
print("W=13 cases involving 63:", len(res), "; excluded by mixed bounds:", sum(1 for r in res if r[2]))
