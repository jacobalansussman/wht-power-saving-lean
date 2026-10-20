#!/usr/bin/env python3
"""fc-reopt: repair an integral class solution (SOL-*.json) whose end of clock violates the suffix-interval condition,
WITHOUT changing the global pattern multiset (so the price and the profile are unchanged): swap the patterns of single
bins between a class open in the violated suffix and a class of the same capacity that is not, choosing swaps that
add slots of the short widths to the suffix. After each swap: every Hall interval row (all A <= B, every width) and every
suffix condition (last LAST ticks, every set of widths) is re-checked; a swap that breaks one is undone.
usage: repair_end.py LINE HUB3815 ORDER SOL.json OUT.json [LAST=120] [MAXSWAPS=400]"""
import sys, os, json
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hall_fc as HF
line, hub, orderp, solp, outp = sys.argv[1:6]
LAST = int(sys.argv[6]) if len(sys.argv) > 6 else 120; MAXS = int(sys.argv[7]) if len(sys.argv) > 7 else 400
order = json.loads(open(orderp).read()); T = 1 + max(r for _, r in order)
if line == '3945':
    sup = [26, 2, 18, 53, 96, 213]; local = {1: 2440, 2: 725, 3: 224, 4: 102, 5: 20, 8: 8}
else:
    prov = json.load(open(os.path.join(os.path.dirname(orderp), 'PROVIDERS.json')))
    inv = Counter(prov['gauge44']['widths']); sup = [inv.get(w, 0) for w in (4, 5, 6, 7, 8, 9)]
    local = {int(k): v for k, v in prov['gauge44']['actual_paid_profile'].items()}
H = HF.load(os.path.join(hub, 'research/staged544-clock-ceiling/construction/hall_lp.py'), T, sup, local, 24, False, line == '3945')
cap, lo, hi = H.classes(order); N = H.N; K = len(cap); SUP = np.array(sup)
sol = json.load(open(solp))
cnt = defaultdict(int)
for k, j, v in sol['cols']: cnt[(int(k), int(j))] += int(round(v))
SUBS = list(H.SUBSETS)
SL = {c: np.array([[sum(int(P[j, w]) for w in range(6) if S >> w & 1) for S in SUBS] for j in range(len(P))]) for c, (P, _) in H.PATS.items()}
lo1 = lo + 1


def demand():
    d = np.zeros((K, 6))
    for (k, j), n in cnt.items():
        if n: d[k] += n * H.PATS[int(cap[k])][0][j]
    return d


def hall_ok(d):
    Bidx = np.arange(N)
    for w in range(6):
        M = np.zeros((N, N)); np.add.at(M, (lo1, hi), d[:, w])
        Dm = np.cumsum(np.cumsum(M[::-1], axis=0)[::-1], axis=1)
        V = Dm - SUP[w] * (Bidx[None, :] - Bidx[:, None] + 1)
        V[np.tril_indices(N, -1)] = -np.inf
        if V.max() > 1e-9: return False
    return True


def suffix_viol():
    out = []
    for A in range(N - LAST, N):
        ov = np.maximum(np.minimum(hi, N - 1) - np.maximum(lo1, A) + 1, 0)
        capS = np.zeros(len(SUBS))
        for (k, j), n in cnt.items():
            if n and ov[k] > 0: capS += n * np.minimum(SL[int(cap[k])][j], ov[k])
        need = H.SUBSET_SUP * (N - A)
        bad = np.where(capS < need)[0]
        if len(bad): out.append((A, bad, capS, need))
    return out

viol = suffix_viol(); print('initial suffix violations:', [(A, [int(SUBS[i]) for i in b[:3]]) for A, b, c, n in viol][:10], flush=True)
assert hall_ok(demand()), 'Hall rows already violated'
swaps = 0; tried = 0
while viol and swaps < MAXS:
    A, bad, capS, need = viol[0]
    ov = np.maximum(np.minimum(hi, N - 1) - np.maximum(lo1, A) + 1, 0)
    i = bad[0]                                       # the short width set
    late = [k for k in range(K) if ov[k] > 0 and k % 9 != 4]          # never touch the head classes (children 8/10)
    early = [k for k in range(K) if ov[k] == 0 and k % 9 != 4]
    bycls = defaultdict(list)
    for (k, j), n in cnt.items():
        if n: bycls[k].append(j)
    cands = []
    for kl in late:
        c = int(cap[kl])
        for j in bycls[kl]:
            for ke in early:
                if int(cap[ke]) != c: continue
                for j2 in bycls[ke]:
                    if j2 == j: continue
                    gain = min(SL[c][j2][i], ov[kl]) - min(SL[c][j][i], ov[kl])
                    if gain > 0: cands.append((-gain, -int(hi[ke]), kl, j, ke, j2))
    cands.sort()
    done = False
    for g, _, kl, j, ke, j2 in cands[:400]:
        cnt[(kl, j)] -= 1; cnt[(kl, j2)] += 1; cnt[(ke, j2)] -= 1; cnt[(ke, j)] += 1
        tried += 1
        if hall_ok(demand()):
            done = True; break
        cnt[(kl, j)] += 1; cnt[(kl, j2)] -= 1; cnt[(ke, j2)] += 1; cnt[(ke, j)] -= 1
    if not done: print('no Hall-safe improving swap for A=%d set %d among %d candidates' % (A, int(SUBS[i]), len(cands))); break
    swaps += 1
    viol = suffix_viol()
    print('swap %d (tried %d): class %d pattern %d <-> class %d pattern %d (gain %d on width set %d at A=%d); violations left %d' % (
        swaps, tried, kl, j, ke, j2, -g, int(SUBS[i]), A, len(viol)), flush=True)
json.dump(dict(cols=[[k, j, n] for (k, j), n in sorted(cnt.items()) if n], repaired_from=solp, swaps=swaps,
               remaining_violations=len(viol)), open(outp, 'w'))
print('written', outp, 'swaps', swaps, 'remaining suffix violations', len(viol))
