#!/usr/bin/env python3
"""cx-medium / price_p112.py (own code, stdlib; python3 -I).  The hub template shrunk to ONE bank pair (P = 1):
three stages GF, GB, GF on one pair of 112-port banks (the three-shear swap: y += x, x -= y, y += x gives X = -y,
Y = x), windows of 8 fresh coordinates + the shared line, block m = 1 + 8*3 = 25.  Local word: the hub's U8 helper
circuit, our lowering E1 (fsp-price/data/words.json): per call profile {1:2548, 2:338, 3:446, 4:162, 5:72, 8:8}
(sources, targets, helpers' in-call blocks, the 8 copies of rank 8) and helper leftover widths
{4:80, 5:32, 7:54, 8:122, 9:216} (504 helpers, 3778 levels).
At P = 1 every data array is active in every stage: X climbs 8 per call (1 -> 25), Y 8 per call (0 -> 24), no idle
levels at all, so no helper can be hosted on a data array.  Every helper piece goes to a BANK: an extra array (in W)
whose 25 levels are the in-call blocks of the helpers it hosts (at most one piece per call, distinct calls) plus one
pad block.  Banks are D-neutral; D = 224 - 3*64 = 32 per replica.
Pricing: an-ours' anatomy_ours.price (the Lean fold hypothesis, s = 40), at m = 25.
Variants: (a) ideal packing (every bank exactly full: no pads; T = 25 replicas so the bank count is integral);
          (b) first-fit decreasing with one piece per call per bank (generous: no frame-dimension chaining imposed);
          (c) as (b) but each bank only within ONE replica (a bank is a register of one run; 3 calls);
          (d) helpers as fresh scratch arrays instead (not in W; every in-call helper block is a real loss of D)."""
import sys, os, json
from collections import Counter
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../anatomy-our-lean-words')))
from anatomy_ours import price, rank, narrow, add

PROF = {1: 2548, 2: 338, 3: 446, 4: 162, 5: 72, 8: 8}
WID = {4: 80, 5: 32, 7: 54, 8: 122, 9: 216}
V, CST, M, CALLS = 112, 64, 25, 3
assert rank(PROF) == 2 * V * 8 + CST + sum(w * n for w, n in WID.items()) == 5634

def ffd(T, per_replica_banks=False):
    items = []
    for c in range(CALLS * T):
        for w, n in WID.items(): items += [(w, c)] * n
    items.sort(key=lambda t: -t[0])
    bins = []                              # [load, set(calls), replica]
    open_by_rep = {}
    for w, c in items:
        rep = c // CALLS
        cand = open_by_rep.setdefault(rep if per_replica_banks else 0, [])
        for b in cand:
            if b[0] + w <= M and c not in b[1]:
                b[0] += w; b[1].add(c); break
        else:
            b = [w, {c}]; cand.append(b)
        # drop full bins from the candidate list to keep it short
        if cand and cand[-1][0] == M: pass
    allb = [b for lst in open_by_rep.values() for b in lst]
    return [M - b[0] for b in allb]

def evaluate(T, pads, label):
    H = Counter({r: CALLS * T * n for r, n in PROF.items()})
    for p in pads:
        if p: H[p] += 1
    W = 2 * V * T + len(pads)
    D = M * W - rank(H)
    assert D == T * (2 * V - CALLS * CST), (D, T)
    a = float(price(dict(H), W, m=M))
    N = narrow(dict(H), M)
    print(f"{label:58s} T={T:3d} W={W:7d} banks={len(pads):6d} pads={sum(pads):6d} D={D:5d} N={N:10.1f} D/N={D/N:.4e}  price a={a:.7e}"
          f"  (x{a/1.6797550e-3:.3f} of the hub)")
    return dict(label=label, T=T, W=W, banks=len(pads), pad_levels=sum(pads), D=D, N=N, a=a, hist=dict(sorted(H.items())))

out = {}
T = 25
per_call = sum(w * n for w, n in WID.items())
nb = CALLS * T * per_call // M; assert CALLS * T * per_call % M == 0
out['ideal'] = evaluate(T, [0] * nb, '(a) ideal packing, no pads')
for T in (1, 8, 25):
    out[f'ffd_T{T}'] = evaluate(T, ffd(T), '(b) first-fit, banks shared across replicas')
    out[f'ffd_rep_T{T}'] = evaluate(T, ffd(T, True), '(c) first-fit, banks within one replica')
# (d) helpers as scratch arrays outside W: their in-call blocks are pure loss
H = Counter({r: CALLS * n for r, n in PROF.items()}); W = 2 * V
D = M * W - rank(H)
print(f"(d) helpers as scratch outside W: D per replica = {D} (gain 224 - centres 192 - helper levels {CALLS * per_call}) -> no saving")
out['scratch'] = dict(D=D)
# first-order decomposition of the ideal case
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'price_p112.json'), 'w'), indent=1)
