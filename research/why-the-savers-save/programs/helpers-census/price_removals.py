#!/usr/bin/env python3
"""u-helpers-hub: price helper removals on #942 (T37 owner packing of the U8 word) with m-anatomy's pricer.

Control: #942 profile {1:848490,2:112554,3:148518,4:53946,5:23976,8:2664,48:11086}, m 73, W 33152, deficit 11840
         -> b 1664727258740678855772/10^24, kappa 1661960547713947/10^18 (both checked exactly).
Model of a removal (WHAT-IF, not a construction): N helpers deleted from the U8 word in all 333 calls
(37 replicas x 9 stages).  Their runs leave the local profile; their width returns to the idle stock.
  optimistic : all freed width is re-hosted so that it is paid as rank-48 idle children (best case; needs a re-pack)
  hole       : no re-pack: each deleted helper leaves a hole of its own width in its owner bin, paid as one child of
               rank w per call (feasible in the accounting, since T37's bins are full apart from 11,086 empty 48-bins)
The rank deficit, m and W do not change (helpers never move D: format-limits Lemma H).
"Typical" helper = the mean U8 helper: runs = (U8 helper run histogram)/504, width 3778/504 = 7.50.
Usage: python3 -I price_removals.py [u8_census.json]   (about a minute; certified on the 1e-24 grid)
"""
import sys, os, json
sys.path.insert(0, os.environ.get('PRICER_DIR', '.'))
from pricer import price_b, kappa_from, C_PR589, dec, b_fast
from fractions import Fraction as F
from collections import Counter
cen = json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'u8_census.json')))
P942 = {1: 848490, 2: 112554, 3: 148518, 4: 53946, 5: 23976, 8: 2664, 48: 11086}; M, W, CALLS = 73, 33152, 333
B942 = F(1664727258740678855772, 10**24); K942 = F(1661960547713947, 10**18)
out = {}
b0 = price_b(P942, M, W); k0 = kappa_from(b0, C_PR589)[0]
assert sum(r*n for r, n in P942.items()) == M*W - 11840
out['control'] = dict(b=dec(b0, 25), kappa=dec(k0, 19), b_matches_942=(b0 == B942), kappa_matches_942=(k0 == K942))
print('control', out['control'])
HR = {int(k): x for k, x in cen['helper_runs_total'].items()}             # per call, all 504 helpers
RW = {int(w): {int(k): x for k, x in d.items()} for w, d in cen['runs_per_helper_by_width'].items()}
WI = {int(w): x for w, x in cen['helper_width_inventory'].items()}
def apply(delta_runs_global, freed, holes_rank=None, holes=0, mode='optimistic'):
    G = Counter(P942)
    for r, n in delta_runs_global.items(): G[r] -= n
    if mode == 'optimistic':   # re-pack relaxation: hosted width fills the idle bins smallest-capacity first (whole bins)
        for r in (48, 1): G[r] -= {48: 11086, 1: 6}[r]          # take out T37's own residual
        hosted = 37*9*3778                                       # the control's hosted width
        hosted -= freed
        res = Counter()
        for cap, nb in ((8, 672*37), (16, 448*37), (32, 448*37), (48, 448*37)):
            full = min(nb, hosted//cap); hosted -= full*cap; left = nb - full
            if left and hosted:
                res[cap-hosted] += 1; left -= 1; hosted = 0
            res[cap] += left
        assert hosted == 0
        for r, x in res.items(): G[r] += x
    else:   # holes: 'holes' children whose ranks sum to freed (ranks floor/ceil of the mean hole)
        lo = freed//holes; hi_n = freed - lo*holes; G[lo] += holes-hi_n; G[lo+1] += hi_n
    assert all(x >= 0 for x in G.values()), G
    G = {r: x for r, x in G.items() if x}
    assert sum(r*x for r, x in G.items()) == M*W - 11840
    return G
def priced(G):
    b = price_b(G, M, W); k = kappa_from(b, C_PR589)[0]
    return b, k
rows = []
def row(name, G):
    b, k = priced(G)
    r = dict(case=name, b=dec(b, 13), b_gain_pct=float((b/b0-1)*100), kappa=dec(k, 16), kappa_gain_pct=float((k/k0-1)*100),
             binding='bit (kappa capped by c)' if b > C_PR589 else 'complex')
    rows.append(r); print(json.dumps(r)); return b
res = {}
for N in (0, 1, 10, 28, 50, 100, 112, 168, 200, 252, 392, 504):
    # typical: N/504 of the helper run histogram, in all 333 calls (global counts rounded to integers)
    d = {r: round(F(CALLS*N*x, 504)) for r, x in HR.items()}
    freed = sum(r*n for r, n in d.items())
    for mode in ('optimistic', 'hole'):
        if N == 0 and mode == "hole": continue
        G = apply(d, freed, holes=CALLS*N, mode=mode)
        res[(N, mode)] = row(f'typical x{N} ({mode})', G)
# by width class: one helper of each width, all calls (mean runs of that width class)
for w in sorted(WI):
    d = {r: round(F(CALLS*x, WI[w])) for r, x in RW[w].items()}; freed = sum(r*n for r, n in d.items())
    for mode in ('optimistic', 'hole'):
        row(f'one width-{w} helper ({mode})', apply(d, freed, holes=CALLS, mode=mode))
# the 112 non-plane helpers (entering 1->9 w8 x56, 4->9 w5 x32, 5->9 w4 x24): what-if all removed
out['rows'] = rows
# --- kappa: what the bit side must reach.  kappa is set by min(b, c); c = PR589's 1.7334270e-3.
cap = kappa_from(F(2, 1000), C_PR589)[0]
out['kappa_cap_with_PR589'] = dec(cap, 16)
# bit side: PR589 profile (m-anatomy levers), with the inherited surcharge; 288 calls, helper width 8
mb, Wb = 64, 28416
HB = {1: 381024, 2: 400416, 3: 22752, 4: 26880, 5: 29376, 6: 17376, 7: 26784, 42: 288}
BLOC = {1: 1323, 2: 1342, 3: 79, 4: 44, 5: 102, 6: 11, 7: 19}
c0 = price_b(HB, mb, Wb, bad=True); assert c0 == C_PR589, (c0, C_PR589)
# a typical bit helper: width 8, runs in proportion to the bit local profile (data + helpers; the bit helper runs are
# not separated here), scaled to mass 8.  Freed width paid at rank 42 (optimistic) - the bit idle cap.
bm = sum(r*x for r, x in BLOC.items())
def c_after(Nb, runs_mode='typical'):
    G = Counter(HB)
    if runs_mode == 'typical':
        d = {r: round(F(288*Nb*8*x, bm)) for r, x in BLOC.items()}
    else:
        d = {1: 288*Nb*8}
    for r, n in d.items(): G[r] -= n
    freed = sum(r*n for r, n in d.items()); q, rem = divmod(freed, 42); G[42] += q
    if rem: G[rem] += 1
    assert all(x >= 0 for x in G.values())
    return b_fast({r: x for r, x in G.items() if x}, mb, Wb, bad=True)
need = []
for N in (50, 100, 200, 392, 504):
    bN = res[(N, 'optimistic')]
    if bN <= C_PR589:
        need.append(dict(complex_removed=N, b=dec(bN, 13), bit_needed='none (complex still binds)')); continue
    found = {}
    for mode in ('typical', 'rank1'):
        lo, hi = 0, (439 if mode == 'typical' else 381024//(288*8))
        if c_after(hi, mode) < bN: found[mode] = f'> {hi} (model limit)'; continue
        while hi - lo > 1:
            mid = (lo+hi)//2
            if c_after(mid, mode) >= bN: hi = mid
            else: lo = mid
        found[mode] = hi
    k = kappa_from(bN, bN*(1+F(1, 10**9)))[0]
    need.append(dict(complex_removed=N, b=dec(bN, 13), bit_helpers_to_remove_per_call=found, kappa_if_both=dec(k, 16)))
    print(json.dumps(need[-1]))
out['bit_side_needed'] = need
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'price_removals.json'), 'w'), indent=1)
