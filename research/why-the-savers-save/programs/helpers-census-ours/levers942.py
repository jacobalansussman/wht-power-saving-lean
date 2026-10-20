"""u-helpers-ours: helper levers priced on #942 (U8 word x 333 calls) with the m-anatomy pricer (hub engine re-implementation).
Every figure is a WHAT-IF (profile change only, no construction) unless it is the control.
Usage: python3 -I levers942.py census942.json levers942.json"""
import sys, json, os, statistics as st
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
from census942 import price_delta, pct, B942, LOC
C = json.load(open(sys.argv[1])); regs = {int(q): v for q, v in C['registers'].items()}
H = {q: v for q, v in regs.items() if v['cls'] in ('factory', 'face', 'M')}
out = {'control_b': str(price_delta({}, certify=True)), 'control_equals_942': price_delta({}, certify=True) == B942}
def rm(qs):
    d = Counter()
    for q in qs:
        for r in H[q]['runs']: d[r] -= 1
    return dict(d)
# 1. one helper removed (every call): exact price per helper, summarised by class and width
per = {}
for q in H:
    b = price_delta(rm([q])); per[q] = pct(b)
tab = {}
for q, p in per.items():
    k = (H[q]['cls'], sum(H[q]['runs'])); tab.setdefault(k, []).append(p)
out['one_helper_removed_pct_b'] = {f'{c} w{w} (n={len(v)})': dict(min=round(min(v), 4), median=round(st.median(v), 4), max=round(max(v), 4)) for (c, w), v in sorted(tab.items())}
out['one_helper_removed_all_mean_pct_b'] = round(st.mean(per.values()), 4)
# 2. classes removed / counts down to the floors (what-ifs; the helpers removed are the ones a design would most likely lose)
order_cheap = sorted(H, key=lambda q: per[q])          # least valuable first
order_dear = sorted(H, key=lambda q: -per[q])          # most valuable first
M = [q for q in H if H[q]['cls'] == 'M']; FACE = [q for q in H if H[q]['cls'] == 'face']; FAC = [q for q in H if H[q]['cls'] == 'factory']
cases = {
 'remove all 112 M copies (504 -> 392 incl. 3 centre-hosts; packet-plane bound with second copies)': M,
 'remove all 168 face copies': FACE,
 'remove M + face (-> 224-ish: packet-plane one-copy bound)': M + FACE,
 'remove 112 helpers, most valuable first': order_dear[:112],
 'remove 112 helpers, least valuable first': order_cheap[:112],
 'down to 280 (targets visit planes, two copies; toy V)': order_cheap[:len(H) - 280 + 3],
 'down to 112 (targets visit planes, one copy; toy V)': order_cheap[:len(H) - 112 + 3],
 'down to 107 (hub sink-only floor)': order_cheap[:len(H) - 107 + 3],
 'all helpers removed (strict floor 0)': list(H),
}
out['count_whatifs'] = {}
for name, qs in cases.items():
    b = price_delta(rm(qs)); out['count_whatifs'][name] = dict(removed=len(qs), b=f'{float(b):.7e}', pct=round(pct(b), 3))
# 3. stops (runs) instead of count: same helpers, fewer runs
def merged(runs, k):
    r = sorted(runs)
    while len(r) > k: a = r.pop(0); b = r.pop(0); r.append(a + b); r.sort()
    return r
def stops_case(k, classes=('factory', 'face', 'M')):
    d = Counter()
    for q, v in H.items():
        if v['cls'] in classes and len(v['runs']) > k:
            for r in v['runs']: d[r] -= 1
            for r in merged(v['runs'], k): d[r] += 1
    return dict(d)
out['stop_whatifs'] = {}
for k in (1, 2, 3):
    b = price_delta(stops_case(k)); out['stop_whatifs'][f'every helper in at most {k} runs (count 504 kept)'] = dict(b=f'{float(b):.7e}', pct=round(pct(b), 3))
d1 = Counter()
for q, v in H.items():
    if len(v['runs']) > 1:
        r = sorted(v['runs']); a, b2 = r[0], r[1]; d1[a] -= 1; d1[b2] -= 1; d1[a + b2] += 1
b = price_delta(dict(d1)); out['stop_whatifs']['one stop fewer on every helper (two smallest runs merged)'] = dict(b=f'{float(b):.7e}', pct=round(pct(b), 3))
# 4. ledger: shares of the gross saving per replica (8*112 = 896 moves), via the F-shares at the root
cl = C['classes']; helpersF = sum(cl[k]['F_share_pct'] for k in ('factory', 'face', 'M'))
dataF = cl['data x']['F_share_pct'] + cl['data y']['F_share_pct']
out['ledger_per_replica'] = dict(gross=896, copies=576, deficit=320, copies_pct_gross=round(100 * 576 / 896, 1),
    helpers_pct_D=round(helpersF, 1), helpers_pct_gross=round(helpersF * 320 / 896, 1), data_pct_gross=round(dataF * 320 / 896, 1),
    centres_pct_gross=round(cl['centres']['F_share_pct'] * 320 / 896, 2), idle_pct_gross=round(C['idle_F_share_pct'] * 320 / 896, 2))
json.dump(out, open(sys.argv[2], 'w'), indent=1); print(json.dumps(out, indent=1))
