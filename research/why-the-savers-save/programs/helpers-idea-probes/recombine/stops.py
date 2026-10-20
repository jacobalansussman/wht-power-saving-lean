#!/usr/bin/env python3
"""g-recombine: why does each helper run stop where it does?  Replays the U8 CANONICAL-WORD with u8_census.py's climb
rule (same profile) and records, for every paid climb of a helper, the phase, the partner register and its class, and
whether the helper was pulled up to the partner's frame, pulled the partner up, or both moved to a larger gate frame.
Also counts climbs made only by the final frame (exit-frame climbs).
Usage: python3 -I stops.py CANONICAL-WORD.json.gz [out.json]"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return bb
B = [basis(f) for f in F]; dims = [len(b) for b in B]
def inside(a, b):
    bb = B[b]
    for x in B[a]:
        for z in bb: x = min(x, x ^ z)
        if x: return False
    return True
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
cur = list(c['start']); log = defaultdict(list)
def climb(r, f, ph, partner, pf):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not inside(cur[r], f): return
    if partner is None: how = 'final'
    elif pf == f: how = 'to-partner'
    else:
        jd = len(basis(B[cur[r]] + B[pf]))
        how = 'join' if jd == dims[f] else 'above-join'
    log[r].append(dict(ph=ph, r=dims[f]-dims[cur[r]], fr=dims[cur[r]], to=dims[f], partner=None if partner is None else cls(partner), how=how))
    cur[r] = f
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f, ph, r, f)
        continue
    for e in c.get(ph) or []:
        if e[0] == 'add':
            a, b, f = e[1], e[2], e[-1]
            fa, fb = cur[a], cur[b]
            climb(a, f, ph, b, fb); climb(b, f, ph, a, fa)
        else:
            climb(e[1], e[-1], ph, e[1], e[-1])
for r in range(n): climb(r, c['final'][r], 'final', None, None)
H = range(2*v, n)
out = {}
out['helper_runs_by_how'] = Counter(x['how'] for r in H for x in log[r])
out['helper_mass_by_how'] = Counter()
for r in H:
    for x in log[r]: out['helper_mass_by_how'][x['how']] += x['r']
out['helper_runs_by_phase'] = Counter(x['ph'] for r in H for x in log[r])
out['helper_rank1_runs_by_partner'] = Counter((x['partner'], x['how']) for r in H for x in log[r] if x['r'] == 1)
out['helper_runs_by_partner'] = Counter((x['partner'], x['r']) for r in H for x in log[r])
out['final_climbs'] = Counter((x['fr'], x['to']) for r in H for x in log[r] if x['how'] == 'final')
for cl, rr in (('x', range(v)), ('y', range(v, 2*v))):
    out[cl+'_runs_by_how'] = Counter(x['how'] for r in rr for x in log[r])
    out[cl+'_final_climbs'] = Counter((x['fr'], x['to']) for r in rr for x in log[r] if x['how'] == 'final')
# chain shapes of helpers: tuple of dims visited
out['helper_chain_shapes'] = Counter(tuple([dims[c['start'][r]]] + [x['to'] for x in log[r]]) for r in H).most_common(25)
o = {k: ({str(a): b for a, b in val.items()} if isinstance(val, Counter) else [[list(a), b] for a, b in val]) for k, val in out.items()}
print(json.dumps(o, indent=1))
if len(sys.argv) > 2: json.dump(o, open(sys.argv[2], 'w'), indent=1)
# --- detail: helper runs by (phase, partner class, how, rank, from->to)
det = Counter((x['ph'], x['partner'], x['how'], x['fr'], x['to']) for r in H for x in log[r])
print('\nDETAIL phase partner how from to : count')
for k, x in sorted(det.items(), key=lambda kv: -kv[1])[:40]: print(k, x)
