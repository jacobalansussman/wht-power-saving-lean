#!/usr/bin/env python3
"""u-helpers-hub: census of the 504 helpers of the U8 local word (#869/#870 lineage, used unchanged by #890 and #942).

Replays the U8 CANONICAL-WORD under the climb rule that reproduces the hub's local profile exactly
(phases A, centres, precondition, cut, compensation, B, tail, cleanup; reads climb too; one paid child of rank
dim(f)-dim(cur) per frame change; then every register climbs to its final frame), as found by
bound-check/k-word/tools/hub_climb_model.py.  Frames are F2-subspaces of F2^9 given by integer bases.
Usage: python3 -I u8_census.py CANONICAL-WORD.json.gz [out.json]
"""
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
cur = list(c['start']); runs = defaultdict(list); visits = defaultdict(list)  # register -> list of frames visited
for r in range(n): visits[r].append(cur[r])
def climb(r, f, ph):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not inside(cur[r], f): return
    runs[r].append(dims[f]-dims[cur[r]]); cur[r] = f; visits[r].append(f)
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f, ph)
        continue
    for e in c.get(ph) or []:
        climb(e[1], e[-1], ph)
        if e[0] == 'add': climb(e[2], e[-1], ph)
for r in range(n): climb(r, c['final'][r], 'final')
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
prof = Counter(); byc = defaultdict(Counter)
for r in range(n):
    for k in runs[r]: prof[k] += 1; byc[cls(r)][k] += 1
hub = {int(k): x for k, x in c['local_profile'].items()}
cen = set(c['centers'])
# centres: their rank-8 climb is the copy charge (cst 64); the hub profile lists it as 8 children of rank 8
width = {r: sum(runs[r]) for r in range(2*v, n)}
start_dim = {r: dims[c['start'][r]] for r in range(n)}
final_dim = {r: dims[c['final'][r]] for r in range(n)}
winv = Counter(width.values())
# helper classes by entrance (start frame dim) and end (final frame dim)
hcls = Counter((start_dim[r], final_dim[r], width[r]) for r in range(2*v, n))
# rank-3 frames visited, by class of visitor
r3 = defaultdict(lambda: Counter())
for r in range(n):
    for f in visits[r]:
        if dims[f] == 3: r3[tuple(sorted(B[f]))][cls(r)] += 1
nplanes = len(r3)
vis_hist = Counter(sum(cn.values()) for cn in r3.values())
helper_runs = byc['h']
out = dict(
  registers=dict(sources=v, targets=v, helpers=R, total=n),
  replayed_profile=dict(sorted(prof.items())), hub_profile=dict(sorted(hub.items())),
  profile_match_hub_minus_centres_plus_8x8=(dict(prof + Counter({8: 8})) == hub) or dict(prof) == hub,
  runs_by_class={k: dict(sorted(x.items())) for k, x in byc.items()},
  mass_by_class={k: sum(a*b for a, b in x.items()) for k, x in byc.items()},
  centres=[[r, cls(r), runs[r]] for r in c['centers']],
  helper_width_inventory=dict(sorted(winv.items())), total_helper_width=sum(width.values()),
  helper_runs_total=dict(sorted(helper_runs.items())), helper_runs_count=sum(helper_runs.values()),
  helper_classes_startdim_finaldim_width={f'{a}->{b} w{w}': k for (a, b, w), k in sorted(hcls.items())},
  runs_per_helper_by_width={w: dict(sorted(Counter(k for r in width if width[r] == w for k in runs[r]).items())) for w in sorted(winv)},
  rank3_frames_visited=nplanes, rank3_visitors_hist=dict(sorted(vis_hist.items())),
  rank3_visits_by_class=dict(sum((cn for cn in r3.values()), Counter())),
)
print(json.dumps(out, indent=1))
if len(sys.argv) > 2: json.dump(out, open(sys.argv[2], 'w'), indent=1)
# --- second pass: which helpers visit a rank-3 (packet) plane, by class; dims visited per helper
vis3 = {r: any(dims[f] == 3 for f in visits[r]) for r in range(n)}
k2 = Counter((start_dim[r], final_dim[r], width[r], vis3[r]) for r in range(2*v, n))
print('helpers (startdim, finaldim, width, visits a rank-3 plane):', dict(sorted(k2.items())))
dimpath = Counter(tuple(dims[f] for f in visits[r]) for r in range(2*v, n))
print('helper dim paths (top 15):', dimpath.most_common(15))
print('source dim paths:', Counter(tuple(dims[f] for f in visits[r]) for r in range(v)).most_common(6))
print('target dim paths:', Counter(tuple(dims[f] for f in visits[r]) for r in range(v, 2*v)).most_common(6))
# per plane: which classes visit; and the 8 sources per block
pc = Counter(tuple(sorted(cn.items())) for cn in r3.values()); print('per-plane visitor classes:', pc)
# --- third pass: run histograms of the 392 plane-visiting helpers and of the 112 non-plane (entering) helpers
np_runs = Counter(k for r in range(2*v, n) if not vis3[r] for k in runs[r])
pl_runs = Counter(k for r in range(2*v, n) if vis3[r] for k in runs[r])
extra = dict(nonplane_helpers=sum(1 for r in range(2*v, n) if not vis3[r]), nonplane_runs=dict(sorted(np_runs.items())),
             nonplane_width=sum(a*b for a, b in np_runs.items()), plane_helpers=sum(1 for r in range(2*v, n) if vis3[r]),
             plane_runs=dict(sorted(pl_runs.items())), plane_width=sum(a*b for a, b in pl_runs.items()))
print(json.dumps(extra))
if len(sys.argv) > 2:
    out.update(extra); json.dump(out, open(sys.argv[2], 'w'), indent=1)
