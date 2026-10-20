#!/usr/bin/env python3
"""lever-reprice: head test on a hub owner FILL (staged544 format: rows [time, helper, owner_replica, child, port, offset]).

Own code. Semantics follow rl-check/tools/chrono.py (the T74 head test), adapted to the FILL format that the hub's
staged544 packages use (window/family semantics as read in the hub's owner/verify.py, re-implemented, not imported):
  - calls: ORDER[t] = (stage, replica); each call runs that replica's template block for that stage (its local sweep,
    then the bridges / idle sweeps / terminal events up to the next local sweep); the two global events come last.
  - family of (owner o, child c, port p) = o*896 + (2*pair(c)+side(c))*112 + p (pairs/sides from PLAN).
  - touches of a family: local sweep of its (replica, pair); bridges with its pair as src or dst; idle sweeps of its
    (replica, pair, side) for children with a remainder; terminal_Y_translations / terminal (every pair);
    hosted helper roles on it (at the hosting call's local-sweep position); translate_every_live_family (free tail).
  - remainder of bin (o, c, p) = rank(c) - cells used by hosted roles; widths from the final canonical provider words
    (dim helper_final - dim start_extra over GF(2), per helper in helper_origin order).
Classes, exactly as rl-check: HEAD-literal (nothing touches the family before the bin's idle sweep), HEAD-net(argued)
(only hosted roles / idle sweeps before), TAIL-literal (only free translations after), neither.
Memory: keeps first/last touch positions per family only.
args: PKGDIR [--control-bridge-first] [--fill FILE] [--T 360]
"""
import json, gzip, sys
from collections import Counter, defaultdict
pkg = sys.argv[1].rstrip('/')
CONTROL = '--control-bridge-first' in sys.argv
fillp = sys.argv[sys.argv.index('--fill') + 1] if '--fill' in sys.argv else pkg + '/FILL.json.gz'
rd = lambda p: json.loads(gzip.decompress(open(p, 'rb').read())) if p.endswith('.gz') else json.load(open(p))
inp = pkg + '/inputs/'
order = [tuple(x) for x in rd(inp + 'ORDER.json')]
T = 1 + max(r for s, r in order); V = 112
disp = {(s, r): tag for s, r, tag in rd(inp + 'DISPATCH.json')}
tmpl = rd(inp + 'TEMPLATE.json.gz'); plan = rd(inp + 'PLAN.json.gz')
src = tmpl['execution_schedule']
if CONTROL:   # negative control (as rl-check): move the bridge reading pair 3 before pair 3's first idle sweep
    i14 = next(j for j, e in enumerate(src) if e['kind'] == 'idle_cover_sweep' and e['pair'] == 3 and e['from_level'] == -1)
    i15 = next(j for j, e in enumerate(src) if e['kind'] == 'bridge' and e['src'] == 3)
    assert i15 == i14 + 1
    src = src[:i14] + [src[i15], src[i14]] + src[i15 + 1:]
prank = {}
for c in tmpl['idle_children']: prank[c['child']] = c['rank']
pair_of = {int(k): v for k, v in plan['pairs'].items()}; side_of = {int(k): v for k, v in plan['sides'].items()}
stage_pair = {e['stage']: e['pair'] for e in tmpl['recipe'] if e['kind'] == 'local'}
blocks, globals_, stage = {}, [], None
for j, e in enumerate(src):
    if e['kind'] in ('bank_padding_cover_sweep', 'translate_every_live_family'): globals_.append((j, e)); continue
    if e['kind'] == 'local_cover_sweep': stage = e['stage']; blocks[stage] = []
    blocks[stage].append((j, e))
clock = {c: t for t, c in enumerate(order)}
assert len(clock) == len(order) == 9 * T
assert all(clock[s, r] < clock[s + 1, r] for r in range(T) for s in range(8)), 'per-replica chronology'
# positions are (t, k) with k the index inside the block; globals at (10**9, j)
INF = 10 ** 9
# widths
def dim(rows):
    piv = {}
    for v in rows:
        while v:
            l = v.bit_length()
            if l not in piv: piv[l] = v; break
            v ^= piv[l]
    return len(piv)
widths = {}
for tag, pid in (('base', 'gauge44'), ('extra', 'extra')):
    w = rd(inp + f'WORD-{tag}.json.gz'); nat = rd(inp + f'old-NATIVE-{tag}.json.gz')
    origin = nat['helper_origin']; full = [1 << i for i in range(9)]
    widths[pid] = [dim(w['helper_final'].get(str(r), full)) - dim(w['start_extra'].get(str(r), [])) for r in origin]
# per (replica, pair): positions of local/bridge/terminal touches; per (replica, pair, side): idle sweeps
first_hard = {}   # (rep, pair) -> first pos of local/bridge/terminal touch
last_hard = {}    # (rep, pair) -> last pos of a non-free touch (local/bridge/terminal)
sweeps = defaultdict(list)  # (rep, pair, side) -> [(pos, child)]
sweep_pos = {}    # (rep, child) -> pos
callpos = {}
def upd(dmin, dmax, key, pos):
    if key not in dmin or pos < dmin[key]: dmin[key] = pos
    if dmax is not None and (key not in dmax or pos > dmax[key]): dmax[key] = pos
for t, (s, r) in enumerate(order):
    for k, (j, e) in enumerate(blocks[s]):
        pos = (t, k); kd = e['kind']
        if kd == 'local_cover_sweep':
            callpos[t] = pos; upd(first_hard, last_hard, (r, stage_pair[s]), pos)
        elif kd == 'bridge':
            upd(first_hard, last_hard, (r, e['src']), pos); upd(first_hard, last_hard, (r, e['dst']), pos)
        elif kd == 'idle_cover_sweep':
            for sd, i in enumerate(e['child_indices']):
                if i not in prank: continue
                assert pair_of[i] == e['pair'] and side_of[i] == sd
                sweeps[r, e['pair'], sd].append((pos, i)); sweep_pos[r, i] = pos
        elif kd in ('terminal_Y_translations', 'terminal'):
            for p in range(4): upd(first_hard, last_hard, (r, p), pos)
        else: raise SystemExit('unknown kind ' + kd)
# hosted roles
fill = rd(fillp); rows = fill['rows']
assert fill['columns'] == ['time', 'helper', 'owner_replica', 'child', 'port', 'offset']
H = len(rows) // len(order); assert H * len(order) == len(rows)
used = Counter(); first_role = {}
for jrow, (t, h, o, c, p, z) in enumerate(rows):
    assert t == jrow // H and h == jrow % H
    wdt = widths[disp[order[t]]][h]
    assert z == used[o, c, p], 'compact prefix'
    used[o, c, p] = z + wdt
    assert used[o, c, p] <= prank[c], 'capacity'
    f = (o, pair_of[c], side_of[c], p); pos = callpos[t]
    if f not in first_role or pos < first_role[f]: first_role[f] = pos
    # windows: role must be strictly before the bin's idle sweep (the hub's 'deadline')
    assert pos < sweep_pos[o, c], 'role after deadline'
del rows, fill
cls = Counter(); idle = Counter()
for o in range(T):
    for c, cap in prank.items():
        pr, sd = pair_of[c], side_of[c]; d = sweep_pos[o, c]
        others_before = [ps for ps, i in sweeps[o, pr, sd] if i != c and ps < d]
        for p in range(V):
            rem = cap - used.get((o, c, p), 0)
            if not rem: continue
            idle[rem] += 1
            f = (o, pr, sd, p)
            fh = first_hard.get((o, pr), (INF, 0)); fr = first_role.get(f, (INF, 0))
            fo = min(others_before) if others_before else (INF, 0)
            lit_head = fh > d and fr > d and fo > d
            net_head = fh > d
            lh = last_hard.get((o, pr), (-1, 0))
            lit_tail = lh < d   # only free translations (globals) after
            k = 'HEAD-literal' if lit_head else ('HEAD-net(argued)' if net_head else ('TAIL-literal' if lit_tail else 'neither'))
            cls[c, cap, rem, pr, sd, k] += 1
loc = {}
for tag, pid in (('base', 'gauge44'), ('extra', 'extra')):
    loc[pid] = {int(r): n for r, n in rd(inp + f'WORD-{tag}.json.gz')['report']['local_profile'].items()}
prof = Counter(idle)
for call in order: prof.update(loc[disp[call]])
out = dict(pkg=pkg, fill=fillp, control_bridge_first=CONTROL, T=T, calls=len(order), roles_per_call=H,
           idle_profile=dict(sorted(idle.items())), full_profile=dict(sorted(prof.items())), rank_mass=sum(r*n for r, n in prof.items()),
           classes=[dict(child=k[0], child_rank=k[1], remainder=k[2], pair=k[3], side=k[4], cls=k[5], count=n) for k, n in sorted(cls.items())])
tot = defaultdict(Counter)
for k, n in cls.items(): tot[k[5]][k[2]] += n
out['remainder_rank_by_class'] = {c: dict(sorted(v.items())) for c, v in sorted(tot.items())}
out['totals_by_class'] = {c: sum(v.values()) for c, v in sorted(tot.items())}
print(json.dumps(out, indent=1))
