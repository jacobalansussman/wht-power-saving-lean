#!/usr/bin/env python3
"""hx-price (wht45/hubx, 20 Oct 2026): independent literal-head recount and paid profile for a #2979-line FILL.

OWN code, written for this check; it imports nothing from heads_fill.py, recount.py, evalfill.py or fc-reopt.
Definition (reopt/BRIEF.md, rl-check semantics): a literal head is a rank-48 bin (owner o, child c, port p) of child 8
or 10 that is unhosted (no FILL row names it), whose idle sweep is the first operation on its own array (the family
(o, pair(c), side(c), p)), acting on that array alone, at most one per array.

Differences from heads_fill.py by design:
  - pairs and sides are derived from the TEMPLATE's idle sweeps (child_indices = [side 0, side 1]) and only
    cross-checked against PLAN.json;
  - role widths come from the hub author's published allocation/PROVIDERS.json (per provider tag), cross-checked by an
    own GF(2) rank of WORD helper_final / start_extra (column-pivot elimination over sets of bits);
  - chronology: one global integer sequence number per executed event, families keyed by tuples, first touch kept
    with the event that made it.
usage: hx_heads.py --fill FILL.json.gz --order ORDER.json --dispatch DISPATCH.json --template TEMPLATE.json.gz
                   --plan PLAN.json.gz --providers PROVIDERS.json --word-base WORD-base.json.gz --word-extra WORD-extra.json.gz
                   [--control-bridge-first] [--profile PROFILE.json]
"""
import argparse, gzip, json, sys
from collections import Counter, defaultdict

ap = argparse.ArgumentParser()
for k in ('fill', 'order', 'dispatch', 'template', 'plan', 'providers', 'word-base', 'word-extra'):
    ap.add_argument('--' + k, required=True)
ap.add_argument('--profile'); ap.add_argument('--control-bridge-first', action='store_true')
a = ap.parse_args()


def load(p):
    raw = open(p, 'rb').read()
    return json.loads(gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw)


def need(c, msg):
    if not c:
        raise SystemExit('FAIL: ' + msg)


order = [tuple(x) for x in load(a.order)]
dispatch = {(s, r): tag for s, r, tag in load(a.dispatch)}
tmpl = load(a.template); plan = load(a.plan); prov = load(a.providers)
T = len({r for s, r in order}); S = len({s for s, r in order})
need(len(order) == S * T and len(set(order)) == len(order), 'ORDER is not every (stage, replica) once')
need(set(dispatch) == set(order), 'DISPATCH covers ORDER')
when = {c: i for i, c in enumerate(order)}
need(all(when[s, r] < when[s + 1, r] for r in range(T) for s in range(S - 1)), 'stages in order per replica')

# --- template: ranks, pairs and sides (from the idle sweeps), stage -> pair of local sweeps
rank = {}
for e in tmpl['idle_children']:
    need(rank.setdefault(e['child'], e['rank']) == e['rank'], 'one rank per child')
ports = sorted({e['port'] for e in tmpl['idle_children']}); V = len(ports)
need(ports == list(range(V)), 'ports 0..V-1')
sched = list(tmpl['execution_schedule'])
if a.control_bridge_first:
    i = next(j for j, e in enumerate(sched) if e['kind'] == 'idle_cover_sweep' and 8 in e['child_indices'])
    need(sched[i + 1]['kind'] == 'bridge' and sched[i + 1]['src'] == sched[i]['pair'], 'control: bridge after the sweep')
    sched[i], sched[i + 1] = sched[i + 1], sched[i]
pair, side = {}, {}
for e in sched:
    if e['kind'] == 'idle_cover_sweep':
        for sd, c in enumerate(e['child_indices']):
            need(c not in pair, 'child swept twice'); pair[c] = e['pair']; side[c] = sd
need(all(pair[c] == plan['pairs'][str(c)] and side[c] == plan['sides'][str(c)] for c in rank), 'PLAN pairs/sides agree')
local_pair = {e['stage']: e['pair'] for e in tmpl['recipe'] if e['kind'] == 'local'}
# blocks: each event belongs to the stage of the most recent local sweep; the two global events run once at the end
GLOBAL = ('bank_padding_cover_sweep', 'translate_every_live_family')
block = defaultdict(list); cur = None
for e in sched:
    if e['kind'] in GLOBAL: continue
    if e['kind'] == 'local_cover_sweep': cur = e['stage']
    block[cur].append(e)
need(sorted(block) == list(range(S)), 'one block per stage')

# --- widths: published PROVIDERS.json, cross-checked by own GF(2) ranks of the words
def gf2rank(vectors):
    rows = [int(v) for v in vectors]; r = 0
    bits = set()
    for v in rows: bits |= {i for i in range(v.bit_length()) if v >> i & 1}
    for b in sorted(bits, reverse=True):
        piv = next((k for k in range(r, len(rows)) if rows[k] >> b & 1), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for k in range(len(rows)):
            if k != r and rows[k] >> b & 1: rows[k] ^= rows[r]
        r += 1
    return r


width = {tag: list(v['widths']) for tag, v in prov.items()}
H = len(next(iter(width.values())))
gf2 = {}
for tag, path in (('gauge44', a.word_base), ('extra', a.word_extra)):
    w = load(path); removed = set(w['removed']); R = w['report']['R']
    # helper identities 224..703 minus the removed ones (the word's 480 helpers in increasing order; R = 480)
    ids = [h for h in range(224, 224 + R) if h not in removed]
    need(len(ids) == H, 'helper count %s' % tag)
    full = [1 << i for i in range(w['report']['h'])]
    gf2[tag] = [gf2rank(w['helper_final'].get(str(h), full)) - gf2rank(w['start_extra'].get(str(h), [])) for h in ids]
    need(Counter(gf2[tag]) == Counter({int(k): v for k, v in w['report']['width_inventory'].items()}), 'width inventory %s' % tag)
widths_agree = {tag: gf2[tag] == width[tag] for tag in width}
need(all(widths_agree.values()), 'PROVIDERS widths == own GF(2) widths')

# --- fill rows
fill = load(a.fill)
need(fill['columns'] == ['time', 'helper', 'owner_replica', 'child', 'port', 'offset'], 'fill columns')
rows = fill['rows']
need(len(rows) == len(order) * H, 'one row per (call, helper)')
used = Counter(); roles_on = defaultdict(list)
for i, (t, h, o, c, p, z) in enumerate(rows):
    need((t, h) == divmod(i, H), 'canonical row order')
    need(c in rank and 0 <= p < V and 0 <= o < T, 'role bounds')
    need(z == used[o, c, p], 'compact prefix (offset = cells used so far)')
    used[o, c, p] += width[dispatch[order[t]]][h]
    need(used[o, c, p] <= rank[c], 'capacity')
    roles_on[o, pair[c], side[c], p].append(t)

# --- chronology: global sequence numbers
first = {}                      # family (o, pair, side, port) -> (seq, what)
sweep_seq = {}                  # (o, c) -> seq of the idle sweep of child c on replica o
local_seq = {}                  # call t -> seq of its local sweep
seq = 0
def touch_pair(r, pr, what):
    for sd in (0, 1):
        for p in range(V):
            first.setdefault((r, pr, sd, p), (seq, what))
for t, (s, r) in enumerate(order):
    for e in block[s]:
        seq += 1; k = e['kind']
        if k == 'local_cover_sweep':
            local_seq[t] = seq; touch_pair(r, local_pair[s], 'local')
        elif k == 'bridge':
            touch_pair(r, e['src'], 'bridge'); touch_pair(r, e['dst'], 'bridge')
        elif k in ('terminal_Y_translations', 'terminal'):
            for pr in range(4): touch_pair(r, pr, k)
        elif k == 'idle_cover_sweep':
            for sd, c in enumerate(e['child_indices']):
                if c not in rank: continue
                sweep_seq[r, c] = seq
                for p in range(V): first.setdefault((r, e['pair'], sd, p), (seq, ('sweep', c)))
        else: need(False, 'unknown event ' + k)
# hosted roles touch the owner's family at the hosting call's local sweep; they must precede the bin's own sweep
first_role = {}
for f, ts in roles_on.items():
    first_role[f] = min(local_seq[t] for t in ts)
for i, (t, h, o, c, p, z) in enumerate(rows):
    need(local_seq[t] < sweep_seq[o, c], 'role before its bin\'s idle sweep (deadline)')
del rows

heads = Counter(); idle = Counter(); fams_with_head = Counter(); nonhead48 = Counter()
for o in range(T):
    for c in rank:
        for p in range(V):
            rem = rank[c] - used.get((o, c, p), 0)
            if rem: idle[rem] += 1
            if rank[c] != 48 or rem != 48: continue
            f = (o, pair[c], side[c], p)
            fs, what = first[f]
            fr = first_role.get(f, 10 ** 12)
            if what == ('sweep', c) and fs == sweep_seq[o, c] and fr > fs:
                heads[c] += 1; fams_with_head[f] += 1
            else:
                nonhead48[c] += 1
need(all(v == 1 for v in fams_with_head.values()), 'at most one head per array')
loc = {tag: {int(k): v for k, v in prov[tag]['actual_paid_profile'].items()} for tag in prov}
profile = Counter(idle)
for call in order: profile.update(loc[dispatch[call]])
out = dict(T=T, calls=len(order), roles_per_call=H, roles=len(order) * H, W=len({(o, pr, sd, p) for (o, pr, sd, p) in first}),
           control_bridge_first=a.control_bridge_first, widths_providers_equal_own_gf2=widths_agree,
           literal_heads=sum(heads.values()), heads_by_child=dict(sorted(heads.items())),
           unhosted_rank48_not_head=dict(sorted(nonhead48.items())),
           idle_profile=dict(sorted(idle.items())), paid_profile=dict(sorted(profile.items())),
           rank_mass=sum(r * n for r, n in profile.items()))
if a.profile:
    claimed = {int(k): v for k, v in load(a.profile).items()}
    out['profile_equals_PROFILE_json'] = claimed == dict(profile)
print(json.dumps(out, indent=1))
