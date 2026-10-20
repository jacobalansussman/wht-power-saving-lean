"""aA-data (wht45/allheads, 19 Oct 2026): x-tables' own_to_host2979.py with ONE change, the FILL pin: design A's
every-head-kept fill (reopt FOUND.md, o-model force-t74, sha256 764395fc...); everything else (ORDER, DISPATCH, PROFILE,
PROVIDERS, words, TEMPLATE, PLAN pins, every check) unchanged.  Original description:
x-tables (wht45/lever2979, 19 Oct 2026): #2979's literal owner fill -> OUR cross-replica host table (first code
path).  A copy of p-own's tools/own_to_host.py (wht45/lever, the #3945 converter) with the instance changed and nothing
else: T = 74 replicas, 666 calls, R = 408 slots, two providers; the same format "p-own host table v1x", the same
window rule, the same hub-rule checks on the raw FILL.  What differs from the #3945 converter, and why:
  * #2979 ships no native word, so helper_origin comes from x-gate's native form (p-lower's transport); it is checked
    to be range(224, 704) minus the word's 72 removed helpers (the hub's borrow numbering), as for #3945;
  * the expected inventory, local profiles, dispatch counts and remainders are NOT literals of #3945: they are read
    from the hub's own #2979 files (PROVIDERS.json widths and actual_paid_profile, DISPATCH.json, PROFILE.json), and
    the table must reproduce them exactly (full profile == PROFILE.json);
  * TEMPLATE and PLAN are #3945's copies, sha-checked equal to #2979's vendor pins (21297b20, d566cfb8);
  * g-heads' G0 file is optional (x-gate's G0 of #2979, when it exists, is compared on cross-replica roles and the
    remainder profile).
The hub's files are read as data; no hub code is run.

usage: python3 -I own_to_host2979.py ALLOC_DIR DATA_DIR VENDOR_DIR RESIDUALS.json OUT_PREFIX [G0.json]
  ALLOC_DIR  = #2979 research/rank7-joint-endframes/allocation (FILL.json.gz, PROFILE.json, ORDER.json,
               DISPATCH.json, PROVIDERS.json);  DATA_DIR = .../data (WORD-base.json.gz, WORD-extra.json.gz);
  VENDOR_DIR = a directory with TEMPLATE.json.gz (21297b20...) and PLAN.json.gz (d566cfb8...);
  RESIDUALS  = x-tables in/residuals-2979.json (mk_resid2979.py, from x-gate's gcert/1w words)."""
import sys, json, gzip, hashlib, collections, time as _time
t0 = _time.time()
alloc, datad, vend, residf, outp = sys.argv[1:6]
g0f = sys.argv[6] if len(sys.argv) > 6 else None
CH = ['A2', 'B2', 'A3', 'B3', 'A4', 'B4', 'F1', 'F2', 'F3']
RANK = [16, 8, 32, 8, 48, 8, 48, 32, 16]
KZ = [0, 2, 0, 4, 0, 6, 3, 5, 7]; KK = [2, 1, 4, 1, 6, 1, 6, 4, 2]
KIND = {(1, 16): 0, (1, 8): 1, (2, 32): 2, (2, 8): 3, (3, 48): 4, (3, 8): 5, (0, 48): 6, (1, 32): 7, (2, 16): 8}
HPAIR = [1, 1, 2, 2, 3, 3, 0, 1, 2]
PROV = {'extra': 0, 'gauge44': 1}; PNAME = {0: 'extra', 1: 'gauge44'}; TAG = {0: 'extra', 1: 'base'}
T, R, V = 74, 408, 112; NC = 9 * T; NH = 2016 * T
def rd(p):
    b = open(p, 'rb').read(); return json.loads(gzip.decompress(b) if p.endswith('.gz') else b)
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def need(c, why):
    if not c: raise SystemExit('REFUSED: ' + why)
files = {'FILL': alloc + '/FILL.json.gz', 'PROFILE': alloc + '/PROFILE.json', 'ORDER': alloc + '/ORDER.json',
         'DISPATCH': alloc + '/DISPATCH.json', 'PROVIDERS': alloc + '/PROVIDERS.json',
         'TEMPLATE': vend + '/TEMPLATE.json.gz', 'PLAN': vend + '/PLAN.json.gz',
         'WORD-base': datad + '/WORD-base.json.gz', 'WORD-extra': datad + '/WORD-extra.json.gz', 'residuals': residf}
if g0f: files['g0'] = g0f
shas = {k: sha(v) for k, v in files.items()}
PIN = {'FILL': '764395fcc2d883e9dc1ecfedce2131c745dc7ad1b9c41d7b52ef54da16d7f71c',
       'ORDER': '519d1d5dbb997a219fcfcad34c0b03e7838bef24199784d7e964fa10090f66e1',
       'DISPATCH': '65e435bd68e97243c003677271b88485a29d4d0fae2dff1354100404a145e74c',
       'PROFILE': '94a819747f3ced41f1544ada0036c983fc0c92847b9644380f877a06cd7931cf',
       'PROVIDERS': 'f2edb7922f3029a5d607437e737775ff3371bca49ebab62e5e2495d16d3bf573',
       'TEMPLATE': '21297b20f1f2dc76eb9dd8124416068c3cb096bab9df17fa4da19bc68dfe082e',
       'PLAN': 'd566cfb87b71639d4138e1f89b543e2dd6a7772070226a2a14c799a6835bf21e',
       'WORD-base': 'abaf67bfcb257528658c906b612f57025617592bd7889828f92b97a5f4f4300d',
       'WORD-extra': 'e0eaf64c4d5d324693d43889bce1b9e53e8da00402105c8c978fd102f208bed7'}
for k, v in PIN.items(): need(shas[k] == v, k + ' sha (the #2979 pin / ALLOCATION-PINS.json)')
order = rd(files['ORDER']); dispatch = rd(files['DISPATCH']); tpl = rd(files['TEMPLATE']); plan = rd(files['PLAN'])
prof_hub = {int(a): b for a, b in rd(files['PROFILE']).items()}
provs = rd(files['PROVIDERS']); resid = rd(residf); g0 = rd(g0f) if g0f else None
need(sorted(provs) == ['extra', 'gauge44'], 'two providers extra, gauge44')
ok = {}
# ---- widths: hub rule (own GF(2) rank) vs our lowering's residuals vs the hub's PROVIDERS.json widths
def rank2(rows):
    piv = {}
    for x in rows:
        while x:
            b = x.bit_length()
            if b in piv: x ^= piv[b]
            else: piv[b] = x; break
    return len(piv)
FULL = [1 << i for i in range(9)]
res, lprof = {}, {}
for p in (0, 1):
    w = rd(files['WORD-' + TAG[p]]); origin = resid[PNAME[p]]['helper_origin']
    need(origin == sorted(set(range(224, 704)) - set(w['removed'])) and len(origin) == R, 'helper origin')
    hubw = [rank2(w['helper_final'].get(str(r), FULL)) - rank2(w['start_extra'].get(str(r), [])) for r in origin]
    ours = resid[PNAME[p]]['residuals']
    need(ours == hubw, 'our residuals == hub widths, index for index (' + PNAME[p] + ')')
    need(ours == provs[PNAME[p]]['widths'], 'our residuals == PROVIDERS.json widths (' + PNAME[p] + ')')
    res[p] = ours
    lp = {int(a): b for a, b in w['report']['local_profile'].items()}
    gh = {int(a): b for a, b in resid[PNAME[p]]['block_hist'].items()}
    pp = {int(a): b for a, b in provs[PNAME[p]]['actual_paid_profile'].items()}
    need(lp == gh == pp, 'local profile == our lowered histogram == PROVIDERS.json actual_paid_profile (' + PNAME[p] + ')')
    lprof[p] = lp
ok['widths_hub_eq_ours'] = True
provider_slot_diffs = [q for q in range(R) if res[0][q] != res[1][q]]
# ---- children -> kinds
pairs = {int(a): b for a, b in plan['pairs'].items()}; sides = {int(a): b for a, b in plan['sides'].items()}
crank = {c['child']: c['rank'] for c in tpl['idle_children']}
need(sorted(crank) == list(range(18)), 'children')
kind_of = {c: KIND[(pairs[c], crank[c])] for c in range(18)}
need(sorted((kind_of[c], sides[c]) for c in range(18)) == [(k, s) for k in range(9) for s in (0, 1)], 'children <-> (kind, side)')
child_of = [[None, None] for _ in range(9)]
for c in range(18): child_of[kind_of[c]][sides[c]] = c
need(all(HPAIR[kind_of[c]] == pairs[c] for c in range(18)), 'hub pair of kinds')
if g0:
    gk = {int(c): v for c, v in g0['child_to_kind'].items()}
    need(all(CH[kind_of[c]] == gk[c]['kind'] and sides[c] == gk[c]['side'] for c in range(18)), 'child_to_kind == G0')
# ---- calls, order, dispatch
need(len(order) == NC and sorted(map(tuple, order)) == sorted((s, r) for s in range(9) for r in range(T)), 'order covers calls')
kof = lambda s, r: 9 * r + s
time = [None] * NC
for t, (s, r) in enumerate(order): time[kof(s, r)] = t
order_k = [kof(s, r) for s, r in order]
need(all(time[9 * r + l] < time[9 * r + l + 1] for r in range(T) for l in range(8)), 'per-replica chronology')
dmap = {(s, r): tag for s, r, tag in dispatch}
need(len(dmap) == len(dispatch) == NC, 'dispatch coverage')
disp = [PROV[dmap[(k % 9, k // 9)]] for k in range(NC)]
hub_dc = collections.Counter(tag for s, r, tag in dispatch)
# ---- windows: ours vs the hub's template windows
win = [None] * (9 * T)
for r in range(T):
    for c in range(9):
        lo = 0 if KZ[c] == 0 else time[9 * r + KZ[c] - 1] + 1
        win[9 * r + c] = [lo, time[9 * r + KZ[c] + KK[c] - 1]]
spair = {e['stage']: e['pair'] for e in tpl['recipe'] if e['kind'] == 'local'}
hubwin = {}
for o in range(T):
    touched = {}; tick = None
    for e in tpl['execution_schedule']:
        kd = e['kind']
        if kd == 'local_cover_sweep': tick = time[kof(e['stage'], o)]; touched[spair[e['stage']]] = tick
        elif kd == 'bridge': touched[e['src']] = tick; touched[e['dst']] = tick
        elif kd == 'idle_cover_sweep':
            for c in e['child_indices']:
                if c in crank:
                    need((o, c) not in hubwin, 'repeated owner child'); hubwin[(o, c)] = (touched.get(e['pair'], -1), tick)
need(len(hubwin) == 18 * T, 'hub windows')
wcmp = collections.Counter()
for (o, c), (lh, hh) in hubwin.items():
    lo, hi = win[9 * o + kind_of[c]]
    wcmp['equal' if (lo, hi) == (lh + 1, hh) else 'differ'] += 1
ok['windows_ours_eq_hub_template'] = wcmp == {'equal': 18 * T}
# ---- FILL -> rows
fill = rd(files['FILL']); rows = fill['rows']
need(fill['columns'] == ['time', 'helper', 'owner_replica', 'child', 'port', 'offset'] and fill.get('dispatch_overrides', []) == [], 'fill schema')
need(len(rows) == NC * R, 'role count')
hr = [[-1] * R for _ in range(NC)]; offs = [[-1] * R for _ in range(NC)]
used = [0] * NH; bad = collections.Counter(); cross = 0; fams = set()
for j, row in enumerate(rows):
    t, h, o, c, p, z = row
    if (t, h) != divmod(j, R): bad['helper_identity'] += 1
    if not (0 <= o < T and 0 <= c < 18 and 0 <= p < V): bad['range'] += 1; continue
    if h == 0: fams = set()
    k = order_k[t]; kd = kind_of[c]; i = 2016 * o + 224 * kd + 112 * sides[c] + p
    lh, hh = hubwin[(o, c)]
    if not (lh < t <= hh): bad['hub_window'] += 1
    lo, hi = win[9 * o + kd]
    if not (lo <= t <= hi): bad['our_window'] += 1
    f = (o, pairs[c], sides[c], p)
    if f in fams: bad['family_repeated'] += 1
    fams.add(f)
    w = res[disp[k]][h]
    if z != used[i]: bad['compact_prefix'] += 1
    if z + w > crank[c]: bad['capacity'] += 1
    used[i] = z + w
    if hr[k][h] != -1: bad['slot_twice'] += 1
    hr[k][h] = i; offs[k][h] = z
    if o != k // 9: cross += 1
del rows, fill
ok['fill_rules'] = not bad
ok['all_slots_hosted'] = all(e >= 0 for x in hr for e in x)
ok['rows_injective'] = all(len(set(x)) == R for x in hr)
rem = [RANK[i % 2016 // 224] - used[i] for i in range(NH)]
ok['capacity'] = min(rem) >= 0
remh = collections.Counter(x for x in rem if x)
# expected remainders from the hub's files alone: PROFILE.json minus the calls' local profiles
exp_rem = collections.Counter(prof_hub)
for k in range(NC): exp_rem.subtract(lprof[disp[k]])
exp_rem = {a: b for a, b in exp_rem.items() if b}
ok['remainders_eq_PROFILE_minus_calls'] = dict(remh) == exp_rem
prof = collections.Counter()
for k in range(NC): prof.update(lprof[disp[k]])
prof.update(remh)
ok['profile_eq_PROFILE_json'] = dict(prof) == prof_hub
ok['rank_total'] = sum(rem) + sum(res[disp[k]][q] for k in range(NC) for q in range(R)) == 224 * sum(RANK) * T
ok['dispatch_counts_eq_DISPATCH_json'] = (disp.count(0), disp.count(1)) == (hub_dc['extra'], hub_dc['gauge44'])
if g0:
    ok['g0_agrees'] = (g0['G0']['cross_replica_roles'] == cross and
                       {int(a): b for a, b in g0['G0']['remainder_profile'].items()} == dict(remh))
summ = dict(producer='aA-data own_to_host2979H.py (x-tables own_to_host2979.py, FILL = design A 764395fc)', T=T, calls=NC, R=R, NH=NH, nb=0,
            roles=NC * R, cross_replica_roles=cross, cross_share=round(cross / (NC * R), 6),
            dispatch={'extra': disp.count(0), 'gauge44': disp.count(1)},
            provider_slot_diffs=provider_slot_diffs, inventory={str(p): dict(sorted(collections.Counter(res[p]).items())) for p in (0, 1)},
            windows=dict(wcmp), remainders=dict(sorted(remh.items())), profile=dict(sorted(prof.items())),
            profile_hub=dict(sorted(prof_hub.items())), rank_mass=sum(a * b for a, b in prof.items()),
            fill_rule_errors=dict(bad), checks=ok, all_ok=all(ok.values()), inputs_sha256=shas,
            seconds=round(_time.time() - t0, 1))
data = dict(format='p-own host table v1x (lg-a2a4 host table v1 + offs, local_profile, child_of, hpair)',
            T=T, R=R, NH=NH, nb=0, kinds=CH, ranks=RANK, kz=KZ, kk=KK, hpair=HPAIR, child_of=child_of,
            providers={'0': 'extra', '1': 'gauge44'}, res={'0': res[0], '1': res[1]},
            local_profile={str(p): {str(a): b for a, b in sorted(lprof[p].items())} for p in (0, 1)},
            disp=disp, time=time, order_k=order_k, windows=win, hrows=hr, offs=offs, hrem=rem,
            source={k: shas[k] for k in ('FILL', 'ORDER', 'DISPATCH', 'PROFILE', 'PROVIDERS', 'TEMPLATE', 'PLAN', 'WORD-base', 'WORD-extra', 'residuals')})
need(summ['all_ok'], 'a conversion check failed: ' + json.dumps(summ))
raw = json.dumps(data, separators=(',', ':')).encode()
with open(outp + '-hosttab.json.gz', 'wb') as f: f.write(gzip.compress(raw, compresslevel=6, mtime=0))
summ['hosttab_sha256'] = sha(outp + '-hosttab.json.gz'); summ['hosttab_content_sha256'] = hashlib.sha256(raw).hexdigest()
json.dump(summ, open(outp + '-audit.json', 'w'), indent=1)
print(json.dumps(summ)); print('wrote', outp + '-hosttab.json.gz', outp + '-audit.json', 'ALL OK')
