"""bB-data (wht45/allheads, 19 Oct 2026; copy of p-own's own_to_host.py, only the FILL pin (all-heads fill 9c06677d) and the expected cross-replica count (1,161,423, = G0 on the new fill) changed): #3945's literal owner fill -> OUR cross-replica host table (first code path).
Own code. The hub's files are read as data; the hub's verify.py was read as a reference, never run.

FORMAT "p-own host table v1x" = lg-a2a4 host table v1 (hub_to_host.py; the format p-a3's mk_kdata.py reads), plus
offsets and bookkeeping.  Generalised to T = 360 replicas, R = 408 slots, a global order and two providers:
  calls k = 9 rho + l (replica rho, stage l), NC = 9 T; slot q = the FILL's helper index h (0..407, the compacted
    borrow numbering of g-heads' gcert/1b, helper_origin[q] in the 480-helper word);
  hrows[k][q] = host index 2016 rho_h + 224 c + 112 s + t of slot q of call k (rho_h = FILL owner_replica, c = our
    window kind of the FILL child, s = side, t = FILL port); every slot is hosted (no banks, nb = 0, NH = 2016 T);
  offs[k][q] = the FILL offset (cells [offs, offs + res) of the bin);
  res[p][q] = residual width of slot q in provider p (p = 0 extra, 1 gauge44): our lowering's (g-heads) residuals,
    checked equal, index for index, to the hub's widths dim(helper_final) - dim(start_extra);
  disp[k] = provider of call k; time[k] = position of call k in ORDER; order_k = calls in ORDER;
  windows[9 rho + c] = [lo, hi] (inclusive) by OUR rule: lo = 0 if KZ[c] = 0 else time[9 rho + KZ[c] - 1] + 1,
    hi = time[9 rho + KZ[c] + KK[c] - 1]; checked equal to the hub's template windows (lo_h, hi_h] for every
    (owner, child);
  hrem[i] = rank(kind of i) - load(i), load(i) = sum of res[disp[k]][q] over the entries hrows[k][q] = i;
  local_profile[p] = the provider's paid call profile (word report; = g-heads' lowered block histogram).
Kinds c = 0..8: A2 B2 A3 B3 A4 B4 F1 F2 F3, ranks 16 8 32 8 48 8 48 32 16; hub pair of a kind: 1 1 2 2 3 3 0 1 2.
Child -> kind: (hub pair, rank) -> kind, side = PLAN side (g-heads child_to_kind; children 8/10 = A4 x/y).

usage: python3 -I own_to_host.py OWNER_DIR RESIDUALS.json G0.json OUT_PREFIX
  OWNER_DIR = #3945 research/.../owner (FILL.json.gz, PROFILE.json, inputs/{ORDER,DISPATCH}.json,
  inputs/{TEMPLATE,PLAN,WORD-base,WORD-extra,old-NATIVE-base,old-NATIVE-extra}.json.gz);
  RESIDUALS.json = g-heads data/residuals-3945.json; G0.json = g-heads out/g0-3945.json (block histograms)."""
import sys, json, gzip, hashlib, collections, time as _time
t0 = _time.time()
own, residf, g0f, outp = sys.argv[1:5]
CH = ['A2', 'B2', 'A3', 'B3', 'A4', 'B4', 'F1', 'F2', 'F3']
RANK = [16, 8, 32, 8, 48, 8, 48, 32, 16]
KZ = [0, 2, 0, 4, 0, 6, 3, 5, 7]; KK = [2, 1, 4, 1, 6, 1, 6, 4, 2]
KIND = {(1, 16): 0, (1, 8): 1, (2, 32): 2, (2, 8): 3, (3, 48): 4, (3, 8): 5, (0, 48): 6, (1, 32): 7, (2, 16): 8}
HPAIR = [1, 1, 2, 2, 3, 3, 0, 1, 2]
PROV = {'extra': 0, 'gauge44': 1}; PNAME = {0: 'extra', 1: 'gauge44'}; TAG = {0: 'extra', 1: 'base'}
T, R, V = 360, 408, 112; NC = 9 * T; NH = 2016 * T
def rd(p):
    b = open(p, 'rb').read(); return json.loads(gzip.decompress(b) if p.endswith('.gz') else b)
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def need(c, why):
    if not c: raise SystemExit('REFUSED: ' + why)
I = own + '/inputs/'
files = {'FILL': own + '/FILL.json.gz', 'PROFILE': own + '/PROFILE.json', 'ORDER': I + 'ORDER.json',
         'DISPATCH': I + 'DISPATCH.json', 'TEMPLATE': I + 'TEMPLATE.json.gz', 'PLAN': I + 'PLAN.json.gz',
         'WORD-base': I + 'WORD-base.json.gz', 'WORD-extra': I + 'WORD-extra.json.gz',
         'old-NATIVE-base': I + 'old-NATIVE-base.json.gz', 'old-NATIVE-extra': I + 'old-NATIVE-extra.json.gz',
         'residuals': residf, 'g0': g0f}
shas = {k: sha(v) for k, v in files.items()}
need(shas['FILL'] == '9c06677d66b04436b02ec8d88036327bdc8fb12e6fd4f71a1b9c7e915c10c2d3', 'FILL sha (design B all-heads fill, reopt FOUND.md)')
need(shas['TEMPLATE'] == '21297b20f1f2dc76eb9dd8124416068c3cb096bab9df17fa4da19bc68dfe082e', 'TEMPLATE sha (T74 template)')
order = rd(files['ORDER']); dispatch = rd(files['DISPATCH']); tpl = rd(files['TEMPLATE']); plan = rd(files['PLAN'])
prof_hub = {int(a): b for a, b in rd(files['PROFILE']).items()}
resid = rd(residf); g0 = rd(g0f)
ok = {}
# ---- widths: hub rule (own GF(2) rank) vs our lowering's residuals
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
    w = rd(files['WORD-' + TAG[p]]); nat = rd(files['old-NATIVE-' + TAG[p]]); origin = nat['helper_origin']
    need(origin == sorted(set(range(224, 704)) - set(w['removed'])) and len(origin) == R, 'helper origin')
    hubw = [rank2(w['helper_final'].get(str(r), FULL)) - rank2(w['start_extra'].get(str(r), [])) for r in origin]
    ours = resid[PNAME[p]]['residuals']
    need(resid[PNAME[p]]['helper_origin'] == origin, 'our lowering uses the same helper order')
    need(ours == hubw, 'our residuals == hub widths, index for index (' + PNAME[p] + ')')
    need(collections.Counter(ours) == {5: 28, 6: 18, 7: 53, 8: 96, 9: 213}, 'inventory')
    res[p] = ours
    lp = {int(a): b for a, b in w['report']['local_profile'].items()}
    gh = {int(a): b for a, b in g0['inputs'][PNAME[p]]['block_hist'].items()}
    need(lp == gh == {1: 2482, 2: 712, 3: 218, 4: 109, 5: 20, 8: 8}, 'local profile == our lowered histogram')
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
gk = {int(c): v for c, v in g0['child_to_kind'].items()}
need(all(CH[kind_of[c]] == gk[c]['kind'] and sides[c] == gk[c]['side'] for c in range(18)), 'child_to_kind == g-heads')
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
ok['remainders'] = remh == {3: 23328, 5: 4320, 48: 137952}
prof = collections.Counter()
for k in range(NC): prof.update(lprof[disp[k]])
prof.update(remh)
ok['profile_eq_PROFILE_json'] = dict(prof) == prof_hub
ok['mass_deficit'] = sum(a * b for a, b in prof.items()) == 65088 * T
ok['cross_replica_roles'] = cross == 1161423   # design B (g-heads' G0 on the new fill: 1,161,423; #3945's own fill: 1,152,628)
ok['dispatch_counts'] = (disp.count(0), disp.count(1)) == (1672, 1568)
ok['g0_agrees'] = (g0['G0']['cross_replica_roles'] == cross and
                   {int(a): b for a, b in g0['G0']['remainder_profile'].items()} == dict(remh))
summ = dict(producer='p-own own_to_host.py', T=T, calls=NC, R=R, NH=NH, nb=0, roles=NC * R,
            cross_replica_roles=cross, cross_share=round(cross / (NC * R), 6),
            dispatch={'extra': disp.count(0), 'gauge44': disp.count(1)},
            provider_slot_diffs=provider_slot_diffs, inventory=dict(sorted(collections.Counter(res[0]).items())),
            windows=dict(wcmp), remainders=dict(sorted(remh.items())), profile=dict(sorted(prof.items())),
            profile_hub=dict(sorted(prof_hub.items())), rank_mass=sum(a * b for a, b in prof.items()),
            fill_rule_errors=dict(bad), checks=ok, all_ok=all(ok.values()), inputs_sha256=shas,
            seconds=round(_time.time() - t0, 1))
data = dict(format='p-own host table v1x (lg-a2a4 host table v1 + offs, local_profile, child_of, hpair)',
            T=T, R=R, NH=NH, nb=0, kinds=CH, ranks=RANK, kz=KZ, kk=KK, hpair=HPAIR, child_of=child_of,
            providers={'0': 'extra', '1': 'gauge44'}, res={'0': res[0], '1': res[1]},
            local_profile={str(p): {str(a): b for a, b in sorted(lprof[p].items())} for p in (0, 1)},
            disp=disp, time=time, order_k=order_k, windows=win, hrows=hr, offs=offs, hrem=rem,
            source={k: shas[k] for k in ('FILL', 'ORDER', 'DISPATCH', 'TEMPLATE', 'PLAN', 'WORD-base', 'WORD-extra', 'residuals')})
need(summ['all_ok'], 'a conversion check failed: ' + json.dumps(summ))
raw = json.dumps(data, separators=(',', ':')).encode()
with open(outp + '-hosttab.json.gz', 'wb') as f: f.write(gzip.compress(raw, compresslevel=6, mtime=0))
summ['hosttab_sha256'] = sha(outp + '-hosttab.json.gz'); summ['hosttab_content_sha256'] = hashlib.sha256(raw).hexdigest()
json.dump(summ, open(outp + '-audit.json', 'w'), indent=1)
print(json.dumps(summ)); print('wrote', outp + '-hosttab.json.gz', outp + '-audit.json', 'ALL OK')
