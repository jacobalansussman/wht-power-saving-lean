#!/usr/bin/env python3
"""g-heads (wht45/lever, 19 Oct 2026): gate G0 of rl-plan section 6 for a staged544 T360 word (#3945, #3815), in OUR
lowering (not the hub's controller). Adapted from wht43/bound-check/lg-g0/tools/g0_ours.py (T74, PASS); the network
part (parsing OUR hosted nine-stage network from the Lean source, the gate read/write sets from Shea's Gates.lean,
windows, segments, the chronology of every remainder block, the controls) is lg-g0's code unchanged; what changes:
  * the per-call words: OUR gcert/1b lowerings (tools/lower_b.py) of the hub's native canonical BORROW words (re-run on
    Colab C1 from the PR's own word files), no copies, exit-route residuals;
  * the call data: the staged544 package's ORDER (call order, T = 360 replicas), DISPATCH (provider per call), PLAN
    (pairs and sides of the 18 children), TEMPLATE (child ranks) and the literal FILL rows
    [time, helper, owner_replica, child, port, offset] (1,321,920 roles), all read as DATA;
  * width of a role = OUR residual of that helper in that call's provider (checked: the FILL's offsets are exactly the
    compact prefix sums of our widths in time order, and no bin overflows its rank).
Questions (G0): (1) which remainder blocks are literal heads in OUR network (first non-identity operation on their
array, on that array alone, nothing hosted in that array before)? (2) remainder profile? (3) call-block ranks <= 8?
(4) does OUR lowering's global profile equal the hub's PROFILE.json digit for digit?
usage: g0_staged.py BC PKG FILL GCERT_BASE GCERT_EXTRA PROFILE OUT.json [--gates Gates.lean] [--label L]
       (BC = .../wht43/bound-check; PKG = a staged544 owner dir with inputs/{ORDER,DISPATCH,PLAN.json.gz,TEMPLATE.json.gz})
"""
import sys, re, json, gzip, hashlib
from collections import Counter, defaultdict

a = sys.argv[1:]
BC, PKG, FILLP, GB, GE, PROFP, OUT = a[:7]
GATES = a[a.index('--gates') + 1] if '--gates' in a else BC + '/lg-g0/colab/Gates.lean'
LABEL = a[a.index('--label') + 1] if '--label' in a else ''
V = 112
rd = lambda p: json.loads(gzip.decompress(open(p, 'rb').read())) if p.endswith('.gz') else json.load(open(p))
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
res = dict(label=LABEL, inputs={})

# ---------------------------------------------------------------- 1. the per-call words (our lowering)
def blocks(c):
    """climb histogram of a gcert/1b word (warps free), lower_b.check's accounting"""
    v, R = c['v'], c['R']; dim = [len(f) for f in c['frames']]
    cur = list(c['start']); H = Counter()
    def regs(g):
        return [g[2]] if g[0] == 'warp' else [g[2]] + [t for t, _, _ in g[3]] + list(g[4])
    def climb(r, f):
        if cur[r] != f:
            H[dim[f] - dim[cur[r]]] += 1; cur[r] = f
    for ph in ('A', 'B'):
        if ph == 'B':
            for _, s, f in c['ret']:
                assert cur[s] == f; H[dim[f]] += 1
        for g in c[ph]:
            if g[0] == 'warp':
                assert cur[g[2]] == g[4]; cur[g[2]] = g[1]
            else:
                for r in regs(g): climb(r, g[1])
    for r in range(2 * v + R): climb(r, c['final'][r])
    return dict(sorted(H.items()))

words = {}
for prov, p in (('gauge44', GB), ('extra', GE)):
    c = rd(p)
    assert c['format'] == 'gcert/1b' and c.get('copies') == 0 and c['v'] == V
    dim = [len(f) for f in c['frames']]
    E = set(c['exits']); v, R, h = c['v'], c['R'], c['h']
    resq = [(dim[c['final'][2 * v + q]] if 2 * v + q in E else h) - dim[c['start'][2 * v + q]] for q in range(R)]
    assert resq == c['residuals']
    hist = blocks(c)
    words[prov] = dict(R=R, res=resq, hist=hist)
    res['inputs'][prov] = dict(gcert=p, gcert_sha256=sha(p), R=R, h=h, block_hist={str(k): x for k, x in hist.items()},
                               max_block_rank=max(hist), res_hist={str(k): x for k, x in sorted(Counter(resq).items())},
                               res_sum=sum(resq), exits=len(E))

# ---------------------------------------------------------------- 2. our network, parsed from the Lean source (lg-g0, unchanged)
XH = a[a.index('--xhost') + 1] if '--xhost' in a else BC + '/h-L3-build/lean/Work/XHost'   # release copy: --xhost added
hst = open(XH + '/HStage.lean').read()
hce = open(XH + '/HCert.lean').read()
hk = open(XH + '/HostK.lean').read()
res['inputs']['lean'] = {f: sha(XH + '/' + f) for f in ('HStage.lean', 'HCert.lean', 'HostK.lean')}
runtxt = re.search(r'lemma GHost401\.run \(ρ : Fin 8\) :.*?have R := (.*?)\n', hce, re.S).group(1)
RUN = re.findall(r'N\.(\w+) ρ', runtxt)
steps = {}
for part in re.split(r'\n(?=/-- )', hst):
    m = re.search(r'lemma GHost401\.(\w+)', part)
    if not m: continue
    name = m.group(1)
    bul = [l.strip() for l in part.split('\n') if l.strip().startswith('· ')][-10:]
    assert len(bul) == 10, name
    tch = {}
    for comp, l in enumerate(bul[:8], start=1):
        if name.startswith('gate'):
            k = None if '(hk rfl).elim' in l else 'gate_write'
        elif name.startswith('adapt'):
            k = None if 'GReach.refl' in l else 'adapter'
            mm = re.search(r'hHP_(\w\d)_\d|Hp (\d) \d', l)
            if k: k = ('adapter', mm.group(1) or 'kind%s' % mm.group(2))
        elif name.startswith('stage'):
            if 'jX' in l or 'jY' in l: k = 'local'
            elif 'hHZn_' in l: k = ('hostuse', re.search(r'hHZn_(\w\d)_', l).group(1))
            else: raise SystemExit('stage case ' + l)
        elif name.startswith('climb'):
            k = None if 'GReach.refl' in l else ('climb', int(re.search(r'hrem (\d)', l).group(1)))
        elif name == 'fin':
            if 'hHE_' in l: k = ('climb', {'F1': 6, 'F2': 7, 'F3': 8}[re.search(r'hHE_(F\d)', l).group(1)])
            elif 'GReach.ofX' in l: k = 'free_final'
            else: k = None if 'GReach.refl' in l else 'unknown'
        else:
            raise SystemExit('step ' + name)
        if k is not None: tch[comp] = k
    steps[name] = tch
assert set(RUN) <= set(steps), set(RUN) - set(steps)
gtxt = open(GATES).read()
res['inputs']['Gates.lean'] = dict(path=GATES, sha256=sha(GATES))
RQ = {'rqX1': 1, 'rqY1': 2, 'rqX2': 3, 'rqY2': 4, 'rqX3': 5, 'rqY3': 6, 'rqX4': 7, 'rqY4': 8}
def gdef(nm):
    return re.search(r'def ' + nm + r' : BRun4 Γ T nb C → [^\n]*:=\s*\n\s*mkR4 (.*)\n', gtxt).group(1)
GREAD = {}
for name in RUN:
    if not name.startswith('gate'): continue
    pn, kn = re.search(r'lemma GHost401\.' + name + r' .*?gscQ (\w+) (\w+)\)', hst).groups()
    pp = [RQ.get(x) for x in re.findall(r'=> (rq\w+) \(', gdef(pn))]
    kk = [int(x) for x in re.findall(r'=> (-?\d+)\)', gdef(kn))]
    assert len(pp) == 10 and len(kk) == 10, (name, pp, kk)
    wr = {i + 1 for i in range(8) if kk[i]}
    rdd = {pp[i] for i in range(8) if kk[i]}
    assert wr == set(steps[name]), (name, wr, steps[name])
    assert not wr & rdd
    for c in rdd: steps[name][c] = 'gate_read'
    GREAD[name] = dict(p=pn, k=kn, written=sorted(wr), read=sorted(rdd), coeff={str(i + 1): kk[i] for i in range(8) if kk[i]})
KN = ['A2', 'B2', 'A3', 'B3', 'A4', 'B4', 'F1', 'F2', 'F3']
for name, t in steps.items():
    for c, k in list(t.items()):
        if isinstance(k, tuple) and k[0] == 'adapter' and k[1].startswith('kind'):
            t[c] = ('adapter', KN[int(k[1][4:])])
inwin_src = re.search(r'def inWin \(c l : Nat\) : Bool :=\s*match c with(.*?)\n\n', hk, re.S).group(1)
rank_src = re.search(r'def rankOf \(c : Nat\) : Nat :=\s*match c with(.*?)\n\n', hk, re.S).group(1)
RANK = {int(x): int(y) for x, y in re.findall(r'\| (\d) => (\d+)', rank_src)}
def inWin(c, l):
    expr = dict(re.findall(r'\| (\d) => (.*)', inwin_src))[str(c)]
    return bool(eval(expr.replace('&&', ' and ').replace('≤', '<='), {'l': l}))
KIND = {}
for name in RUN:
    for comp, k in steps[name].items():
        if isinstance(k, tuple) and k[0] == 'climb':
            KIND.setdefault(k[1], [name, []])[1].append(comp)
assert sorted(KIND) == list(range(9)) and all(len(c) == 2 for _, c in KIND.values())
kind_pair = {c: (comps[0] + 1) // 2 for c, (_, comps) in KIND.items()}
res['our_network'] = dict(gates=GREAD, run_order=RUN,
                          kinds={KN[c]: dict(climb_step=s, comps=comps, pair=kind_pair[c], rank=RANK[c]) for c, (s, comps) in KIND.items()})
SEG = {}
cur = None
for name in RUN:
    if name.startswith('adapt'): cur = int(name[5:]); SEG[cur] = []
    SEG[cur].append(name)
assert sorted(SEG) == list(range(9)) and all(SEG[l][1] == 'stage%d' % l for l in range(9))
END = [('padR', {}), ('shifts', {2: 'free_shift', 4: 'free_shift', 6: 'free_shift', 8: 'free_shift'}),
       ('brec401', {c: 'relabel' for c in range(1, 9)})]
assert 'mkD4 (fun _ _ => 0) (fun g t => N.u g t) (fun _ _ => 0) (fun g t => N.u g t) (fun _ _ => 0) (fun g t => N.u g t) (fun _ _ => 0) (fun g t => N.u g t)' in hce

# ---------------------------------------------------------------- 3. the staged544 call data
inp = PKG.rstrip('/') + '/inputs/'
order = [tuple(x) for x in rd(inp + 'ORDER.json')]
T = 1 + max(r for _, r in order)
disp = rd(inp + 'DISPATCH.json')
prov = {(s, r): p for s, r, p in disp}
assert len(prov) == len(order) == 9 * T and set(prov.values()) == {'gauge44', 'extra'}
plan = rd(inp + 'PLAN.json.gz'); tmpl = rd(inp + 'TEMPLATE.json.gz')
for f in ('ORDER.json', 'DISPATCH.json', 'PLAN.json.gz', 'TEMPLATE.json.gz'):
    res['inputs'][f] = sha(inp + f)
res['inputs']['FILL'] = dict(path=FILLP, sha256=sha(FILLP))
res['inputs']['PROFILE'] = dict(path=PROFP, sha256=sha(PROFP))
per = defaultdict(list)
for t, (l, r) in enumerate(order): per[r].append(l)
assert all(x == list(range(9)) for x in per.values()) and len(per) == T
hubpair = {int(k): x for k, x in plan['pairs'].items()}
hubside = {int(k): x for k, x in plan['sides'].items()}
hubrank = {}
for e in tmpl['idle_children']:
    hubrank.setdefault(e['child'], set()).add(e['rank'])
assert all(len(s) == 1 for s in hubrank.values())
hubrank = {ch: s.pop() for ch, s in hubrank.items()}
ch2kind = {}
for ch in hubpair:
    ks = [c for c in range(9) if kind_pair[c] == hubpair[ch] + 1 and RANK[c] == hubrank.get(ch)]
    assert len(ks) == 1, (ch, ks)
    ch2kind[ch] = ks[0]
res['child_to_kind'] = {str(ch): dict(kind=KN[k], hub_pair=hubpair[ch], side=hubside[ch], rank=hubrank[ch]) for ch, k in sorted(ch2kind.items())}

def build(segs, end):
    pos = 0; callpos = {}; stepos = {}
    for t, (l, r) in enumerate(order):
        for n in segs[l]:
            stepos[r, n] = pos
            if n == 'stage%d' % l: callpos[t] = pos
            pos += 1
    endpos = {}
    for n, _ in end: endpos[n] = pos; pos += 1
    return callpos, stepos, endpos

fill = rd(FILLP)
assert fill['columns'] == ['time', 'helper', 'owner_replica', 'child', 'port', 'offset']
ROWS = fill['rows']; del fill
HPC = len(ROWS) // len(order)
assert HPC * len(order) == len(ROWS) and all(words[p]['R'] == HPC for p in words)
res['fill'] = dict(rows=len(ROWS), roles_per_call=HPC, calls=len(order), replicas=T)

def analyse(segs, end, label):
    callpos, stepos, endpos = build(segs, end)
    load = Counter(); bad = Counter(); onecall = Counter(); win = {}
    uses = defaultdict(list)
    cross = 0
    for j, (t, hlp, rep, ch, q, off) in enumerate(ROWS):
        assert t == j // HPC and hlp == j % HPC
        l, rcall = order[t]
        c = ch2kind[ch]; hp, sd = hubpair[ch], hubside[ch]
        comp = 2 * (hp + 1) - 1 + sd
        assert comp in KIND[c][1]
        w = words[prov[(l, rcall)]]['res'][hlp]
        key = (rep, c, sd, q)
        if off != load[key]: bad['offset_not_compact_prefix_of_our_widths'] += 1
        load[key] += w
        onecall[(t, key)] += 1
        if rcall != rep: cross += 1
        p = callpos[t]
        if (rep, c, comp) not in win:
            climb = stepos[rep, KIND[c][0]]
            lastown = max([stepos[rep, n] for l2 in range(9) for n in segs[l2]
                           if stepos[rep, n] < climb and comp in steps[n]
                           and not (isinstance(steps[n][comp], tuple) and steps[n][comp][0] in ('adapter', 'hostuse') and steps[n][comp][1] == KN[c])],
                          default=-1)
            win[rep, c, comp] = (lastown, climb)
        lastown, climb = win[rep, c, comp]
        if not (lastown < p < climb): bad['outside_our_window'] += 1
        if rcall == rep and not inWin(c, l): bad['samerun_not_inWin'] += 1
        uses[(rep, comp, q)].append((p, KN[c], t))
    over = [k for k, x in load.items() if x > RANK[k[1]]]
    rem = {}
    for rep in range(T):
        for c in range(9):
            for sd in (0, 1):
                for q in range(V):
                    r = RANK[c] - load[(rep, c, sd, q)]
                    if r: rem[(rep, c, sd, q)] = r
    prof = Counter(rem.values())
    cls = Counter(); ex = {}
    for (rep, c, sd, q), r in rem.items():
        comp = KIND[c][1][sd]
        climb = stepos[rep, KIND[c][0]]
        ev = []
        for l2 in range(9):
            for n in segs[l2]:
                k = steps[n].get(comp)
                if k is None or n == KIND[c][0]: continue
                p = stepos[rep, n]
                if isinstance(k, tuple) and k[0] in ('adapter', 'hostuse'):
                    kind_used = any(u[1] == k[1] for u in uses[(rep, comp, q)])
                    ev.append((p, n + ':' + str(k[0]) + ':' + k[1], not kind_used))
                else:
                    ev.append((p, n + ':' + (k if isinstance(k, str) else k[0]), k in ('free_final',)))
        for p, kname, t in uses[(rep, comp, q)]:
            ev.append((p, 'hosted:' + kname, False))
        for n, tch in end:
            if comp in tch: ev.append((endpos[n], n + ':' + tch[comp], tch[comp] in ('free_shift', 'relabel')))
        before = [e for e in ev if e[0] < climb and not e[2]]
        after = [e for e in ev if e[0] > climb and not e[2]]
        k = 'HEAD-literal' if not before else 'TAIL-literal' if not after else 'neither'
        hub_child = next(ch for ch, kk in ch2kind.items() if kk == c and hubside[ch] == sd)
        cls[(hub_child, KN[c], 'xy'[sd], r, k)] += 1
        exk = (hub_child, r, k)
        if exk not in ex:
            ex[exk] = dict(replica=rep, port=q, comp=comp, climb_pos=climb,
                           first_nonidentity_before=[list(map(str, e)) for e in sorted(before)[:3]],
                           first_nonfree_after=[list(map(str, e)) for e in sorted(after)[:3]],
                           identity_touches_before=len([e for e in ev if e[0] < climb and e[2]]),
                           free_touches_after=[e[1] for e in sorted(ev) if e[0] > climb and e[2]])
    tot = Counter(); tot48 = Counter()
    for k, n in cls.items():
        tot[k[4]] += n
        if k[3] == 48: tot48[k[4]] += n
    return dict(label=label, roles=len(ROWS), cross_replica_roles=cross, violations=dict(bad), loads_over_rank=len(over),
                one_use_per_call_and_role=max(onecall.values()) == 1,
                remainders=len(rem), remainder_profile={str(k): x for k, x in sorted(prof.items())},
                classes=[dict(hub_child=k[0], kind=k[1], side=k[2], remainder=k[3], cls=k[4], count=n) for k, n in sorted(cls.items())],
                totals_by_class=dict(sorted(tot.items())), rank48_by_class=dict(sorted(tot48.items())),
                heads_rank48_unhosted=tot48.get('HEAD-literal', 0),
                examples={f'{x}:{y}:{z}': e for (x, y, z), e in sorted(ex.items())})

main = analyse(SEG, END, 'our network, %s order' % LABEL)
res['G0'] = main
seg_n1 = {l: list(s) for l, s in SEG.items()}
seg_n1[5] = ['adapt5', 'stage5', 'gateA4', 'climbA4']
n1 = analyse(seg_n1, END, 'NEGATIVE: gateA4 moved before climbA4 (the A4 heads must vanish)')
res['control_gateA4_before_climbA4'] = dict(totals_by_class=n1['totals_by_class'], rank48_by_class=n1['rank48_by_class'], violations=n1['violations'])
del n1
seg_n2 = {l: list(s) for l, s in SEG.items()}; seg_n2[8] = [n for n in seg_n2[8] if n != 'gateC2']
n2 = analyse(seg_n2, END, 'POSITIVE (tail detector): gate C2 deleted; pair-1 F1 blocks must become tails')
res['control_without_gateC2'] = dict(totals_by_class=n2['totals_by_class'], rank48_by_class=n2['rank48_by_class'],
                                     classes=[x for x in n2['classes'] if x['kind'] == 'F1'])
del n2

# ---------------------------------------------------------------- 4. call blocks and the global profile
calls = Counter(prov[(l, r)] for l, r in order)
prof = Counter()
for p, n in calls.items():
    for r, k in words[p]['hist'].items(): prof[r] += n * k
for r, k in main['remainder_profile'].items(): prof[int(r)] += k
hub = {int(k): x for k, x in rd(PROFP).items()}
res['call_blocks'] = dict(calls_per_provider=dict(calls), max_rank_over_both_words=max(max(words[p]['hist']) for p in words),
                          all_call_blocks_rank_le_8=all(max(words[p]['hist']) <= 8 for p in words),
                          global_profile=dict(sorted(prof.items())), hub_profile=dict(sorted(hub.items())),
                          global_profile_equals_hub=dict(sorted(prof.items())) == dict(sorted(hub.items())),
                          rank_mass=sum(r * k for r, k in prof.items()))
res['ok'] = main['heads_rank48_unhosted'] >= 1 and not main['violations'] and main['loads_over_rank'] == 0
json.dump(res, open(OUT, 'w'), indent=1, default=str)
print(json.dumps({k: res[k] for k in ('label', 'call_blocks', 'child_to_kind', 'fill', 'ok')}, indent=1))
print(json.dumps({k: x for k, x in main.items() if k not in ('examples',)}, indent=1))
print('controls:', json.dumps(res['control_gateA4_before_climbA4']), json.dumps(res['control_without_gateC2']))
