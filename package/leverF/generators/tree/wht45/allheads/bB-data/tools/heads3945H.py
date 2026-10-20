"""bB-data (wht45/allheads, 19 Oct 2026; copy of p-own's heads.py with the counts of design B: 80,640 heads =
40,320 child 8 + 40,320 child 10, and 57,312 non-head F1 bins, checked against g-heads' G0 classes on the new fill;
nothing else changed).  Original header: the exact list of the literal rank-48 heads of #3945, from OUR host table,
cross-checked three ways: (a) against the hub's FILL read directly (bins no FILL row names; own recount),
(b) against the template facts that make an unhosted child-8/10 bin a literal head (heads_fill.py's facts, re-checked
here), (c) against g-heads' G0 classes (counts per child, remainder 48, chronology in our network).

Head = a bin i = 2016 o + 224*4 + 112 s + t of kind A4 (hub children 8 (s = 0, x) and 10 (s = 1, y), hub pair 3,
rank 48) with NO hosted entry in the table.  Its remainder is 48, so the host chain (o, A4, (s, t)) of our network
pays exactly ONE block of rank 48 on it: the climb of that chain over its window, whose line set is the WHOLE pending
set D of the A4 climb at (o, s, t) (cells 0..47 of the bin, none hosted; in Lean the `D` of `host_framesK` for host
j = (s, t) of chain (o, 4), with |D| = 48 and load 0).  Non-heads: the unhosted F1 bins (children 12/13, kind 6,
80,640), not heads because gate C2 reads pair 1 after the F1 climb (g-heads).
usage: python3 -I heads.py HOSTTAB.json.gz FILL.json.gz TEMPLATE.json.gz G0.json OUT.json.gz OUT-summary.json"""
import sys, json, gzip, hashlib, collections
tabf, fillf, tplf, g0f, outf, sumf = sys.argv[1:7]
rd = lambda p: json.loads(gzip.decompress(open(p, 'rb').read())) if p.endswith('.gz') else json.load(open(p))
d = rd(tabf); T, R, NH = d['T'], d['R'], d['NH']; NC = 9 * T
hosted = set()
for row in d['hrows']: hosted.update(row)
child_of = d['child_of']; win = d['windows']; rem = d['hrem']
def unhosted(kind):
    out = []
    for o in range(T):
        for s in (0, 1):
            for t in range(112):
                i = 2016 * o + 224 * kind + 112 * s + t
                if i not in hosted:
                    lo, hi = win[9 * o + kind]
                    out.append([i, o, child_of[kind][s], s, t, lo, hi])
    return out
heads = unhosted(4); non = unhosted(6)
ok = {}
ok['children'] = child_of[4] == [8, 10] and child_of[6] == [12, 13]
byc = collections.Counter(h[2] for h in heads); byn = collections.Counter(h[2] for h in non)
ok['heads_80640'] = len(heads) == 80640 and byc == {8: 40320, 10: 40320}
ok['nonheads_57312'] = len(non) == 57312 and set(byn) <= {12, 13}
ok['heads_rem_48'] = all(rem[h[0]] == 48 for h in heads) and all(rem[h[0]] == 48 for h in non)
r48 = {i for i, x in enumerate(rem) if x == 48}
ok['rem48_bins_eq_heads_plus_nonheads'] = r48 == {h[0] for h in heads} | {h[0] for h in non}
ok['one_head_per_array'] = len({h[0] for h in heads}) == len(heads)
# (a) the hub's FILL, directly
fill = rd(fillf); named = {(o, c, p) for t, h, o, c, p, z in fill['rows']}; del fill
fh = sorted((o, c, p) for o in range(T) for c in (8, 10) for p in range(112) if (o, c, p) not in named)
fn = sorted((o, c, p) for o in range(T) for c in (12, 13) for p in range(112) if (o, c, p) not in named)
ok['fill_recount_heads_equal'] = fh == sorted((h[1], h[2], h[4]) for h in heads)
ok['fill_recount_nonheads_equal'] = fn == sorted((h[1], h[2], h[4]) for h in non)
# (b) template facts (T37-envelope template, sha 21297b20; the facts heads_fill.py relies on)
tpl = rd(tplf); S = tpl['execution_schedule']
cap = {c['child']: c['rank'] for c in tpl['idle_children']}
sp = {e['stage']: e['pair'] for e in tpl['recipe'] if e['kind'] == 'local'}
i14 = next(j for j, e in enumerate(S) if e['kind'] == 'idle_cover_sweep' and e.get('child_indices') == [8, 10])
pre = S[:i14]
ok['template_heads_first'] = (all(cap[c] == 48 for c in (8, 10, 12, 13)) and S[i14]['pair'] == 3 and S[i14]['from_level'] == -1
    and not any(e['kind'] == 'local_cover_sweep' and sp[e['stage']] == 3 for e in pre)
    and not any(e['kind'] == 'bridge' and 3 in (e['src'], e['dst']) for e in pre)
    and not any(e['kind'] == 'idle_cover_sweep' and e['pair'] == 3 for e in pre))
i21 = next(j for j, e in enumerate(S) if e['kind'] == 'idle_cover_sweep' and e.get('child_indices') == [12, 13])
ok['template_F1_read_after'] = any(e['kind'] == 'bridge' and e['src'] == 0 for e in S[i21 + 1:])
# (c) g-heads G0
g0 = rd(g0f)['G0']; cls = {(c['hub_child'], c['remainder'], c['cls']): c['count'] for c in g0['classes']}
ok['g_heads_agrees'] = (cls.get((8, 48, 'HEAD-literal')) == 40320 and cls.get((10, 48, 'HEAD-literal')) == 40320
    and all(cls.get((c, 48, 'neither'), 0) == byn.get(c, 0) for c in (12, 13))
    and g0['heads_rank48_unhosted'] == 80640 and g0['rank48_by_class'] == {'HEAD-literal': 80640, 'neither': 57312})
data = dict(format='p-own heads v1', T=T, columns=['bin', 'owner', 'child', 'side', 'port', 'win_lo', 'win_hi'],
            kind='A4 (our kind 4, pair 4 = hub pair 3), rank 48', count=len(heads), by_child={str(a): b for a, b in sorted(byc.items())},
            lines='every head: the whole pending set D of the A4 climb of its bin, |D| = 48, cells 0..47, load 0 (remainder 48)',
            heads=heads,
            nonheads=dict(kind='F1 (our kind 6, pair 1 = hub pair 0), rank 48', count=len(non),
                          by_child={str(a): b for a, b in sorted(byn.items())}, rows=non),
            source=dict(hosttab_sha256=hashlib.sha256(open(tabf, 'rb').read()).hexdigest()))
raw = json.dumps(data, separators=(',', ':')).encode()
open(outf, 'wb').write(gzip.compress(raw, compresslevel=6, mtime=0))
summ = dict(checks=ok, all_ok=all(ok.values()), heads=len(heads), by_child=dict(byc), nonheads=len(non), nonheads_by_child=dict(byn),
            heads_sha256=hashlib.sha256(open(outf, 'rb').read()).hexdigest(), heads_content_sha256=hashlib.sha256(raw).hexdigest(),
            head_bins_sha256=hashlib.sha256(json.dumps([h[0] for h in heads]).encode()).hexdigest(),
            first_heads=heads[:3], last_heads=heads[-3:])
json.dump(summ, open(sumf, 'w'), indent=1); print(json.dumps(summ)); print('ALL OK' if summ['all_ok'] else 'SOME CHECK FAILED')
