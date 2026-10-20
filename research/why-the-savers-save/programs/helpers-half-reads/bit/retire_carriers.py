"""d-halfread BIT SIDE: retire single-source carriers of #589's tuned bit word by direct reads.  For a carrier h that
copies source s at s's line frame in phase A and is then read by mini-node helpers: every gate that reads h reads s
instead (at the same gate frame, so s climbs there, or for carrier 470 at s's current frame [33,77]), h's own gates
go, h leaves D, C is rebuilt as the literal reverse (wlib.rebuild_C), h is deleted and registers renumbered.
Checked by the hub's exact replay (wlib.paid) and then by xcheck_p.py.  Usage: python3 -I retire_carriers.py [h ...]"""
import sys, os, json, copy, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import wlib
w0 = wlib.load(os.path.join(HERE, 'bitE8mut-p111-tuned.json.gz'))
V = w0['v']; base = wlib.paid(w0)
def retire(w, hs):
    w = copy.deepcopy(w)
    for h in hs:
        own = [g for g in w['A'] if g[1] == h]
        assert len(own) == 1 and len(own[0][2]) == 1 and own[0][2][0] < V, ('not a single-source carrier', h, own)
        s = own[0][2][0]
        w['A'] = [g for g in w['A'] if g[1] != h]
        for sec in ('A', 'B'):
            new = []
            for f, t, ss in w[sec]:
                if h in ss:
                    if h == 470: f = w['frames'].index([33, 77])
                    ss = sorted(set(ss) ^ {h, s})
                new.append([f, t, ss])
            w[sec] = new
        for g in w['D']:
            if h in g[2]: g[2].remove(h)
        w['D'] = [g for g in w['D'] if g[2]]
    w['C'] = wlib.rebuild_C(w)
    # delete the registers and renumber
    dead = sorted(hs)
    f = lambda q: q - sum(1 for d in dead if d < q)
    for sec in 'DABC':
        assert all(g[1] not in dead and not (set(g[2]) & set(dead)) for g in w[sec]), 'carrier still used'
        w[sec] = [[g[0], f(g[1]), [f(x) for x in g[2]]] for g in w[sec]]
    w['ret'] = [f(q) for q in w['ret']]
    for d in reversed(dead): del w['start'][d]; del w['final'][d]
    w['R'] -= len(dead); w['stats'] = dict(w['stats']); w['stats']['extra_carriers'] -= len(dead)
    return w
cands = [462, 463, 464, 467, 469, 470, 471]
res = {}
for h in cands:
    try:
        w = retire(w0, [h]); p = wlib.paid(w); res[h] = dict(ok=True, profile=p, delta={r: p.get(r, 0) - base.get(r, 0) for r in set(p) | set(base) if p.get(r, 0) != base.get(r, 0)})
    except AssertionError as e:
        res[h] = dict(ok=False, why=str(e)[:120])
    print(h, res[h])
ok = [h for h in cands if res[h]['ok']]
best = None
for k in range(len(ok), 0, -1):
    for combo in itertools.combinations(ok, k):
        try:
            w = retire(w0, list(combo)); p = wlib.paid(w); best = (combo, p, w); break
        except AssertionError: pass
    if best: break
if best:
    combo, p, w = best
    print('largest jointly retirable set', combo, 'R', w['R'], 'profile', p, 'mass', sum(r*n for r, n in p.items()))
    sha = wlib.save(w, os.path.join(HERE, 'bitE8mut-p111-tuned-retire%d.json.gz' % len(combo)))
    json.dump(dict(base=base, single=res, joint=list(combo), joint_profile=p, R=w['R'], sha256=sha), open(os.path.join(HERE, 'retire_carriers.json'), 'w'), indent=1, default=str)
    print('saved', 'bitE8mut-p111-tuned-retire%d.json.gz' % len(combo), sha)
