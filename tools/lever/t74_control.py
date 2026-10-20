#!/usr/bin/env python3
"""g-heads: regression control for tools/g0_staged.py. Re-packages the T74 data of lg-g0 (bh-gap's gcert/1 lowerings of
the #1839 words, the #1839 PLAN order / dispatch / OWNERS, the T37 template that #3945 pins byte for byte) in the
staged544 package format (ORDER, DISPATCH, PLAN, TEMPLATE, FILL rows, PROFILE) and the words in gcert/1b (no copies,
warps = [], levels = {}, residuals on the exit route). g0_staged.py on this package must reproduce lg-g0's PASS:
heads 11,010 (child 8: 6,150, child 10: 4,860), remainders {48: 22,748}, global profile == the hub's.
usage: t74_control.py BC OUTDIR"""
import sys, os, json, gzip, copy, importlib.util
BC, OUT = sys.argv[1], sys.argv[2]
rd = lambda p: json.loads(gzip.decompress(open(p, 'rb').read())) if p.endswith('.gz') else json.load(open(p))
wz = lambda p, x: open(p, 'wb').write(gzip.compress(json.dumps(x, sort_keys=True, separators=(',', ':')).encode(), mtime=0))
spec = importlib.util.spec_from_file_location('make_cert', BC + '/k-combine/tools/make_cert.py')
mc = importlib.util.module_from_spec(spec); spec.loader.exec_module(mc)
os.makedirs(OUT + '/inputs', exist_ok=False)
plan = rd(BC + '/bh-gap/dl1839/PLAN.json')
order = [list(x) for x in plan['order']]
json.dump(order, open(OUT + '/inputs/ORDER.json', 'w'))
json.dump(rd(BC + '/bh-gap/dl1839/DISPATCH.json'), open(OUT + '/inputs/DISPATCH.json', 'w'))
wz(OUT + '/inputs/PLAN.json.gz', dict(pairs=plan['pairs'], sides=plan['sides']))
open(OUT + '/inputs/TEMPLATE.json.gz', 'wb').write(open(BC + '/rl-check/hubsrc/inputs-GLOBAL-TEMPLATE.json.gz', 'rb').read())
json.dump(plan['profile'], open(OUT + '/PROFILE.json', 'w'))
own = rd(BC + '/bh-gap/dl1839/OWNERS.json.gz')
own.sort(key=lambda o: (o['time'], o['helper']))
rows = [[o['time'], o['helper'], o['owner_replica'], o['child'], o['port'], o['offset']] for o in own]
wz(OUT + '/FILL.json.gz', dict(columns=['time', 'helper', 'owner_replica', 'child', 'port', 'offset'], rows=rows))
for prov, w in (('gauge44', 'T74-gauge44-566'), ('extra', 'T74-extra-566')):
    c = rd(BC + f'/bh-gap/data/{w}.gcert.json.gz')
    d = mc.partial_uncopy(c, set())
    v, R, h = d['v'], d['R'], d['h']
    dim = [len(f) for f in d['frames']]
    E = set(d['exits'])
    d['residuals'] = [(dim[d['final'][2 * v + q]] if 2 * v + q in E else h) - dim[d['start'][2 * v + q]] for q in range(R)]
    d.update(format='gcert/1b', copies=0, warps=[], levels={})
    wz(OUT + f'/{prov}.gcert1b.json.gz', d)
print('T74 control package written to', OUT, 'rows', len(rows))
