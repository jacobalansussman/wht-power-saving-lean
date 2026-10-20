#!/usr/bin/env python3
"""d-hpe-multigate: build the explicit U8 word with the chosen rank-7 end groups (frames only), refit, recount, save.
Usage: python3 -I make_word.py FULL-WORD.json.gz lazy_cones.json OUT.json.gz q1,q2,..."""
import sys, json, gzip, hashlib
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from build_ends import build
from refit import refit, recount, S
from cones_lib import basis
w = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read())); lz = json.load(open(sys.argv[2]))['helpers']
qs = [int(x) for x in sys.argv[4].split(',')]
w0, _ = build(w, {}); p0, wd0 = recount(w0)
w2, info = build(w, {q: tuple(lz[str(q)]['Ulazy']) for q in qs}); moves = refit(w2); p, wd = recount(w2)
raw = json.dumps(w2, separators=(',', ':')).encode(); open(sys.argv[3], 'wb').write(gzip.compress(raw, mtime=0))
diff = []
for sec in ('A', 'B'):
    for i, (a, b) in enumerate(zip(w[sec], w2[sec])):
        if basis(a[-1]) != basis(b[-1]): diff.append([sec, i, list(basis(a[-1])), list(basis(b[-1]))])
for k, (a, b) in enumerate(zip(w['early_cleanup'], w2['early_cleanup'])):
    if basis(a[1]) != basis(b[1]): diff.append(['early_cleanup', k, list(basis(a[1])), list(basis(b[1]))])
for q in qs: diff.append(['helper_final', q, list(basis(w['helper_final'][str(q)])), list(basis(w2['helper_final'][str(q)]))])
scal = all(a[:-1] == b[:-1] for sec in ('A', 'B', 'cut', 'precondition', 'cleanup') for a, b in zip(w[sec], w2[sec]))
out = dict(helpers=qs, refit_moves=moves, local_profile=p, dprofile={r: p.get(r, 0)-p0.get(r, 0) for r in sorted(set(p) | set(p0)) if p.get(r, 0) != p0.get(r, 0)},
           width_total=wd, dwidth=wd-wd0, dS_local=S(p)-S(p0), annotation_changes=diff, scalar_program_unchanged=scal,
           content_sha256=hashlib.sha256(raw).hexdigest())
print(json.dumps(out))
json.dump(out, open(sys.argv[3].replace('.json.gz', '.recipe.json'), 'w'), indent=1)
