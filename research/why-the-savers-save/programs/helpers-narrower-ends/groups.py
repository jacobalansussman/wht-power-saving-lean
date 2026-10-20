#!/usr/bin/env python3
"""d-hpe-multigate: every rank-8-ended U8 helper lowered to its least (lazy-cone) end, one group of helpers sharing that
end at a time, then all groups together; exact local profiles after frame refitting.
Usage: python3 -I groups.py FULL-WORD.json.gz lazy_cones.json out.json"""
import sys, json, gzip
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from collections import defaultdict
from build_ends import build
from refit import refit, recount, S
w = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read())); lz = json.load(open(sys.argv[2]))['helpers']
w0, _ = build(w, {}); p0, wd0 = recount(w0); S0 = S(p0)
grp = defaultdict(list)
for q, U in w['helper_final'].items():
    if len(U) == 8: grp[tuple(lz[q]['Ulazy'])].append(int(q))
res = []
for gi, (U, qs) in enumerate(sorted(grp.items(), key=lambda kv: min(kv[1]))):
    w2, info = build(w, {q: U for q in qs})
    if w2 is None: res.append(dict(group=gi, helpers=qs, infeasible=info)); continue
    refit(w2); p, wd = recount(w2)
    d = {r: p.get(r, 0)-p0.get(r, 0) for r in sorted(set(p) | set(p0)) if p.get(r, 0) != p0.get(r, 0)}
    res.append(dict(group=gi, helpers=qs, U=list(U), cold=[q for q in qs if str(q) not in w['start_extra']],
                    dprofile=d, dwidth=wd-wd0, dS=round(S(p)-S0, 4)))
    print(res[-1])
allq = {q: tuple(lz[str(q)]['Ulazy']) for qs in grp.values() for q in qs}
w2, info = build(w, allq); refit(w2); p, wd = recount(w2)
print('ALL 70:', {r: p.get(r, 0)-p0.get(r, 0) for r in sorted(set(p) | set(p0))}, 'dwidth', wd-wd0, 'dS', round(S(p)-S0, 3))
json.dump(dict(groups=res, all=dict(dprofile={r: p.get(r, 0)-p0.get(r, 0) for r in sorted(set(p) | set(p0))}, dwidth=wd-wd0, dS=S(p)-S0)), open(sys.argv[3], 'w'), indent=1)
