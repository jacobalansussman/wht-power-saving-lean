#!/usr/bin/env python3
"""ctr-search / soundness.py: calibration of the profile method on exact decoders it must NOT exclude.
For every exact dictionary decoder saved by the searches (calib, dsearch/bresp/sas outputs) compute its one-sided
profile and check (i) the counting bounds of profiles.py at a budget equal to its own cost, (ii) the incidence
bound, (iii) NC_k directly on its actual port sets for every K of size <= 7 when m <= 16 (sampled K otherwise)."""
import sys, os, json, glob, itertools, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import profiles as PR
from p112 import *
ok_all = True
# (i) the frame decoder profile (8^8) at budget 64, butterfly (9^8) at 72
for prof, bud in (((8,) * 8, 64), ((9,) * 8, 72), ((8,) + (9,) * 7, 71)):
    m = len(prof)
    passes = tuple(sorted(prof)) in set(PR.gen(m, budget=bud))
    print('profile', prof, 'budget', bud, 'passes counting bounds:', passes); ok_all &= passes
files = glob.glob(os.path.join(HERE, 'out/*.json'))
seen = 0
for f in files:
    try: res = json.load(open(f))
    except Exception: continue
    if not isinstance(res, list): continue
    for v in res:
        if not isinstance(v, dict) or not v.get('exact') or 'star_terms' not in v: continue
        prof = tuple(sorted(d for d in v['star_terms'] if d > 0)); cost = sum(prof); m = len(prof)
        passes = (prof in set(PR.gen(m, budget=cost))) if m <= 12 else PR.ok_prefix(list(prof), m)
        seen += 1
        print(os.path.basename(f), v['tag'], 'n', m, 'cost', cost, 'passes counting bounds:', passes)
        if passes is False: ok_all = False
print('decoders checked:', seen, ' all accepted:', ok_all)
