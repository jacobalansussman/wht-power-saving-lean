#!/usr/bin/env python3
"""cx-mirror / chaincost2.py: chaincost.py's first-order cleaning cost for the best twisted mirror of the EXHAUSTIVE
run (job-twistall/remote-out/*/twistall.json; A = I + M, M given by its column images)."""
import os, json, glob, importlib.util
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('cc', os.path.join(here, 'chaincost.py'))
src = open(os.path.join(here, 'chaincost.py')).read().split("out = [cost_case(I9, 'natural mirror')]")[0]
ns = {'__file__': os.path.join(here, 'chaincost.py'), '__name__': 'cc'}
exec(compile(src, 'chaincost.py', 'exec'), ns)
f = sorted(glob.glob(os.path.join(here, 'all-twists', 'twistall.json')))[0]
r = json.load(open(f)); M = r['best']['M']
A = [(1 << i) ^ M[i] for i in range(9)]
res = ns['cost_case'](A, 'best twisted mirror (all)')
res['source'] = os.path.relpath(f, here)
json.dump(res, open(os.path.join(here, 'chaincost2.json'), 'w'), indent=1)
