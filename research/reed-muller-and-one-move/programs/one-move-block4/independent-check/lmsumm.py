#!/usr/bin/env python3
"""lmsumm.py: summary of my LM trial files (count, lowest costs, histogram of floors)."""
import json, sys, glob, collections
for fn in sys.argv[1:]:
    rs = [json.loads(l) for l in open(fn) if l.strip()]
    c = sorted(r['cost'] for r in rs)
    h = collections.Counter(round(x, 1) for x in c)
    print('%-28s trials %3d  lowest %s  floors %s  max|x| %.1f' % (fn.split('/')[-1], len(rs), [round(x, 4) for x in c[:4]],
          h.most_common(4), max([r['norm'] for r in rs] or [0])))
