#!/usr/bin/env python3
"""seeds15.py -- clock over seeds 1-15 (fpgaF.jsonl seeds 1-5 + fpgaF10.jsonl seeds 6-15) per tag."""
import json, statistics
v = {}
for f in ('fpgaF.jsonl', 'fpgaF10.jsonl'):
    for l in open(f):
        d = json.loads(l)
        if d.get('summary') or 'fmax_mhz' not in d: continue
        v.setdefault(d['tag'], {})[d['seed']] = (d['fmax_mhz'], d['dut_lc'])
for t, s in v.items():
    f = [x[0] for x in s.values()]
    print(t, 'seeds', sorted(s), 'cells', sorted({x[1] for x in s.values()}),
          'clock med %.1f [%.1f - %.1f]' % (statistics.median(f), min(f), max(f)))
