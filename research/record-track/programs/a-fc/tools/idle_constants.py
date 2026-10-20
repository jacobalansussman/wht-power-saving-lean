#!/usr/bin/env python3
"""fc-reopt: the idle histograms of a #3945-line design at the #3815-stage widths (old-FILL, old words) and the final
widths (FILL, final words), and the number of bins whose used width changes, computed with our own code (widths =
dim helper_final - dim start_extra over GF(2), as in o-model's offsets3945.py). usage: idle_constants.py OWNER_DIR"""
import json, gzip, sys
from collections import Counter
from pathlib import Path
R = Path(sys.argv[1]); I = R / 'inputs'
rd = lambda p: json.loads(gzip.decompress(p.read_bytes()) if p.suffix == '.gz' else p.read_bytes())
def rank(rows):
    piv = {}
    for x in rows:
        while x:
            b = x.bit_length()
            if b in piv: x ^= piv[b]
            else: piv[b] = x; break
    return len(piv)
order = rd(I / 'ORDER.json'); dispatch = {(s, r): p for s, r, p in rd(I / 'DISPATCH.json')}
full = [1 << i for i in range(9)]
def widths(prefix):
    out = {}
    for tag, pid in (('base', 'gauge44'), ('extra', 'extra')):
        word = rd(I / f'{prefix}WORD-{tag}.json.gz'); origin = rd(I / f'old-NATIVE-{tag}.json.gz')['helper_origin']
        out[pid] = [rank(word['helper_final'].get(str(r), full)) - rank(word['start_extra'].get(str(r), [])) for r in origin]
    return out
caps = {c['child']: c['rank'] for c in rd(I / 'TEMPLATE.json.gz')['idle_children']}
def walk(fill, wd):
    used = Counter()
    for t, h, o, c, p, z in fill['rows']:
        used[o, c, p] += wd[dispatch[tuple(order[t])]][h]
    idle = Counter()
    for o in range(360):
        for c, cap in caps.items():
            for p in range(112):
                l = cap - used.get((o, c, p), 0)
                assert l >= 0
                if l: idle[l] += 1
    return used, idle
ou, oi = walk(rd(I / 'old-FILL.json.gz'), widths('old-')); nu, ni = walk(rd(R / 'FILL.json.gz'), widths(''))
changed = sum(1 for k in set(ou) | set(nu) if ou.get(k, 0) != nu.get(k, 0))
print(json.dumps(dict(old_idle=dict(sorted(oi.items())), new_idle=dict(sorted(ni.items())), changed=changed)))
