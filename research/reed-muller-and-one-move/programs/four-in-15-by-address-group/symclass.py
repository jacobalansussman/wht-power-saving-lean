#!/usr/bin/env python3
"""symclass.py: which placement classes of cyclic H (order 3, 5, 15) contain a SYMMETRIC generator (g = g^T; then
every virtual target h h^T = h^2 lies in H itself, e.g. the Singer cycle written in a self-dual basis of F_16)."""
import re
exec(open('placements.py').read().split("want = ")[0])
classes = []; cur = None
for line in open('out_placements.txt'):
    mo = re.match(r'order (\d+)', line)
    if mo: cur = int(mo.group(1)); continue
    mo = re.search(r'class size\s+(\d+)\s+gen rows \[([\d, ]+)\]\s+(.*)', line)
    if mo: classes.append((cur, int(mo.group(1)), [int(x) for x in mo.group(2).split(',')], mo.group(3)))
for (n, sz, g, fx) in classes:
    if n not in (3, 5, 15): continue
    cs = {cyc(mm(mm(o, g), inv[key(o)])) for o in O4}
    nsym = 0; ex = None
    for s in cs:
        for h in s:
            if list(h) != I and list(h) == tr(list(h)): nsym += 1; ex = list(h); break
    print('order %2d class size %2d gen %s %-18s subgroups with a symmetric generator: %d %s' % (n, sz, g, fx, nsym, ex if ex else ''))
