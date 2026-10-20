#!/usr/bin/env python3
"""whichclass.py: which placement class (out_placements.txt) contains e3-rmbuild's Singer powers (x^4 + x + 1,
polynomial basis) and the six fixed-point-free order-3 survivors."""
import re, itertools
exec(open('placements.py').read().split("want = ")[0])
classes = []; cur = None
for line in open('out_placements.txt'):
    mo = re.match(r'order (\d+)', line)
    if mo: cur = int(mo.group(1)); continue
    mo = re.search(r'class size\s+(\d+)\s+gen rows \[([\d, ]+)\]', line)
    if mo: classes.append((cur, int(mo.group(1)), [int(x) for x in mo.group(2).split(',')]))
def conjset(g):
    return {cyc(mm(mm(o, g), inv[key(o)])) for o in O4}
def singer_rows(p):
    cols = []
    for j in range(m):
        v = 1 << j
        for _ in range(p):
            v <<= 1
            if v >> m: v ^= 0b10011
        cols.append(v)
    return tr(cols)
tests = {'singer^1 (order 15)': singer_rows(1), 'singer^3 (order 5)': singer_rows(3), 'singer^5 (order 3)': singer_rows(5)}
for i, g in enumerate([[3, 1, 12, 4], [3, 1, 8, 12], [5, 10, 1, 2], [5, 8, 1, 10], [9, 6, 2, 1], [9, 4, 6, 1]]): tests['Z3fpf%d' % i] = g
for name, g in tests.items():
    cs = conjset(g)
    hit = [c for c in classes if cyc(c[2]) in cs]
    print(name, 'rows', g, '-> class', hit)
