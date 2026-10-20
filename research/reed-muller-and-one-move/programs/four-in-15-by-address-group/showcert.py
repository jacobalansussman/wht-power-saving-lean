#!/usr/bin/env python3
"""showcert.py: print the nonzero dual rows of an mtlp.py dump in words (rationalised multipliers)."""
import sys, json
from fractions import Fraction as Fr
d = json.load(open(sys.argv[1])); m = int(sys.argv[2])
idx = [(i, j) for i in range(m) for j in range(i, m)]
def formstr(F):
    t = []
    for b, (i, j) in enumerate(idx):
        if (F >> b) & 1: t.append("x%d" % (i + 1) if i == j else "x%dx%d" % (i + 1, j + 1))
    return "+".join(t)
def dname(i): return format(i + 1, '0%db' % m)[::-1]
for (row, rhs, lab), y in zip(d['rows'], d['dual']):
    l = Fr(-y).limit_denominator(10 ** 6)
    if l == 0: continue
    k = lab[0]
    if k == 'cover': s = "cover row of F = %s" % formstr(int(lab[1]))
    elif k == 'del': s = "deletion row of G = %s" % formstr(int(lab[1]))
    elif k == 'piece<=slack': s = "cell (%s,%s) <= slack(%s)" % (dname(int(lab[1])), dname(int(lab[2])), formstr(int(lab[3])))
    else: s = "cell (%s,%s) <= c_%s" % (dname(int(lab[1])), dname(int(lab[2])), dname(int(lab[1])))
    print("  multiplier %-8s %s  (a = %s)" % (l, s, d['aF'][lab[1]] if k in ('cover', 'del') else (d['aF'][lab[3]] if k == 'piece<=slack' else '-')))
