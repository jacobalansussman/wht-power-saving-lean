#!/usr/bin/env python3
"""fixc_classes.py: fixc.py (deletion + rank cover, exact, all orbit splits) on one generator per placement class of the
fixed-point-free cyclic groups of order 3 and 5 (out_placements.txt), shapes given on the command line (W:M,...)."""
import re, subprocess, sys
sys.path.insert(0, '.')
from forms import setup
GL = setup(4)['GL']
shapes = [tuple(map(int, s.split(':'))) for s in sys.argv[1].split(',')]
cur = None
for line in open('out_placements.txt'):
    mo = re.match(r'order (\d+)', line)
    if mo: cur = int(mo.group(1)); continue
    mo = re.search(r'class size\s+(\d+)\s+gen rows \[([\d, ]+)\]\s+(.*)', line)
    if not mo or 'fixed-point-free' not in mo.group(3) or cur not in (3, 5): continue
    g = [int(x) for x in mo.group(2).split(',')]; idx = GL.index(g)
    for W, M in shapes:
        r = subprocess.run(['python3', 'fixc.py', '--H', 'elt:%d' % idx, '--W', str(W), '--M', str(M), '--all', '--allpart', '4'],
                           capture_output=True, text=True).stdout.strip().splitlines()
        print('order %d class size %2d gen %s  %d in %d: %s' % (cur, int(mo.group(1)), g, W, M, r[-1] if r else '?'), flush=True)
