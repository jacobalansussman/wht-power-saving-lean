#!/usr/bin/env python3
"""chain_classes.py (w8-rmexact): for one generator per fixed-point-free placement class of order 3 and 5 and a shape
W in M, run the three exact tests on EVERY split of the M moves over the direction orbits:
fixc.py (deletion + rank cover with all set partitions of 4-point supports), then epat.py (forced zeros + matching)
on the fixc survivors, then rcover.py (rectangle partitions) on the epat survivors.  Prints the survivors per class."""
import re, subprocess, sys, ast
sys.path.insert(0, '.')
from forms import setup
GL = setup(4)['GL']
W, M = map(int, sys.argv[1].split(':'))
orders = [int(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [3, 5]
def run(cmd): return subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines()
def comps_of(lines):
    out = []
    for l in lines:
        mo = re.match(r'^(\([\d, ]+\)) \(\'pass\'', l)
        if mo: out.append(ast.literal_eval(mo.group(1)))
    return out
cur = None
for line in open('out_placements.txt'):
    mo = re.match(r'order (\d+)', line)
    if mo: cur = int(mo.group(1)); continue
    mo = re.search(r'class size\s+(\d+)\s+gen rows \[([\d, ]+)\]\s+(.*)', line)
    if not mo or 'fixed-point-free' not in mo.group(3) or cur not in orders: continue
    g = [int(x) for x in mo.group(2).split(',')]; H = 'elt:%d' % GL.index(g)
    base = ['--H', H, '--W', str(W), '--M', str(M)]
    L1 = run(['python3', 'fixc.py'] + base + ['--all', '--allpart', '4'])
    c1 = comps_of(L1); tot = re.search(r'of (\d+)', L1[-1]).group(1)
    c2 = []
    if c1:
        fn = 'tmp_comps_%d_%s.txt' % (cur, H.replace(':', ''))
        open(fn, 'w').write(''.join(','.join(map(str, c)) + '\n' for c in c1))
        c2 = comps_of(run(['python3', 'epat.py'] + base + ['--comps', fn]))
    c3 = []
    if c2:
        c3 = comps_of(run(['python3', 'rcover.py'] + base + ['--comps', ';'.join(','.join(map(str, c)) for c in c2)]))
    print('order %d class size %2d gen %s  %d in %d: splits %s; after rank cover %d; after E-pattern %d; after rectangles %d: %s'
          % (cur, int(mo.group(1)), g, W, M, tot, len(c1), len(c2), len(c3), c3), flush=True)
