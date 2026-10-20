#!/usr/bin/env python3
"""soundness_classes.py: the three exact tests must PASS every known unit.  For one generator per fixed-point-free
class of order 3 and 5, run fixc.py, epat.py and rcover.py on the split of the plain schedules
1 in 4 (directions 1,2,4,8), 2 in 8 and 4 in 16 (the same, repeated) and on the slice-type schedule {7,11,13,14} x 4
(sum of z z^T = I)."""
import re, subprocess, sys
sys.path.insert(0, '.')
from forms import setup
D = setup(4); GL = D['GL']; ap = D['apply']
def run(cmd): return subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines()
cur = None
for line in open('out_placements.txt'):
    mo = re.match(r'order (\d+)', line)
    if mo: cur = int(mo.group(1)); continue
    mo = re.search(r'class size\s+(\d+)\s+gen rows \[([\d, ]+)\]\s+(.*)', line)
    if not mo or 'fixed-point-free' not in mo.group(3) or cur not in (3, 5): continue
    g = [int(x) for x in mo.group(2).split(',')]; H = 'elt:%d' % GL.index(g)
    grp = [[1, 2, 4, 8]]; h = g
    while h != [1, 2, 4, 8]: grp.append(h); h = D['matmul'](g, h)
    orb = {}; orbs = []
    for x in range(1, 16):
        if x in orb: continue
        o = sorted({ap(hh, x) for hh in grp})
        for y in o: orb[y] = len(orbs)
        orbs.append(o)
    out = []
    for W, base in ((1, [1, 2, 4, 8]), (2, [1, 2, 4, 8]), (4, [1, 2, 4, 8]), (4, [7, 11, 13, 14])):
        dl = base * W; comp = [0] * len(orbs)
        for z in dl: comp[orb[z]] += 1
        cs = ','.join(map(str, comp)); M = len(dl); args = ['--H', H, '--W', str(W), '--M', str(M)]
        r1 = run(['python3', 'fixc.py'] + args + ['--comp', cs, '--allpart', '4'])
        r2 = run(['python3', 'epat.py'] + args + ['--dirs', ','.join(map(str, dl))])
        r3 = run(['python3', 'rcover.py'] + args + ['--comps', cs])
        ok = ('passing: 1 of 1' in r1[-1]) and ('E-pattern: 1 of 1' in r2[-1]) and ('pass: 1 of 1' in r3[-1])
        out.append('%d in %d %s (%s): %s' % (W, M, base, cs, 'pass' if ok else 'FAIL ' + r1[-1] + ' | ' + r2[-1] + ' | ' + r3[-1]))
    print('order %d class %s:' % (cur, g), '; '.join(out), flush=True)
