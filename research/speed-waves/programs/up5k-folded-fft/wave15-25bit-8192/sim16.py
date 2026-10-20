#!/usr/bin/env python3
"""sim16.py -- av-5: check the scaled 16-bit floor baseline (f16_dsp6rr_8k.v, as-6's ports) against gen24.model with
W=16, pre=0 (as-6's arithmetic at any N; checked equal to as-6's own ffold.model at 4,096 first)."""
import sys, os, subprocess, numpy as np
sys.path.insert(0, '.')
sys.path.insert(1, '..')
import gen24, ffold
rng = np.random.default_rng(3)
for t in range(3):
    xr = rng.integers(-32768, 32768, 4096); xi = rng.integers(-32768, 32768, 4096)
    a = ffold.model(xr, xi); b = gen24.model(xr, xi, 12, pre=0, W=16)
    assert (a[0] == b[0]).all() and (a[1] == b[1]).all()
print('gen24.model(W=16, pre=0) == as-6 ffold.model at 4,096 on 3 random frames')
L = 13; N = 1 << L; nf = int(sys.argv[1]); skip = int(sys.argv[2]); work = sys.argv[3]; os.makedirs(work, exist_ok=True)
F = gen24.frames(nf + skip, L)[skip:]
with open(os.path.join(work, 'in.hex'), 'w') as f:
    for xr, xi in F:
        for p, q in zip(xr, xi):
            f.write('%04x%04x\n' % (int(p) & 0xFFFF, int(q) & 0xFFFF))
tb = ffold.TB.replace('reg [11:0] ad', 'reg [12:0] ad').replace('4096', str(N)) % dict(top='ffold', nw=nf * N - 1, nf=nf)
open(os.path.join(work, 'tb.v'), 'w').write(tb)
exe = os.path.join(work, 'a.out')
r = subprocess.run([gen24.BIN + '/iverilog', '-g2012', '-o', exe, '-DNO_ICE40_DEFAULT_ASSIGNMENTS', os.path.abspath('f16_dsp6rr_8k.v'), os.path.join(work, 'tb.v'), gen24.CELLS], capture_output=True, text=True, cwd=work)
if r.returncode: print(r.stderr[-2000:]); sys.exit(1)
subprocess.run([gen24.BIN + '/vvp', '-n', exe], capture_output=True, text=True, cwd=work)
lines = open(os.path.join(work, 'out.txt')).read().split('\n'); pos = 0; bad = 0; cyc = []
for xr, xi in F:
    cyc.append(int(lines[pos].split()[1])); pos += 1
    mr, mi = gen24.model(xr, xi, L, pre=0, W=16)
    for u in range(N):
        v = lines[pos + u]
        if 'x' in v or int(v[:4], 16) != (int(mr[u]) & 0xFFFF) or int(v[4:], 16) != (int(mi[u]) & 0xFFFF): bad += 1
    pos += N
print('f16_dsp6rr_8k: frames %d (skip %d) samples %d mismatching %d clocks %s' % (nf, skip, nf * N, bad, sorted(set(cyc))))
print('RESULT', 'PASS' if bad == 0 else 'FAIL')
