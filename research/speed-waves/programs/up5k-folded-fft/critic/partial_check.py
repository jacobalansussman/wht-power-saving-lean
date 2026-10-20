"""partial_check.py WORKDIR SKIP NF -- au-3: check the COMPLETE frames already written to WORKDIR/out.txt by a running
ffold.py simulation (same comparison as ffold.sim: every word against the bit-true model), without waiting for the end."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ffold
work, skip, nf = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
F = ffold.frames(nf + skip)[skip:]
lines = open(os.path.join(work, 'out.txt')).read().split('\n')
pos = 0; done = 0; bad = 0
for fi, (xr, xi) in enumerate(F):
    if pos + 1 + ffold.N > len(lines) - 1: break
    cyc = int(lines[pos].split()[1]); pos += 1
    mr, mi = ffold.model(xr, xi); got = lines[pos:pos + ffold.N]; pos += ffold.N
    nb = sum(1 for u in range(ffold.N) if ('x' in got[u] or int(got[u][:4], 16) != (int(mr[u]) & 0xFFFF) or int(got[u][4:], 16) != (int(mi[u]) & 0xFFFF)))
    bad += nb; done += 1
    print('frame %d (global %d): clocks %d mismatches %d' % (fi, fi + skip, cyc, nb))
print('complete frames %d of %d, samples %d, mismatching %d' % (done, nf, done * ffold.N, bad))
