# au-2: per-frame mismatch counts (and how many outputs are x) for a sim.py work folder
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen import model
rtl, work = sys.argv[1], sys.argv[2]
meta = json.load(open(rtl + '.json')); n = meta['n']; N = 1 << n; lat = meta['lat']; P = meta['P']; name = meta['name']
inp = [l.strip() for l in open(os.path.join(work, name + '_in.hex'))]
def s16(v): return v - 65536 if v & 0x8000 else v
sr = np.array([s16(int(h[:4], 16)) for h in inp]); si = np.array([s16(int(h[4:], 16)) for h in inp])
def _p(v):
    try: return int(v)
    except ValueError: return 999999
out = np.array([[_p(v) for v in ln.split()] for ln in open(os.path.join(work, name + '_out.txt')) if len(ln.split()) == 2], dtype=np.int64)
posmap = np.array([sum(1 << P[t] for t in range(n) if (q >> t) & 1) for q in range(N)])
f = 0
while (f * N + lat - 1 + N) <= len(out):
    mr, mi = model(sr[f * N:(f + 1) * N], si[f * N:(f + 1) * N], n)
    base = f * N + lat - 1
    hr = out[base:base + N, 0][np.argsort(posmap)]; hi = out[base:base + N, 1][np.argsort(posmap)]
    bad = (hr != np.array(mr)) | (hi != np.array(mi)); nx = int(((hr == 999999) | (hi == 999999)).sum())
    ex = ''
    if bad.any():
        k = int(np.argmax(bad)); ex = 'first pos %d got (%d,%d) want (%d,%d)' % (k, hr[k], hi[k], mr[k], mi[k])
    print('frame', f, 'mismatches', int(bad.sum()), 'x', nx, ex); f += 1
