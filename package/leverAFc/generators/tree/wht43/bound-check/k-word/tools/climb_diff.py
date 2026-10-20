"""Where does the gcert-lowering histogram differ from the hub's local profile?  Replays the climbs of the lowered
word (A, cut, P+B, final) and of the same with the hub's cleanup and tail appended (as charged by the hub)."""
import sys, os, json, gzip, importlib.util
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lower_exit as LE
spec = importlib.util.spec_from_file_location('gx', sys.argv[2]); gx = importlib.util.module_from_spec(spec); spec.loader.exec_module(gx)
c0 = LE.read(sys.argv[1])
cert, ok, total = LE.lower(c0, gx)
print('lowered (copies) ', dict(sorted(total.items())))
# hub-style: same word, with tail and cleanup gates (remapped) climbing too, before the final climb
c = json.loads(json.dumps(c0)); v, R = c['v'], c['R']; base = 2 * v + R
remap = lambda r: base + r if r < v else r
dims = [len(f) for f in c['frames']]
cur = c['start'] + [0] * v; fin = c['final'] + [1] * v
H = Counter(); per = {}
def climb(r, f, tag):
    if cur[r] != f and f is not None:
        d = dims[f] - dims[cur[r]]; assert d > 0, (r, tag)
        H[d] += 1; per.setdefault(r, []).append((tag, d)); cur[r] = f
for r in range(v): climb(base + r, c['start'][r], 'copy')
order = [('A', c['A']), ('centres', [['out', f, remap(r)] for r, f in zip(c['centers'], c['center_frames'])])]
order += [('cut', c['cut']), ('P', c['precondition']), ('B', c['B'])]
if len(sys.argv) > 3: order += [('cleanup', c['cleanup']), ('tail', c.get('tail') or [])]
for tag, gates in order:
    for e in gates:
        if tag == 'centres': climb(e[2], e[1], tag); continue
        regs = [remap(e[1])] + ([remap(e[2])] if e[0] == 'add' else [])
        if tag in ('cut',): continue
        for r in regs: climb(r, e[-1], tag)
for r in range(base + v): climb(r, fin[r], 'final')
print('replay', 'with cleanup+tail' if len(sys.argv) > 3 else 'without', dict(sorted(H.items())))
prof = Counter({int(k): n for k, n in c0['local_profile'].items()}); prof[1] += v; prof[8] += v
print('hub profile+copies', dict(sorted(prof.items())))
