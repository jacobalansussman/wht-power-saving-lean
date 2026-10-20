"""x-exit-words: in the bank route (helpers are bank scratch, never restored: Shea drops the cleanup; we also drop
the tail), let every helper (word helpers AND source copies) EXIT at its NATURAL frame = its frame after its last
gate in the lowered word (A, centres, P, B), instead of climbing to the hub's final frame.  Then the local histogram
loses those final climbs and the residual of helper q is dim(natural exit) - dim(entrance).  Data roles keep their
interface ends (x full, y dim h-1), centres keep their frames.  Prints the histogram, the residual multiset, checks
the deficit identity, writes a gcert/1 file with the natural exits and runs our gx checks (E5, E7, check1x)."""
import sys, os, json, gzip, importlib.util
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lower_exit as LE
wp, gxp, out = sys.argv[1], sys.argv[2], sys.argv[3]
spec = importlib.util.spec_from_file_location('gx', gxp); gx = importlib.util.module_from_spec(spec); spec.loader.exec_module(gx)
c0 = LE.read(wp)
cert, ok, total = LE.lower(c0, gx)
v, R, h = cert['v'], cert['R'], cert['h']; n = 2 * v + R
fr = [tuple(f) for f in cert['frames']]; dim = [len(f) for f in fr]
cur = list(cert['start'])
for ph in ('A', 'B'):
    if ph == 'B':
        pass
    for g in cert[ph]:
        if g[0] == 'out' and g[2] >= 2 * v and len(g[3]) == 0 and len(g[4]) == 0 and g in cert['B'][-n:]:
            continue
        for r in gx.gate_regs(g)[0]:
            if cur[r] != g[1]: cur[r] = g[1]
# the lowered B ends with one frame-only 'out' gate per register (its final frame): undo those for slots
nat = list(cur)
cur = list(cert['start'])
Bcore = cert['B'][:-n]   # strip the per-register final 'out' gates appended by the lowering
assert all(g[0] == 'out' and not g[3] and not g[4] for g in cert['B'][-n:])
for g in cert['A'] + Bcore:
    for r in gx.gate_regs(g)[0]: cur[r] = g[1]
for _, r, f in cert['ret']: pass
newfinal = list(cert['final'])
for r in range(2 * v, n): newfinal[r] = cur[r]
c2 = dict(cert); c2['final'] = newfinal; c2['B'] = Bcore + [['out', f, r, [], []] for r, f in enumerate(newfinal)]
# recompute blocks / N with the same climb replay as the lowering
H = {k: Counter() for k in ('x', 'y', 's', 'c')}; cur = list(c2['start'])
cls = lambda r: 'x' if r < v else 'y' if r < 2 * v else 's'
def climb(r, f):
    if cur[r] != f:
        assert dim[f] > dim[cur[r]] and gx.inside(fr[cur[r]], fr[f]); H[cls(r)][dim[f] - dim[cur[r]]] += 1; cur[r] = f
for ph, gates in (('A', c2['A']), ('B', c2['B'])):
    if ph == 'B':
        for _, r, f in c2['ret']: assert cur[r] == f; H['c'][dim[f]] += 1
    for g in gates:
        for r in gx.gate_regs(g)[0]: climb(r, g[1])
for r in range(n): climb(r, c2['final'][r])
c2['blocks'] = {k: {str(a): b for a, b in sorted(p.items())} for k, p in H.items()}
tot = Counter()
for p in H.values(): tot.update(p)
c2['N'] = sum(a * b for a, b in tot.items())
res_old = [dim[cert['final'][r]] - dim[cert['start'][r]] for r in range(2 * v, n)]
res_new = [dim[newfinal[r]] - dim[cert['start'][r]] for r in range(2 * v, n)]
print('slots', n - 2 * v, 'residual sum: hub exits', sum(res_old), '-> natural exits', sum(res_new))
print('residuals natural', dict(sorted(Counter(res_new).items())))
print('histogram natural (incl. copies, centres)', dict(sorted(tot.items())), 'N', c2['N'])
# deficit identity of the route: 72 (N - residual) must equal the reference
print('N - residual: hub-exit lowering', cert['N'] - sum(res_old), ' natural', c2['N'] - sum(res_new))
chk = {}
for lab, f in (('scal_E5', lambda: gx.scal(c2)), ('q2_E7', lambda: gx.q2(c2)),
               ('check1x_without_E4', lambda: LE.check1x_from(gx, drop_e4=True)(c2, scalar=False))):
    try: f(); chk[lab] = True
    except AssertionError as e: chk[lab] = str(e)[:200]
print('checks', chk)
open(out, 'wb').write(gzip.compress(json.dumps(c2, sort_keys=True).encode(), mtime=0))
json.dump(dict(histogram=dict(sorted(tot.items())), residuals=res_new, checks=chk), open(out + '.info.json', 'w'))
