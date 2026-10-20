"""k-lower: design the zero-check batches so that every packed replay passes (exact mirror of GS.zCheck on the y-free
twin; the binding limit is the 64-bit tag of Sussman's packed vectors, sum_T digit_T tau(T) < 2^63):
  * exits split into LIGHT (bit batches over the light list, prefix of k x registers) and HEAVY (one batch each:
    prefix of k x registers then the exit), every heavy exit separately;
  * non-exit sources (x and non-exit slots) packed greedily (first fit, heaviest first) into chunks of <= v registers.
Rows: all slots for exit batches, the exit list for chunks.
usage: zdesign.py <twin.gcert.json.gz> <out.json> [sw]"""
import sys, json, gzip, os
from fractions import Fraction as Q
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zmirror as zm
from zplan import zCheck
c0 = json.loads(gzip.open(sys.argv[1], 'rt').read())
sw = int(sys.argv[3]) if len(sys.argv) > 3 else 36
c = zm.lean_data(c0)
v, R = c0['v'], c0['R']; nr = 2 * v + R
P = (10, sw, 32, 32, 1)
el = sorted(c0['exits']); E = set(el)
slots = list(range(2 * v, nr))
ok = lambda gl, rows: zm.evaluate(zCheck, c, P, gl, rows)['value']
# column weights (max accumulation x 32), as zheavy
Pp = [dict() for _ in range(nr)]; Nn = [dict() for _ in range(nr)]
for r in range(nr):
    if not v <= r < 2 * v: Pp[r][r] = Q(1)
mx = defaultdict(lambda: Q(0))
def add(t, s, a):
    pos, neg = (Pp[s], Nn[s]) if a >= 0 else (Nn[s], Pp[s])
    for dst, sr in ((Pp[t], pos), (Nn[t], neg)):
        for key, value in sr.items():
            z = dst.get(key, Q(0)) + abs(a) * value; dst[key] = z
            if z > mx[key]: mx[key] = z
for g in c0['A'] + c0['B']:
    for t, a, b in g[3]: add(t, g[2], Q(a, b))
wt = {k: float(m * 32) for k, m in mx.items()}
def bitB(l, b, val): return [e for e in l if ((l.index(e) >> b) & 1) == val]
def with_prefix(gl, rows, kmax=None):
    for k in range(0, (kmax if kmax is not None else v - len(gl)) + 1):
        g2 = list(range(v, v + k)) + gl
        if len(g2) > v: break
        if ok(g2, rows): return k
    return None
# exits: raise the heavy threshold until all light bit batches pass
thr_list = [2 ** e for e in (40, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20)]
# NOTE: fillers are y registers v .. v+k-1 (never read or written in the y-free twin; never rows here)
design = None
for thr in thr_list:
    heavy = sorted(e for e in el if wt.get(e, 0) > thr)
    light = [e for e in el if e not in heavy]
    nb = max(1, (len(light) - 1).bit_length())
    bits = []
    good = True
    for b in range(nb):
        for val in (0, 1):
            gl = bitB(light, b, val)
            k = with_prefix(gl, slots)
            if k is None: good = False; break
            bits.append([b, val, k])
        if not good: break
    print('thr 2^%d heavy %d light %d bits ok %s' % (thr.bit_length() - 1, len(heavy), len(light), good), flush=True)
    if good:
        design = dict(thr=thr, heavy=heavy, light=light, nb=nb, bits=bits); break
assert design, 'no threshold works'
hb = []
for e in design['heavy']:
    k = with_prefix([e], slots)
    assert k is not None, ('heavy exit', e)
    hb.append([e, k])
design['heavy_batches'] = hb
print('heavy prefixes', hb, flush=True)
# chunks: first fit, heaviest first
nel = sorted((r for r in range(nr) if (r < v or 2 * v <= r) and r not in E), key=lambda r: -wt.get(r, 0))
chunks = []
for r in nel:
    placed = False
    for ch in chunks:
        if len(ch) < v and ok(sorted(ch + [r]), el):
            ch.append(r); placed = True; break
    if not placed:
        k = with_prefix([r], el)
        assert k is not None, ('single chunk fails', r)
        chunks.append(list(range(v, v + k)) + [r])
chunks = [sorted(ch) for ch in chunks]
design['chunks'] = chunks
print('chunks', len(chunks), [len(ch) for ch in chunks], flush=True)
# final verification and margin
allb = [(list(range(v, v + k)) + bitB(design['light'], b, val), slots) for b, val, k in design['bits']] + \
       [(list(range(v, v + k)) + [e], slots) for e, k in hb] + [(ch, el) for ch in chunks]
design['verified'] = all(ok(g, r) for g, r in allb)
design['spar'] = list(P)
design['n_batches'] = len(allb)
json.dump(design, open(sys.argv[2], 'w'))
print('verified', design['verified'], 'batches', len(allb))
