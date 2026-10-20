#!/usr/bin/env python3
"""d-bit-exits / anatomy.py: where do PR589's bit helpers start and end, and what forces them to FULL?

Reads the tuned #589 word (DreamingOfClouds, on eumemic #555/#543). Replays D, A, scatter, B, C with nested frames
exactly as xcheck_p.py does, and records for every register its frame chain. Reports:
  - per helper: start dim, frame at the end of B (pre-cleanup), number of C gates writing it / reading it,
    the dims of its runs, the rank of its final climb into FULL;
  - which C gates on a helper could run below FULL (partners' frames at that moment);
  - a 'least-frame cleanup' replay: every C gate at the closure of the join of its operands' current frames,
    then sources forced to FULL, helpers left where they end (an exit at its own frame).
usage: python3 -I anatomy.py WORD.json.gz [out.json]
"""
import gzip, json, sys
from collections import Counter, defaultdict
from fractions import Fraction as Q

w = json.load(gzip.open(sys.argv[1], 'rt'))
V, R = w['v'], w['R']; N = 2 * V + R
L = [tuple(r) for r in w['lines']]
FR = w['frames']

def rank(vecs):
    B = {}
    for rr in vecs:
        row = list(map(Q, rr))
        for c in sorted(B):
            q = row[c]
            if q: row = [a - q * b for a, b in zip(row, B[c])]
        c = next((i for i, x in enumerate(row) if x), None)
        if c is not None:
            q = row[c]; B[c] = [x / q for x in row]
    return B

def inspan(B, vec):
    row = list(map(Q, vec))
    for c in sorted(B):
        q = row[c]
        if q: row = [a - q * b for a, b in zip(row, B[c])]
    return not any(row)

mask = [sum(1 << i for i in f) for f in FR]
dim = [len(rank([L[i] for i in f])) for f in FR]
FULL = next(k for k, d in enumerate(dim) if d == 8); ZERO = next(k for k, d in enumerate(dim) if d == 0)
key2frame = {m: k for k, m in enumerate(mask)}
newframes = []                       # frames not in the catalogue (closures of joins)

def closure_of(labels):
    B = rank([L[i] for i in labels])
    inside = [i for i in range(V) if inspan(B, L[i])]
    return sum(1 << i for i in inside), len(B)

def join(fs):
    m = 0
    for f in fs: m |= mask[f]
    if m in key2frame:
        return key2frame[m]
    cm, d = closure_of([i for i in range(V) if m >> i & 1])
    if cm in key2frame: return key2frame[cm]
    k = len(mask); mask.append(cm); dim.append(d); key2frame[cm] = k; newframes.append(k)
    return k

def le(a, b): return mask[a] & mask[b] == mask[a]

def kind(r): return 'x' if r < V else ('y' if r < 2 * V else 'h')

def replay(Cmode='full', force_src_full=False):
    cur = list(w['start']); val = [1 << i for i in range(N)]
    runs = defaultdict(list)          # register -> list of (section, from_dim, to_dim)
    H = Counter()
    endB = None
    def climb(r, f, sec):
        a = cur[r]
        if a == f: return
        assert le(a, f) and dim[f] > dim[a], (r, a, f)
        H[dim[f] - dim[a]] += 1; runs[r].append((sec, dim[a], dim[f])); cur[r] = f
    def gate(sec, f, t, ss):
        for r in [t, *ss]: climb(r, f, sec)
        for s in ss: val[t] ^= val[s]
    for f, t, ss in w['D']: gate('D', f, t, ss)
    for f, t, ss in w['A']: gate('A', f, t, ss)
    tops = [dim[cur[q]] for q in w['ret']]; snap = [val[q] for q in w['ret']]
    for t, ks in enumerate(w['scatter']):
        for k in ks: val[V + t] ^= snap[k]
    for f, t, ss in w['B']: gate('B', f, t, ss)
    endB = list(cur)
    for f, t, ss in w['C']:
        if Cmode == 'full': g = FULL
        else:
            g = join([cur[r] for r in [t, *ss]])
            if force_src_full and any(r < V for r in [t, *ss]): g = FULL
        gate('C', g, t, ss)
    for r in range(N):
        fin = w['final'][r]
        if r < V or V <= r < 2 * V: climb(r, fin, 'final')     # sources to FULL, targets to t-perp
        elif Cmode == 'full': climb(r, fin, 'final')
    for r in range(N):
        assert val[r] == ((1 << r) ^ (1 << (r - V)) if V <= r < 2 * V else 1 << r), r
    paid = H + Counter(tops)
    return dict(paid=dict(sorted(paid.items())), mass=sum(k * n for k, n in paid.items()), runs=runs, endB=endB,
                final=list(cur))

base = replay('full')
assert base['paid'] == {1: 1323, 2: 1342, 3: 79, 4: 44, 5: 102, 6: 11, 7: 19}, base['paid']
assert base['mass'] == 8 * R + 63 + 14 * V
out = dict(control=dict(paid=base['paid'], mass=base['mass']))

# per-helper anatomy
Cw = Counter(); Cr = Counter(); partners = defaultdict(set)
for f, t, ss in w['C']:
    Cw[t] += 1
    for s in ss: Cr[s] += 1
    for r in [t, *ss]:
        for q in [t, *ss]:
            if q != r: partners[r].add(q)
helpers = range(2 * V, N)
endBdim = Counter(dim[base['endB'][h]] for h in helpers)
lastclimb = Counter()
for h in helpers:
    rr = [x for x in base['runs'][h] if x[0] in ('C', 'final')]
    assert len(rr) <= 1
    lastclimb[(rr[0][2] - rr[0][1]) if rr else 0] += 1
startdim = Counter(dim[w['start'][h]] for h in helpers)
out['helpers'] = dict(
    count=R, start_dims=dict(startdim), end_of_B_dims=dict(sorted(endBdim.items())),
    rank_of_climb_into_FULL=dict(sorted(lastclimb.items())),
    C_writes_per_helper=dict(sorted(Counter(Cw[h] for h in helpers).items())),
    C_reads_of_helper=dict(sorted(Counter(Cr[h] for h in helpers).items())),
    helpers_with_C_partner_kinds=dict(Counter('+'.join(sorted({kind(q) for q in partners[h]})) for h in helpers)),
)
out['sources'] = dict(end_of_B_dims=dict(sorted(Counter(dim[base['endB'][s]] for s in range(V)).items())),
                      climbs_in_C=dict(sorted(Counter(sum(1 for x in base['runs'][s] if x[0] in ('C', 'final')) for s in range(V)).items())))
out['runs_by_class'] = {}
for cls, regs in (('x', range(V)), ('y', range(V, 2 * V)), ('h', helpers)):
    c = Counter(b - a for r in regs for (_, a, b) in base['runs'][r])
    out['runs_by_class'][cls] = dict(sorted(c.items()))

# least-frame cleanup variants
for name, fs in (('leastframe_all', False), ('leastframe_sources_full', True)):
    r = replay('least', fs)
    fd = Counter(dim[r['final'][h]] for h in helpers)
    out[name] = dict(paid=r['paid'], mass=r['mass'], mass_saving=base['mass'] - r['mass'],
                     helper_final_dims=dict(sorted(fd.items())),
                     source_C_climbs=dict(sorted(Counter(sum(1 for x in r['runs'][s] if x[0] in ('C', 'final')) for s in range(V)).items())),
                     helper_width_saving=sum(8 - dim[r['final'][h]] for h in helpers),
                     new_frames=len(newframes))
print(json.dumps(out, indent=1))
if len(sys.argv) > 2: json.dump(out, open(sys.argv[2], 'w'), indent=1)
