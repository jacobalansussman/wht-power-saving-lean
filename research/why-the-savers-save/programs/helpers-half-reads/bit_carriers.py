"""d-halfread BIT SIDE, part 2: the 7 copy-like helpers of the bit word are single-source carriers (extra carriers
of #555: a helper that reads a source of block B at its line frame in phase A for a D4 mini-node).  The rank-1
analogue of a half read: the reading mini-node helper takes the row straight from the source register while that
source still sits at its line frame, i.e. the reader adds a rank-1 stop at the source's line frame.  Tests, for every
read of each carrier: reader's frame before inside the source's line frame, reader's next frame containing it, the
reader idle in between, the source still at its line frame and holding exactly its own row.
Usage: python3 -I bit_carriers.py WORD.json.gz bit_copies.json"""
import sys, json, gzip
from collections import Counter, defaultdict
w = json.load(gzip.open(sys.argv[1], 'rt')); flagged = json.load(open(sys.argv[2]))['flagged']
V = w['v']; R = w['R']; N = 2*V + R
lines = w['lines']
mask = [sum(1 << i for i in f) for f in w['frames']]
def le(a, b): return mask[a] & mask[b] == mask[a]
from fractions import Fraction as Q
def rank(f):
    rows = [list(map(Q, lines[i])) for i in w['frames'][f]]; r = 0; m = [row[:] for row in rows]
    cols = len(lines[0]); piv = 0
    for c in range(cols):
        p = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if p is None: continue
        m[r], m[p] = m[p], m[r]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                q = m[i][c]/m[r][c]; m[i] = [a - q*b for a, b in zip(m[i], m[r])]
        r += 1
    return r
dim = [rank(f) for f in range(len(w['frames']))]
cur = list(w['start']); chain = defaultdict(list); val = [1 << i for i in range(N)]
for q in range(N): chain[q].append((-1, cur[q]))
T = 0; events = defaultdict(list); valat = {}
def gates(sec):
    global T
    for f, t, ss in w[sec]:
        for r in [t, *ss]:
            if cur[r] != f: cur[r] = f; chain[r].append((T, f))
        for r in [t, *ss]: events[r].append(T)
        for s in ss: valat[(T, s)] = val[s]
        for s in ss: val[t] ^= val[s]
        T += 1
gates('D'); gates('A')
snap = [val[q] for q in w['ret']]
for t, ks in enumerate(w['scatter']):
    for k in ks: val[V+t] ^= snap[k]
gates('B'); gates('C')
def frame_at(q, t):
    f = chain[q][0][1]
    for tt, ff in chain[q]:
        if tt <= t: f = ff
    return f
def next_after(q, t):
    for tt, ff in chain[q]:
        if tt > t: return tt, ff
    return None, None
stat = Counter(); rows = []
for F in flagged:
    h = F['helper']
    for sec, t, reader, cl, nterm, others in F['reads']:
        srcs = [o for o in others if o < V]
        if not srcs: stat['duplicate held by a non-source'] += 1; continue
        s = srcs[0]; Ls = w['start'][s]
        # the carrier read the source at some earlier time; the reader would read s directly before time t
        f0 = frame_at(reader, t - 1); tn, fn = next_after(reader, t)
        last_busy = max([e for e in events[reader] if e < t] or [-1])
        src_at_line = frame_at(s, t - 1) == Ls
        src_clean = val is not None
        Sc = frame_at(s, t - 1)
        conds = [("reader frame before inside the source's current frame", le(f0, Sc)), ("reader next frame contains the source's current frame", fn is not None and le(Sc, fn)),
                 ("reader frame before inside the source line frame", le(f0, Ls)),
                 ('reader next frame contains the line frame', fn is not None and le(Ls, fn)),
                 ('source still at its line frame at the read', src_at_line)]
        first = next((nm for nm, ok in conds if not ok), 'all frame conditions pass')
        stat[first] += 1
        rows.append(dict(carrier=h, reader=reader, source=s, time=t, reader_frame_dim_before=dim[f0],
                         reader_next_dim=dim[fn] if fn is not None else None, source_frame_dim=dim[frame_at(s, t-1)], verdict=first))
for x, y in stat.items(): print(y, x)
for r in rows: print(r)
# ---- option A: the reader reads the source at the carrier-read gate's own frame G (the source climbs to G)
print('\noption A: source read directly at the gate frame (source climbs to G)')
gate_frame = {}
T2 = 0
for sec in ('D', 'A', 'B', 'C'):
    for f, t, ss in w[sec]:
        gate_frame[T2] = f; T2 += 1
for row in rows:
    s, t = row['source'], row['time']; G = gate_frame[t]
    sf = frame_at(s, t - 1)
    later = [ff for tt, ff in chain[s] if tt > t]
    ok1 = le(sf, G); ok2 = all(le(G, ff) for ff in later)
    print(row['carrier'], 'reader', row['reader'], 'source', s, 'G dim', dim[G], 'source frame inside G', ok1, "source's later frames contain G", ok2,
          '(source later dims %s)' % [dim[ff] for ff in later][:4])
