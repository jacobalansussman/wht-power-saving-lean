"""d-halfread BIT SIDE: does any helper of #589's bit word (PR #555's punctured R439 on #543's R420 L3 skeleton, with
#589's two transports; file research/e8-bit-p111-tuned/tune/bitE8mut-p111-tuned.json.gz) exist only to give a second
reader the same row?  Replays the word over F2 exactly as the hub's wlib.replay does (sections D, A, scatter, B, C) and,
for every helper read, records the clean row (source columns) it delivers. A helper is 'copy-like' if at every one of
its reads some OTHER register holds the identical clean row at that moment.
Usage: python3 -I bit_copies.py WORD.json.gz out.json"""
import sys, json, gzip
from collections import Counter, defaultdict
w = json.load(gzip.open(sys.argv[1], 'rt'))
V = w['v']; R = w['R']; N = 2*V + R
SRC = (1 << V) - 1
val = [1 << i for i in range(N)]
deliveries = defaultdict(list)  # helper -> [(time, reader, clean row)]
holders_at = {}
T = 0
cls = lambda q: 'x' if q < V else ('y' if q < 2*V else 'h')
def snapshot_rows():
    d = defaultdict(set)
    for q in range(N):
        c = val[q] & SRC
        if c: d[c].add(q)
    return d
copy_like = Counter(); examples = []
pending = []
def gates(sec):
    global T
    for f, t, ss in w[sec]:
        rows = None
        for s in ss:
            if 2*V <= s:
                c = val[s] & SRC
                if c:
                    if rows is None: rows = snapshot_rows()
                    others = rows[c] - {s}
                    deliveries[s].append((sec, T, t, cls(t), bin(c).count('1'), sorted(others)[:5]))
        for s in ss: val[t] ^= val[s]
        T += 1
gates('D'); gates('A')
snap = [val[q] for q in w['ret']]
for t, ks in enumerate(w['scatter']):
    for k in ks: val[V+t] ^= snap[k]
gates('B'); gates('C')
ok = all(val[r] == ((1 << r) ^ (1 << (r-V)) if V <= r < 2*V else 1 << r) for r in range(N))
summary = Counter(); flagged = []
for h in range(2*V, N):
    d = [x for x in deliveries[h] if x[0] != 'C']
    if not d: summary['no clean reads outside C'] += 1; continue
    every = all(x[5] for x in d)
    some = any(x[5] for x in d)
    tg = sorted({x[2] for x in d})
    summary['every clean read duplicated elsewhere' if every else ('some duplicated' if some else 'no read duplicated')] += 1
    if every: flagged.append(dict(helper=h, reads=d[:8]))
print('endpoint exact:', ok)
print(dict(summary))
for f in flagged[:12]: print(f)
json.dump(dict(endpoint=ok, summary=dict(summary), flagged=flagged), open(sys.argv[2], 'w'), indent=0)
