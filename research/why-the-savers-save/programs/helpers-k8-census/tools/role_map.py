#!/usr/bin/env python3
"""k8-pairs: role map U8 -> k8 for the in-block helpers of given U8 blocks: a helper's role = (its in-block line and
plane visits as source sets, entrance dim, final dim, width) under the census climb rule.  Sources/targets keep their
canonical ids (ports identical).  Prints, per block, how many roles match one-to-one, and writes {block: {u8_id: k8_id}}.
Usage: python3 -I role_map.py U8CANON K8CANON U8_BLOCKS blocks(comma) out.json"""
import sys, json, gzip
from collections import Counter, defaultdict
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return sorted(bb, reverse=True)
def inspan(x, bb):
    for z in bb: x = min(x, x ^ z)
    return x == 0
def roles(path, blocks):
    c = json.loads(gzip.decompress(open(path, 'rb').read()))
    F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
    B = [basis(f) for f in F]; dims = [len(b) for b in B]
    spans = [basis([ports[s] for s in x]) for x in blocks]
    cur = list(c['start']); vis = defaultdict(list); runs = defaultdict(int)
    def climb(r, f):
        if f is None or cur[r] == f: return
        if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
        runs[r] += dims[f] - dims[cur[r]]; cur[r] = f; vis[r].append(f)
    for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
        if ph == 'centres':
            for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
            continue
        for e in c.get(ph) or []:
            climb(e[1], e[-1])
            if e[0] == 'add': climb(e[2], e[-1])
    for r in range(n): climb(r, c['final'][r])
    out = defaultdict(dict)
    for r in range(2*v, n):
        for bi, sp in enumerate(spans):
            fr = [tuple(s for s in range(v) if inspan(ports[s], B[f])) for f in vis[r] if dims[f] in (2, 3) and all(inspan(x, sp) for x in B[f])]
            if fr: out[bi][r] = (tuple(fr), dims[c['start'][r]], dims[c['final'][r]], runs[r])
    return out
blocks = json.load(open(sys.argv[3]))['blocks']; want = [int(x) for x in sys.argv[4].split(',')]
U = roles(sys.argv[1], blocks); K = roles(sys.argv[2], blocks)
res = {}
for b in want:
    ku = defaultdict(list); kk = defaultdict(list)
    for r, role in U[b].items(): ku[role].append(r)
    for r, role in K[b].items(): kk[role].append(r)
    m = {}; um = []; km = []
    for role in set(ku) | set(kk):
        a, z = sorted(ku.get(role, [])), sorted(kk.get(role, []))
        for x, y in zip(a, z): m[x] = y
        um += a[len(z):]; km += z[len(a):]
    # second pass: match leftovers on frames only (width/ends may differ by layers)
    fu = defaultdict(list); fk = defaultdict(list)
    for r in um: fu[U[b][r][0]].append(r)
    for r in km: fk[K[b][r][0]].append(r)
    m2 = {}
    for fr in set(fu) & set(fk):
        for x, y in zip(sorted(fu[fr]), sorted(fk[fr])): m2[x] = y
    print(f'block {b}: U8 in-block helpers {len(U[b])}, k8 {len(K[b])}; exact role matches {len(m)}, frame-only matches {len(m2)}, '
          f'unmatched U8 {len(um) - len(m2)}, unmatched k8 {len(km) - len(m2)}; ids equal in {sum(1 for x, y in {**m, **m2}.items() if x == y)}')
    res[b] = dict(exact={str(x): y for x, y in m.items()}, frames_only={str(x): y for x, y in m2.items()})
json.dump(res, open(sys.argv[5], 'w'), indent=0)
