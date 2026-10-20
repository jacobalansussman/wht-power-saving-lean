#!/usr/bin/env python3
"""k8-pairs: compare the rank-3 frames (as source 4-sets) visited by two words under the census climb rule, with
visitors by class; relative to the U8 block partition (d-pair-market u8_blocks.json).
Usage: python3 -I planes_cmp.py WORD_U8 WORD_K8 U8_BLOCKS out.json"""
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
def visits(path):
    c = json.loads(gzip.decompress(open(path, 'rb').read()))
    F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
    B = [basis(f) for f in F]; dims = [len(b) for b in B]
    cur = list(c['start']); vis = defaultdict(list)
    for r in range(n): vis[r].append(cur[r])
    def climb(r, f):
        if f is None or cur[r] == f: return
        if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
        cur[r] = f; vis[r].append(f)
    for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
        if ph == 'centres':
            for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
            continue
        for e in c.get(ph) or []:
            climb(e[1], e[-1])
            if e[0] == 'add': climb(e[2], e[-1])
    for r in range(n): climb(r, c['final'][r])
    cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
    P = defaultdict(list)
    for r in range(n):
        for f in vis[r]:
            if dims[f] == 3:
                S = tuple(s for s in range(v) if inspan(ports[s], B[f]))
                P[S].append((cls(r), r))
    return P, c
U, cu = visits(sys.argv[1]); K, ck = visits(sys.argv[2])
blocks = json.load(open(sys.argv[3]))['blocks']; blk = {s: i for i, x in enumerate(blocks) for s in x}
def kind(S): return 'in-block' if len({blk[s] for s in S}) == 1 else 'cross'
print('U8 planes', len(U), Counter(kind(S) for S in U), ' k8 planes', len(K), Counter(kind(S) for S in K))
only_u = [S for S in U if S not in K]; only_k = [S for S in K if S not in U]
print('planes only in U8:', len(only_u), ' only in k8:', len(only_k))
for S in only_k: print('  k8-only plane', S, 'blocks', [blk[s] for s in S], 'visitors', K[S])
for S in only_u: print('  U8-only plane', S, 'blocks', [blk[s] for s in S], 'visitors', U[S])
per = Counter(); 
for S, vv in K.items():
    if len(vv) != len(U.get(S, [])): per[(kind(S), blk[S[0]] if kind(S)=='in-block' else None, len(U.get(S, [])), len(vv))] += 1
print('k8 planes whose visitor count differs from U8 (kind, block, U8 count, k8 count):', sorted(per.items(), key=str))
json.dump(dict(u8={str(S): v for S, v in U.items()}, k8={str(S): v for S, v in K.items()}), open(sys.argv[4], 'w'))
