#!/usr/bin/env python3
"""u-helpers-hub toy: the packet-plane visit model, solved exactly (max-flow), on one Hamming block or all 14.

Rules (the hub's, structural-gain/factory README, Scholz on icekylinx's PR531 factory):
  * frames are F2-subspaces of F2^9; a register's frames form a nested chain (monotone), so it visits at most ONE
    rank-3 frame;  [variant 'descents=d': up to d+1 rank-3 frames per register]
  * each packet plane P (rank 3, spanned by the 4 source labels of a face packet) must be visited by k distinct
    registers: k = 2 to create one copy (a two-register addition at P), k = 3 for two retained copies;
  * a source register s starts at <u_s>, so it can visit P only if u_s in P;
  * a target register t ends at u_t^perp, so it can visit P only if P is inside u_t^perp  [variant 'targets'];
    the hub's interface forbids free use of output registers  [default: targets excluded];
  * helpers start at 0 and can visit any plane.
Least number of helpers = ceil((sum_P k - max plane visits the data registers can supply) / planes per register).
Note 'block0 + targets' lets one block use all 112 targets; it is shown only as a control of the matching code.
The planes and port labels are read from the U8 CANONICAL-WORD (the 168 rank-3 frames its registers visit).
Usage: python3 -I toy_visits.py CANONICAL-WORD.json.gz [--block i]
"""
import sys, json, gzip
from collections import Counter, defaultdict, deque
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return sorted(bb, reverse=True)
def inspan(x, bb):
    for z in bb: x = min(x, x ^ z)
    return x == 0
dot = lambda a, b: bin(a & b).count('1') & 1
ports = c['ports']
# the planes: every rank-3 frame visited by some register in the replay (same rule as u8_census.py)
B = [basis(f) for f in F]; dims = [len(b) for b in B]
cur = list(c['start']); visits = defaultdict(set)
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
    cur[r] = f
    if dims[f] == 3: visits[r].add(tuple(B[f]))
for r in range(n):
    if dims[cur[r]] == 3: visits[r].add(tuple(B[cur[r]]))
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        climb(e[1], e[-1])
        if e[0] == 'add': climb(e[2], e[-1])
for r in range(n): climb(r, c['final'][r])
planes = sorted({p for s in visits.values() for p in s})
# source register s carries label ports[s]; check that its start frame is <u_s> and the target end frame is u_t^perp
src_ok = all(B[c['start'][s]] == basis([ports[s]]) for s in range(v))
perp = lambda u: [x for x in range(1, 512) if not dot(x, u)]
tgt_ok = all(B[c['final'][v+t]] == basis(perp(ports[t])) for t in range(v))
# blocks: group planes by the set of sources they contain (connected components of the source/plane incidence)
srcs_in = {P: [s for s in range(v) if inspan(ports[s], list(P))] for P in planes}
parent = list(range(v))
def find(a):
    while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
    return a
for P in planes:
    ss = srcs_in[P]
    for a in ss[1:]: parent[find(a)] = find(ss[0])
blocks = defaultdict(list)
for P in planes: blocks[find(srcs_in[P][0])].append(P)
blocks = [blocks[k] for k in sorted(blocks)]
def maxflow(planes_sel, k, use_targets, per_reg):
    """bipartite b-matching: data registers (cap per_reg) -> planes (cap k); returns max data visits"""
    regs = [('x', s) for s in range(v)] + ([('y', t) for t in range(v)] if use_targets else [])
    ok = lambda reg, P: (inspan(ports[reg[1]], list(P)) if reg[0] == 'x' else all(not dot(x, ports[reg[1]]) for x in P))
    nodes = ['S', 'T'] + regs + list(planes_sel); idx = {x: i for i, x in enumerate(nodes)}
    cap = defaultdict(int); adj = defaultdict(set)
    def edge(a, b, w): cap[(a, b)] += w; adj[a].add(b); adj[b].add(a)
    for reg in regs:
        edge(idx['S'], idx[reg], per_reg)
        for P in planes_sel:
            if ok(reg, P): edge(idx[reg], idx[P], 1)
    for P in planes_sel: edge(idx[P], idx['T'], k)
    flow = 0
    while True:
        prev = {idx['S']: None}; q = deque([idx['S']])
        while q and idx['T'] not in prev:
            a = q.popleft()
            for b in adj[a]:
                if b not in prev and cap[(a, b)] > 0: prev[b] = a; q.append(b)
        if idx['T'] not in prev: return flow
        b = idx['T']
        while prev[b] is not None: a = prev[b]; cap[(a, b)] -= 1; cap[(b, a)] += 1; b = a
        flow += 1
def solve(planes_sel, k, use_targets=False, per_reg=1):
    data = maxflow(planes_sel, k, use_targets, per_reg)
    need = k*len(planes_sel) - data              # plane visits left for helpers
    # a helper can visit per_reg planes (distinct); with per_reg = 1 this is exact; with per_reg = 2 it is exact as
    # long as the leftover visits are spread over >= 2 planes (true in every case printed)
    return -(-need // per_reg)
out = dict(source_frames_are_labels=src_ok, target_ends_are_perps=tgt_ok, planes=len(planes), blocks=len(blocks),
           planes_per_block=sorted(Counter(len(b) for b in blocks).items()),
           sources_per_plane=sorted(Counter(len(srcs_in[P]) for P in planes).items()),
           planes_containing_a_source=sorted(Counter(sum(1 for P in planes if s in srcs_in[P]) for s in range(v)).items()),
           targets_whose_perp_contains_plane=sorted(Counter(sum(1 for t in range(v) if all(not dot(x, ports[t]) for x in P)) for P in planes).items()))
b0 = blocks[0]
small = {}
for name, args in [('one copy (k=2)', (2,)), ('two copies (k=3)', (3,)), ('two copies + targets', (3, True)),
                   ('one copy + targets', (2, True)), ('two copies, one descent (2 planes/register)', (3, False, 2)),
                   ('two copies + targets + one descent', (3, True, 2))]:
    small[name] = dict(block0=solve(b0, *args), all14=solve(planes, *args))
out['least_helpers'] = small
out['U8_actual'] = dict(plane_visitors_per_plane=3, helper_plane_visitors=392, helpers_not_visiting_a_plane=112, total_helpers=504)
print(json.dumps(out, indent=1))
json.dump(out, open('toy_visits.json', 'w'), indent=1)
