#!/usr/bin/env python3
"""g-break probes on the U8 CANONICAL word (exact; seconds).
 P1 block span: rank of each block's 8 labels; hyperplanes of the span and their label counts (H1/H2)
 P2 entering helpers (start dim>0) as extra plane visitors: planes P with start frame inside P (H10/H6)
 P3 liveness: value replay (sources symbolic, every helper with its own dirt symbol); at the A->compensation boundary,
    how many registers per block hold a nonzero clean row that a target later reads (H3/H4)
 P4 copy redundancy: at each target read of a plane helper, does another register hold the same clean row
    (up to scale) at that moment?  (H3)
Usage: python3 -I probes.py WORD"""
import sys, json, gzip, itertools
from fractions import Fraction as Fr
from collections import Counter, defaultdict, deque
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
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
B = [basis(f) for f in F]; dims = [len(b) for b in B]
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
ORDER = ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']
# --- census climb (as u8_census) to find the plane of each register
cur = list(c['start']); P_of = {}; path = defaultdict(list)
for r in range(n):
    path[r].append(cur[r])
    if dims[cur[r]] == 3: P_of[r] = tuple(B[cur[r]])
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
    cur[r] = f; path[r].append(f)
    if dims[f] == 3: P_of[r] = tuple(B[f])
for ph in ORDER:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        climb(e[1], e[-1])
        if e[0] == 'add': climb(e[2], e[-1])
for r in range(n): climb(r, c['final'][r])
planes = sorted(set(P_of.values()))
srcs_in = {P: [s for s in range(v) if inspan(ports[s], list(P))] for P in planes}
# blocks
parent = list(range(v))
def find(a):
    while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
    return a
for P in planes:
    for a in srcs_in[P][1:]: parent[find(a)] = find(srcs_in[P][0])
blk_src = defaultdict(list)
for s in range(v): blk_src[find(s)].append(s)
blocks = sorted(blk_src.values())
bid = {s: i for i, S in enumerate(blocks) for s in S}
pblk = {P: bid[srcs_in[P][0]] for P in planes}
out = {}
# P1
def subspaces_hyper(span):
    # all hyperplanes of the span: kernels of nonzero functionals on the span basis
    k = len(span); vecs = [sum(span[i] for i in range(k) if m >> i & 1) and 0 or 0 for m in range(1)]
    elems = []
    for m in range(1, 1 << k):
        x = 0
        for i in range(k):
            if m >> i & 1: x ^= span[i]
        elems.append((m, x))
    H = []
    for f in range(1, 1 << k):
        H.append([x for m, x in elems if bin(m & f).count('1') % 2 == 0])
    return H
p1 = []
for i, S in enumerate(blocks):
    span = basis([ports[s] for s in S]); labs = set(ports[s] for s in S)
    hc = Counter()
    for Hp in subspaces_hyper(span):
        hc[sum(1 for x in Hp if x in labs)] += 1
    # targets (labels ports[t]) whose perp contains the span
    tperp = sum(1 for t in range(v) if all(not dot(x, ports[t]) for x in span))
    p1.append(dict(rank=len(span), hyperplane_label_counts=dict(sorted(hc.items())), targets_perp_containing_span=tperp,
                   used_planes=sum(1 for P in planes if pblk[P] == i)))
out['P1_block_spans'] = dict(Counter(json.dumps(x, sort_keys=True) for x in p1))
# P2: entering helpers
ent = [r for r in range(2*v, n) if dims[c['start'][r]] > 0]
lab = {ports[s]: s for s in range(v)}
p2 = Counter(); compat = defaultdict(list)
for r in ent:
    st = B[c['start'][r]]; isplane = r in P_of
    cp = [P for P in planes if all(inspan(x, list(P)) for x in st)]
    # planes compatible with its existing chain: P contains start and P inside some later frame of its chain
    later = [f for f in path[r][1:]]
    cp_chain = [P for P in cp if any(all(inspan(x, B[f]) for x in P) for f in later)]
    lineislabel = (len(st) == 1 and st[0] in lab)
    p2[(len(st), isplane, lineislabel, len(cp) > 0, len(cp_chain) > 0)] += 1
    if not isplane: compat[r] = cp_chain
out['P2_entering_helpers(startdim,visits_plane,line_is_label,has_compatible_plane,has_plane_inside_its_chain)'] = {str(k): x for k, x in sorted(p2.items())}
# max-flow: replace plane helpers by non-plane entering helpers on chain-compatible planes (each one plane)
# count: per plane need 3 visitors; sources supply 112, current helpers 392; how many plane visits can the 112 non-plane
# entering helpers supply (bipartite matching, plane capacity 3 - sources visiting it - 1 (keep >=1 plane helper))?
def bmatch(left, right_cap, adj):
    # simple augmenting-path b-matching (left cap 1)
    match = defaultdict(list); owner = {}
    def try_(u, seen):
        for P in adj[u]:
            if P in seen: continue
            seen.add(P)
            if len(match[P]) < right_cap[P]: match[P].append(u); owner[u] = P; return True
            for w in list(match[P]):
                if try_(w, seen):
                    match[P].remove(w); match[P].append(u); owner[u] = P; return True
        return False
    tot = 0
    for u in left:
        if try_(u, set()): tot += 1
    return tot
srcvis = Counter(P_of[s] for s in range(v) if s in P_of)
capP = {P: 3 - srcvis[P] for P in planes}
out['P2_nonplane_entering_supply_chain_compatible'] = bmatch([r for r in compat], capP, compat)
out['P2_nonplane_entering_supply_any_plane_containing_start'] = bmatch([r for r in compat], capP,
    {r: [P for P in planes if all(inspan(x, list(P)) for x in B[c['start'][r]])] for r in compat})
# P3/P4: value replay. symbols: ('x',s) sources; ('d',r) helper dirt for every helper; ('t',t) target initial
val = [dict() for _ in range(n)]
for s in range(v): val[s] = {('x', s): Fr(1)}
for t in range(v, 2*v): val[t] = {('t', t): Fr(1)}
for r in range(2*v, n): val[r] = {('d', r): Fr(1)}
clean = lambda r: {k: x for k, x in val[r].items() if k[0] == 'x'}
def normkey(d):
    if not d: return None
    ks = sorted(d); a = d[ks[0]]
    return tuple((k, d[k]/a) for k in ks)
plane_helpers = [r for r in range(2*v, n) if r in P_of]
treads = defaultdict(list)    # helper -> list of (phase, idx)
redundant = Counter(); total_reads = Counter(); live_at = {}
step = 0
for ph in ORDER:
    if ph == 'centres': continue
    if ph == 'compensation':
        # P3 snapshot at A/centres/precondition/cut -> compensation boundary
        live_at['before_compensation'] = {r: normkey(clean(r)) for r in range(n)}
    E = c.get(ph) or []
    for e in E:
        step += 1
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den)
            if cls(a) == 'y' and b in P_of and cls(b) == 'h':
                total_reads[ph] += 1
                kb = normkey(clean(b))
                if kb is not None:
                    others = [r for r in range(n) if r != b and r != a and normkey(clean(r)) == kb]
                    if others: redundant[(ph, tuple(sorted(Counter(cls(r) for r in others).items())))] += 1
                treads[b].append(step)
            for key, x in list(val[b].items()):
                y = val[a].get(key, 0) + k*x
                if y: val[a][key] = y
                else: val[a].pop(key, None)
        elif e[0] == 'scale':
            _, a, num, den = e[:4]; k = Fr(num, den); val[a] = {key: x*k for key, x in val[a].items()}
out['P4_target_reads_of_plane_helpers_by_phase'] = dict(total_reads)
out['P4_reads_where_another_register_holds_same_clean_row'] = {str(k): x for k, x in redundant.items()}
# P3: per block, registers holding a nonzero clean row (support within block) before compensation that targets read later
snap = live_at['before_compensation']
p3 = Counter(); p3d = Counter()
for i, S in enumerate(blocks):
    Sset = set(S)
    holders = [r for r in range(2*v, n) if r in treads and snap[r] is not None and all(k[1] in Sset for k, _ in snap[r])]
    rows = Counter(snap[r] for r in holders)
    p3[(len(holders), len(rows))] += 1
    p3d[tuple(sorted(Counter(len(k) for k in rows).items()))] += 1
out['P3_per_block(helpers holding a clean in-block row before compensation and read by targets later, distinct rows)'] = {str(k): x for k, x in p3.items()}
out['P3_support_sizes_of_those_rows'] = {str(k): x for k, x in p3d.items()}
print(json.dumps(out, indent=1))
json.dump(out, open(sys.argv[2] if len(sys.argv) > 2 else '/dev/null', 'w'), indent=1)
