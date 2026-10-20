#!/usr/bin/env python3
"""g-recombine: which cold helpers of the U8 word could be PROMOTED to enter at their first gate frame?
Hub rule (complex-t21-owners/vendor/pr870/independent/audit_local.py, PROOF-NOTE-PROMO): an entering helper may not be
read by a compensation gate ('nonzero entrance response') nor touched by precondition/compensation; frames must stay
monotone. A cold helper (start frame 0) is promotable to E = its first gate frame if no precondition or compensation
gate touches it.  Reports the count and the run that would disappear.
Usage: python3 -I promo_scan.py CANONICAL-WORD.json.gz"""
import sys, json, gzip
from collections import Counter
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return bb
dims = [len(basis(f)) for f in F]
touched_comp = Counter(); touched_pre = Counter()
for e in c['compensation']:
    for q in (e[1:3] if e[0] == 'add' else e[1:2]): touched_comp[q] += 1
for e in c['precondition']:
    for q in (e[1:3] if e[0] == 'add' else e[1:2]): touched_pre[q] += 1
cold = [r for r in range(2*v, n) if dims[c['start'][r]] == 0]
ent = [r for r in range(2*v, n) if dims[c['start'][r]] > 0]
print('cold helpers', len(cold), 'entering', len(ent))
print('cold touched by compensation', sum(1 for r in cold if touched_comp[r]), ' by precondition', sum(1 for r in cold if touched_pre[r]))
print('entering touched by compensation', sum(1 for r in ent if touched_comp[r]), ' by precondition', sum(1 for r in ent if touched_pre[r]))
# first gate frame of each cold helper in the main program (A, cut, B ...) skipping precondition/compensation
first = {}
for ph in ['A', 'cut', 'B', 'tail', 'cleanup']:
    for e in c.get(ph) or []:
        for q in (e[1:3] if e[0] == 'add' else e[1:2]):
            if q in cold and q not in first and e[-1] is not None: first[q] = dims[e[-1]]
free = [r for r in cold if not touched_comp[r] and not touched_pre[r]]
print('promotable cold helpers (untouched by the wrapper):', len(free), ' first-gate dims', Counter(first.get(r) for r in free))
print('non-promotable cold helpers first-gate dims', Counter(first.get(r) for r in cold if r not in free))
print('compensation gates per touched cold helper', Counter(touched_comp[r] for r in cold if touched_comp[r]).most_common(8))
# --- group promotion: cold helpers linked by frame-0 gates (any phase); a group can enter at E = intersection of
# its members' first non-zero gate frames, if every member is untouched by the wrapper and no frame-0 gate links the
# group to a register outside it (sources/targets/other helpers).
B = [basis(f) for f in F]
def inter(a, b):
    # intersection of two subspaces given by bases (ints over F2^9): brute force over the span of a
    span = {0}
    for x in a: span |= {y ^ x for y in span}
    bb = basis(b)
    def ins(x):
        for z in sorted(bb, reverse=True): x = min(x, x ^ z)
        return x == 0
    return basis([x for x in span if x and ins(x)])
parent = {r: r for r in cold}
def find(a):
    while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
    return a
outside = Counter(); firstnz = {}
for ph in ['A', 'cut', 'B', 'tail', 'cleanup']:
    for e in c.get(ph) or []:
        regs = e[1:3] if e[0] == 'add' else e[1:2]; f = e[-1]
        if f is None: continue
        if dims[f] == 0:
            cs = [q for q in regs if q in parent]
            if len(cs) == 2: parent[find(cs[0])] = find(cs[1])
            for q in cs:
                if len(cs) < len(regs): outside[q] += 1
        else:
            for q in regs:
                if q in parent and q not in firstnz: firstnz[q] = f
groups = {}
for r in cold: groups.setdefault(find(r), []).append(r)
gain = Counter(); summary = Counter()
for g in groups.values():
    ok = all(not touched_comp[r] and not touched_pre[r] and not outside[r] for r in g)
    E = None
    for r in g:
        E = list(B[firstnz[r]]) if E is None else inter(E, B[firstnz[r]])
    d = len(E) if E is not None else 0
    summary[(len(g), ok, d)] += 1
    if ok and d: gain[d] += len(g)
print('groups (size, promotable, dim E):', sorted(summary.items()))
print('helpers that could enter at dim d (d: count):', dict(gain), ' width saved per call:', sum(d*k for d, k in gain.items()))
