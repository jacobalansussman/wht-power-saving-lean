#!/usr/bin/env python3
"""g-small idea 'data registers as helpers inside the call': for each U8 helper B, find data registers D (source or target)
that have no gate during B's whole lifetime and whose frame chain can absorb B's chain monotonically:
  frame(D just before B's first gate) <= B's frame after its first gate, and
  B's last frame <= frame(D at its next gate after B's lifetime, or D's final frame).
Linearity: B's endpoint holds for arbitrary dirt, so D's content (any X/Y-dependent value) can serve as B's dirt.
Also helper-helper chaining (B2 starts after B1 ends, frames nest) for comparison with d-G06.
Usage: python3 -I u8_idle_merge.py CANONICAL-WORD.json.gz"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return bb
B_ = [basis(f) for f in F]; dims = [len(b) for b in B_]
def sub(a, b):  # frame index a inside frame index b
    for x in B_[a]:
        y = x
        for z in B_[b]: y = min(y, y ^ z)
        if y: return False
    return True
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
cur = list(c['start']); T = 0
gates = defaultdict(list)   # r -> [(T, frame after)]
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not sub(cur[r], f): return
    cur[r] = f
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        T += 1
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f); gates[r].append((T, cur[r]))
        continue
    for e in c.get(ph) or []:
        T += 1
        regs = [e[1]] + ([e[2]] if e[0] == 'add' else [])
        for r in regs: climb(r, e[-1])
        for r in regs: gates[r].append((T, cur[r]))
Tend = T + 1
fin = c['final']
def frame_before(r, t):  # frame of r just before time t
    f = c['start'][r]
    for tt, ff in gates[r]:
        if tt < t: f = ff
        else: break
    return f
def next_frame_after(r, t):
    for tt, ff in gates[r]:
        if tt > t: return ff, tt
    return fin[r], Tend
helpers = [r for r in range(2*v, n) if gates[r]]
life = {r: (gates[r][0][0], gates[r][-1][0]) for r in helpers}
data = list(range(2*v))
ok = defaultdict(list)
for b in helpers:
    t0, t1 = life[b]; f_first = gates[b][0][1]; f_last = gates[b][-1][1]
    for d in data:
        if any(t0 <= tt <= t1 for tt, _ in gates[d]): continue
        fd = frame_before(d, t0)
        if not sub(fd, f_first): continue
        fn, tn = next_frame_after(d, t1)
        if not sub(f_last, fn): continue
        ok[b].append(d)
print('helpers with a gate', len(helpers))
print('helpers absorbable by some idle data register (frames nest):', len(ok))
print('by data class:', Counter(cls(d) for b in ok for d in ok[b]))
print('lifetime lengths of helpers (gate count) for absorbable ones:', Counter(len(gates[b]) for b in ok).most_common(10))
# greedy disjoint assignment: a data register can host several helpers if their lifetimes are disjoint and chained frames nest
assign = {}; used = defaultdict(list)
for b in sorted(ok, key=lambda b: life[b][0]):
    for d in ok[b]:
        if all(life[b][0] > life[o][1] or life[b][1] < life[o][0] for o in used[d]):
            # chained frames: previous hosted helper's last frame <= this one's first frame if earlier
            prev = [o for o in used[d] if life[o][1] < life[b][0]]
            if prev:
                o = max(prev, key=lambda o: life[o][1])
                if not sub(gates[o][-1][1], gates[b][0][1]): continue
            assign[b] = d; used[d].append(b); break
print('greedy hosted helpers:', len(assign), ' data registers used:', len(used))
wid = lambda r: dims[fin[r]] - dims[c['start'][r]]
print('widths of greedily absorbed helpers:', Counter(dims[gates[b][-1][1]] - dims[c['start'][b]] for b in assign))
ex = sorted(assign)[:8]
for b in ex:
    d = assign[b]
    print(f' helper {b} life {life[b]} frames {[dims[f] for _, f in gates[b]][:12]} -> data {d} ({cls(d)}) frame before {dims[frame_before(d, life[b][0])]} next {dims[next_frame_after(d, life[b][1])[0]]}')
json.dump(dict(absorbable=len(ok), greedy=len(assign), pairs={str(b): ok[b][:10] for b in ok}), open('u8_idle_merge.json', 'w'), indent=1)
# diagnostics: which condition fails
idle = Counter(); fr1 = Counter(); both = 0
for b in helpers:
    t0, t1 = life[b]
    ids = [d for d in data if not any(t0 <= tt <= t1 for tt, _ in gates[d])]
    idle[len(ids) > 0] += 1
    f1 = [d for d in ids if sub(frame_before(d, t0), gates[b][0][1])]
    fr1[len(f1) > 0] += 1
print('helpers with >=1 idle data register during their life:', idle[True], ' and start-frame nesting:', fr1[True])
print('helper lifetime as fraction of call, histogram (tenths):', sorted(Counter(int(10*(life[b][1]-life[b][0])/Tend) for b in helpers).items()))
