#!/usr/bin/env python3
"""d-bit-exits / kernel_t.py: time-shifted response kernels for #589 helper entrances.

A pivot helper p starts at its first gate frame E = F1(p) (or not at all), and D no longer reads it. To cancel p's
target response K_p, insert 'd += p' gates at times tau when p still holds only its dirt s_p (no write to p yet) and
the donor d sits at exactly the frame of p at that time (E), so no register climbs anywhere it did not climb before.
The inserted gate adds R_d(tau) s_p to the targets, where R_d(tau) is the target response of an injection into d right
after gate tau (adjoint replay of A, scatter, B over F2). The cleanup is the literal reverse at FULL (it undoes the
inserted gate too), so every z register is restored. p's response vanishes iff sum of the chosen R_d(tau) = K_p.

Reports, per helper: whether K_p lies in the F2 span of its zero-cost injections; then a greedy disjoint selection
(pivots never donors), filtered also by the pipeline's residual-isometry family (entrance frames spanned by the
coordinate lines 0..7: research/e8-r740-replica32 bank_routes.validate_isometries accepts only coordinate swaps).
usage: python3 -I kernel_t.py WORD.json.gz [out.json]"""
import gzip, json, os, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames_lib import Frames

w = json.load(gzip.open(sys.argv[1], 'rt'))
V, R = w['v'], w['R']; N = 2 * V + R
F = Frames(w)
helpers = list(range(2 * V, N))
isz = lambda r: not (V <= r < 2 * V)
# forward pass: frames over time for A (index t = position in A; time 'before gate t')
A, B = w['A'], w['B']
cur = list(w['start'])
for f, t, ss in w['D']:
    pass                                               # D at frame 0, no climbs
frame_before = []                                      # frame_before[t][r] would be huge; store events instead
first_write = {}; first_touch = {}
history = defaultdict(list)                            # r -> list of (t, frame) climbs during A
for k, (f, t, ss) in enumerate(A):
    for r in [t, *ss]:
        if cur[r] != f:
            history[r].append((k, f)); cur[r] = f
        first_touch.setdefault(r, k)
    first_write.setdefault(t, k)
def frame_at(r, k):
    """frame of r just before A gate k (after D)"""
    fr = w['start'][r]
    for (kk, f) in history[r]:
        if kk < k: fr = f
        else: break
    return fr
# adjoint pass: response (bitmask of targets) of an injection into register r just before A gate k
resp = [0] * N
for t in range(V): resp[V + t] = 1 << t
for f, t, ss in reversed(B):
    for s in ss: resp[s] ^= resp[t]
for t, ks in enumerate(w['scatter']):
    for k in ks: resp[w['ret'][k]] ^= 1 << t
respA = [None] * (len(A) + 1)
respA[len(A)] = list(resp)
inj = defaultdict(list)                                # pivot candidate p -> list of (k, d, R_d)
# for the selection we need R_d(k) for donors at frame F1(p) at time k <= first_write(p); collect lazily
need = {}
for p in helpers:
    E = A[first_touch[p]][0]
    need[p] = (E, first_write.get(p, len(A)))
for k in range(len(A) - 1, -1, -1):
    f, t, ss = A[k]
    for s in ss: resp[s] ^= resp[t]
    respA[k] = list(resp)
# K from D reads, check K = adjoint response at the start of A
K = defaultdict(int)
for f, t, ss in w['D']:
    for s in ss: K[s] ^= 1 << (t - V)
assert all(K[h] == respA[0][h] for h in helpers), 'K != adjoint response'
# candidate injections per pivot: times k in [0, first_write(p)] (inject just before gate k), donor d != p, z register,
# d at frame exactly E at time k, and p at E too (p's frame before its first touch is its entrance E; after a read of p
# at E it is still E until it climbs)
def p_frame_ok(p, k, E):
    fr = frame_at(p, k)
    return fr == E or fr == w['start'][p]
cands = {}
for p in helpers:
    E, fw = need[p]
    lst = []
    for k in range(0, fw + 1):
        if not p_frame_ok(p, k, E): break
        for d in range(N):
            if d == p or not isz(d): continue
            if frame_at(d, k) == E:
                lst.append((k, d, respA[k][d]))
    cands[p] = lst
def in_span(vecs, target):
    basis = {}
    for v, tag in vecs:
        x = v; s = {tag}
        while x:
            b = x.bit_length() - 1
            if b in basis: x ^= basis[b][0]; s ^= basis[b][1]
            else: basis[b] = (x, s); break
    x = target; s = set()
    while x:
        b = x.bit_length() - 1
        if b not in basis: return None
        x ^= basis[b][0]; s ^= basis[b][1]
    return s
coord = set(range(8))
def coordinate_frame(E):
    labs = F.labels(E)
    sub = [i for i in labs if i in coord]
    return len(sub) == F.dim[E] and F.dim[E] > 0
ok = {}; stats = Counter()
for p in helpers:
    E, fw = need[p]
    vecs = [(r, (k, d)) for (k, d, r) in cands[p]]
    s = in_span(vecs, K[p])
    stats['candidates>0' if cands[p] else 'no zero-cost injection'] += 1
    if s is not None:
        ok[p] = sorted(s); stats['K_p in span'] += 1
        if coordinate_frame(E): stats['K_p in span and E coordinate'] += 1
# greedy disjoint selection: pivots never donors; process pivots by entrance dim (2 first), then register
chosen = {}; donors = set()
for p in sorted(ok, key=lambda p: (-F.dim[need[p][0]], p)):
    if p in donors: continue
    E, fw = need[p]
    vecs = [(r, (k, d)) for (k, d, r) in cands[p] if d not in chosen]
    s = in_span(vecs, K[p])
    if s is None: continue
    chosen[p] = sorted(s); donors |= {d for (k, d) in s}
out = dict(stats=dict(stats), chosen=len(chosen),
           chosen_entrance_dims=dict(sorted(Counter(F.dim[need[p][0]] for p in chosen).items())),
           chosen_coordinate=sum(1 for p in chosen if coordinate_frame(need[p][0])),
           chosen_coordinate_dims=dict(sorted(Counter(F.dim[need[p][0]] for p in chosen if coordinate_frame(need[p][0])).items())),
           donor_count=len(donors), mean_injections=sum(len(v) for v in chosen.values()) / max(1, len(chosen)))
print(json.dumps(out, indent=1))
if len(sys.argv) > 2:
    out['pivots'] = {str(p): dict(E=need[p][0], injections=v) for p, v in chosen.items()}
    json.dump(out, open(sys.argv[2], 'w'), indent=1)
