#!/usr/bin/env python3
"""d-bit-exits / kernel_z.py: zero-cost response-kernel entrances for #589 helpers, with re-timed climbs.

Pivot p starts at its first gate frame E_p = F1(p) (saving its first run), D stops reading p, and gates 'd += p' are
inserted at times k <= first write of p, at a frame G that costs nothing: for both p and d, G is either the register's
current frame c_r(k) or the frame of its next gate n_r(k) (then the register climbs early along the same run). The
cleanup is the literal reverse at FULL. p's response vanishes iff the chosen R_d(k) sum to K_p (F2).
Then a greedy selection (pivots never donors, donors never pivots), the modified word, an exact replay on all 661 columns
(symbolic helper dirt, nested frames, D at frame 0) and the local paid profile.
usage: python3 -I kernel_z.py WORD.json.gz OUTDIR"""
import gzip, json, os, sys, copy
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames_lib import Frames

w = json.load(gzip.open(sys.argv[1], 'rt')); outdir = sys.argv[2]; os.makedirs(outdir, exist_ok=True)
V, R = w['v'], w['R']; N = 2 * V + R; F = Frames(w)
A, B = w['A'], w['B']; helpers = list(range(2 * V, N))
isz = lambda r: not (V <= r < 2 * V)
touches = defaultdict(list)                 # r -> list of (k, frame) for every A gate touching r
first_write = {}
for k, (f, t, ss) in enumerate(A):
    for r in [t, *ss]: touches[r].append((k, f))
    first_write.setdefault(t, k)
def cn(r, k, start):
    """current frame before A gate k and frame of the next touch at index >= k (None if no more A touches)"""
    c = start; n = None
    for kk, f in touches[r]:
        if kk < k: c = f
        else: n = f; break
    return c, n
resp = [0] * N
for t in range(V): resp[V + t] = 1 << t
for f, t, ss in reversed(B):
    for s in ss: resp[s] ^= resp[t]
for t, ks in enumerate(w['scatter']):
    for k in ks: resp[w['ret'][k]] ^= 1 << t
respA = [None] * (len(A) + 1); respA[len(A)] = list(resp)
for k in range(len(A) - 1, -1, -1):
    f, t, ss = A[k]
    for s in ss: resp[s] ^= resp[t]
    respA[k] = list(resp)
K = {h: respA[0][h] for h in helpers}
F1 = {h: touches[h][0][1] for h in helpers}

def options(p, banned):
    E = F1[p]; fw = first_write.get(p, len(A)); out = []
    for k in range(0, fw + 1):
        cp, np_ = cn(p, k, E)
        Gs = {cp} | ({np_} if np_ is not None else set())
        for d in range(N):
            if d == p or not isz(d) or d in banned: continue
            cd, nd = cn(d, k, w['start'][d])
            for G in Gs & ({cd} | ({nd} if nd is not None else set())):
                if F.dim[G] == 0: continue
                out.append((respA[k][d], (k, d, G)))
    return out

def solve(vecs, target):
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
    return sorted(s)

feasible = {}
for p in helpers:
    s = solve(options(p, set()), K[p])
    if s is not None: feasible[p] = s
chosen = {}; donors = set()
for p in sorted(feasible, key=lambda p: (-F.dim[F1[p]], p)):
    if p in donors: continue
    s = solve(options(p, set(chosen) | {p}), K[p])
    if s is None: continue
    if any(d in chosen for _, d, _ in s): continue
    chosen[p] = s; donors |= {d for _, d, _ in s}

def build(chosen):
    w2 = copy.deepcopy(w); P = set(chosen)
    w2['D'] = [[f, t, [s for s in ss if s not in P]] for f, t, ss in w['D']]
    w2['D'] = [g for g in w2['D'] if g[2]]
    ins = defaultdict(list)
    for p, lst in chosen.items():
        w2['start'][p] = F1[p]
        for k, d, G in lst: ins[k].append([G, d, [p]])
    A2 = []
    for k, g in enumerate(A):
        A2 += ins.get(k, []); A2.append(g)
    A2 += ins.get(len(A), [])
    w2['A'] = A2
    z = [g for g in w2['A'] + w2['B'] if isz(g[1])]
    w2['C'] = [[F.FULL, t, list(ss)] for f, t, ss in reversed(z)]
    return w2

def replay(w):
    cur = list(w['start']); val = [1 << i for i in range(N)]; H = Counter()
    def climb(r, f):
        a = cur[r]
        if a == f: return
        assert F.le(a, f) and F.dim[f] > F.dim[a], ('climb', r, a, f)
        H[F.dim[f] - F.dim[a]] += 1; cur[r] = f
    for f, t, ss in w['D']:
        assert F.dim[f] == 0 and all(F.dim[cur[r]] == 0 for r in [t, *ss])
    for sec in ('D', 'A'):
        for f, t, ss in w[sec]:
            for r in [t, *ss]: climb(r, f)
            for s in ss: val[t] ^= val[s]
    tops = [F.dim[cur[q]] for q in w['ret']]; snap = [val[q] for q in w['ret']]
    for t, ks in enumerate(w['scatter']):
        for k in ks: val[V + t] ^= snap[k]
    for sec in ('B', 'C'):
        for f, t, ss in w[sec]:
            for r in [t, *ss]: climb(r, f)
            for s in ss: val[t] ^= val[s]
    for r in range(N): climb(r, w['final'][r])
    for r in range(N):
        assert val[r] == ((1 << r) ^ (1 << (r - V)) if V <= r < 2 * V else 1 << r), ('endpoint', r)
    return dict(sorted((H + Counter(tops)).items()))

# keep only pivots whose combined construction replays (drop greedily on failure)
order = sorted(chosen, key=lambda p: (-F.dim[F1[p]], p)); kept = {}
for p in order:
    trial = dict(kept); trial[p] = chosen[p]
    try:
        replay(build(trial)); kept = trial
    except AssertionError as e:
        pass
w2 = build(kept); prof = replay(w2)
coord = lambda E: all(i < 8 for i in F.labels(E)) and F.dim[E] > 0 and len(F.labels(E)) == F.dim[E]
out = dict(feasible_alone=len(feasible), feasible_dims=dict(sorted(Counter(F.dim[F1[p]] for p in feasible).items())),
           greedy=len(chosen), replay_kept=len(kept), kept_entrance_dims=dict(sorted(Counter(F.dim[F1[p]] for p in kept).items())),
           kept_coordinate_entrances=sum(1 for p in kept if coord(F1[p])),
           inserted_gates=sum(len(v) for v in kept.values()), local_profile=prof, mass=sum(r * n for r, n in prof.items()),
           ENT=sum(F.dim[F1[p]] for p in kept))
print(json.dumps(out, indent=1))
json.dump(dict(out, pivots={str(p): v for p, v in kept.items()}), open(os.path.join(outdir, 'kernel_z.json'), 'w'), indent=1)
open(os.path.join(outdir, 'bitE8mut-p111-kz.json.gz'), 'wb').write(
    gzip.compress(json.dumps(w2, sort_keys=True, separators=(',', ':')).encode(), compresslevel=9, mtime=0))
