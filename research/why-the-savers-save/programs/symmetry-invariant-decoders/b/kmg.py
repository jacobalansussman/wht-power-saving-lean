"""kmg.py (t-km): Kramer-Mesner search for star decoders on P112 invariant under a prescribed CYCLIC group <g> of any
order n (g = port permutation 'perm{n}' in gens.npz; it comes from a linear map of F_2^9, so it permutes flats).
Frames are flats of the label matroid (closed port sets; cost = dim).  A G-invariant family is a union of G-orbits of
flats.  Shape: one BIG orbit type (size m, dim d), enumerated over all its orbit representatives (vectorised), plus a
multiset of SMALL orbits (any size and dim) with total cost <= COST - m*d and K >= KMIN frames in all.
Exact tests (any payload field), as km7b.py: stage 0 count (8 distinct frames meet every octad), clique Hall on all
1996 maximal cliques, pinning.  usage: python -I kmg.py n m d COST [KMIN] [MAXSMALLDIM]  -> kmg_n{n}_m{m}_d{d}_c{COST}.json"""
import sys, json, time, itertools
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
import numpy as np
from collections import Counter, defaultdict
n = int(sys.argv[1]); m = int(sys.argv[2]); d = int(sys.argv[3]); COST = int(sys.argv[4])
KMIN = int(sys.argv[5]) if len(sys.argv) > 5 else 9
budget = COST - m * d
import os as _os0
SEC = _os0.environ.get('SECOND')                       # "m2xd2": a second enumerated orbit type
m2, d2 = (tuple(int(a) for a in SEC.split('x')) if SEC else (0, 0))
budget -= m2 * d2
perm = np.load('gens.npz')['perm%d' % n].astype(np.int64); inv = np.argsort(perm)
t0 = time.time()
# flats by true dim (a mask first produced at dim k has true dim k)
import os
CACHE = 'flats_cache.npz'
if not os.path.exists(CACHE):
    seen = set(); allfl = {}
    for k in range(1, 9):
        _, M = pool(k, minports=1); M = np.unique(M, axis=0)
        Pk = np.packbits(M, axis=1); del M
        keys = [r.tobytes() for r in Pk]
        new_ = [i for i, kk in enumerate(keys) if kk not in seen]
        seen.update(keys[i] for i in new_)
        allfl['d%d' % k] = Pk[new_]
        print('  flats of true dim', k, len(new_), flush=True)
    np.savez(CACHE, **allfl)
cz = np.load(CACHE)
_ST = _os0.environ.get('SMALLTYPES')
need_dims = sorted(set([d] + ([d2] if SEC else []) + ([int(x.split('x')[1]) for x in _ST.split(',') if int(x.split('x')[1]) < 9] if _ST else [k for k in range(1, 9) if k <= budget])))
FL = {k: np.unpackbits(cz['d%d' % k], axis=1)[:, :NP].astype(bool) for k in need_dims}
print('flats by dim:', {k: len(v) for k, v in FL.items()}, '%.0fs' % (time.time() - t0), flush=True)
def orbit_info(M):
    """orbit size and representative flag, in chunks (memory)"""
    size = np.zeros(len(M), dtype=np.int64); isrep = np.zeros(len(M), dtype=bool)
    for a in range(0, len(M), 100000):
        A = M[a:a + 100000]
        P0 = np.packbits(A, axis=1)
        k0 = np.concatenate([P0, np.zeros((len(P0), 2), np.uint8)], axis=1).view('>u8').astype(np.uint64)
        keys = [k0]
        for k in range(1, n):
            A = A[:, inv]
            Pk = np.packbits(A, axis=1)
            keys.append(np.concatenate([Pk, np.zeros((len(Pk), 2), np.uint8)], axis=1).view('>u8').astype(np.uint64))
        K = np.stack(keys)                                   # n x c x 2
        same = [(K[k] == K[0]).all(1) for k in range(n)]
        # orbit size = n / |stabiliser| ; stabiliser = {k : g^k F = F}
        stab = np.sum(same, axis=0)
        size[a:a + len(P0)] = n // stab
        # representative = lexicographically least key among the images
        lt = np.ones(len(P0), bool)
        for k in range(1, n):
            hi0, lo0, hik, lok = K[0][:, 0], K[0][:, 1], K[k][:, 0], K[k][:, 1]
            lt &= (hi0 < hik) | ((hi0 == hik) & (lo0 <= lok))
        isrep[a:a + len(P0)] = lt
    return size, isrep
REPS = {}
for k, M in FL.items():
    size, isrep = orbit_info(M)
    for s in sorted(set(size.tolist())):
        sel = isrep & (size == s)
        REPS[(s, k)] = M[sel]
if _os0.environ.get('ALLOW9') == '1':
    REPS[(1, 9)] = np.ones((1, NP), bool)            # the whole space: one frame reading every port, cost 9
print('orbit representatives (size, dim) -> count:', {str(a): len(b) for a, b in sorted(REPS.items())}, '%.0fs' % (time.time() - t0), flush=True)
def orbit_frames(mask, s):
    out = [mask]; A = mask
    for k in range(1, s): A = A[inv]; out.append(A)
    return out
# octad orbits
OC = [tuple(sorted(o)) for o in OCTP]; sn = set(); oreps = []
for o in OC:
    if o in sn: continue
    q = o
    for k in range(n): sn.add(q); q = tuple(sorted(perm[list(q)]))
    oreps.append(o)
OR = np.array(oreps)
# maximal cliques
OM = [sum(1 << T for T in range(NP) if ORTH[S, T]) for S in range(NP)]
MC = []
def bk(Rr, Pm, Xm):
    if not Pm and not Xm: MC.append(Rr); return
    u = (Pm | Xm).bit_length() - 1; c = Pm & ~OM[u]
    while c:
        low = c & -c; v = low.bit_length() - 1; c ^= low
        bk(Rr + [v], Pm & OM[v], Xm & OM[v]); Pm &= ~low; Xm |= low
bk([], (1 << NP) - 1, 0)
def full_tests(frames):
    S = np.array(frames)
    J = [frozenset(np.nonzero(S[:, T])[0].tolist()) for T in range(NP)]
    if any(not j for j in J): return 'port unread'
    def hall_ok(cl):
        match = {}
        def aug(T, sn_):
            for j in J[T]:
                if j in sn_: continue
                sn_.add(j)
                if j not in match or aug(match[j], sn_): match[j] = T; return True
            return False
        return all(aug(T, set()) for T in cl)
    nf = sum(1 for c in MC if not hall_ok(c))
    if nf: return 'clique Hall fails'
    single = {T: next(iter(J[T])) for T in range(NP) if len(J[T]) == 1}
    for a, b in itertools.combinations(single, 2):
        if ORTH[a, b] and single[a] == single[b]: return 'orthogonal pair one frame'
    bad = [Sx for Sx in range(NP) if all(any(single.get(T) == j for T in np.nonzero(ORTH[Sx])[0]) for j in J[Sx])]
    if bad: return 'pinned'
    return 'PASS'
# big piece: per rep, count of its frames meeting each octad rep
BIG = REPS.get((m, d), np.zeros((0, NP), bool))
POP16 = np.array([bin(i).count('1') for i in range(1 << 16)], dtype=np.int64)
def meet_counts(R, s):
    cnt = np.zeros((len(R), len(oreps)), dtype=np.int64); A = R.copy()
    for k in range(s):
        cnt += A[:, OR].any(2); A = A[:, inv]
    return cnt
cb = meet_counts(BIG, m)
# small pieces: all orbit reps with cost <= budget (excluding nothing; the big type may also appear as small)
small = []
import os as _os
ST = _os.environ.get('SMALLTYPES')                 # e.g. "1x6,1x7,2x3" restricts the small orbit types
STs = None if not ST else {tuple(int(a) for a in x.split('x')) for x in ST.split(',')}
MAXP = int(_os.environ.get('MAXPIECES', '99'))
for (s, k), R in REPS.items():
    if s * k <= budget and (STs is None or (s, k) in STs):
        for i in range(len(R)): small.append((s, k, R[i]))
print('big reps', len(BIG), 'small orbit pieces within budget %d:' % budget, len(small), Counter((s, k) for s, k, _ in small), flush=True)
smc = np.array([meet_counts(r[None], s)[0] for s, k, r in small]) if small else np.zeros((0, len(oreps)), np.int64)
combos = []
def rec(start, cur, cost, nfr):
    if nfr + m + m2 >= KMIN: combos.append(tuple(cur))
    if len(cur) >= MAXP: return
    for j in range(start, len(small)):
        s, k, _ = small[j]
        if cost + s * k <= budget: rec(j, cur + [j], cost + s * k, nfr + s)
rec(0, [], 0, 0)
print('small combos:', len(combos), flush=True)
stats = Counter(); surv = []
SECR = REPS.get((m2, d2), np.zeros((0, NP), bool)) if SEC else [None]
c2all = meet_counts(SECR, m2) if SEC else None
print('second orbit type', SEC, 'reps', len(SECR), flush=True)
ADD = np.array([(smc[list(cbn)].sum(0) if cbn else np.zeros(len(oreps), np.int64)) for cbn in combos])
addmax = ADD.max(0) if len(ADD) else np.zeros(len(oreps), np.int64)
cb8 = cb.astype(np.int16)
for i2 in range(len(SECR)):
    c2 = c2all[i2] if SEC else 0
    pre = np.nonzero(((cb8 + (c2 + addmax)[None, :].astype(np.int16)) >= 8).all(1))[0]   # pair prefilter
    stats['pair prefilter pass'] += len(pre)
    if len(pre) == 0: continue
    for ci, cbn in enumerate(combos):
        add = ADD[ci] + c2
        ok = pre[((cb[pre] + add[None, :]) >= 8).all(1)]
        stats['families'] += len(BIG); stats['stage0 pass'] += len(ok)
        if len(ok):
            # stage 0b: octad Hall on the octad-orbit representatives, vectorised (port bitmasks over frames)
            Wb = np.zeros((len(ok), NP), dtype=np.int64); bit = 0
            A_ = BIG[ok].copy()
            for k in range(m): Wb |= A_.astype(np.int64) << bit; A_ = A_[:, inv]; bit += 1
            extra = []
            if SEC: extra += orbit_frames(SECR[i2], m2)
            for j in cbn: extra += orbit_frames(small[j][2], small[j][0])
            for fr in extra: Wb |= fr.astype(np.int64)[None, :] << bit; bit += 1
            if bit <= 62:
                pm_ = Wb[:, OR]; Nq = np.zeros((256,) + pm_.shape[:2], dtype=np.int64); worst = np.full(pm_.shape[:2], -8)
                for q in range(1, 256):
                    low = (q & -q).bit_length() - 1
                    Nq[q] = Nq[q & (q - 1)] | pm_[:, :, low]
                    cntq = POP16[Nq[q] & 0xFFFF] + POP16[(Nq[q] >> 16) & 0xFFFF] + POP16[(Nq[q] >> 32) & 0xFFFF] + POP16[(Nq[q] >> 48) & 0xFFFF]
                    worst = np.maximum(worst, bin(q).count('1') - cntq)
                keep = (worst <= 0).all(1)
                stats['octad Hall fails'] += int((~keep).sum()); ok = ok[keep]
        for i in ok:
            frames = orbit_frames(BIG[i], m)
            if SEC: frames += orbit_frames(SECR[i2], m2)
            for j in cbn: frames += orbit_frames(small[j][2], small[j][0])
            v = full_tests(frames); stats[v] += 1
            if v == 'PASS':
                surv.append(dict(big=np.nonzero(BIG[i])[0].tolist(), second=(None if not SEC else np.nonzero(SECR[i2])[0].tolist()),
                                 small=[(small[j][0], small[j][1], np.nonzero(small[j][2])[0].tolist()) for j in cbn],
                                 cost=m * d + m2 * d2 + sum(small[j][0] * small[j][1] for j in cbn)))
print('n %d big (size %d, dim %d) second %s cost <= %d K >= %d: %s survivors %d  %.0fs' % (n, m, d, SEC, COST, KMIN, dict(stats), len(surv), time.time() - t0), flush=True)
json.dump(dict(n=n, m=m, d=d, SECOND=SEC, COST=COST, KMIN=KMIN, SMALLTYPES=ST, MAXPIECES=MAXP, reps={str(a): len(b) for a, b in REPS.items()}, small=len(small),
               combos=len(combos), stats=dict(stats), survivors=surv[:200], n_survivors=len(surv)),
          open('kmg_n%d_m%d_d%d_c%d%s%s.json' % (n, m, d, COST, ('_' + SEC) if SEC else '', '_with9' if _os0.environ.get('ALLOW9') == '1' else ''), 'w'), indent=1)
