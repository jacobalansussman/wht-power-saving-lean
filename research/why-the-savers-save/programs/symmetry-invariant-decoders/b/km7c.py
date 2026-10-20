"""km7b.py (t-km): Kramer-Mesner search, G = <g> of prime order p (gens.npz), for G-invariant star decoders on P112
with cost <= COST: one G-orbit of p frames of dim d (orbit representative of a non-fixed flat) + a multiset of
G-fixed flats (fixed_flats_g{p}.json).  Exact necessary tests, any payload field:
  stage 0  count: every octad meets >= 8 distinct frames (one octad per G-orbit suffices);
  stage 1  Hall on every maximal clique of the orthogonality graph on P112 (1260 octads + 512 7-cliques + 224 6-cliques);
  stage 2  pinning: a port T read by one frame j has rho_T = c e_j; port S fails if every frame of S is pinned by
           a port orthogonal to S (then rho_S lies in span rho_{N(S)}); also two orthogonal ports with the same single frame.
km7c variant: a SECOND G-orbit of p frames of small dim d2 (1 or 2) is added (all its orbit representatives).
usage: python -I km7c.py p d d2 COST [KMIN]   -> kmc_p{p}_d{d}_e{d2}_c{COST}.json"""
import sys, json, time, itertools
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
import numpy as np
from collections import Counter
p = int(sys.argv[1]); d = int(sys.argv[2]); d2 = int(sys.argv[3]); COST = int(sys.argv[4]); KMIN = int(sys.argv[5]) if len(sys.argv) > 5 else 9
z = np.load('gens.npz'); perm = z['perm%d' % p].astype(np.int64); inv = np.argsort(perm)
fixed = [(np.array(sorted(f)), dd) for f, dd in json.load(open('fixed_flats_g%d.json' % p)) if dd < 9 or __import__('os').environ.get('ALLOW9') == '1']
OM = [sum(1 << T for T in range(NP) if ORTH[S, T]) for S in range(NP)]
MC = []
def bk(Rr, Pm, Xm):
    if not Pm and not Xm: MC.append(Rr); return
    u = (Pm | Xm).bit_length() - 1; c = Pm & ~OM[u]
    while c:
        low = c & -c; v = low.bit_length() - 1; c ^= low
        bk(Rr + [v], Pm & OM[v], Xm & OM[v]); Pm &= ~low; Xm |= low
bk([], (1 << NP) - 1, 0)
OC = [tuple(sorted(o)) for o in OCTP]; seen = set(); oreps = []
for o in OC:
    if o in seen: continue
    q = o
    for k in range(p): seen.add(q); q = tuple(sorted(perm[list(q)]))
    oreps.append(o)
OR = np.array(oreps)
t0 = time.time()
ds, M = pool(d, minports=1); M = np.unique(M, axis=0)
dims = np.array([f2dim([LAB[i] for i in np.nonzero(r)[0]]) for r in M]); M = M[dims == d]
imgs = [M]
for k in range(1, p): imgs.append(imgs[-1][:, inv])
keyv = [np.array([int.from_bytes(bytes(r), 'big') for r in np.packbits(A, axis=1)], dtype=object) for A in imgs]
canon = [min(keyv[k][i] for k in range(p)) for i in range(len(M))]
isfix = np.array([all(keyv[k][i] == keyv[0][i] for k in range(p)) for i in range(len(M))])
rep = np.array([keyv[0][i] == canon[i] and not isfix[i] for i in range(len(M))])
R = M[rep]
print('p', p, 'd', d, 'flats', len(M), 'orbit reps', len(R), 'octad orbits', len(oreps), 'maximal cliques', len(MC),
      Counter(len(c) for c in MC), '%.0fs' % (time.time() - t0), flush=True)
def reps_of(dd):
    _, Mx = pool(dd, minports=1); Mx = np.unique(Mx, axis=0)
    dx = np.array([f2dim([LAB[i] for i in np.nonzero(r)[0]]) for r in Mx]); Mx = Mx[dx == dd]
    im = [Mx]
    for k in range(1, p): im.append(im[-1][:, inv])
    kv = [[int.from_bytes(bytes(r), 'big') for r in np.packbits(A_, axis=1)] for A_ in im]
    keep = [i for i in range(len(Mx)) if kv[0][i] == min(kv[k][i] for k in range(p)) and any(kv[k][i] != kv[0][i] for k in range(p))]
    return Mx[keep]
R2 = reps_of(d2)
print('secondary dim', d2, 'orbit reps', len(R2), flush=True)
W = np.zeros((len(R), NP), dtype=np.int64); A = R.copy()
for k in range(p): W |= A.astype(np.int64) << k; A = A[:, inv]
POP16 = np.array([bin(i).count('1') for i in range(1 << 16)], dtype=np.int64)
pc = lambda x: POP16[x & 0xFFFF] + POP16[(x >> 16) & 0xFFFF]
ORm = np.zeros((len(R), len(oreps)), dtype=np.int64)
for j in range(8): ORm |= W[:, OR[:, j]]
cnt_orbit = pc(ORm)                                                        # n x octad-orbits
fmeet = np.array([[bool(set(f.tolist()) & set(o)) for o in oreps] for f, dd in fixed], dtype=np.int64)  # nf x orbits
budget = COST - p * d - p * d2; combos = []
def rec(start, cur, cost):
    if len(cur) + 2 * p >= KMIN: combos.append(tuple(cur))
    for j in range(start, len(fixed)):
        if cost + fixed[j][1] <= budget: rec(j, cur + [j], cost + fixed[j][1])
rec(0, [], 0)
print('fixed combos within budget', budget, ':', len(combos), flush=True)
def full_tests(frames):
    S = np.array(frames); K = len(S)
    J = [frozenset(np.nonzero(S[:, T])[0].tolist()) for T in range(NP)]
    if any(not j for j in J): return 'port unread'
    def hall_ok(cl):
        match = {}
        def aug(T, sn):
            for j in J[T]:
                if j in sn: continue
                sn.add(j)
                if j not in match or aug(match[j], sn): match[j] = T; return True
            return False
        return all(aug(T, set()) for T in cl)
    nf = sum(1 for c in MC if not hall_ok(c))
    if nf: return 'clique Hall fails on %d maximal cliques' % nf
    single = {T: next(iter(J[T])) for T in range(NP) if len(J[T]) == 1}
    for a, b in itertools.combinations(single, 2):
        if ORTH[a, b] and single[a] == single[b]: return 'orthogonal pair with one shared frame'
    bad = [Sx for Sx in range(NP) if all(any(single.get(T) == j for T in np.nonzero(ORTH[Sx])[0]) for j in J[Sx])]
    if bad: return 'pinned ports %d' % len(bad)
    return 'PASS'
stats = Counter(); survivors = []
for i2 in range(len(R2)):
    W2 = np.zeros(NP, dtype=np.int64); A2 = R2[i2].copy()
    for k in range(p): W2 |= A2.astype(np.int64) << k; A2 = A2[inv]
    OR2 = np.zeros(len(oreps), dtype=np.int64)
    for j in range(8): OR2 |= W2[OR[:, j]]
    c2 = pc(OR2)
    for cb in combos:
        cf = fmeet[list(cb)].sum(0) if cb else np.zeros(len(oreps), dtype=np.int64)
        ok = np.nonzero(((cnt_orbit + (c2 + cf)[None, :]) >= 8).all(1))[0]
        stats['stage0 pass'] += len(ok); stats['candidates'] += len(R)
        for i in ok:
            frames = []
            m = R[i].copy()
            for k in range(p): frames.append(m.copy()); m = m[inv]
            m = R2[i2].copy()
            for k in range(p): frames.append(m.copy()); m = m[inv]
            for fi in cb:
                mm = np.zeros(NP, bool); mm[fixed[fi][0]] = True; frames.append(mm)
            verdict = full_tests(frames)
            stats[verdict.split(' on ')[0] if 'Hall' in verdict else verdict.split(' ')[0]] += 1
            if verdict == 'PASS':
                survivors.append(dict(rep=np.nonzero(R[i])[0].tolist(), rep2=np.nonzero(R2[i2])[0].tolist(), fixed_idx=list(cb)))
print('cost <=', COST, 'stats', dict(stats), 'survivors', len(survivors), '%.0fs' % (time.time() - t0), flush=True)
json.dump(dict(p=p, d=d, d2=d2, secondary_reps=int(len(R2)), COST=COST, flats=int(len(M)), reps=int(len(R)), combos=len(combos), stats=dict(stats),
               survivors=survivors), open('kmc_p%d_d%d_e%d_c%d%s.json' % (p, d, d2, COST, '_with9' if __import__('os').environ.get('ALLOW9') == '1' else ''), 'w'), indent=1)
