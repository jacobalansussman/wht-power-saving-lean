"""Label-space search with repair: a seed relabel (merge a register's two consecutive frames, or one gate to a frame of
its chain) is closed under the realisability rules by propagation: a gate forced above (DAG/section successor) is
raised to the join; a gate on a shared register whose frame becomes incomparable is raised to the join if it is at
least as high, else lowered to the meet; bounds are checked.  The closed relabel set is costed exactly."""
import sys, os, time, random
from collections import Counter
sys.path.insert(0, '.')
import fw, search
from search import le, join, meet
from lsearch import LS

def close(L, seed, limit=400):
    R = dict(seed); stack = list(seed.items())
    fn = lambda g: R.get(g, L.f[g])
    while stack:
        g, U = stack.pop()
        if R.get(g) != U and g in R: continue
        if not (le(L.lo[g], U) and le(U, L.hi[g])): return None
        for h in L.after[g]:
            Vh = fn(h)
            if not le(U, Vh):
                V2 = join(Vh, U); R[h] = V2; stack.append((h, V2))
        for h in L.before[g]:
            Vh = fn(h)
            if not le(Vh, U):
                V2 = meet(Vh, U); R[h] = V2; stack.append((h, V2))
        for x in L.regs[g]:
            for h in L.G[x]:
                if h == g or h in L.after[g] or h in L.before[g]: continue
                Vh = fn(h)
                if le(Vh, U) or le(U, Vh): continue
                V2 = join(Vh, U) if len(Vh) >= len(U) else meet(Vh, U)
                R[h] = V2; stack.append((h, V2))
        if len(R) > limit: return None
    return R

if __name__ == '__main__':
    W = search.Word(fw.load(os.environ.get('WORD', 'K8-FULL-WORD.json.gz')))
    L = LS(W); c0 = L.total(); print('cost0', c0, flush=True)
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0; minutes = float(sys.argv[2]) if len(sys.argv) > 2 else 5
    rng = random.Random(seed)
    stats = Counter(); hist = Counter()
    for x in range(W.n):
        for S in L.moves(x):
            R = close(L, S)
            if R is None: stats['closure failed'] += 1; continue
            d = L.check(R)
            if d is None: stats['illegal after closure'] += 1; continue
            stats['legal'] += 1; hist[round(d, 2)] += 1
            if d < -1e-9: L.apply(R); stats['improving'] += 1; print('improve', x, round(d, 4), len(R), flush=True)
    print('closure pass', dict(stats), 'cost', L.total(), flush=True)
    print('deltas', sorted(hist.items())[:25], flush=True)
    t0 = time.time(); acc = 0
    while time.time() - t0 < minutes * 60:
        x = rng.randrange(W.n); mv = L.moves(x)
        if not mv: continue
        R = close(L, rng.choice(mv))
        if R is None: continue
        d = L.check(R)
        if d is None or d > 1e-9: continue
        L.apply(R); acc += 1
        if d < -1e-9: print(f'{time.time()-t0:.0f}s NEW BEST {L.total():.4f} gain {c0-L.total():.4f}', flush=True)
    print('walk accepted', acc, 'final', L.total(), 'gain', c0 - L.total(), flush=True)
    import pickle; pickle.dump(L.f, open(f'ls2_final_{seed}.pkl', 'wb'))
    w2 = L.to_word(); fw.save(w2, f'W-ls2-{seed}.json.gz'); print('saved', w2['report']['local_profile'], flush=True)
