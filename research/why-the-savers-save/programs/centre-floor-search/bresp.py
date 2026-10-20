#!/usr/bin/env python3
"""ctr-search / bresp.py: best-response dynamics with kicks for dictionary (one-sided normal form) decoders.
Each port in turn takes the representation of least marginal one-sided cost (number of its directions whose current
user-label span does not already contain the port's label); kicks force random ports to random representations;
plateau moves accepted.  Final decoders are rebuilt explicitly and checked exactly (dsearch.verify)."""
import sys, os, json, time, random, pickle, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dsearch import *

def span_elems(labs):
    basis = []
    for v in labs:
        for x in basis: v = min(v, v ^ x)
        if v: basis.append(v)
    els = {0}
    for x in basis: els |= {e ^ x for e in els}
    return els, len(basis)

class BR:
    def __init__(s, P, reps, seed=0):
        s.P, s.reps = P, reps; s.p = len(P)
        s.rng = random.Random(seed)
        s.users = [dict() for _ in range(s.p)]
        s.span = [({0}, 0) for _ in range(s.p)]
        s.choice = [None] * N
    def _recalc(s, k):
        s.span[k] = span_elems(list(s.users[k].keys()))
    def remove(s, b):
        F = s.reps[b][s.choice[b]]
        for k in F:
            u = LABS[b]; s.users[k][u] -= 1
            if not s.users[k][u]: del s.users[k][u]
            s._recalc(k)
    def add(s, b, c):
        s.choice[b] = c
        for k in s.reps[b][c]:
            u = LABS[b]; s.users[k][u] = s.users[k].get(u, 0) + 1
            if u not in s.span[k][0]: s._recalc(k)
    def best(s, b):
        u = LABS[b]
        bestv, bl = None, []
        for c, F in enumerate(s.reps[b]):
            v = 0
            for k in F:
                if u not in s.span[k][0]: v += 1
            # tie-break: prefer representations whose directions already have users
            if bestv is None or v < bestv: bestv, bl = v, [c]
            elif v == bestv: bl.append(c)
        return s.rng.choice(bl), bestv
    def total(s): return sum(sp[1] for sp in s.span)
    def sweep(s):
        order = list(range(N)); s.rng.shuffle(order)
        for b in order:
            s.remove(b)
            c, _ = s.best(b)
            s.add(b, c)
        return s.total()
    def init(s, choice=None):
        s.users = [dict() for _ in range(s.p)]; s.span = [({0}, 0) for _ in range(s.p)]
        for b in range(N):
            s.add(b, choice[b] if choice is not None else s.rng.randrange(len(s.reps[b])))
        return s.total()

def run(P, reps, seed, rounds, kick, init=None, log=print):
    br = BR(P, reps, seed)
    br.init(init)
    cur = br.total()
    for _ in range(3): cur = br.sweep()
    best, bestc = cur, list(br.choice)
    for r in range(rounds):
        saved = list(br.choice)
        ks = br.rng.sample(range(N), kick)
        for b in ks:
            br.remove(b); br.add(b, br.rng.randrange(len(reps[b])))
        prev = None
        for _ in range(6):
            cur = br.sweep()
            if cur == prev: break
            prev = cur
        if cur <= best:
            if cur < best: log(f'   seed {seed} round {r}: {cur}')
            best, bestc = cur, list(br.choice)
        else:
            br.init(saved)
    return best, bestc

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--fams', default='FDR')
    ap.add_argument('--maxsize', type=int, default=3)
    ap.add_argument('--seeds', type=int, default=4)
    ap.add_argument('--rounds', type=int, default=200)
    ap.add_argument('--kick', type=int, default=6)
    ap.add_argument('--fromcoord', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    P = pool(a.fams)
    cache = os.path.join(HERE, f'out/reps_{a.fams}_{a.maxsize}.pkl')
    if os.path.exists(cache):
        reps = pickle.load(open(cache, 'rb'))
    else:
        t0 = time.time(); reps = reps_for(P, maxsize=a.maxsize)
        pickle.dump(reps, open(cache, 'wb')); print(f'reps built ({time.time()-t0:.0f} s)', flush=True)
    print(f'pool {a.fams}: {len(P)} directions; reps per port {min(map(len,reps))}-{max(map(len,reps))}', flush=True)
    init = None
    if a.fromcoord:
        frame = [P.index(tuple(1 if j == i else 0 for j in range(8))) for i in range(8)]
        init = []
        for b in range(N):
            T = TET[TETOF[b]]
            want = tuple(sorted(frame[i] for i in T))
            init.append(next(c for c, F in enumerate(reps[b]) if tuple(sorted(F)) == want))
    res = []
    t0 = time.time()
    for sd in range(a.seeds):
        best, bc = run(P, reps, sd, a.rounds, a.kick, init=init)
        v = verify(P, reps, bc, f'{a.fams} seed {sd}')
        print(f'seed {sd}: best one-sided {best}; verified exact {v["exact"]} n {v["n"]} pair {v["pair"]} star {v["star_gather"]} '
              f'{v.get("REFUSED","")} terms {sorted(v["star_terms"], reverse=True)[:24]} ({time.time()-t0:.0f} s)', flush=True)
        res.append(v)
    if a.out: json.dump(res, open(a.out, 'w'), indent=1, default=str)
