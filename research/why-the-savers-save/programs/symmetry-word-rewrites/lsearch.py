"""t-egraph: label-space search -- frames without a fixed order.

Observation (relist.py): any order sorted by the frames gives the same exact optimum, so the real object is the
frame labelling f(g).  A labelling is realisable by some legal order iff
  (i)  for every register x the frames of its gates are pairwise comparable (a chain) inside [start(x), end(x)];
  (ii) f(g) <= f(h) for every non-commuting pair g -> h on a common register (kept in their original order), and for
       gates on a common register in different sections (A < B < tail; sections kept);
  (iii) static bounds: low-end cones (A/B gates), centre frames (A gates <= centre frame, B gates >= it, and the last
       A frame of a centre equals it).
Then the order 'sort by dim, ties topologically' realises it (each chain is sorted by dim; DAG edges go up in dim or
stay equal).  The cost is the sum over registers of the climbs of its chain of distinct frames.
Moves: relabel one gate to a frame used on its registers; merge a register's two consecutive frames (all its gates at
one frame move to the other).  Every winner is rebuilt as a word (order by frames) and audited by the hub's checker.
"""
import sys, os, time, random, json
from collections import defaultdict, Counter
sys.path.insert(0, '.')
import fw, search, reorder2
from search import le, join, meet, FC

class LS:
    def __init__(self, W):
        self.W = W; O = reorder2.OW2(W); self.O = O
        self.N = O.N; self.regs = W.regs; self.sec = list(W.sec)
        self.f = [W.f[g] for g in range(self.N)]
        self.G = defaultdict(list)
        for g in range(self.N):
            for x in self.regs[g]: self.G[x].append(g)
        rank = {'A': 0, 'B': 1, 'T': 2}
        self.before = defaultdict(set); self.after = defaultdict(set)   # g -> gates that must be <= / >= in frame
        for x, gs in self.G.items():
            for i in range(len(gs)):
                for j in range(i + 1, len(gs)):
                    a, b = gs[i], gs[j]
                    if rank[self.sec[a]] < rank[self.sec[b]] or not O.commute(a, b):
                        self.after[a].add(b); self.before[b].add(a)
        self.lo = []; self.hi = []
        for g in range(self.N):
            lo = (); hi = W.full
            for x in self.regs[g]:
                lo = join(lo, W.start[x]); hi = meet(hi, W.end[x])
                if self.sec[g] != 'T': hi = meet(hi, W.ub[x])
                if x in W.cframe:
                    if self.sec[g] == 'A': hi = meet(hi, W.cframe[x])
                    else: lo = join(lo, W.cframe[x])
            self.lo.append(lo); self.hi.append(hi)
        self.centreA = {x: [g for g in self.G[x] if self.sec[g] == 'A'] for x in W.cframe}

    def reg_cost(self, x, f):
        W = self.W
        S = {f(g) for g in self.G.get(x, ())} | {W.start[x], W.end[x]}
        S = sorted(S, key=len)
        c = 0.0
        for a, b in zip(S, S[1:]):
            if len(a) == len(b) or not le(a, b): return None
            c += FC[len(b) - len(a)]
        return c

    def reg_profile(self, x):
        W = self.W
        S = sorted({self.f[g] for g in self.G.get(x, ())} | {W.start[x], W.end[x]}, key=len)
        return Counter(len(b) - len(a) for a, b in zip(S, S[1:]))

    def total(self):
        return sum(self.reg_cost(x, lambda g: self.f[g]) for x in range(self.W.n))

    def check(self, R):
        """R: {g: U}.  Returns cost delta or None if illegal."""
        fn = lambda g: R.get(g, self.f[g])
        W = self.W
        for g, U in R.items():
            if not (le(self.lo[g], U) and le(U, self.hi[g])): return None
            for h in self.before[g]:
                if not le(fn(h), U): return None
            for h in self.after[g]:
                if not le(U, fn(h)): return None
        touched = {x for g in R for x in self.regs[g]}
        for x in touched:
            if x in W.cframe and self.centreA[x]:
                if max((len(fn(g)) for g in self.centreA[x])) != len(W.cframe[x]): return None
        d = 0.0
        for x in touched:
            c1 = self.reg_cost(x, fn)
            if c1 is None: return None
            d += c1 - self.reg_cost(x, lambda g: self.f[g])
        return d

    def apply(self, R):
        for g, U in R.items(): self.f[g] = U

    def moves(self, x):
        W = self.W
        by = defaultdict(list)
        for g in self.G.get(x, ()): by[self.f[g]].append(g)
        S = sorted(set(by) | {W.start[x], W.end[x]}, key=len)
        mv = []
        for i, U in enumerate(S):
            if U not in by: continue
            if i + 1 < len(S): mv.append({g: S[i + 1] for g in by[U]})
            if i > 0: mv.append({g: S[i - 1] for g in by[U]})
        # single-gate relabels to any frame of x's chain
        for g in self.G.get(x, ()):
            for U in S:
                if U != self.f[g]: mv.append({g: U})
        return mv

    def to_word(self):
        """Order: sections kept; inside, sort by frame dim with ties broken topologically (original index)."""
        import copy, heapq
        W = self.W; w = copy.deepcopy(W.w)
        out = {}
        for s in ('A', 'B', 'T'):
            gs = [g for g in range(self.N) if self.sec[g] == s]
            gset = set(gs)
            indeg = {g: 0 for g in gs}
            for g in gs:
                for h in self.after[g]:
                    if h in gset: indeg[h] += 1
            heap = [(len(self.f[g]), g) for g in gs if indeg[g] == 0]; heapq.heapify(heap); order = []
            while heap:
                _, g = heapq.heappop(heap); order.append(g)
                for h in self.after[g]:
                    if h in gset:
                        indeg[h] -= 1
                        if indeg[h] == 0: heapq.heappush(heap, (len(self.f[h]), h))
            assert len(order) == len(gs)
            out[s] = order
        A, B, T = [], [], []
        for g in out['A']: A.append(W.w['A'][W.seq[g][1]][:-1] + [list(self.f[g])])
        for g in out['B']: B.append(W.w['B'][W.seq[g][1]][:-1] + [list(self.f[g])])
        for g in out['T']: T.append([W.w['early_cleanup'][W.seq[g][1]][0], list(self.f[g])])
        w['A'], w['B'], w['early_cleanup'] = A, B, T
        p = fw.paid_profile(w)
        w['report']['local_profile'] = {str(d): c for d, c in sorted(p.items())}
        w['report']['local_rank_mass'] = sum(d * c for d, c in p.items())
        return w

if __name__ == '__main__':
    W = search.Word(fw.load(os.environ.get('WORD', 'K8-FULL-WORD.json.gz')))
    L = LS(W)
    c0 = L.total(); print('cost0', c0, flush=True)
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0; minutes = float(sys.argv[2]) if len(sys.argv) > 2 else 5
    rng = random.Random(seed)
    # greedy full pass
    stats = Counter()
    for x in range(W.n):
        for R in L.moves(x):
            d = L.check(R)
            if d is None: stats['illegal'] += 1; continue
            stats['legal'] += 1
            if d < -1e-9: L.apply(R); stats['improving'] += 1; print('improve', x, round(d, 4), flush=True)
            elif abs(d) < 1e-9: stats['neutral'] += 1
            else: stats['worse'] += 1
    print('greedy pass', dict(stats), 'cost', L.total(), flush=True)
    # plateau walk
    t0 = time.time(); acc = 0; best = L.total()
    while time.time() - t0 < minutes * 60:
        x = rng.randrange(W.n); mv = L.moves(x)
        if not mv: continue
        R = rng.choice(mv); d = L.check(R)
        if d is None or d > 1e-9: continue
        L.apply(R); acc += 1
        if d < -1e-9:
            best = L.total(); print(f'{time.time()-t0:.0f}s NEW BEST {best:.4f} gain {c0-best:.4f}', flush=True)
    print('walk accepted', acc, 'final', L.total(), flush=True)
    import pickle; pickle.dump(L.f, open(f'ls_final_{seed}.pkl', 'wb'))
    w2 = L.to_word(); fw.save(w2, f'W-ls-{seed}.json.gz'); print('saved', w2['report']['local_profile'])
