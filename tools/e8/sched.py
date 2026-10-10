#!/usr/bin/env python3
"""sched.py -- re-schedule a gcert/1 E8 certificate: same single additions, other frames.

The certificate is cut into its single additions  tgt += (num/den) * src.  For every register the additions that
name it form runs of reads (the register is the source) and runs of writes (it is the target).  Inside one run the
additions commute, so their order is free; between runs the order is kept.  A frame assignment is legal exactly
when, for every register, the frames of each run form a chain under inclusion and every run lies at or above the
run before it (and above the start frame, below the final frame).  The cost of a register is the sum over its
blocks of rank * ln(45 / rank): the first-order weight (figure ~ D / total weight).

Local search with two moves:  push(add, larger frame)  and  pull(add, smaller frame), each dragging along the
additions that would otherwise become illegal (joins / meets).  Nothing here is trusted: the result is rebuilt as
an e8.Design and goes through the two published checkers.
usage: sched.py <in> <out> [seed] [annealing sweeps] [start temperature]
NOTE: the annealing loop stops after budget = 400 seconds of wall-clock time (the recorded run needs about 45)."""
import sys, math, random, time, heapq
sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import e8

M = 45
FW = [0.0] + [r * math.log(M / r) for r in range(1, 10)]


class Frames:
    def __init__(self):
        self.basis, self.mask, self.dim = [], [], []
        self.by_basis, self.by_mask, self._join = {}, {}, {}

    def add(self, basis):
        basis = e8.echelon(basis)
        i = self.by_basis.get(basis)
        if i is None:
            vecs = [0]
            for b in basis:
                vecs += [x ^ b for x in vecs]
            m = 0
            for x in vecs:
                m |= 1 << x
            i = len(self.basis)
            self.basis.append(basis); self.mask.append(m); self.dim.append(len(basis))
            self.by_basis[basis] = i; self.by_mask[m] = i
        return i

    def leq(self, a, b):
        return self.mask[a] & self.mask[b] == self.mask[a]

    def meet(self, a, b):
        m = self.mask[a] & self.mask[b]
        i = self.by_mask.get(m)
        if i is None:
            i = self.add([x for x in range(1, 512) if m >> x & 1])
        return i

    def join(self, a, b):
        if a == b:
            return a
        k = (a, b) if a < b else (b, a)
        i = self._join.get(k)
        if i is None:
            i = self._join[k] = self.add(self.basis[a] + self.basis[b])
        return i


class Sched:
    def __init__(self, c, split_totals_fixed=True):
        self.c = c
        self.v, self.R, self.h = c["v"], c["R"], c["h"]
        v = self.v
        self.nreg = 2 * v + self.R
        self.F = Frames()
        F = self.F
        fid = [F.add(tuple(f)) for f in c["frames"]]
        self.ZERO, self.FULL = fid[0], fid[1]
        self.start = [fid[f] for f in c["start"]]
        self.final = [fid[f] for f in c["final"]]
        self.tot = set(s for _, s, _ in c["ret"])
        self.tgt, self.src, self.coef, self.fr, self.fixed = [], [], [], [], []
        for ph in "AB":
            for g in c[ph]:
                if g[0] == "out":
                    self.dropped_stops = getattr(self, "dropped_stops", 0) + len(g[4]) + (0 if g[3] else 1)
                    for t, n, d in g[3]:
                        self._add(t, g[2], n, d, fid[g[1]])
                elif g[0] == "in":
                    for s, n, d in g[3]:
                        self._add(g[2], s, n, d, fid[g[1]])
                else:
                    raise ValueError(g[0])
        self.n = len(self.tgt)
        # runs per register
        self.segs = [[] for _ in range(self.nreg)]
        self.segidx = [dict() for _ in range(self.nreg)]
        last = [None] * self.nreg
        for g in range(self.n):
            for p, role in ((self.tgt[g], "w"), (self.src[g], "r")):
                if last[p] != role:
                    self.segs[p].append([]); last[p] = role
                self.segs[p][-1].append(g)
                self.segidx[p][g] = len(self.segs[p]) - 1
        # ---- phase rule.  B0 = the additions that, in the input, must come after an addition into a target
        # (order between runs + chain order inside runs).  Rule kept by every move: inside one run an addition
        # outside B0 never stands at a frame strictly above an addition of B0.  Then nothing outside B0 ever has
        # to wait for a target, so all of it (the totals included) can be placed before the scatter.
        self.mixed = [[False] * len(sg) for sg in self.segs]
        self.inB = [False] * self.n
        self.cost = [self.reg_eval(r, {}) for r in range(self.nreg)]
        assert all(x is not None for x in self.cost), "the input certificate is not legal in my model"
        succ = self.arcs()
        st = [g for g in range(self.n) if v <= self.tgt[g] < 2 * v]
        while st:
            a = st.pop()
            if self.inB[a]:
                continue
            self.inB[a] = True; st.extend(succ[a])
        assert not any(self.inB[g] and self.fixed[g] for g in range(self.n)), "a total waits for a target in the input"
        for r in range(self.nreg):
            for i, sg in enumerate(self.segs[r]):
                k = sum(1 for g in sg if self.inB[g])
                self.mixed[r][i] = 0 < k < len(sg)
        assert all(self.reg_eval(r, {}) is not None for r in range(self.nreg))

    def arcs(self):
        """succ[a] = additions that must come after a: next group of the same register (runs in order, frames
        growing inside a run)"""
        dim = self.F.dim
        succ = [[] for _ in range(self.n)]
        for r in range(self.nreg):
            groups = []
            for seg in self.segs[r]:
                byf = {}
                for g in seg:
                    byf.setdefault(self.fr[g], []).append(g)
                for f in sorted(byf, key=dim.__getitem__):
                    groups.append(byf[f])
            for A, B in zip(groups, groups[1:]):
                for a in A:
                    succ[a].extend(B)
        return succ

    def _add(self, t, s, n, d, f):
        self.tgt.append(t); self.src.append(s); self.coef.append((n, d)); self.fr.append(f)
        self.fixed.append(t in self.tot)          # the additions into the totals stay at the full frame

    # ---- cost / legality of one register under tentative frames
    def chain(self, r, ov=None):
        """the stops of register r (frame ids, start excluded, final included), or None if illegal"""
        F = self.F; dim = F.dim; mask = F.mask; fr = self.fr
        cur = self.start[r]; out = []
        mixed = self.mixed[r]
        for si, seg in enumerate(self.segs[r]):
            if ov:
                fs = set(ov[g] if g in ov else fr[g] for g in seg)
            else:
                fs = set(fr[g] for g in seg)
            if mixed[si] and len(fs) > 1:
                inB = self.inB
                if ov:
                    hiA = max(dim[ov[g] if g in ov else fr[g]] for g in seg if not inB[g])
                    loB = min(dim[ov[g] if g in ov else fr[g]] for g in seg if inB[g])
                else:
                    hiA = max(dim[fr[g]] for g in seg if not inB[g])
                    loB = min(dim[fr[g]] for g in seg if inB[g])
                if hiA > loB:
                    return None
            for f in (sorted(fs, key=dim.__getitem__) if len(fs) > 1 else fs):
                if f == cur:
                    continue
                if dim[f] <= dim[cur] or mask[cur] & mask[f] != mask[cur]:
                    return None
                out.append(f); cur = f
        f = self.final[r]
        if f != cur:
            if dim[f] <= dim[cur] or mask[cur] & mask[f] != mask[cur]:
                return None
            out.append(f)
        return out

    def reg_eval(self, r, ov):
        ch = self.chain(r, ov)
        if ch is None:
            return None
        dim = self.F.dim
        d0 = dim[self.start[r]]; w = 0.0
        for f in ch:
            w += FW[dim[f] - d0]; d0 = dim[f]
        return w

    def total(self):
        return sum(self.cost)

    def hist(self):
        """block histogram {rank: count} of one invocation, the 8 scratch copies included"""
        dim = self.F.dim; H = {}
        for r in range(self.nreg):
            d0 = dim[self.start[r]]
            for f in self.chain(r):
                H[dim[f] - d0] = H.get(dim[f] - d0, 0) + 1; d0 = dim[f]
        for k, s, f in self.c["ret"]:
            H[self.h] = H.get(self.h, 0) + 1
        return H

    # ---- moves
    def push(self, seeds, cap=400):
        F = self.F; fr = self.fr
        new = dict(seeds); stack = list(seeds)
        while stack:
            a = stack.pop(); Fa = new[a]
            for p in (self.tgt[a], self.src[a]):
                if not F.leq(Fa, self.final[p]):
                    return None
                sa = self.segidx[p][a]; sg = self.segs[p]
                for si in range(sa, len(sg)):
                    for b in sg[si]:
                        if b == a:
                            continue
                        Fb = new[b] if b in new else fr[b]
                        if F.leq(Fa, Fb):
                            continue
                        if si == sa and F.leq(Fb, Fa):
                            continue
                        if self.fixed[b]:
                            return None
                        new[b] = F.join(Fb, Fa); stack.append(b)
                        if len(new) > cap:
                            return None
        return new

    def pull(self, seeds, cap=400):
        F = self.F; fr = self.fr
        new = dict(seeds); stack = list(seeds)
        while stack:
            a = stack.pop(); Fa = new[a]
            for p in (self.tgt[a], self.src[a]):
                if not F.leq(self.start[p], Fa):
                    return None
                sa = self.segidx[p][a]; sg = self.segs[p]
                for si in range(0, sa + 1):
                    for b in sg[si]:
                        if b == a:
                            continue
                        Fb = new[b] if b in new else fr[b]
                        if F.leq(Fb, Fa):
                            continue
                        if si == sa and F.leq(Fa, Fb):
                            continue
                        if self.fixed[b]:
                            return None
                        new[b] = F.meet(Fb, Fa); stack.append(b)
                        if len(new) > cap:
                            return None
        return new

    def delta(self, new):
        regs = set()
        for g in new:
            regs.add(self.tgt[g]); regs.add(self.src[g])
        d = 0.0; nc = {}
        for r in regs:
            w = self.reg_eval(r, new)
            if w is None:
                return None, None
            nc[r] = w; d += w - self.cost[r]
        return d, nc

    def apply(self, new, nc):
        for g, f in new.items():
            self.fr[g] = f
        for r, w in nc.items():
            self.cost[r] = w

    def stops(self, r):
        return self.chain(r)

    def candidates(self, g):
        """frames worth trying for addition g: the stops of its two registers, and meets / joins of pairs of them"""
        F = self.F; f0 = self.fr[g]
        t, s = self.tgt[g], self.src[g]
        ct = [f for f in self.chain(t) + [self.start[t]] if F.dim[f] > 0]
        cs = [f for f in self.chain(s) + [self.start[s]] if F.dim[f] > 0]
        up, dn = set(), set()
        for f in ct + cs:
            if f == f0:
                continue
            if F.leq(f0, f): up.add(f)
            elif F.leq(f, f0): dn.add(f)
        for a in ct:
            for b in cs:
                if a == b:
                    continue
                m = F.meet(a, b)
                if m != f0 and F.dim[m] > 0:
                    if F.leq(f0, m): up.add(m)
                    elif F.leq(m, f0): dn.add(m)
                j = F.join(a, b)
                if j != f0:
                    if F.leq(f0, j): up.add(j)
                    elif F.leq(j, f0): dn.add(j)
        return up, dn

    def try_move(self, seeds, direction, T=0.0, rng=random):
        new = self.push(seeds) if direction > 0 else self.pull(seeds)
        if new is None:
            return 0
        d, nc = self.delta(new)
        if d is None:
            return 0
        if d < -1e-9 or (T > 0 and rng.random() < math.exp(-d / T)) or (T == 0 and abs(d) <= 1e-9 and rng.random() < 0.3):
            self.apply(new, nc)
            return 1 if d < -1e-9 else 2
        return 0

    def sweep(self, T=0.0, rng=random, stopmoves=True):
        """one pass over all additions (single-add moves) and over all stops of all registers (stop moves)"""
        order = [g for g in range(self.n) if not self.fixed[g]]
        rng.shuffle(order)
        acc = 0
        for g in order:
            up, dn = self.candidates(g)
            cands = [(f, 1) for f in up] + [(f, -1) for f in dn]
            rng.shuffle(cands)
            for f, di in cands:
                if self.fr[g] == f:
                    continue
                ok = (self.F.leq(self.fr[g], f) if di > 0 else self.F.leq(f, self.fr[g]))
                if not ok:
                    continue
                if self.try_move({g: f}, di, T, rng) == 1:
                    acc += 1
        if stopmoves:
            regs = list(range(self.nreg)); rng.shuffle(regs)
            for r in regs:
                ch = self.chain(r)
                full = [self.start[r]] + ch
                for i, f in enumerate(ch):
                    mine = [g for seg in self.segs[r] for g in seg if self.fr[g] == f and not self.fixed[g]]
                    if not mine:
                        continue
                    opts = []
                    if i + 1 < len(ch): opts.append((ch[i + 1], 1))
                    if self.F.dim[full[i]] > 0: opts.append((full[i], -1))
                    rng.shuffle(opts)
                    for f2, di in opts:
                        if self.try_move({g: f2 for g in mine}, di, T, rng) == 1:
                            acc += 1
                            break
        return acc

    def ysweep(self, rng=random):
        """compound move: a target reads, at one lower frame G, all the helpers whose previous stop is G"""
        v = self.v; acc = 0
        for y in range(v, 2 * v):
            groups = {}
            for seg in self.segs[y]:
                for g in seg:
                    q = self.src[g]
                    if q < 2 * v:
                        continue
                    ch = [self.start[q]] + self.chain(q)
                    if self.fr[g] in ch:
                        i = ch.index(self.fr[g])
                        if i > 0 and self.F.dim[ch[i - 1]] > 0:
                            groups.setdefault(ch[i - 1], []).append(g)
            for G, gs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
                if len(gs) < 2:
                    continue
                if self.try_move({g: G for g in gs if self.fr[g] != G}, -1, 0.0, rng) == 1:
                    acc += 1
        return acc

    def snapshot(self):
        return list(self.fr), list(self.cost)

    def restore(self, snap):
        self.fr, self.cost = list(snap[0]), list(snap[1])

    # ---- back to a design
    def order(self):
        """a legal order of the additions: the kept order between runs, chain order inside runs"""
        dim = self.F.dim; n = self.n
        succ = self.arcs(); indeg = [0] * n
        for a in range(n):
            for b in succ[a]:
                indeg[b] += 1
        v = self.v
        inB = self.inB

        def key(g):
            t, s, f = self.tgt[g], self.src[g], self.fr[g]
            if v <= t < 2 * v:
                return (3, dim[f], f, t, s, g)
            if inB[g]:
                return (2, dim[f], f, g, 0, 0)
            if s < v:
                return (0, dim[f], f, t, s, g) if t in self.tot else (0, dim[f], f, s, t, g)
            return (1, dim[f], f, g, 0, 0)
        heap = [(key(g), g) for g in range(n) if indeg[g] == 0]
        heapq.heapify(heap)
        out = []
        while heap:
            _, g = heapq.heappop(heap)
            out.append(g)
            for b in succ[g]:
                indeg[b] -= 1
                if indeg[b] == 0:
                    heapq.heappush(heap, (key(b), b))
        assert len(out) == n, "cycle in the order: should be impossible"
        return out

    def design(self, name="sched.py", merge=True):
        from fractions import Fraction as Fr
        v = self.v; F = self.F
        od = self.order()
        lastA = max(i for i, g in enumerate(od) if self.fixed[g])
        assert all(not (v <= self.tgt[g] < 2 * v) for g in od[:lastA + 1]), "a target is named before the scatter"
        d = e8.Design(name=name)
        d.R = self.R
        i = 0
        while i < len(od):
            g = od[i]; f = self.fr[g]; ph = "A" if i <= lastA else "B"
            lim = lastA + 1 if i <= lastA else len(od)
            j = i + 1
            if merge:
                # run of additions into the same register
                used = {self.tgt[g], self.src[g]}
                while j < lim and self.fr[od[j]] == f and self.tgt[od[j]] == self.tgt[g] and self.src[od[j]] not in used:
                    used.add(self.src[od[j]]); j += 1
                if j == i + 1:
                    used = {self.tgt[g], self.src[g]}
                    while j < lim and self.fr[od[j]] == f and self.src[od[j]] == self.src[g] and self.tgt[od[j]] not in used:
                        used.add(self.tgt[od[j]]); j += 1
                    if j > i + 1:
                        d.out(ph, F.basis[f], self.src[g], [(self.tgt[a], Fr(*self.coef[a])) for a in od[i:j]])
                        i = j
                        continue
            d.inn(ph, F.basis[f], self.tgt[g], [(self.src[a], Fr(*self.coef[a])) for a in od[i:j]])
            i = j
        for k, s, f in self.c["ret"]:
            d.retain(s)
        d.scat = ("table", self.c["scat"]["table"])
        return d


def polish(S, rng, verbose=True, t0=None, maxsweeps=40):
    t0 = t0 or time.time()
    for it in range(maxsweeps):
        acc = S.sweep(0.0, rng) + S.ysweep(rng)
        if verbose:
            print("  greedy sweep %d: improving moves %d weight %.2f  %.0fs" % (it, acc, S.total(), time.time() - t0)); sys.stdout.flush()
        if acc == 0:
            break


def optimise(c, seed=0, sweeps=40, T0=0.0, verbose=True, stopmoves=True, budget=400):
    """greedy to a local optimum; if T0 > 0, then `sweeps` annealing sweeps from T0 down to 0, each followed by nothing,
    and a final greedy polish; the best state seen is returned"""
    rng = random.Random(seed)
    S = c if isinstance(c, Sched) else Sched(c)
    t0 = time.time()
    if verbose:
        print("start: adds %d weight %.2f (dropped pure stops / extras: %d)" % (S.n, S.total(), getattr(S, "dropped_stops", 0))); sys.stdout.flush()
    polish(S, rng, verbose, t0)
    best = S.total(); snap = S.snapshot()
    if T0 > 0:
        for it in range(sweeps):
            T = T0 * (1 - it / float(sweeps))
            S.sweep(T, rng, stopmoves)
            if verbose and it % 5 == 0:
                print("anneal %d T %.3f: weight %.2f best %.2f frames %d  %.0fs" % (it, T, S.total(), best, len(S.F.basis), time.time() - t0)); sys.stdout.flush()
            if it % 5 == 4 or it == sweeps - 1:
                keep = S.snapshot()
                polish(S, rng, False, t0, 6)
                if S.total() < best - 1e-9:
                    best = S.total(); snap = S.snapshot()
                    if verbose: print("   new best after polish: %.2f" % best); sys.stdout.flush()
                S.restore(keep)
            if time.time() - t0 > budget:
                break
        S.restore(snap)
        polish(S, rng, verbose, t0)
    return S


if __name__ == "__main__":
    src = sys.argv[1]; dst = sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    sweeps = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    T0 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    c = e8.load(src)
    S = optimise(c, seed, sweeps, T0)
    H = S.hist()
    f = e8.figure(H, 9, S.v, S.R, c["cst"])
    print("model: weight %.2f predicted figure %d hist %s" % (S.total(), f["whole_block"], dict(sorted(H.items()))))
    # This text becomes the field `derived_from` of the certificate.  The published certificate was written with the
    # working name of this pass ("i-stops") and of its input file in it; both are kept, because any other text gives
    # the same circuit with another checksum.
    d = S.design(name="i-stops: re-scheduled %s (seed %d)" % (src.split("/")[-1], seed))
    cc = e8.build_cert(d)
    q = e8.quick(cc)
    print({k: q.get(k) for k in ("accepted_ref", "accepted_mirror", "ref_msg", "mirror_msg", "R", "cst", "whole_block", "adds", "most_regs_one_gate")})
    print("gates A %d B %d blocks %s" % (len(cc["A"]), len(cc["B"]), cc["blocks"]))
    if q["accepted_ref"] and q["accepted_mirror"]:
        e8.dump(cc, dst)
        print("wrote", dst)
