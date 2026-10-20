"""t-egraph: commuting-gate reorders of the k8 word, each re-framed by an (local) exact MILP.

Rewrite rules (scalar-exact in both orientations, since the dual applies the same gate list transposed):
  R1 add(t1,s1) and add(t2,s2) commute iff t1 != s2 and t2 != s1   (two adds into one accumulator commute)
  R2 scale(t) and add(t2,s2) commute iff t not in {t2,s2};  scale/scale always commute
Moves stay inside a section (A, B or tail), so the cut, compensation and every dirt/cone map are unchanged.
After a reorder the frames are re-optimised: positions whose legal interval [P,Q] changed, plus the moved gates,
are freed together with the chains through them, and a MILP over the pool candidates picks the cheapest frames.
"""
import sys, time, math, random, json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
from collections import defaultdict, Counter
sys.path.insert(0, '.')
import fw, search
from search import le, join, meet, FC

class OW:
    def __init__(self, W):
        self.W = W
        self.N = len(W.seq)
        self.kind = [e[0] for _, _, e in W.seq]
        self.regs = W.regs
        self.sec = W.sec
        self.order = list(range(self.N))          # position -> gid
        self.f = {g: W.f[g] for g in range(self.N)}  # gid -> frame
        self.refresh()

    def commute(self, g, h):
        a, b = self.kind[g], self.kind[h]
        rg, rh = self.regs[g], self.regs[h]
        if a == 'scale' and b == 'scale': return True
        if a == 'scale': return rg[0] not in rh
        if b == 'scale': return rh[0] not in rg
        return rg[0] != rh[1] and rh[0] != rg[1]

    def refresh(self):
        W = self.W
        self.posof = {g: i for i, g in enumerate(self.order)}
        self.chain = defaultdict(list)            # register -> gids in order
        for g in self.order:
            for r in self.regs[g]: self.chain[r].append(g)
        self.cidx = {}
        for r, gs in self.chain.items():
            for i, g in enumerate(gs): self.cidx[(r, g)] = i
        self.frozen = {}
        for q, U in W.cframe.items():
            ga = [g for g in self.chain[q] if self.sec[g] == 'A']
            if ga: self.frozen[ga[-1]] = U
        self.PQ()

    def PQ(self):
        W = self.W
        cur = list(W.start); P = {}
        for g in self.order:
            rr = self.regs[g]; U = ()
            for r in rr: U = join(U, cur[r])
            if g in self.frozen: U = join(U, self.frozen[g])
            P[g] = U
            for r in rr: cur[r] = U
        cur = list(W.end); Q = {}
        for g in reversed(self.order):
            rr = self.regs[g]; U = W.full
            for r in rr:
                U = meet(U, cur[r])
                if self.sec[g] != 'T': U = meet(U, W.ub[r])
                if self.sec[g] == 'A' and r in W.cframe: U = meet(U, W.cframe[r])
            if g in self.frozen: U = meet(U, self.frozen[g])
            Q[g] = U
            for r in rr: cur[r] = U
        self.P, self.Q = P, Q
        self.feasible = all(le(P[g], Q[g]) for g in self.order)

    def ok_at(self, g, U):
        W = self.W
        for r in self.regs[g]:
            if not le(U, W.end[r]) or not le(W.start[r], U): return False
            if self.sec[g] != 'T' and not le(U, W.ub[r]): return False
            if self.sec[g] == 'A' and r in W.cframe and not le(U, W.cframe[r]): return False
        if g in self.frozen and U != self.frozen[g]: return False
        return True

    def reg_cost(self, r, fr):
        W = self.W; c = 0.0; last = W.start[r]
        for g in self.chain.get(r, ()):
            U = fr[g]
            if U != last: c += FC[len(U) - len(last)]; last = U
        if W.end[r] != last: c += FC[len(W.end[r]) - len(last)]
        return c

    def total_cost(self):
        return sum(self.reg_cost(r, self.f) for r in range(self.W.n))

    def profile(self):
        W = self.W; p = Counter()
        for r in range(W.n):
            last = W.start[r]
            for g in self.chain.get(r, ()):
                U = self.f[g]
                if U != last: p[len(U) - len(last)] += 1; last = U
            if W.end[r] != last: p[len(W.end[r]) - len(last)] += 1
        return p

    def local_milp(self, free, time_limit=60):
        """Re-optimise frames of gids in `free` (others fixed). Returns (new frames dict for free, delta) or None."""
        W = self.W
        free = set(free)
        # fixed positions must lie in [P,Q] and be legal; otherwise free them
        while True:
            regsF = {r for g in free for r in self.regs[g]}
            add = set()
            for r in regsF:
                for g in self.chain[r]:
                    if g not in free:
                        U = self.f[g]
                        if not (le(self.P[g], U) and le(U, self.Q[g]) and self.ok_at(g, U)): add.add(g)
            if not add: break
            free |= add
            if len(free) > 3000: return None
        # also fixed positions with frame outside [P,Q] anywhere else make the order infeasible for them; check globally
        for g in self.order:
            if g not in free:
                U = self.f[g]
                if not (le(self.P[g], U) and le(U, self.Q[g])):
                    free.add(g)
        regsF = {r for g in free for r in self.regs[g]}
        # candidates
        pool = {}
        for r in regsF:
            s = {W.start[r], W.end[r]}
            for g in self.chain[r]: s.update((self.P[g], self.Q[g], self.f[g]))
            pool[r] = s
        cand = {}
        for g in free:
            s = set()
            for r in self.regs[g]: s |= pool.get(r, set())
            c = [U for U in s if le(self.P[g], U) and le(U, self.Q[g]) and self.ok_at(g, U)]
            if not c: return None
            cand[g] = sorted(c, key=lambda U: (len(U), U))
        xid = {}; nv = 0
        for g in sorted(free):
            for U in cand[g]: xid[(g, U)] = nv; nv += 1
        rows, cols, vals, lo, hi, ycost = [], [], [], [], [], []
        def addrow(ent, l, h):
            rr = len(lo)
            for c_, v_ in ent: rows.append(rr); cols.append(c_); vals.append(v_)
            lo.append(l); hi.append(h)
        for g in free: addrow([(xid[(g, U)], 1) for U in cand[g]], 1, 1)
        for r in regsF:
            seqp = [('S', None)] + [(('k' if g in free else 'F'), g) for g in self.chain[r]] + [('E', None)]
            for (ta, a), (tb, b) in zip(seqp, seqp[1:]):
                CA = [W.start[r]] if ta == 'S' else ([self.f[a]] if ta == 'F' else cand[a])
                CB = [W.end[r]] if tb == 'E' else ([self.f[b]] if tb == 'F' else cand[b])
                ys = [(U, V) for U in CA for V in CB if le(U, V)]
                if not ys: return None
                base = nv + len(ycost)
                for (U, V) in ys: ycost.append(FC[len(V) - len(U)] if U != V else 0.0)
                idx = {uv: base + i for i, uv in enumerate(ys)}
                for U in CA:
                    ent = [(idx[(U, V)], 1) for V in CB if (U, V) in idx]
                    if ta == 'k': addrow(ent + [(xid[(a, U)], -1)], 0, 0)
                    else: addrow(ent, 1, 1)
                for V in CB:
                    ent = [(idx[(U, V)], 1) for U in CA if (U, V) in idx]
                    if tb == 'k': addrow(ent + [(xid[(b, V)], -1)], 0, 0)
                    else: addrow(ent, 1, 1)
        ntot = nv + len(ycost)
        cvec = np.zeros(ntot); cvec[nv:] = ycost
        A = coo_matrix((vals, (rows, cols)), shape=(len(lo), ntot)).tocsr()
        integ = np.zeros(ntot); integ[:nv] = 1
        res = milp(cvec, constraints=LinearConstraint(A, lo, hi), integrality=integ, bounds=Bounds(0, 1),
                   options=dict(time_limit=time_limit))
        if res.x is None or res.status not in (0,): return None
        newf = {g: max(cand[g], key=lambda U: res.x[xid[(g, U)]]) for g in free}
        return newf, res.fun, regsF

    def move(self, g, newpos):
        """Move gid g to position newpos (in the order list without g). Returns the old order for undo."""
        old = list(self.order)
        self.order.remove(g); self.order.insert(newpos, g)
        return old

    def try_move(self, g, newpos, accept=lambda d: d < -1e-9):
        oldorder = list(self.order); oldP, oldQ = self.P, self.Q; oldfroz = self.frozen
        oldcost_regs = None
        self.order.remove(g); self.order.insert(newpos, g)
        self.refresh()
        if not self.feasible:
            self.order = oldorder; self.refresh(); return None
        changed = {h for h in self.order if self.P[h] != oldP.get(h) or self.Q[h] != oldQ.get(h)}
        changed.add(g)
        r = self.local_milp(changed)
        if r is None:
            self.order = oldorder; self.refresh(); return None
        newf, obj, regsF = r
        # delta: cost over regsF with old order/frames vs new
        newfr = dict(self.f); newfr.update(newf)
        newc = sum(self.reg_cost(q, newfr) for q in regsF)
        # old cost of regsF under old order
        curorder = self.order; self.order = oldorder; self.refresh_chains_only()
        oldc = sum(self.reg_cost(q, self.f) for q in regsF)
        self.order = curorder; self.refresh_chains_only()
        d = newc - oldc
        if accept(d):
            self.f = newfr; self.refresh(); return d
        self.order = oldorder; self.refresh(); return None

    def refresh_chains_only(self):
        self.chain = defaultdict(list)
        for g in self.order:
            for r in self.regs[g]: self.chain[r].append(g)

    def candidate_moves(self):
        """For each register, consecutive gates g1<g2 on its chain in one section that commute: move g2 just
        before g1 (if it commutes with everything in between), or g1 just after g2."""
        mv = []
        for r, gs in self.chain.items():
            for g1, g2 in zip(gs, gs[1:]):
                if self.sec[g1] != self.sec[g2] or not self.commute(g1, g2): continue
                i, j = self.posof[g1], self.posof[g2]
                between = self.order[i + 1:j]
                if all(self.commute(g2, h) for h in between): mv.append((g2, i))           # g2 to position i
                if all(self.commute(g1, h) for h in between): mv.append((g1, j))  # g1 just after g2
        return mv

    def to_word(self):
        import copy
        W = self.W; w = copy.deepcopy(W.w)
        A, B, T = [], [], []
        for g in self.order:
            s, i, e = W.seq[g]
            U = list(self.f[g])
            if s == 'A': A.append(W.w['A'][i][:-1] + [U])
            elif s == 'B': B.append(W.w['B'][i][:-1] + [U])
            else: T.append([W.w['early_cleanup'][i][0], U])
        w['A'], w['B'] = A, B
        # tail order: early_cleanup lists (cleanup index, U); its climbed order is the list order
        w['early_cleanup'] = T
        p = fw.paid_profile(w)
        w['report']['local_profile'] = {str(d): c for d, c in sorted(p.items())}
        w['report']['local_rank_mass'] = sum(d * c for d, c in p.items())
        return w
