"""t-egraph: reorder moves v2 -- long-range window moves and moves across the A/B cut.

The cut (targets += c*centre) is a fixed block between A and B, represented by a marker gid CUT.  A gate may cross it
iff it commutes with every cut gate: it must not write a centre (the cut reads centres) and must not touch a target.
A gate's section (A or B) is its side of the marker; the tail (early cleanup) stays a separate block.
Every move keeps the scalar map of A+cut+B (so compensation, dirt, cones are unchanged); frames re-optimised by MILP.
"""
import sys, random
sys.path.insert(0, '.')
import reorder
from reorder import OW
from collections import defaultdict

class OW2(OW):
    def __init__(self, W):
        self.W = W
        self.N = len(W.seq)
        self.CUT = self.N
        self.kind = [e[0] for _, _, e in W.seq] + ['cut']
        self.regs = list(W.regs) + [()]
        self.sec0 = list(W.sec) + ['C']
        self.sec = list(self.sec0)
        nA = W.nA
        self.order = list(range(nA)) + [self.CUT] + list(range(nA, self.N))
        self.f = {g: W.f[g] for g in range(self.N)}
        self.f[self.CUT] = ()
        v = W.v; centres = set(W.cframe)
        cutw = W.w['cut']
        cut_t = {e[1] for e in cutw}; cut_s = {e[2] for e in cutw}
        self.cross_ok = []
        for g in range(self.N):
            rr = self.regs[g]
            ok = self.sec0[g] in ('A', 'B') and not any(v <= r < 2 * v for r in rr)
            if ok:
                if self.kind[g] == 'add':
                    t, s = rr
                    ok = t not in cut_s and t not in cut_t and s not in cut_t
                else:
                    ok = rr[0] not in cut_s and rr[0] not in cut_t
            self.cross_ok.append(ok)
        self.cross_ok.append(False)
        self.refresh()

    def commute(self, g, h):
        if g == self.CUT: return self.cross_ok[h]
        if h == self.CUT: return self.cross_ok[g]
        return OW.commute(self, g, h)

    def refresh(self):
        # dynamic sections
        side = 'A'
        for g in self.order:
            if g == self.CUT: side = 'B'; continue
            if self.sec0[g] == 'T': self.sec[g] = 'T'
            else: self.sec[g] = side
        OW.refresh(self)

    def refresh_chains_only(self):
        self.chain = defaultdict(list)
        for g in self.order:
            for r in self.regs[g]: self.chain[r].append(g)

    def window(self, g):
        p = self.posof[g]; lo = p; hi = p
        same = lambda h: (self.sec0[h] == 'T') == (self.sec0[g] == 'T') or h == self.CUT
        while lo - 1 >= 0 and same(self.order[lo - 1]) and self.commute(g, self.order[lo - 1]): lo -= 1
        while hi + 1 < len(self.order) and same(self.order[hi + 1]) and self.commute(g, self.order[hi + 1]): hi += 1
        return lo, hi

    def random_window_move(self, rng):
        g = rng.randrange(self.N)
        lo, hi = self.window(g)
        if hi == lo: return None
        p = self.posof[g]
        # meaningful targets: positions that pass at least one gate sharing a register (or the cut)
        rs = set(self.regs[g])
        marks = [i for i in range(lo, hi + 1) if i != p and (self.order[i] == self.CUT or rs & set(self.regs[self.order[i]]))]
        if not marks: return None
        i = rng.choice(marks)
        # new index in the list without g: passing position i
        newpos = i if i < p else i
        return g, newpos

    def to_word(self):
        import copy, fw
        W = self.W; w = copy.deepcopy(W.w)
        A, B, T = [], [], []
        side = 'A'
        for g in self.order:
            if g == self.CUT: side = 'B'; continue
            s, i, e = W.seq[g]
            U = list(self.f[g])
            if s == 'T': T.append([W.w['early_cleanup'][i][0], U]); continue
            src = W.w['A'][i] if s == 'A' else W.w['B'][i]
            (A if side == 'A' else B).append(src[:-1] + [U])
        w['A'], w['B'], w['early_cleanup'] = A, B, T
        p = fw.paid_profile(w)
        w['report']['local_profile'] = {str(d): c for d, c in sorted(p.items())}
        w['report']['local_rank_mass'] = sum(d * c for d, c in p.items())
        return w
