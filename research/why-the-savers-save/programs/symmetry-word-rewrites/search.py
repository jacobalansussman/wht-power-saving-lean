"""t-egraph: frame re-assignment search on the k8 word (gate order fixed, later with commuting reorders).

State = one frame per climbed position (A, B, tail).  Legality enforced here (and re-checked by audit_local):
  * every register's chain start <= f(p1) <= f(p2) <= ... <= end  (subspace containment)
  * targets untouched in A (true in the word; frames of B/T only), centres end A exactly at their centre frame
  * low-end helper cones: endB[k] <= end[q] for k in the inverse-row / dirt-column cone of each low-end q
Cost = sum over climbs of d*ln(73/d) (monotone proxy of b at fixed mass); exact b from the pricer at the end.
Moves: raise a run of equal frames of one register to its next frame (forward join propagation), or lower a
run to its previous frame (backward meet propagation); single-gate retargets to a partner's neighbour frame.
"""
import sys, math, random, json, time
from collections import Counter, defaultdict
from fractions import Fraction as Fr
sys.path.insert(0, '.')
import fw

_B = {}
def mask(U):
    m = _B.get(U)
    if m is None:
        els = {0}
        for b in U: els |= {x ^ b for x in els}
        m = 0
        for x in els: m |= 1 << x
        _B[U] = m
    return m
_J = {}
def join(U, V):
    k = (U, V)
    r = _J.get(k)
    if r is None: r = _J[k] = fw.basis(list(U) + list(V))
    return r
_M = {}
def meet(U, V):
    k = (U, V)
    r = _M.get(k)
    if r is None:
        m = mask(U) & mask(V); els = [x for x in range(1, 512) if m >> x & 1]
        r = _M[k] = fw.basis(els)
    return r
def le(U, V): return mask(U) & ~mask(V) == 0

FC = {d: (d * math.log(73 / d) if d else 0.0) for d in range(0, 20)}

class Word:
    def __init__(self, w):
        self.w = w
        v, h, R, n, full, start, end = fw.geometry(w)
        self.v, self.n, self.full, self.start, self.end = v, n, full, start, end
        self.seq = fw.climbed(w)
        self.regs = [tuple(fw.registers(e)) for _, _, e in self.seq]
        self.sec = [s for s, _, _ in self.seq]
        self.f = [fw.basis(e[-1]) for _, _, e in self.seq]
        self.nA = len(w['A'])
        self.pos = defaultdict(list)
        for k, rr in enumerate(self.regs):
            for r in rr: self.pos[r].append(k)
        self.idx = {}
        for r, ps in self.pos.items():
            for i, k in enumerate(ps): self.idx[(r, k)] = i
        # bounds for A/B positions: cones of low ends
        self.ub = [full] * n
        self._cones()
        # centres: A positions <= centre frame; the last A position of a centre is frozen
        self.frozen = set()
        self.cframe = {}
        for q, U in zip(w['centers'], w['center_frames']):
            U = fw.basis(U); self.cframe[q] = U
            pa = [k for k in self.pos[q] if k < self.nA]
            if pa: self.frozen.add(pa[-1])
            else: assert start[q] == U
        for k, rr in enumerate(self.regs):
            if k < self.nA:
                assert not any(v <= r < 2 * v for r in rr), 'target in A'

    def _cones(self):
        from fractions import Fraction as F
        w = self.w; v = self.v; n = self.n
        lowends = set(map(int, w.get('helper_final', {})))
        if not lowends: return
        # prefix = precondition + compensation + main; compute p (forward) and pinv like audit_local
        import importlib.util
        spec = importlib.util.spec_from_file_location('al', 'audit_local_py39.py'); al = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(al)
        ident = lambda: [{k: F(1)} for k in range(n)]
        main = w['A'] + w['cut'] + w['B']
        raw = ident()
        for e in main: al.apply(e, raw)
        dirt = [{} for _ in range(n)]
        for y in range(v, 2 * v): dirt[y] = {k: a for k, a in raw[y].items() if k >= 2 * v}
        for e in w['precondition']: al.apply(e, dirt)
        comp = [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2 * v) for k, a in sorted(dirt[y].items())]
        prefix = w['precondition'] + comp + main
        p = ident(); pinv = ident()
        for e in prefix: al.apply(e, p)
        for e in reversed(prefix): al.apply(al.invert(e), pinv)
        self.conekeys = {}
        for q in lowends:
            U = self.end[q]
            ks = set(pinv[q]) | {k for k in range(n) if k != q and p[k].get(q)}
            self.conekeys[q] = ks
            for k in ks: self.ub[k] = meet(self.ub[k], U)

    # ---- cost
    def chain_cost(self, r, f):
        c = 0.0; last = self.start[r]
        for k in self.pos.get(r, ()):
            U = f(k)
            if U != last:
                c += FC[len(U) - len(last)]; last = U
        if self.end[r] != last: c += FC[len(self.end[r]) - len(last)]
        return c
    def chain_profile(self, r, f):
        p = Counter(); last = self.start[r]
        for k in self.pos.get(r, ()):
            U = f(k)
            if U != last: p[len(U) - len(last)] += 1; last = U
        if self.end[r] != last: p[len(self.end[r]) - len(last)] += 1
        return p
    def total_cost(self):
        g = lambda k: self.f[k]
        return sum(self.chain_cost(r, g) for r in range(self.n))
    def profile(self):
        g = lambda k: self.f[k]; p = Counter()
        for r in range(self.n): p.update(self.chain_profile(r, g))
        return p

    # ---- legality of a single frame at a position
    def ok_at(self, k, U):
        for r in self.regs[k]:
            if not le(U, self.end[r]): return False
            if not le(self.start[r], U): return False
            if self.sec[k] != 'T' and not le(U, self.ub[r]): return False
            if k < self.nA and r in self.cframe and not le(U, self.cframe[r]): return False
        if k in self.frozen and U != self.f[k]: return False
        return True

    def try_set(self, assigns):
        """assigns: {k: U}. Apply with propagation; return (changes dict) or None if illegal."""
        ch = {}
        cur = lambda k: ch.get(k, self.f[k])
        stack = []
        for k, U in assigns.items():
            old = cur(k)
            if le(old, U): stack.append(('up', k, U))
            elif le(U, old): stack.append(('dn', k, U))
            else: stack.append(('up', k, U)); stack.append(('dn', k, U))
        guard = 0
        while stack:
            guard += 1
            if guard > 20000: return None
            kind, k, X = stack.pop()
            old = cur(k)
            new = join(old, X) if kind == 'up' else meet(old, X)
            if kind == 'up' and k in assigns and X == assigns[k] and not le(old, X):
                new = X  # direct assignment (handled by both directions)
            if new == old: continue
            if not self.ok_at(k, new): return None
            ch[k] = new
            for r in self.regs[k]:
                ps = self.pos[r]; i = self.idx[(r, k)]
                if i + 1 < len(ps):
                    k2 = ps[i + 1]
                    if not le(new, cur(k2)): stack.append(('up', k2, new))
                if i > 0:
                    k0 = ps[i - 1]
                    if not le(cur(k0), new): stack.append(('dn', k0, new))
        # final consistency of touched registers
        touched = {r for k in ch for r in self.regs[k]}
        for r in touched:
            last = self.start[r]
            for k in self.pos[r]:
                U = cur(k)
                if not le(last, U): return None
                last = U
            if not le(last, self.end[r]): return None
        return ch

    def delta(self, ch):
        touched = {r for k in ch for r in self.regs[k]}
        new = lambda k: ch.get(k, self.f[k]); old = lambda k: self.f[k]
        return sum(self.chain_cost(r, new) - self.chain_cost(r, old) for r in touched)

    def apply(self, ch):
        for k, U in ch.items(): self.f[k] = U

    def runs(self, r):
        """[(frame, [positions])] distinct consecutive frames of r (excluding start/end pseudo-runs)."""
        out = []
        for k in self.pos.get(r, ()):
            U = self.f[k]
            if out and out[-1][0] == U: out[-1][1].append(k)
            else: out.append((U, [k]))
        return out

    def moves_for(self, r):
        rs = self.runs(r)
        seqf = [self.start[r]] + [U for U, _ in rs] + [self.end[r]]
        mv = []
        for i, (U, ks) in enumerate(rs):
            prevU = seqf[i]; nextU = seqf[i + 2]
            if nextU != U: mv.append({k: nextU for k in ks})       # raise run to next frame
            if prevU != U: mv.append({k: prevU for k in ks})       # lower run to previous frame
        return mv

    def to_word(self):
        import copy
        w = copy.deepcopy(self.w)
        for k, (s, i, e) in enumerate(self.seq):
            U = list(self.f[k])
            if s == 'A': w['A'][i][-1] = U
            elif s == 'B': w['B'][i][-1] = U
            else: w['early_cleanup'][i][1] = U
        p = fw.paid_profile(w)
        w['report']['local_profile'] = {str(d): c for d, c in sorted(p.items())}
        w['report']['local_rank_mass'] = sum(d * c for d, c in p.items())
        return w

def greedy(W, rounds=10, verbose=True):
    for it in range(rounds):
        gained = 0.0; acc = 0
        for r in range(W.n):
            for mv in W.moves_for(r):
                ch = W.try_set(mv)
                if ch is None: continue
                d = W.delta(ch)
                if d < -1e-9:
                    W.apply(ch); gained += d; acc += 1
        if verbose: print('round', it, 'accepted', acc, 'gain', round(gained, 4), 'cost', round(W.total_cost(), 4), flush=True)
        if acc == 0: break

if __name__ == '__main__':
    w = fw.load(sys.argv[1] if len(sys.argv) > 1 else 'K8-FULL-WORD.json.gz')
    W = Word(w)
    c0 = W.total_cost(); print('cost0', c0, dict(sorted(W.profile().items())))
    # count legal moves and their deltas
    hist = Counter(); legal = 0; tried = 0; zero = 0
    for r in range(W.n):
        for mv in W.moves_for(r):
            tried += 1
            ch = W.try_set(mv)
            if ch is None: continue
            legal += 1; d = W.delta(ch); hist[round(d, 3)] += 1
    print('tried', tried, 'legal', legal); print(sorted(hist.items())[:40])
    greedy(W)
    print('final profile', dict(sorted(W.profile().items())))
