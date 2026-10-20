"""r-riders / replay.py -- exact two-stage scalar-plus-frame replay of the toy (formal sums, exact Z4 frames).

Every register is a formal sum {(input, physical frame q): Fraction}.  Letters:
  move(r, z)        paid, 1 move: every term q += chi_z                                   (both rules)
  shift(r, w)       free: every term q += 2 chi_w                                          (both rules)
  add(d, s, c)      free: d += c * s.   WIDE: no condition.  EQUAL (hub's local rule, global_k.py op 1):
                    both registers in one-local-frame form and at the SAME local frame (an empty d takes s's frame)
  copy(C, s)        free: open the private centre C as a copy of s (C must be empty); the hub then pays the
                    moves made on C (its op 2: the paid copy is 'copy + moves').  The source and the moves are recorded.
  boundary()        EQUAL2 only: the hub's call boundary (new per-port routes, see Machine.boundary)
  erase(C)          free: WIDE: zero C (engine localM).  EQUAL (hub op 3): allowed only if C == (recorded moves and
                    shifts applied to the recorded source's CURRENT content), i.e. the copy is intact.
Local frame of a term (input l, q) = q + s_l (the toy's per-port route; s_x_T = chi_a, s_y_T = 0, s-wide).
Rules: WIDE = engine game; EQUAL = hub's local rule with ONE identification for both stages (s-wide's toy);
EQUAL2 = hub's local rule with two calls, per-port routes changing at the boundary as in global_k.py / GEOMETRY.
End check: every target register equals its target exactly; every other register (centre) is zero.
usage: python3 -B replay.py   (runs every program below under both rules and prints a table)"""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import toy
chi, add_, F, PSI, E = toy.chi, toy.add, toy.F, toy.PSI, toy.E

class Illegal(Exception): pass

class Machine:
    def __init__(self, rule, regs_live):
        self.rule = rule; self.R = {}; self.paid = 0; self.log = []; self.copyrec = {}
        self.start = dict(toy.START); self.call = 'i'
        for name, inp in regs_live.items(): self.R[name] = {(inp, toy.Z0): Fr(1)}
    def get(self, r): return self.R.setdefault(r, {})
    def lframes(self, r):
        return {add_(q, self.start[inp]) for (inp, q), c in self.get(r).items() if c}
    def lframe(self, r):
        fs = self.lframes(r)
        if len(fs) > 1: raise Illegal('%s holds terms at %d local frames' % (r, len(fs)))
        return next(iter(fs)) if fs else None
    def move(self, r, z):
        self.R[r] = {(inp, add_(q, chi(z))): c for (inp, q), c in self.get(r).items()}
        self.paid += 1; self.log.append(('move', r, z))
        if r in self.copyrec: self.copyrec[r][1].append(('m', z))
    def shift(self, r, w):
        self.R[r] = {(inp, add_(q, tuple(2 * v for v in chi(w)))): c for (inp, q), c in self.get(r).items()}
        self.log.append(('shift', r, w))
        if r in self.copyrec: self.copyrec[r][1].append(('s', w))
    def add(self, d, s, c=1):
        if self.rule in ('EQUAL', 'EQUAL2'):
            fd, fs = self.lframe(d), self.lframe(s)
            if fd is not None and fs is not None and fd != fs:
                raise Illegal('add %s += %s: local frames differ' % (d, s))
        D = self.get(d)
        for k, v in self.get(s).items():
            D[k] = D.get(k, 0) + c * v
            if D[k] == 0: del D[k]
        self.log.append(('add', d, s, c))
        if self.rule in ('EQUAL', 'EQUAL2'): self.lframe(d)
    def setlin(self, d, combo):
        """WIDE only (engine localM): d := sum c * register (overwrites a live register; the hub's local format has no
        such step: its op 1 needs two different registers at one frame)."""
        if self.rule != 'WIDE': raise Illegal('setlin %s: overwrite of a register is not a hub op' % d)
        new = {}
        for r, c in combo:
            for k, v in self.get(r).items():
                new[k] = new.get(k, 0) + c * v
        self.R[d] = {k: v for k, v in new.items() if v}; self.log.append(('setlin', d, combo))
    def copy(self, C, s):
        if self.get(C): raise Illegal('centre %s not empty' % C)
        self.R[C] = dict(self.get(s)); self.copyrec[C] = (s, []); self.log.append(('copy', C, s))
    def erase(self, C):
        if self.rule in ('EQUAL', 'EQUAL2'):
            src, word = self.copyrec[C]
            img = dict(self.get(src))
            for kind, z in word:
                g = chi(z) if kind == 'm' else tuple(2 * v for v in chi(z))
                img = {(inp, add_(q, g)): c for (inp, q), c in img.items()}
            if img != self.get(C): raise Illegal('erase %s: copy no longer equals its source image' % C)
        self.R[C] = {}; self.copyrec.pop(C, None); self.log.append(('erase', C))
    def boundary(self):
        """EQUAL2 only: the call boundary i -> i+1.  New per-port routes from the hub's boundary identities
        H_i = P_(i+1), K_i = 0_(i+1) (d-g-crossstage (2), s-fix recount): X_T's stage-i end frame F becomes its new
        start frame chi_a, Y_T's end frame psi_T becomes the new frame 0.  So s'_x_T = 2 chi_a - F, s'_y_T = -psi_T
        (port dependent).  Every register must be in one-frame form under the new routes."""
        if self.rule != 'EQUAL2': return
        for a in E:
            self.start[('x', a)] = toy.sub(toy.smul(2, chi(a)), F)
            self.start[('y', a)] = toy.sub(toy.Z0, PSI[a])
        self.call = 'i+1'
        for r in self.R:
            try: self.lframe(r)
            except Illegal as e: raise Illegal('at the call boundary: ' + str(e))
    def check(self, T):
        bad = []
        for r, want in T.items():
            have = {k: v for k, v in self.get(r).items() if v}
            if have != {k: Fr(v) for k, v in want.items()}: bad.append(r)
        for r, cont in self.R.items():
            if r not in T and any(cont.values()): bad.append(r + '(not clean)')
        return bad

LIVE = {'X%d' % a: ('x', a) for a in E}
LIVE.update({'Y%d' % a: ('y', a) for a in E})

def climbs(M, riders=False):
    """stage i data climbs (12 paid moves).  riders=True: v-s-wide's three riders, put on at the frames he found:
    y_4 into X_1 at X_1's start (from Y_4 after its first move), y_1 into X_2 at its start (from Y_1 after its first
    move), y_2 into X_4 at X_4's mid point (from Y_2 at its end, before Y_2's scatter)."""
    P = toy.PORTS
    if not riders:
        for a in E:
            for z in P[a]: M.move('X%d' % a, z)
            for z in P[a]: M.move('Y%d' % a, z)
        return
    M.move('Y4', P[4][0]); M.move('Y1', P[1][0])
    M.add('X1', 'Y4'); M.add('X2', 'Y1')                       # riders 1 and 2
    for z in P[1]: M.move('X1', z)
    for z in P[2]: M.move('X2', z)
    M.move('Y4', P[4][1]); M.move('Y1', P[1][1])
    for z in P[2]: M.move('Y2', z)
    M.move('X4', P[4][0]); M.add('X4', 'Y2'); M.move('X4', P[4][1])   # rider 3 at X_4's mid point

def scatter(M):
    """stage i scatter Y_T += x_T @ psi_T through the centre: copy X_T, dir(a_T), shift(a_T) (1 paid move each)."""
    for a in E:
        M.copy('C', 'X%d' % a); M.move('C', a); M.shift('C', a); M.add('Y%d' % a, 'C'); M.erase('C')

def total_from_Y(M):
    """coded total Z = sum_T (y_T + x_T)@F through the centre: copy Y_T (after scatter), dir(a_T) (1 paid move each)."""
    for a in E:
        M.copy('C', 'Y%d' % a); M.move('C', a); M.add('Z', 'C'); M.erase('C')

def prog_A_plain(M):
    climbs(M); scatter(M); M.boundary()
    for a in E:   # stage i+1 (reverse): X_S -= Y_S through the centre
        M.copy('C', 'Y%d' % a); M.move('C', a); M.add('X%d' % a, 'C', -1); M.erase('C')

def prog_B_plain(M):
    climbs(M); scatter(M); total_from_Y(M)

RID = {1: 4, 2: 1, 4: 2}     # X_T carries y_RID[T]
def prog_B_riders(M):
    climbs(M, riders=True)
    for a in E: M.add('Z', 'X%d' % a)                       # the coded total, 0 paid moves (v-s-wide)
    for a in E:                                             # cleanup: y_S @ F through the centre, 1 paid move each
        s = RID[a]
        M.copy('C', 'Y%d' % s); M.move('C', s); M.add('X%d' % a, 'C', -1); M.erase('C')
    scatter(M)

def prog_B_riders_noclean(M):
    climbs(M, riders=True)
    for a in E: M.add('Z', 'X%d' % a)
    scatter(M)

def prog_A_riders(M):
    """the riders in the full two-stage program A: riders, cleanup, scatter, stage i+1."""
    climbs(M, riders=True)
    for a in E:
        s = RID[a]; M.copy('C', 'Y%d' % s); M.move('C', s); M.add('X%d' % a, 'C', -1); M.erase('C')
    scatter(M); M.boundary()
    for a in E:
        M.copy('C', 'Y%d' % a); M.move('C', a); M.add('X%d' % a, 'C', -1); M.erase('C')

def prog_A_riders_late(M):
    """riders kept across the call boundary; their cleanup is done in stage i+1 (scatter of the clean x_T
    is taken from a fresh copy of X_T after its cleanup, so the cleanup comes first in stage i+1 order)."""
    climbs(M, riders=True)
    M.boundary()
    for a in E:
        s = RID[a]; M.copy('C', 'Y%d' % s); M.move('C', s); M.add('X%d' % a, 'C', -1); M.erase('C')
    scatter(M)
    for a in E:
        M.copy('C', 'Y%d' % a); M.move('C', a); M.add('X%d' % a, 'C', -1); M.erase('C')

def prog_A_riders_wide(M):
    """riders; the three cleanup moves y_S@F are KEPT (registers Q_S) and reused as the stage-(i+1) moves:
    X_S := -Q_S needs the engine's free overwrite (WIDE only)."""
    climbs(M, riders=True)
    for a in E:
        s = RID[a]; Q = 'Q%d' % s
        M.copy(Q, 'Y%d' % s); M.move(Q, s); M.add('X%d' % a, Q, -1)
    scatter(M); M.boundary()
    for a in E:
        M.setlin('X%d' % a, [('Q%d' % a, -1)]); M.erase('Q%d' % a)

def prog_B_total_across_boundary(M):
    """plain B, then the coded total Z (formed at local F of call i) is carried into call i+1."""
    prog_B_plain(M); M.boundary()

PROGS = [('A', 'plain (no riders)', prog_A_plain), ('A', 'riders + cleanup', prog_A_riders),
         ('B', 'plain (no riders)', prog_B_plain), ('B', 'riders + cleanup', prog_B_riders),
         ('A', 'riders, cleanup after boundary', prog_A_riders_late),
         ('A', 'riders, cleanup reused (overwrite)', prog_A_riders_wide),
         ('B', 'plain, Z carried into call i+1', prog_B_total_across_boundary),
         ('B', 'riders, NO cleanup (what-if)', prog_B_riders_noclean)]

if __name__ == '__main__':
    for v, name, fn in PROGS:
        T = toy.targets(v)
        for rule in ('EQUAL', 'EQUAL2', 'WIDE'):
            M = Machine(rule, LIVE)
            try:
                fn(M); bad = M.check(T); res = 'END STATE EXACT' if not bad else 'WRONG: ' + ', '.join(bad)
            except Illegal as e:
                res = 'ILLEGAL: %s' % e
            print('%s | %-34s | %-6s | paid moves %2d | %s' % (v, name, rule, M.paid, res))
