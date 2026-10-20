"""r-swaps/swapdoor.py -- exact two-stage replay of s-fix's permutation door on the tiny odd unit, in the style of the hub's
global_k.py (B_K recipe with K = 1: stage 0 forward, stage 1 reverse, stage 2 forward, then the terminal exchange), with
per-port routes.  Standard library only, exact rationals.   usage: python3 -I swapdoor.py

THE UNIT.  h = 4, eight ports: e1..e4 and four weight-3 vectors.  Over Q the weight-3 vectors are signed as the rows of
an antisymmetric conference matrix, so the Q-orthogonality graph equals the F_2 one of f2h4-odd (16 unordered pairs,
32 ordered); their F_2 reductions are exactly the four weight-3 vectors of F_2^4.  All labels are anisotropic.

THE GEOMETRY (the hub's, scaled down: GEOMETRY.json / geometry527.py identities).  Global space U (+) E1 (+) E2,
U = Q^4, E_s = Q^3 (a = h - 1).  Local frames X of a call are idempotents on U: P_t = q q^T/N(q), H = I, 0, K_t = I - P_t
(plus meets and joins used by the variants).  The local frame lifts to the window as
    Lam_0(X) = X (stage 0: window U),   Lam_s(X) = I_{E<s} (+) X (+) I_{E_s}  (s >= 1).
Port route rho_{s,t} (s >= 1): fixes q_t, exchanges K_t and E_s (a GL chart: any linear bijection; the hub: "General GL
charts use gPg^-1"), rho_{0,t} = I.  Global frame of port t's local frame X at stage s: rho Lam_s(X) rho^-1.

THE PER-PORT COVER, at term level.  A register (owner) is at physical frame Lam (its contract at cover slot g is
D_{g Lam g^-1}).  It holds a formal sum of TERMS (input u, instance route r): slot g carries instance g.r of input u,
whose own (class-level) frame is Gamma = r Lam r^-1.  A standard owner of port t at stage s holds terms with
r = rho_{s,t}.  Rules (each is the physical meaning of a hub rule):
  * gate x += c y: legal iff Lam_x == Lam_y (hub: cur[x]==cur[y]==f; then D_{g Lam g^-1} agree at every slot g).
    Terms merge or cancel only when (input, r) agree: the "every input, every slot" reading of the all-column theorem.
  * frame move Lam -> Lam' on a register: paid rank |rank Lam' - rank Lam| for nested frames; acts on all its terms.
  * boundary re-route of owner (bank, port O) from stage s to s+1: every term r -> r M_O, Lam -> M_O^-1 Lam M_O,
    M_O = rho_{s,O}^-1 rho_{s+1,O} (hub: "a data route is right multiplication by rho_{i,S}"; charged address
    movement, rank 0).
  * copy: a scratch register with the same terms and frame; moved, read, then discarded.
The word is right iff, at the end, every owner holds exactly its scalar content with every term on its own port's
route, at the owner's declared frame (then every input carries the kernel D_I: checked as Gamma_end = I - Gamma_start).

VARIANTS of stage 1 (reverse: P-bank target X_S -= 0-bank source z_S), matching e1-e2, e3-e4, w1-w2, w3-w4:
  STD  hub-style per-port star copies (x = 1 on this unit: one rank-1 child per port).
  B1   swap 0-banks of each pair at local 0 at the stage start, deliver through the partner (free gates at local P_S,
       inside K_partner), swap back at the end of the call (local K_S vs K_T).
  B2   as B1, swap back at the start of stage 2 at local 0 (after the boundary re-route).
  B2f  B2 plus a rank-0 re-route putting every term back on its own port's route; measure the frame error left.
  B3   as B1, but the pair descends to the meet K_S cap K_T, swaps back there, and climbs back (repaired door).
  B5   as B1, swap back inside the call at the join P_S + P_T (each owner climbs one), then each content goes home.
  B0   no swap: each source z_S makes its own excursion 0 -> P_S (delivery) -> 0 -> K_S.
  B6   the P-bank door: swap P-banks of each pair at local H at the end of stage 1 (P-bank global frames are
       port-independent), deliver in stage 2 through the partner (free gate at local P_T inside K_S), swap back at H.
  J    no copies in any stage: Y_t climbs to its own line once (0 -> P_t), rides X_t's frames for the whole word (every
       shear a free gate), drops its line at the end (H -> K_t): the cheapest copy-free word; it costs exactly plain.
"""
import sys, json, itertools
from fractions import Fraction as F
sys.dont_write_bytecode = True

h, a = 4, 3
m = h + 2 * a
U = list(range(h)); E = {1: list(range(h, h + a)), 2: list(range(h + a, h + 2 * a))}
LABELS = {'e1': (1, 0, 0, 0), 'e2': (0, 1, 0, 0), 'e3': (0, 0, 1, 0), 'e4': (0, 0, 0, 1),
          'w1': (0, 1, 1, 1), 'w2': (-1, 0, 1, -1), 'w3': (-1, -1, 0, 1), 'w4': (-1, 1, -1, 0)}
PORTS = list(LABELS)
MATCH = [('e1', 'e2'), ('e3', 'e4'), ('w1', 'w2'), ('w3', 'w4')]
PART = {}
for s_, t_ in MATCH: PART[s_] = t_; PART[t_] = s_
# a swap is (x, y) -> (y, -x) (three free gates); after it the partner owner holds SIGN[S] * (content of S)
SIGN = {}
for s_, t_ in MATCH: SIGN[s_] = -1; SIGN[t_] = 1

# ---------------- exact matrix helpers ----------------
def Z(n, k=None): return [[F(0)] * (k or n) for _ in range(n)]
def I(n): return [[F(int(i == j)) for j in range(n)] for i in range(n)]
def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def madd(A, B, s=1): return [[x + s * y for x, y in zip(r, q)] for r, q in zip(A, B)]
def key(A): return tuple(tuple(r) for r in A)
def inv(A):
    n = len(A); M = [r[:] + [F(int(i == j)) for j in range(n)] for i, r in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        iv = 1 / M[c][c]; M[c] = [x * iv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [r[n:] for r in M]
def rank(A):
    M = [r[:] for r in A]; rk = 0; cols = len(M[0]) if M else 0
    for c in range(cols):
        p = next((r for r in range(rk, len(M)) if M[r][c] != 0), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        for r in range(len(M)):
            if r != rk and M[r][c] != 0:
                f = M[r][c] / M[rk][c]; M[r] = [x - f * y for x, y in zip(M[r], M[rk])]
        rk += 1
    return rk
def T(A): return [list(r) for r in zip(*A)]
def image(A): return T(A)                                    # columns span the image
def span_dim(vs): return rank(vs) if vs else 0
def meet_dim(A, B):                                          # dim(im A cap im B)
    ca, cb = image(A), image(B)
    return rank(ca) + rank(cb) - span_dim(ca + cb)
def frame_dist(A, B): return rank(A) + rank(B) - 2 * meet_dim(A, B)
def nested(A, B): return mm(A, B) == A and mm(B, A) == A    # im A inside im B, compatible idempotents
def dot(u, v): return sum(F(x) * F(y) for x, y in zip(u, v))

# ---------------- local frames on U ----------------
def Pq(t):
    q = LABELS[t]; n = dot(q, q); return [[F(q[i]) * F(q[j]) / n for j in range(h)] for i in range(h)]
HU = I(h); ZU = Z(h)
def K(t): return madd(I(h), Pq(t), -1)
def LOC(name):
    """local frame by name: 'P:t', 'H', '0', 'K:t', 'meetK:s:t' (K_s cap K_t), 'joinP:s:t' (P_s + P_t)."""
    p = name.split(':')
    if p[0] == 'P': return Pq(p[1])
    if p[0] == 'H': return HU
    if p[0] == '0': return ZU
    if p[0] == 'K': return K(p[1])
    if p[0] == 'meetK': return madd(madd(I(h), Pq(p[1]), -1), Pq(p[2]), -1)
    if p[0] == 'joinP': return madd(Pq(p[1]), Pq(p[2]))
    raise ValueError(name)
def lift(s, X):
    L = Z(m)
    for i in range(h):
        for j in range(h): L[i][j] = X[i][j]
    for r in range(1, s + 1):
        for i in E[r]: L[i][i] = F(1)
    return L

# ---------------- routes ----------------
def perp_basis(q):
    j = next(k for k in range(h) if q[k]); out = []
    for k in range(h):
        if k == j: continue
        v = [F(0)] * h; v[k] = F(1); v[j] = F(-q[k], q[j]); out.append(v)
    return out
def route(s, t):
    if s == 0: return I(m)
    q = [F(x) for x in LABELS[t]]; b = perp_basis(q); win = U + E[s]
    def emb(u, coords):
        v = [F(0)] * m
        for c, x in zip(coords, u): v[c] = x
        return v
    eps = [emb([F(int(i == j)) for j in range(a)], E[s]) for i in range(a)]
    src = [emb(q, U)] + [emb(x, U) for x in b] + eps          # columns: q, b1..b3, eps1..eps3
    dst = [emb(q, U)] + eps + [emb(x, U) for x in b]          # q -> q, b_i -> eps_i, eps_i -> b_i
    S7 = [[src[c][r] for c in range(len(src))] for r in win]; D7 = [[dst[c][r] for c in range(len(dst))] for r in win]
    R7 = mm(D7, inv(S7)); R = I(m)
    for i, r in enumerate(win):
        for j, c in enumerate(win): R[r][c] = R7[i][j]
    return R
RHO = {(s, t): route(s, t) for s in range(3) for t in PORTS}
RHOI = {k: inv(v) for k, v in RHO.items()}
RKEY = {key(v): k for k, v in RHO.items()}
def glob(s, t, X): return mm(mm(RHO[s, t], lift(s, X)), RHOI[s, t])

# ---------------- geometry checks (the hub's identities, exactly) ----------------
def geometry_report():
    out = []
    q = {t: LABELS[t] for t in PORTS}
    orth = [(s_, t_) for s_ in PORTS for t_ in PORTS if s_ != t_ and dot(q[s_], q[t_]) == 0]
    out.append('labels: %s; norms %s; ordered orthogonal pairs %d (unordered %d); matching %s all orthogonal: %s' % (
        {t: ''.join(str(abs(x) % 2) for x in q[t]) for t in PORTS}, sorted({int(dot(q[t], q[t])) for t in PORTS}),
        len(orth), len(orth) // 2, MATCH, all(dot(q[s_], q[t_]) == 0 for s_, t_ in MATCH)))
    bad = 0; n = 0
    def blockI(blocks, Ublock=None):
        L = Z(m)
        for r in blocks:
            for i in E[r]: L[i][i] = F(1)
        if Ublock is not None:
            for i in range(h):
                for j in range(h): L[i][j] = Ublock[i][j]
        return L
    for t in PORTS:
        for s in (1, 2):
            R = RHO[s, t]; n += 1; bad += mm(R, R) != I(m)
            below = list(range(1, s))
            for X, want in [('P:' + t, blockI(below, HU)), ('H', blockI(below + [s], HU)), ('0', blockI(below, K(t))),
                            ('K:' + t, blockI(below + [s], K(t)))]:
                n += 1; bad += glob(s, t, LOC(X)) != want
        for s in (0, 1):
            n += 2
            bad += glob(s, t, HU) != glob(s + 1, t, Pq(t))
            bad += glob(s, t, K(t)) != glob(s + 1, t, ZU)
    out.append('hub identities rho^2 = I, rho Lam(P|H|0|K) rho = (E<s+U | E<=s+U | E<s+K | E<=s+K), and the boundary '
               'identities H_s = P_(s+1), K_s = 0_(s+1) per port: %d failures of %d' % (bad, n))
    sameP = all(glob(1, s_, Pq(s_)) == glob(1, t_, Pq(t_)) and glob(1, s_, HU) == glob(1, t_, HU) for s_, t_ in MATCH)
    same0 = sum(glob(1, s_, ZU) == glob(1, t_, ZU) for s_, t_ in MATCH)
    out.append('stage 1: P-bank global frames (P and H) port-independent on the matched pairs: %s; 0-bank global start '
               'frames E<1 + K_t equal on matched pairs: %d of %d' % (sameP, same0, len(MATCH)))
    d = [frame_dist(glob(1, s_, K(s_)), glob(1, s_, K(t_))) for s_, t_ in MATCH]
    out.append('own-route distance between the K_S end frame and the partner end frame K_T (S-coordinates), matched '
               'pairs: %s' % d)
    # boundary-transfer condition: a content of port S that crosses the stage boundary s -> s+1 inside port T's owner is
    # re-routed by T's route change M_T = rho(s,T)^-1 rho(s+1,T); it lands on its own route iff rho(s,S) M_T = rho(s+1,S)
    for s in (0, 1):
        ok = sum(mm(RHO[s, S_], mm(RHOI[s, T_], RHO[s + 1, T_])) == RHO[s + 1, S_]
                 for S_ in PORTS for T_ in PORTS if S_ != T_ and dot(LABELS[S_], LABELS[T_]) == 0)
        out.append('boundary-transfer condition rho(%d,S) M_T == rho(%d,S) (content of S crossing the %d->%d boundary in '
                   'T\'s owner lands on its own route), orthogonal ordered pairs: %d of 32' % (s, s + 1, s, s + 1, ok))
    return out

# ---------------- the term machine ----------------
class Fail(Exception): pass

def _gslot():
    import random
    rnd = random.Random(20261013)
    while True:
        g = [[F(rnd.randint(-3, 3)) for _ in range(m)] for _ in range(m)]
        if rank(g) == m: return g, inv(g)
GSLOT, GSLOTI = _gslot()
def bitD(Lam):
    """the hub's bit contract D_P = [[I-P, P], [P, I-P]] at the sample cover slot g, P = g Lam g^-1 (exact)."""
    P = mm(mm(GSLOT, Lam), GSLOTI); Q = madd(I(m), P, -1)
    return [Q[i] + P[i] for i in range(m)] + [P[i] + Q[i] for i in range(m)]

class Machine:
    def __init__(self, label):
        self.label = label; self.reg = {}; self.cost = {}; self.events = []; self.adds = 0; self.issues = []
        self.stage = 0; self.scratch = 0
        self.audit = dict(gates=0, bit_contracts_agree_at_slot_g=0, gates_touching_off_route_terms=0, swap_gates=0,
                          swap_gates_declared_contracts_differ=0, deliveries=0, deliveries_same_route_and_frame=0)
    def pay(self, cat, r):
        if r: self.cost[cat] = self.cost.get(cat, 0) + r; self.cost.setdefault('_profile', {})
    def new(self, name, owner, Lam, terms):
        self.reg[name] = dict(owner=owner, Lam=Lam, terms=dict(terms))
    def move(self, name, X, cat):
        R = self.reg[name]; L2 = lift(self.stage, LOC(X)) if isinstance(X, str) else X
        L1 = R['Lam']
        if nested(L1, L2): r = rank(L2) - rank(L1); kind = 'climb'
        elif nested(L2, L1): r = rank(L1) - rank(L2); kind = 'descent'
        else: raise Fail('non-nested move on %s' % name)
        R['Lam'] = L2
        if r:
            self.cost.setdefault(cat, {}); self.cost[cat][r] = self.cost[cat].get(r, 0) + 1
            if kind == 'descent' and cat != 'copy': self.events.append(('descent', name, r))
    def gate(self, x, y, c, where):
        X, Y = self.reg[x], self.reg[y]
        if X['Lam'] != Y['Lam']:
            raise Fail('ILLEGAL GATE %s += %s*%s (%s): the two registers are at different frames, no common local frame'
                       % (x, c, y, where))
        self.adds += 1; A = self.audit; A['gates'] += 1
        # (i) the physical bit-contract matrices D_(g Lam g^-1) of the two registers at a sample cover slot g
        A['bit_contracts_agree_at_slot_g'] += bitD(X['Lam']) == bitD(Y['Lam'])
        # (ii) terms that are not on their owner's per-port route (declared contract != the content's own contract)
        s_ = self.stage
        own = lambda R: key(RHO[s_, R['owner'][1]])
        if any(r != own(R) for R in (X, Y) for (u, r) in R['terms']): A['gates_touching_off_route_terms'] += 1
        if 'swap' in where:
            A['swap_gates'] += 1
            gx = mm(mm(RHO[s_, X['owner'][1]], X['Lam']), RHOI[s_, X['owner'][1]])
            gy = mm(mm(RHO[s_, Y['owner'][1]], Y['Lam']), RHOI[s_, Y['owner'][1]])
            A['swap_gates_declared_contracts_differ'] += gx != gy
        if 'delivery' in where:
            # same instance and same own frame: every delivered term's route occurs among the target's terms of that port
            A['deliveries'] += 1
            port = where.split('delivery ')[1].split()[0]
            tr = {r for (u, r) in X['terms'] if u.endswith('_' + port)}
            yr = {r for (u, r) in Y['terms'] if u.endswith('_' + port)}
            A['deliveries_same_route_and_frame'] += bool(yr) and yr <= tr
        for k, v in Y['terms'].items():
            X['terms'][k] = X['terms'].get(k, 0) + c * v
            if X['terms'][k] == 0: del X['terms'][k]
        # mixed-route bookkeeping (informational): same input on two different routes in one register
        inputs = {}
        for (u, r) in X['terms']: inputs.setdefault(u, set()).add(r)
        for u, rs in inputs.items():
            if len(rs) > 1: self.issues.append('after %s: register %s holds input %s on %d different routes' % (where, x, u, len(rs)))
    def swap(self, x, y, where, back=False):
        s = -1 if back else 1
        self.gate(x, y, s, where); self.gate(y, x, -s, where); self.gate(x, y, s, where)
    def copy(self, x):
        self.scratch += 1; n = 'c%d' % self.scratch; R = self.reg[x]
        self.reg[n] = dict(owner=R['owner'], Lam=R['Lam'], terms=dict(R['terms'])); return n
    def drop(self, n): del self.reg[n]
    def reroute(self, s):
        for name, R in self.reg.items():
            bank, t = R['owner']; M = mm(RHOI[s, t], RHO[s + 1, t]); Mi = inv(M)
            R['terms'] = {(u, key(mm([list(rr) for rr in r], M))): v for (u, r), v in R['terms'].items()}
            R['Lam'] = mm(mm(Mi, R['Lam']), M)
        self.stage = s + 1

def start_machine(label):
    M = Machine(label)
    for t in PORTS:
        M.new('X_' + t, ('P', t), lift(0, Pq(t)), {('x_' + t, key(RHO[0, t])): F(1)})
        M.new('Y_' + t, ('0', t), lift(0, ZU), {('y_' + t, key(RHO[0, t])): F(1)})
    return M

# ---------------- the calls ----------------
def forward_std(M, s):
    """Y_t += X_t: per-port star copy (X_t copied at P_t, dropped to 0: rank dim P_t = 1), then the climbs."""
    for t in PORTS:
        c = M.copy('X_' + t); M.move(c, '0', 'copy'); M.gate('Y_' + t, c, 1, 'stage %d fwd delivery %s' % (s, t)); M.drop(c)
    for t in PORTS:
        M.move('X_' + t, 'H', 'climb'); M.move('Y_' + t, 'K:' + t, 'climb')

def reverse_std(M, s):
    """X_t -= Y_t (hub reverse lowering: copy the source at K_t, climb the copy to H, signed scatter at H)."""
    for t in PORTS:
        M.move('Y_' + t, 'K:' + t, 'climb'); M.move('X_' + t, 'H', 'climb')
    for t in PORTS:
        c = M.copy('Y_' + t); M.move(c, 'H', 'copy'); M.gate('X_' + t, c, -1, 'stage %d rev delivery %s' % (s, t)); M.drop(c)

def reverse_door(M, s, variant):
    for S_, T_ in MATCH: M.swap('Y_' + S_, 'Y_' + T_, 'stage %d swap at local 0 (%s,%s)' % (s, S_, T_))
    for S_ in PORTS:                                  # owner Y_partner holds z_S; deliver at local P_S (inside K_partner)
        o = 'Y_' + PART[S_]; M.move(o, 'P:' + S_, 'climb')
        M.gate('X_' + S_, o, -SIGN[S_], 'stage %d delivery %s through %s at local P_%s' % (s, S_, o, S_))
    if variant == 'B5':
        for S_, T_ in MATCH:
            M.move('Y_' + S_, 'joinP:%s:%s' % (S_, T_), 'excursion'); M.move('Y_' + T_, 'joinP:%s:%s' % (S_, T_), 'excursion')
            M.swap('Y_' + S_, 'Y_' + T_, 'stage %d swap back at P_S+P_T (%s,%s)' % (s, S_, T_), back=True)
            for A_, B_ in ((S_, T_), (T_, S_)):
                M.move('Y_' + A_, 'P:' + B_, 'excursion')          # drop own line (descent 1)
                M.move('Y_' + A_, 'K:' + A_, 'climb')              # P_B inside K_A: climb 2
        for t in PORTS: M.move('X_' + t, 'H', 'climb')
        return
    for S_ in PORTS: M.move('Y_' + PART[S_], 'K:' + PART[S_], 'climb')
    for t in PORTS: M.move('X_' + t, 'H', 'climb')
    if variant == 'B1':
        for S_, T_ in MATCH: M.swap('Y_' + S_, 'Y_' + T_, 'stage %d swap back at the call end (%s,%s)' % (s, S_, T_), back=True)
    elif variant == 'B3':
        for S_, T_ in MATCH:
            mk = 'meetK:%s:%s' % (S_, T_)
            M.move('Y_' + S_, mk, 'excursion'); M.move('Y_' + T_, mk, 'excursion')
            M.swap('Y_' + S_, 'Y_' + T_, 'stage %d swap back at K_S cap K_T (%s,%s)' % (s, S_, T_), back=True)
            M.move('Y_' + S_, 'K:' + S_, 'excursion'); M.move('Y_' + T_, 'K:' + T_, 'excursion')

def reverse_own_excursion(M, s):
    for t in PORTS:
        M.move('Y_' + t, 'P:' + t, 'excursion'); M.move('X_' + t, 'P:' + t, 'climb')   # X_t already at P_t: no-op
        M.gate('X_' + t, 'Y_' + t, -1, 'stage %d own-excursion delivery %s at local P_%s' % (s, t, t))
        M.move('Y_' + t, '0', 'excursion'); M.move('Y_' + t, 'K:' + t, 'climb'); M.move('X_' + t, 'H', 'climb')

def joint_path_word(M):
    """J: no copies anywhere. Y_t first climbs to its own line (0 -> P_t), then rides X_t's frames for the whole word
    (every shear is a free gate at a common frame), and drops its line at the very end (H -> K_t)."""
    for t in PORTS:
        M.move('Y_' + t, 'P:' + t, 'excursion'); M.gate('Y_' + t, 'X_' + t, 1, 'stage 0 joint delivery %s' % t)
        M.move('X_' + t, 'H', 'climb'); M.move('Y_' + t, 'H', 'climb')
    M.reroute(0)
    for t in PORTS:
        M.gate('X_' + t, 'Y_' + t, -1, 'stage 1 joint delivery %s' % t)
        M.move('X_' + t, 'H', 'climb'); M.move('Y_' + t, 'H', 'climb')
    M.reroute(1)
    for t in PORTS:
        M.gate('Y_' + t, 'X_' + t, 1, 'stage 2 joint delivery %s' % t)
        M.move('X_' + t, 'H', 'climb'); M.move('Y_' + t, 'H', 'climb'); M.move('Y_' + t, 'K:' + t, 'excursion')

def forward_Pdoor(M, s):
    """stage-2 forward with P-banks swapped (owner X_partner holds x_S): deliver Y_S += X_partner at local P_partner."""
    for S_ in PORTS:
        o = 'X_' + PART[S_]; M.move('Y_' + S_, 'P:' + PART[S_], 'climb')          # P_partner inside K_S
        M.gate('Y_' + S_, o, SIGN[S_], 'stage %d P-door delivery %s through %s at local P_%s' % (s, S_, o, PART[S_]))
    for t in PORTS:
        M.move('X_' + t, 'H', 'climb'); M.move('Y_' + t, 'K:' + t, 'climb')

def move_Y_to(M, X):
    for t in PORTS: M.move('Y_' + t, X, 'climb')

# ---------------- expected result and final check ----------------
def expected():
    """K = 1 recipe on one bank pair: fwd (x, y+x); rev (x-(y+x), y+x) = (-y, x+y); fwd (-y, x).  Then exchange."""
    return {('P', t): {'y_' + t: F(-1)} for t in PORTS} | {('0', t): {'x_' + t: F(1)} for t in PORTS}

def final_check(M):
    exp = expected(); rep = dict(scalar_ok=True, routes_ok=True, frames_ok=True, kernel_ok=True, off_route_terms=0,
                                 junk_terms=0, details=[])
    for name, R in M.reg.items():
        bank, t = R['owner']
        sc = {}
        for (u, r), v in R['terms'].items(): sc[u] = sc.get(u, 0) + v
        sc = {u: v for u, v in sc.items() if v}
        if sc != exp[bank, t]: rep['scalar_ok'] = False; rep['details'].append('%s scalar %s, want %s' % (name, sc, exp[bank, t]))
        want_r = key(RHO[2, t]); declared = lift(2, HU if bank == 'P' else K(t))
        if R['Lam'] != declared:
            rep['frames_ok'] = False; rep['details'].append('%s not at its declared end frame' % name)
        for (u, r), v in R['terms'].items():
            if r != want_r:
                rep['off_route_terms'] += 1; rep['routes_ok'] = False
            if u not in exp[bank, t]: rep['junk_terms'] += 1
        # kernel: for the surviving own-route term of the expected input, Gamma_end = I - Gamma_start
        for u in exp[bank, t]:
            v = R['terms'].get((u, want_r))
            if v != exp[bank, t][u]: rep['kernel_ok'] = False; continue
            src_t = u.split('_')[1]; start = lift(0, Pq(src_t) if u.startswith('x') else ZU)
            Gend = mm(mm(RHO[2, t], R['Lam']), RHOI[2, t])
            if madd(Gend, start) != I(m) or mm(Gend, start) != Z(m): rep['kernel_ok'] = False
    return rep

def sample_terms(M):
    for t in PORTS:
        for b in ('Y_', 'X_'):
            R = M.reg[b + t]; want = key(RHO[2, t]); off = [(u, r, v) for (u, r), v in R['terms'].items() if r != want]
            if off:
                u, r, v = off[0]; rr = [list(x) for x in r]
                name = RKEY.get(r)
                cands = {'rho(1,%s).M(%s)' % (p1, p2): key(mm(RHO[1, p1], mm(RHOI[1, p2], RHO[2, p2]))) for p1 in PORTS for p2 in PORTS}
                lab = [k for k, v_ in cands.items() if v_ == r]
                G = mm(mm(rr, R['Lam']), inv(rr)); Gwant = mm(mm(RHO[2, t], R['Lam']), RHOI[2, t])
                return ('example at the end: owner %s%s holds %s*%s on route %s (its own route is rho(2,%s)); that term\'s '
                        'own global frame is %d rank units from the owner\'s declared one; owner holds %d terms in all'
                        % (b, t, v, u, lab or name or 'other', t, frame_dist(G, Gwant), len(R['terms'])))
    return ''

def run(variant):
    M = start_machine(variant); log = {}
    if variant == 'J':
        joint_path_word(M); res = final_check(M); res['stage1_cost'] = {}; res['stage1_percall_contract'] = True
        res['total_cost'] = M.cost; res['additions'] = M.adds; res['issues'] = M.issues[:6]; res['issue_count'] = len(M.issues)
        res['descents'] = len(M.events); res['audit'] = M.audit; res['sample'] = sample_terms(M)
        return res
    try:
        forward_std(M, 0); M.reroute(0)
        if variant == 'STD' or variant == 'B6': reverse_std(M, 1)
        elif variant == 'B0': reverse_own_excursion(M, 1)
        else: reverse_door(M, 1, 'B2' if variant == 'B2f' else variant)
        if variant == 'B6':
            for S_, T_ in MATCH: M.swap('X_' + S_, 'X_' + T_, 'stage 1 P-bank swap at local H (%s,%s)' % (S_, T_))
        st1 = {k: dict(v) for k, v in M.cost.items()}
        # per-call contract of stage 1 (global_k.py line 152: every column as before except target += source)
        pc = all({u for (u, r) in M.reg['Y_' + t]['terms']} == {'x_' + t, 'y_' + t} for t in PORTS)
        M.reroute(1)
        if variant in ('B2', 'B2f'):
            for S_, T_ in MATCH: M.swap('Y_' + S_, 'Y_' + T_, 'stage 2 swap back at local 0 (%s,%s)' % (S_, T_), back=True)
        if variant == 'B2f':
            errs = []
            for t in PORTS:
                R = M.reg['Y_' + t]; want = RHO[2, t]
                bad = {r for (u, r) in R['terms'] if r != key(want)}
                if bad:
                    r = [list(x) for x in next(iter(bad))]
                    G = mm(mm(r, R['Lam']), inv(r))                 # the content's own global frame (all its terms share r)
                    Lfix = mm(mm(RHOI[2, t], G), RHO[2, t])         # its physical frame after a rank-0 re-route to rho_(2,t)
                    errs.append(frame_dist(Lfix, lift(2, ZU)))
            log['B2f_frame_error_after_route_fix'] = errs
            raise Fail('B2f: after the rank-0 route fix every swapped 0-bank content sits %s rank units from its own '
                       'stage-2 start frame (local 0); reaching it costs that many extra paid moves' % sorted(set(errs)))
        if variant == 'B6':
            forward_Pdoor(M, 2)
            for S_, T_ in MATCH: M.swap('X_' + S_, 'X_' + T_, 'stage 2 P-bank swap back at local H (%s,%s)' % (S_, T_), back=True)
        else:
            forward_std(M, 2)
        res = final_check(M); res['stage1_cost'] = st1; res['stage1_percall_contract'] = pc
        res['total_cost'] = M.cost; res['additions'] = M.adds; res['issues'] = M.issues[:6] + (['...'] if len(M.issues) > 6 else [])
        res['issue_count'] = len(M.issues); res['descents'] = len(M.events); res['audit'] = M.audit
        res['sample'] = sample_terms(M)
        return res
    except Fail as e:
        return dict(failed=str(e), stage=M.stage, cost_so_far=M.cost, additions=M.adds, audit=M.audit, **log)

def mass(c): return sum(r * n for cat, d in c.items() if isinstance(d, dict) and not cat.startswith('_') for r, n in d.items())

def main():
    print('r-swaps / swapdoor.py: exact two-stage replay of the permutation door on the tiny odd unit (h = 4, v = 8)')
    for line in geometry_report(): print('  ' + line)
    print()
    rows = []
    for v in ['STD', 'B1', 'B2', 'B2f', 'B3', 'B5', 'B0', 'B6', 'J']:
        r = run(v); rows.append((v, r))
        print('== %s' % v)
        if 'failed' in r:
            print('   FAILS: %s' % r['failed'])
            print('   (failed during stage %d; paid moves until then %s; additions %d)' % (r['stage'], json.dumps(r['cost_so_far']), r['additions']))
            print('   gate audit up to the failure: %s' % json.dumps(r['audit']))
            continue
        c1 = r['stage1_cost']; c1m = {k: v_ for k, v_ in c1.items() if not k.startswith('_')}
        print('   whole-word scalar identity (what global_k.py\'s F2 replay sees): %s' % r['scalar_ok'])
        print('   every term on its own port route at the end (per-port cover): %s  (off-route terms %d, junk terms %d)'
              % (r['routes_ok'], r['off_route_terms'], r['junk_terms']))
        print('   every owner at its declared end frame: %s; kernel D_I on every input (Gamma_end = I - Gamma_start): %s'
              % (r['frames_ok'], r['kernel_ok']))
        print('   stage-1 per-call contract (each 0-bank still holds its own port\'s columns): %s' % r['stage1_percall_contract'])
        print('   paid moves through stage 1 (stages 0+1), by category {rank: count}: %s; rank mass %d' % (json.dumps(c1m), mass(c1m)))
        print('   additions in the whole word: %d; descents (non-copy) %d' % (r['additions'], r['descents']))
        if r['issue_count']: print('   mixed-route warnings: %d, e.g. %s' % (r['issue_count'], r['issues'][0]))
        print('   gate audit: %s' % json.dumps(r['audit']))
        if r['sample']: print('   %s' % r['sample'])
        for d in r.get('details', [])[:4]: print('     ' + d)
    # stage-1 increments relative to STD
    print()
    base = None
    for v, r in rows:
        if 'failed' in r: continue
        c1m = {k: v_ for k, v_ in r['stage1_cost'].items() if not k.startswith('_')}
        if v == 'STD': base = mass(c1m)
    for v, r in rows:
        if 'failed' in r: print('  %-4s stage-1 rank mass: n/a (fails: %s)' % (v, r['failed'].split(':')[0])); continue
        c1m = {k: v_ for k, v_ in r['stage1_cost'].items() if not k.startswith('_')}
        tw = {k: v_ for k, v_ in r['total_cost'].items() if not k.startswith('_')}
        print('  %-4s stage 0+1 rank mass %s (STD %d, delta %s); whole word %d, plain %d (2v arrays x rank %d), deficit %+d; '
              'passes all checks: %s' % (v, mass(c1m) if c1m else 'n/a', base, ('%+d' % (mass(c1m) - base)) if c1m else 'n/a',
              mass(tw), 2 * len(PORTS) * m, m, 2 * len(PORTS) * m - mass(tw),
              r['scalar_ok'] and r['routes_ok'] and r['frames_ok'] and r['kernel_ok']))

if __name__ == '__main__':
    main()
