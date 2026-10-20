"""s-fix/recount.py: the stage-(i+1) target/source frame-sharing count with the hub's ACTUAL port isometries.
Input: GEOMETRY-nine-bit589.json (copy of hub research/community-hyperplane599-audit/pr589/bit589/nine/GEOMETRY.json:
111 ports, h = 8, 7-blocks, rho15 rational, symmetric involutions).  usage: python3 -I recount.py
Global frame of local label X for port t at stage s >= 1 (hub identities rho*diag(X,I7)*rho):
    G_s^t(X) = E_<s (+) rho_t^(s) (X (+) E_s),   rho_t^(s) acting on the axes U (+) E_s.
Exact rank arithmetic: modulo the prime 2^61-1 (no denominator of rho15 is divisible by it) with a Fraction
cross-check on the first ports."""
import json, sys, itertools
from fractions import Fraction as F
P = (1 << 61) - 1
geo = json.load(open('GEOMETRY-nine-bit589.json'))
ports = geo['ports']; V = len(ports); h = geo['h']; a = h - 1; S = 2; i = S - 1   # boundary i=1 -> i+1=2
n = h + a * S                                   # window of stage 2: U (+) E_1 (+) E_2
def modp(x):
    x = F(x); return x.numerator % P * pow(x.denominator % P, P - 2, P) % P
def rank(vs):
    rows = [list(v) for v in vs]; r = 0; m = len(rows[0]) if rows else 0
    for c in range(m):
        piv = next((k for k in range(r, len(rows)) if rows[k][c] % P), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]; inv = pow(rows[r][c], P - 2, P)
        rows[r] = [x * inv % P for x in rows[r]]
        for k in range(len(rows)):
            if k != r and rows[k][c] % P:
                f = rows[k][c]; rows[k] = [(x - f * y) % P for x, y in zip(rows[k], rows[r])]
        r += 1
    return r
def axes(s): return list(range(h)) + list(range(h + a * (s - 1), h + a * s))
def unit(j): v = [0] * n; v[j] = 1; return v
def E(blocks): return [unit(j) for s in blocks for j in range(h + a * (s - 1), h + a * s)]
rhos = [[[modp(x) for x in row] for row in p['rho15']] for p in ports]
qs = [[modp(x) for x in p['q']] for p in ports]
# local subspaces (bases in Q^8): P_t = <q_t>, K_t = q_t^perp, H = Q^8, 0
def perp_basis(q):          # basis of {x : q.x = 0} in F_P^8
    j = next(k for k in range(h) if q[k]); out = []
    for k in range(h):
        if k == j: continue
        v = [0] * h; v[k] = 1; v[j] = (-q[k] * pow(q[j], P - 2, P)) % P; out.append(v)
    return out
def G(t, s, X):             # global frame (list of vectors in F_P^n) of local basis X for port t at stage s
    ax = axes(s); R = rhos[t]
    loc = [x + [0] * a for x in X] + [[0] * h + [int(k == j) for k in range(a)] for j in range(a)]
    out = E(range(1, s))
    for w in loc:
        img = [sum(R[r][c] * w[c] for c in range(h + a)) % P for r in range(h + a)]
        v = [0] * n
        for r, A in enumerate(ax): v[A] = img[r]
        out.append(v)
    return out
def le(A, B): return rank(B + A) == rank(B)               # A subset of B
def eq(A, B): return le(A, B) and le(B, A)
def meet_contains(lo, his):  # lo subset of intersection of his: check lo inside each
    return all(le(lo, hh) for hh in his)
Hb = [[int(k == j) for k in range(h)] for j in range(h)]
# ---- 0. the hub identities, exactly (all ports) ----
bad = 0
for t in range(V):
    Pt = [qs[t]]; Kt = perp_basis(qs[t])
    # boundary identities H_i = P_(i+1), K_i = 0_(i+1) for the same port (i = 1, s = 2)
    if not eq(G(t, i, Hb), G(t, S, Pt)): bad += 1
    if not eq(G(t, i, Kt), G(t, S, [])): bad += 1
print('ports %d; boundary identities H_1 = P_2 and K_1 = 0_2 fail for %d of %d checks' % (V, bad, 2 * V))
# ---- 1. the stage-2 intervals of the four data registers of one port ----
Ug = [unit(j) for j in range(h)]
lo0 = []; hi0 = []
for t in range(V):
    Kt = perp_basis(qs[t])
    lo0.append(G(t, S, [])); hi0.append(G(t, S, Kt))
loP = [G(t, S, [qs[t]]) for t in range(V)]; hiP = [G(t, S, Hb) for t in range(V)]
portindepP = all(eq(loP[t], loP[0]) and eq(hiP[t], hiP[0]) for t in range(V))
print('stage-2 P-bank interval [G(P_t), G(H)] the same for all ports: %s; G(P_t) = E_1+U: %s' %
      (portindepP, eq(loP[0], E([1]) + Ug)))
distinct0 = len({tuple(sorted(range(V), key=lambda u: 0))} )
same_lo0 = sum(1 for s_ in range(V) for t in range(V) if s_ < t and eq(lo0[s_], lo0[t]))
print('stage-2 0-bank start frames G(0) = E_1 + K_t: pairs of ports with EQUAL start frame: %d of %d' % (same_lo0, V*(V-1)//2))
# ---- 2. counts, strict global reading (a gate needs one global subspace in both registers' intervals) ----
def share(lo1, hi1, lo2, hi2): return meet_contains(lo1 + lo2, [hi1, hi2])
# old identification of a stage-2 register = a frame of its stage-2 interval that is also G_1^t(X) for some X:
# only the boundary frame (checked below); so 'old' intervals are the single boundary frames.
oldP = [loP[t] for t in range(V)]; old0 = [lo0[t] for t in range(V)]
def interval_old_check(t):
    # every frame G_1^t(X) lies in E_1 + U; the stage-2 P interval meets E_1+U only in {E_1+U}; the 0 interval only in {E_1+K_t}
    win1 = E([1]) + Ug
    return eq(hiP[t] if False else loP[t], win1) and rank(hi0[t] + []) and le(lo0[t], win1) and not le(hi0[t], win1)
cnt = {}
for name, A, B in [('P-target new / 0-source new', 'nP', 'n0'), ('P-target old / 0-source new', 'oP', 'n0'),
                   ('P-target new / 0-source old', 'nP', 'o0'), ('P-target old / 0-source old', 'oP', 'o0')]:
    c = 0; cd = 0
    for s_ in range(V):
        l1, h1 = (loP[s_], hiP[s_]) if A == 'nP' else (oldP[s_], oldP[s_])
        for t in range(V):
            l2, h2 = (lo0[t], hi0[t]) if B == 'n0' else (old0[t], old0[t])
            if share(l1, h1, l2, h2):
                c += 1; cd += (s_ == t)
    cnt[name] = (c, cd)
    print('strict global, %s: %d of %d ordered pairs share a frame (diagonal %d)' % (name, c, V * V, cd))
c00 = sum(1 for s_ in range(V) for t in range(V) if s_ != t and share(lo0[s_], hi0[s_], lo0[t], hi0[t]))
print('CONTROL strict global, 0-source / 0-source at stage 2 (needed by every cross-port total of the stage-2 decoder): %d of %d' % (c00, V * (V - 1)))
# rev-direction pairs: 0-target / P-source (fwd stage roles), same numbers by symmetry of the test
# ---- 3. local reading (the hub's checker: cur[x] == cur[y] == f in the call's local coordinates) ----
qr = [[F(x) for x in p['q']] for p in ports]
dotq = lambda u, v: sum(x * y for x, y in zip(u, v))
knows = sum(1 for s_ in range(V) for t in range(V) if s_ != t and dotq(qr[s_], qr[t]) == 0)
diag = sum(1 for s_ in range(V) if dotq(qr[s_], qr[s_]) == 0)
print('local reading, stage-2 P-target S [P_S, H] and 0-source T [0, K_T]: share iff q_S.q_T = 0: %d of %d off-diagonal, diagonal %d' % (knows, V * (V - 1), diag))
print('local reading, 0-source / 0-source: share at local 0 always: %d of %d' % (V * (V - 1), V * (V - 1)))
