#!/usr/bin/env python3
"""lagframes.py (t-family) -- the enlarged knows relation and the copy cost on actual Lagrangian frames (h = 3, 4).
A twisted port is the hub's port for an odd label u carried by a symplectic map g (a fixed twist, used by every
total): source register frames X' = g I(L_<u>, L_F), target register frames Y' = g I(L_0, L_{u^perp}),
p' = g (u|u), q' = g (0|u).  Exact knows: S knows T iff X'_T and Y'_S share a frame.  Exact cost of a total that
reads R at one frame and delivers to D at one frame: min d(r0, r*) over r0 in the intersection of the X'_T (T in R)
and r* in the intersection of the Y'_S (S in D) (infinite if either is empty).  Checked against the vector model of
lagfloor.py: knows => omega(p'_T, q'_S) = 0, and cost >= rank omega(span p'_R, span q'_D).
Then: floors of random twisted 4-port families with the EXACT costs (lagmitm, payload F_3).
usage: python3 -I lagframes.py h nfam [v=4] [rand|stabF|stab0]   (uses ../lib/frames.py)"""
import os, sys, random, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
TG = os.path.normpath(os.path.join(HERE, '../lib'))
sys.path.insert(0, TG)
from frames import lag_graph, L_of, span, perp, popc, omega
sys.path.insert(0, HERE)
from lagfloor import rank2
import lagmitm

h = int(sys.argv[1]) if len(sys.argv) > 1 else 3
nfam = int(sys.argv[2]) if len(sys.argv) > 2 else 30
rng = random.Random(2026)
t0 = time.time()
Gr = lag_graph(h); n = Gr['n']; frames = Gr['frames']; idx = Gr['idx']; i0, iF = Gr['i0'], Gr['iF']
masks = [sum(1 << v for v in L) for L in frames]
D = [[h - (popc(masks[i] & masks[j]).bit_length() - 1) for j in range(n)] for i in range(n)]
def interval(a, b):
    return frozenset(z for z in range(n) if D[a][z] + D[z][b] == D[a][b])
odd = [u for u in range(1, 1 << h) if popc(u) & 1]
def std_port(u):
    a = idx[L_of(span([u]), h)]; e = idx[L_of(perp(span([u]), h), h)]
    return interval(a, iF), interval(i0, e), u | (u << h), u << h
LFv = [z | (z << h) for z in range(1, 1 << h)]; L0v = [z << h for z in range(1, 1 << h)]
def rand_sp(k=40):
    """MODE rand: any symplectic map; stabF: transvections along L_F (fix L_F pointwise, so every source interval
    still ends at L_F and any read set has a common frame); stab0: along L_0 (targets keep L_0)."""
    pool_ = {'rand': None, 'stabF': LFv, 'stab0': L0v}[MODE]
    ws = [rng.randrange(1, 1 << (2 * h)) if pool_ is None else rng.choice(pool_) for _ in range(k)]
    def g(v):
        for w in ws:
            if omega(v, w, h): v ^= w
        return v
    return g
def twist(port, g):
    X, Y, p, q = port
    fm = [idx[frozenset(g(v) for v in frames[f])] for f in range(n)]
    return frozenset(fm[f] for f in X), frozenset(fm[f] for f in Y), g(p), g(q)
def knows_exact(P, S, T):
    return S != T and bool(P[T][0] & P[S][1])
def cost_exact(P, R, Dd):
    A = frozenset(range(n)); B = frozenset(range(n))
    for T in R: A &= P[T][0]
    for S in Dd: B &= P[S][1]
    if not A or not B: return None
    return min(D[a][b] for a in A for b in B)
def pair_vec(P, R, Dd):
    return rank2([sum(omega(P[T][2], P[S][3], h) << i for i, S in enumerate(Dd)) for T in R])
out = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)
say('h=%d: %d frames, graph built %.1fs' % (h, n, time.time() - t0))
# control: the hub's ports (identity twist)
P0 = [std_port(u) for u in odd]
ok = all((knows_exact(P0, S, T)) == (popc(odd[S] & odd[T]) % 2 == 0) for S in range(len(odd)) for T in range(len(odd)) if S != T)
subsets = [c for k in range(1, len(odd) + 1) for c in itertools.combinations(range(len(odd)), k)]
okc = all(cost_exact(P0, R, Dd) == pair_vec(P0, R, Dd) for R in subsets for Dd in subsets)
say('control, hub ports (all %d odd labels, no twist): knows == (u_S.u_T = 0): %s; exact cost == pairing rank on all %d (R,D): %s'
    % (len(odd), ok, len(subsets) ** 2, okc))
# random twisted families of 4 ports
v = int(sys.argv[3]) if len(sys.argv) > 3 else 4
MODE = sys.argv[4] if len(sys.argv) > 4 else 'rand'
st = dict(pairs=0, knows_exact=0, knows_vec=0, mismatch_dir=0, cost_cmp=0, cost_eq=0, cost_lt=0, cost_inf=0, floors={})
old_knows = None
for fam in range(nfam):
    P = [twist(std_port(rng.choice(odd)), rand_sp()) for _ in range(v)]
    for S in range(v):
        for T in range(v):
            if S == T: continue
            st['pairs'] += 1
            ke = knows_exact(P, S, T); kv = omega(P[T][2], P[S][3], h) == 0
            st['knows_exact'] += ke; st['knows_vec'] += kv
            if ke and not kv: st['mismatch_dir'] += 1
    subs = [c for k in range(1, v + 1) for c in itertools.combinations(range(v), k)]
    for R in subs:
        for Dd in subs:
            ce = cost_exact(P, R, Dd); pv = pair_vec(P, R, Dd)
            st['cost_cmp'] += 1
            if ce is None: st['cost_inf'] += 1
            elif ce == pv: st['cost_eq'] += 1
            elif ce < pv: st['cost_lt'] += 1
    # floor with EXACT costs and EXACT knows: monkeypatch the cost inside lagmitm through a pattern + cost table
    Nk = [[0 if knows_exact(P, S, T) else 1 for S in range(v)] for T in range(v)]   # Nk[T][S] = 0 iff S knows T
    for T in range(v): Nk[T][T] = 1
    table = {}
    for R in subs:
        for Dd in subs:
            table[(R, Dd)] = cost_exact(P, R, Dd)
    def rank_patch(rows, _t=table):
        raise RuntimeError
    orig = lagmitm.terms
    def terms_exact(N, p, _t=table):
        best, noff = orig(N, p)
        return best, noff
    # exact costs: recompute terms with the table (infinite cost = term not available)
    import itertools as it
    def terms_tab(N, p):
        vv = len(N)
        knows = [[S != T and N[T][S] == 0 for T in range(vv)] for S in range(vv)]
        off = [(S, T) for S in range(vv) for T in range(vv) if S != T and not knows[S][T]]
        oi = {e: i for i, e in enumerate(off)}
        best = {}
        for r in it.product(range(p), repeat=vv):
            if not any(r) or r[next(i for i in range(vv) if r[i])] != 1: continue
            Rs = tuple(T for T in range(vv) if r[T])
            for g in it.product(range(p), repeat=vv):
                if not any(g): continue
                Ds = tuple(S for S in range(vv) if g[S])
                c = table[(Rs, Ds)]
                if c is None: continue
                od = [0] * len(off)
                for S in Ds:
                    for T in Rs:
                        i = oi.get((S, T))
                        if i is not None: od[i] = g[S] * r[T] % p
                key = (tuple(od), tuple(g[S] * r[S] % p for S in range(vv)))
                if key not in best or c < best[key]: best[key] = c
        return best, len(off)
    lagmitm.terms = terms_tab
    lt, c = lagmitm.floor_lt_v(Nk, 3)
    lagmitm.terms = orig
    f = c if lt else v
    st['floors'][f] = st['floors'].get(f, 0) + 1
say('random twisted %d-port families at h=%d (%d families, random odd labels, fixed twists of type %s):' % (v, h, nfam, MODE))
say('  ordered pairs %d: exact knows %d, vector condition omega(p_T,q_S)=0 %d, exact-but-not-vector %d (must be 0)'
    % (st['pairs'], st['knows_exact'], st['knows_vec'], st['mismatch_dir']))
say('  (R,D) totals %d: no common read or delivery frame %d; exact cost == pairing rank %d; exact < rank %d (must be 0)'
    % (st['cost_cmp'], st['cost_inf'], st['cost_eq'], st['cost_lt']))
say('  floors (exact frames, payload F_3): %s  (v = %d)' % (st['floors'], v))
say('time %.1fs' % (time.time() - t0))
open(os.path.join(HERE, 'out', 'lagframes_h%d_v%d_%s.txt' % (h, v, MODE)), 'w').write('\n'.join(out) + '\n')
