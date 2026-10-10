"""xmilp.py -- the integer programme for the WHOLE hierarchy over orbits of the order-64 group.

Class of designs (exactly this):
  * a pool of signed sums, closed under the group (given as a file);
  * the set of sums built is a union of group orbits; the set of uses (target, sum) is a union of group orbits;
  * every needed (target, source) is covered by exactly one use;
  * every built sum is formed exactly once, as the sum or difference of two built sums (or single sources) of the pool; the pairs are taken
    by whole orbits; node types: ORDINARY (one of u + v, u - v is formed from the pair {u, v}) and IN-PLACE PAIR (both
    are formed from the pair, on the two helpers that held u and v);
  * if a sum's own stabiliser moves the pair (the orbit holds k > 1 pairs for one sum) one of them is used and the
    sum is counted as ordinary.
Objective ("helpers in the copy-on-use count"): P = uses + ordinary sums  (an in-place pair costs nothing; a helper
ends in a target or dies in an ordinary sum).  With lam < 1 an ordinary sum that has a single source as a part costs
lam (the x array can be the part that dies).
Used through solve.py.  Needs numpy and scipy (scipy.optimize.milp, which calls the HiGHS solver shipped inside scipy)."""
import sys, time
sys.dont_write_bytecode = True
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
from xlib import *
import hio
from collections import Counter, defaultdict

def build(pool, A, log=print):
    t0 = time.time()
    src = [(1 << T, 0) for T in range(V)]
    syms = src + sorted(s for s in pool if s[0] & (s[0] - 1))
    sid = {s: i for i, s in enumerate(syms)}
    n = len(syms)
    IMG = [[sid[x] for x in A.images(s)] for s in syms]          # closed pool: every image exists
    PG = [A.els[g][0] for g in range(A.n)]
    # ---- orbits of sums
    so = [-1] * n; so_size = []
    for i in range(n):
        if so[i] < 0:
            o = set(IMG[i]);
            for j in o: so[j] = len(so_size)
            so_size.append(len(o))
    # ---- orbits of incidences
    inc = {}; inc_size = []
    for S in range(V):
        for T in bits(NM[S]):
            if (S, T) not in inc:
                o = set((PG[g][S], PG[g][T]) for g in range(A.n))
                for z in o: inc[z] = len(inc_size)
                inc_size.append(len(o))
    # ---- orbits of uses (flags)
    fo = {}; F = []       # F: [size, sum orbit, set of incidence orbits] ; only usable ones
    nbad = 0
    for i in range(n):
        for S in range(V):
            if match(S, syms[i]) and (S, i) not in fo:
                o = set((PG[g][S], IMG[i][g]) for g in range(A.n))
                acc = {}; ok = True; cov = set()
                for S2, j in o:
                    m = syms[j][0]
                    if acc.get(S2, 0) & m: ok = False; break
                    acc[S2] = acc.get(S2, 0) | m
                if ok:
                    for S2, j in o:
                        for T in bits(syms[j][0]): cov.add(inc[(S2, T)])
                    fid = len(F); F.append([len(o), so[i], cov, (S, i)])
                else:
                    fid = -1; nbad += 1
                for z in o: fo[z] = fid
    # ---- decompositions
    bylow = defaultdict(list)
    for i in range(n): bylow[low(syms[i][0])].append(i)
    trip = {}; D = []    # D: [size, w orbit, u orbit, v orbit, k, rep (w,u,v), pair orbit id, single-source part?]
    po = {}; Q = []      # pair orbits: [size, [triple orbit ids]]
    ntr = 0
    for w in range(V, n):
        m, ng = syms[w]
        for u in bylow[low(m)]:
            mu, nu = syms[u]
            if mu == m or mu & ~m or nu != (ng & mu): continue
            mv = m ^ mu
            v = sid.get(norm(mv, ng & mv))
            if v is None: continue
            ntr += 1
            key = (w, min(u, v), max(u, v))
            if key in trip: continue
            o = set()
            for g in range(A.n):
                a, b = IMG[u][g], IMG[v][g]
                o.add((IMG[w][g], min(a, b), max(a, b)))
            pk = (min(u, v), max(u, v))
            if pk not in po:
                oq = set((min(IMG[u][g], IMG[v][g]), max(IMG[u][g], IMG[v][g])) for g in range(A.n))
                for z in oq: po[z] = len(Q)
                Q.append([len(oq), []])
            did = len(D)
            for z in o: trip[z] = did
            k, r = divmod(len(o), so_size[so[w]]); assert r == 0
            D.append([len(o), so[w], so[u], so[v], k, key, po[pk], (u < V or v < V)])
            Q[po[pk]][1].append(did)
    log("sums %d (orbits %d), incidence orbits %d, use orbits usable %d (unusable %d), decompositions %d (orbits %d, with k>1: %d), pair orbits %d  %.0fs" % (
        n - V, len(so_size) - len(set(so[:V])), len(inc_size), len(F), nbad, ntr, len(D), sum(1 for d in D if d[4] > 1), len(Q), time.time() - t0))
    return dict(syms=syms, sid=sid, IMG=IMG, PG=PG, so=so, so_size=so_size, inc=inc, inc_size=inc_size, F=F, fo=fo, D=D, trip=trip, Q=Q, po=po)

def solve(M, lam=1.0, tlim=300, forbid=(), log=print, extra_cost=None, gap=0.0):
    syms, so, so_size, F, D, Q = M["syms"], M["so"], M["so_size"], M["F"], M["D"], M["Q"]
    srcorb = set(so[:V])
    nF, nD, nQ, nO = len(F), len(D), len(Q), len(so_size)
    # variables: x (nF), a (nD), y (nO), c (nQ)
    oX, oA, oY, oC = 0, nF, nF + nD, nF + nD + nO
    nv = oC + nQ
    cost = np.zeros(nv); lb = np.zeros(nv); ub = np.ones(nv)
    for f in range(nF): cost[oX + f] = F[f][0]
    for o in srcorb: lb[oY + o] = 1
    for o in forbid: ub[oY + o] = 0
    rows = []      # (dict col->coef, lo, hi)
    for q in range(nQ):
        ds = [d for d in Q[q][1] if D[d][4] == 1]
        for d in Q[q][1]:
            if D[d][4] > 1:
                cost[oA + d] = so_size[D[d][1]] * (lam if D[d][7] else 1.0)
        if not ds:
            continue
        wq = Q[q][0] * (lam if D[ds[0]][7] else 1.0)
        if len(ds) == 1:
            if D[ds[0]][0] == 2 * Q[q][0]:
                pass                                   # the pair's stabiliser swaps the two sums: always an in-place pair
            else:
                cost[oA + ds[0]] = wq                  # the partner is not available: ordinary
        else:
            assert len(ds) == 2
            cost[oC + q] = wq
            rows.append(({oC + q: 1, oA + ds[0]: -1, oA + ds[1]: 1}, 0, np.inf))
            rows.append(({oC + q: 1, oA + ds[0]: 1, oA + ds[1]: -1}, 0, np.inf))
    if extra_cost: extra_cost(cost, oX, oA, oY, oC)
    # cover
    cov = defaultdict(dict)
    for f in range(nF):
        for rho in F[f][2]: cov[rho][oX + f] = 1
    for rho in range(len(M["inc_size"])):
        rows.append((cov[rho], 1, 1))
    for f in range(nF):
        if F[f][1] not in srcorb: rows.append(({oX + f: 1, oY + F[f][1]: -1}, -np.inf, 0))
    bd = defaultdict(dict)
    for d in range(nD):
        bd[D[d][1]][oA + d] = -1
        for ch in (D[d][2], D[d][3]):
            if ch not in srcorb: rows.append(({oA + d: 1, oY + ch: -1}, -np.inf, 0))
    for o in range(nO):
        if o in srcorb: continue
        r = dict(bd[o]); r[oY + o] = 1
        rows.append((r, 0, 0))                        # a built sum is formed exactly once (copies are made from it)
    # a built sum is used at least once (else the second sum of an in-place pair would be a free, dead helper)
    used = defaultdict(dict)
    for f in range(nF): used[F[f][1]][oX + f] = -1
    for d in range(nD):
        for ch in (D[d][2], D[d][3]): used[ch][oA + d] = -1
    for o in range(nO):
        if o in srcorb: continue
        r = dict(used[o]); r[oY + o] = 1
        rows.append((r, -np.inf, 0))
    R_, C_, V_ = [], [], []; lo = []; hi = []
    for i, (r, l, h) in enumerate(rows):
        for c, x in r.items(): R_.append(i); C_.append(c); V_.append(float(x))
        lo.append(l); hi.append(h)
    Amat = coo_matrix((V_, (R_, C_)), shape=(len(rows), nv)).tocsr()
    t0 = time.time()
    res = milp(cost, constraints=[LinearConstraint(Amat, np.array(lo), np.array(hi))], integrality=np.ones(nv), bounds=Bounds(lb, ub),
               options=dict(time_limit=tlim, mip_rel_gap=gap, disp=False))
    info = dict(status=res.status, message=res.message, fun=None if res.x is None else float(res.fun),
                bound=getattr(res, "mip_dual_bound", None), gap=getattr(res, "mip_gap", None), nodes=getattr(res, "mip_node_count", None),
                nvar=nv, nrow=len(rows), secs=time.time() - t0)
    log("milp: %d variables, %d rows; status %s (%s); objective %s; bound %s; gap %s; nodes %s; %.0fs" % (
        nv, len(rows), res.status, str(res.message)[:50], info["fun"], info["bound"], info["gap"], info["nodes"], info["secs"]))
    if res.x is None: return None, info
    x = res.x
    sol = dict(F=[f for f in range(nF) if x[oX + f] > 0.5], D=[d for d in range(nD) if x[oA + d] > 0.5],
               O=[o for o in range(nO) if x[oY + o] > 0.5 and o not in srcorb])
    return sol, info

def hierarchy(M, sol, A, clean_first=False):
    """the chosen orbits -> a hierarchy (content, frame, kids, need, twin).
    clean_first: the compiler leaves the SECOND sum of an in-place pair (the one with the larger index) on a helper
    that holds halves, and it does not form an in-place pair on such a helper.  With clean_first the sums that are
    parts of an in-place pair get the smaller indices among sums of equal size, so they are the first of their pair."""
    syms, IMG, PG, D, F, trip = M["syms"], M["IMG"], M["PG"], M["D"], M["F"], M["trip"]
    chosenD = set(sol["D"])
    dec = defaultdict(list)                   # w -> [(u, v, triple orbit)]
    for (w, u, v), d in trip.items():
        if d in chosenD: dec[w].append((u, v, d))
    uses = []
    for f in sol["F"]:
        S, i = F[f][3]
        for z in set((PG[g][S], IMG[i][g]) for g in range(A.n)): uses.append(z)
    need_syms = set(); todo = [i for _, i in uses]
    pick = {}
    while todo:
        w = todo.pop()
        if w in need_syms: continue
        need_syms.add(w)
        if w < V: continue
        assert dec[w], "a built sum without a decomposition"
        u, v, d = min(dec[w])
        pick[w] = (u, v, d); todo += [u, v]
    cnt2 = Counter((pick[w][0], pick[w][1]) for w in pick)
    feeds = set()
    for (u, v), c in cnt2.items():
        if c == 2: feeds.add(u); feeds.add(v)
    order = sorted((w for w in need_syms if w >= V), key=lambda w: (bin(syms[w][0]).count("1"), (0 if w in feeds else 1) if clean_first else 0, syms[w]))
    idx = {T: T for T in range(V)}
    h = H(); h.content = [{T: Fr(1)} for T in range(V)]; h.frame = [e8.line(T) for T in range(V)]; h.kids = [None] * V; h.twin = {}
    bypair = defaultdict(list)
    for w in order:
        u, v, d = pick[w]
        cu, cv = h.content[idx[u]], h.content[idx[v]]
        m, ng = syms[w]
        def sg(T): return -1 if ng >> T & 1 else 1
        Tu = next(iter(cu)); Tv = next(iter(cv))
        fu = sg(Tu) * cu[Tu]; fv = sg(Tv) * cv[Tv]
        r = fv * fu
        c = dict(cu)
        for T, x in cv.items(): c[T] = r * x
        assert sym_of(c) == syms[w]
        idx[w] = len(h.content)
        h.content.append(c); h.frame.append(e8.join(h.frame[idx[u]], h.frame[idx[v]])); h.kids.append((idx[u], idx[v], r))
        bypair[(u, v)].append(w)
    npair = 0
    for (u, v), ws in bypair.items():
        if len(ws) == 2:
            a, b = idx[ws[0]], idx[ws[1]]
            assert h.kids[a][:2] == h.kids[b][:2] and h.kids[a][2] != h.kids[b][2]
            h.twin[a] = b; h.twin[b] = a; npair += 1
    h.need = [dict() for _ in range(V)]
    for S, i in uses:
        c = h.content[idx[i]]; T = next(iter(c))
        h.need[S][idx[i]] = Fr(G[S][T], 2) / c[T]
    # exactness
    for S in range(V):
        tot = defaultdict(Fr)
        for u, cf in h.need[S].items():
            for T, x in h.content[u].items(): tot[T] += cf * x
        assert {T: x for T, x in tot.items() if x} == need(S), S
    return h

def save_h(h, path):
    hio.dump_h(dict(content=h.content, frame=h.frame, kids=h.kids, need=h.need, twin=h.twin), path)
