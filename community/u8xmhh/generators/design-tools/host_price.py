"""x-hosting: price idle-axis hosting (DaysSky's hub #442 idea) on Ryan Shea's complex588 Lean route.

Model of the Lean layout (Work/Bridge401/GNetStage.lean, Work/Bridge377/Bank.lean of shea256/fourier-transform-below-nlogn,
jacklightChen #377 / Dugongue #401 layout, our engine): 8 runs (replicas) one after the other, 9 stages each, call k = 9*run+stage,
local pair of stage s = PAIR[s].  Per (replica, side X/Y, port t) the idle climbs ("children") are
  A2 16 after stage 1, B2 8 in the window of stage 2, A3 32 after stage 3, B3 8 (stage 4), A4 48 after stage 5, B4 8 (stage 6),
  finals F1 48 (pair 1, last touch stage 2), F2 32 (pair 2, last touch stage 4), F3 16 (pair 3, last touch stage 6).
A helper residual of width w of call k may be hosted on a child of rank r >= w of a data role that is untouched from its last touch
up to the child, provided call k lies in that window (hub PROOF.md of research/idle-axis-hosting, sections 1-2); the child keeps
rank r - (sum of hosted widths).  One role hosts at most one residual per call (role embedding of one invocation is injective):
we impose  #residuals of a child <= #calls in its window  and check an explicit assignment afterwards.
Unhosted residuals go to banks (<= one residual per call per bank, load <= 73, one pad block per bank).

Integer programme (scipy milp; idea and variable layout follow DaysSky's hosting_core.solve, Apache-2.0, hub
research/idle-axis-hosting): maximise the first-order saving  sum f(r) - f(r - hosted),  f(r) = r ln(73/r).
Then the exact price: Shea's/Sussman's rational moment root and the finite-fill kernel grid point (s = 40), via route_price.price.

usage: python3 -I host_price.py WORDS.json NAME [--copies 0|1] [--defer 0|1] [--out FILE]
"""
import sys, os, json, math, time
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# as run, a local folder holding scipy 1.13.1 and numpy 2.0.2 was added to sys.path here
from route_price import price, unit, widths_from

M, RUNS, STAGES, V, SIDES = 73, 8, 9, 112, 2
INST = SIDES * V                       # child instances per (replica, child type), per orthogonal class
PAIR = [1, 1, 1, 2, 2, 3, 3, 4, 4]
CH = [('A2', 16, 2, 1, 'A'), ('B2', 8, 2, 2, 'B'), ('A3', 32, 3, 3, 'A'), ('B3', 8, 3, 4, 'B'),
      ('A4', 48, 4, 5, 'A'), ('B4', 8, 4, 6, 'B'), ('F1', 48, 1, 2, 'F'), ('F2', 32, 2, 4, 'F'), ('F3', 16, 3, 6, 'F')]


def f(r):
    return r * math.log(M / r) if r > 0 else 0.0


def windows(defer):
    win = {}
    for rp in range(RUNS):
        for name, r, pair, pos, kind in CH:
            if kind == 'A':     # untouched since the start of the unit; earlier runs act on other replicas
                ks = ([] if SAMERUN[0] else [9 * rho + s for rho in range(rp) for s in range(9)]) + [9 * rp + s for s in range(pos + 1)]
            elif kind == 'B':   # between gate layer A and the second climb: the one stage `pos`
                ks = [9 * rp + pos]
            else:               # finals: after the last touch `pos` of the pair; deferred finals also see later runs
                ks = [9 * rp + s for s in range(pos + 1, 9)]
                if defer:
                    ks += [9 * rho + s for rho in range(rp + 1, RUNS) for s in range(9)]
            for k in ks:
                assert PAIR[k % 9] != pair or k // 9 != rp
            win[(rp, name)] = (r, ks)
    return win


def patterns(rank, widths, maxk):
    ws = sorted(set(widths), reverse=True); res = []
    def rec(i, cur, s):
        if i == len(ws):
            res.append(tuple(cur)); return
        w = ws[i]; k = 0
        while s + k * w <= rank and len(cur) + k <= maxk:
            rec(i + 1, cur + [w] * k, s + k * w); k += 1
    rec(0, [], 0)
    return res


SHORT = 6
SAMERUN = [0]   # h-L3-design: 1 = windows restricted to the host's own run (no cross-replica hosting, no deferral)
CACHE = [None]
DROPPED = [0, 0]                       # residuals the expansion could not place (count, width): they go to banks


def combos(r, widths, ks):
    """one residual (or none) per call of a short window; only sums that lower f (or zero)"""
    out = []
    def rec(i, cur, s):
        if i == len(ks):
            if s == 0 or f(r - s) < f(r) - 1e-12:
                out.append(tuple(sorted(cur, reverse=True)))
            return
        rec(i + 1, cur, s)
        for w in widths:
            if s + w <= r:
                rec(i + 1, cur + [(w, ks[i])], s + w)
    rec(0, [], 0)
    return sorted(set(out))


def solve(win, supply, time_limit=900):
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import lil_matrix
    widths = sorted(supply)
    var = []; idx = {}
    for key, (r, ks) in win.items():
        if len(ks) <= SHORT:
            for cb in combos(r, widths, ks):
                idx[('x', key, cb)] = len(var); var.append(('x', key, cb))
            continue
        for pat in patterns(r, widths, min(12, len(ks))):
            idx[('y', key, pat)] = len(var); var.append(('y', key, pat))
        for w in widths:
            for k in ks:
                idx[('z', key, w, k)] = len(var); var.append(('z', key, w, k))
    n = len(var)
    c = np.zeros(n)
    for v, j in idx.items():
        if v[0] == 'y':
            r = win[v[1]][0]; c[j] = -(f(r) - f(r - sum(v[2])))
        elif v[0] == 'x':
            r = win[v[1]][0]; c[j] = -(f(r) - f(r - sum(w for w, k in v[2])))
    rows = []; lo = []; hi = []
    ylist = defaultdict(list); zlist = defaultdict(list); zcall = defaultdict(list)
    xlist = defaultdict(list)
    for v, j in idx.items():
        if v[0] == 'y':
            ylist[v[1]].append((j, v[2]))
        elif v[0] == 'x':
            xlist[v[1]].append(j)
            for w, k in v[2]:
                zcall[(k, w)].append(j)
        else:
            zlist[(v[1], v[2])].append(j); zcall[(v[3], v[2])].append(j)
    zkc = defaultdict(list)                        # (child, call) -> z vars: <= INST residuals of one call per child type
    for v, j in idx.items():
        if v[0] == 'z':
            zkc[(v[1], v[3])].append(j)
    A = lil_matrix((len(win) * (1 + len(widths)) + len(zcall) + len(zkc) + len(xlist), n))
    i = 0
    for key, js in xlist.items():
        for j in js:
            A[i, j] = 1
        lo.append(INST); hi.append(INST); i += 1
    for kk, js in zkc.items():
        for j in js:
            A[i, j] = 1
        lo.append(0); hi.append(INST); i += 1
    for key in win:
        if key in xlist:
            continue
        for j, pat in ylist[key]:
            A[i, j] = 1
        lo.append(INST); hi.append(INST); i += 1
        for w in widths:
            for j, pat in ylist[key]:
                if pat.count(w):
                    A[i, j] = pat.count(w)
            for j in zlist[(key, w)]:
                A[i, j] = -1
            lo.append(0); hi.append(0); i += 1
    for (k, w), js in zcall.items():
        for j in js:
            A[i, j] = 1
        lo.append(0); hi.append(supply[w]); i += 1
    A = A[:i].tocsr()
    t0 = time.time()
    res = milp(c, constraints=LinearConstraint(A, lo, hi), integrality=np.ones(n), bounds=Bounds(np.zeros(n), np.full(n, np.inf)),
               options=dict(time_limit=time_limit, disp=False, mip_rel_gap=1e-6))
    assert res.x is not None, res.message
    x = [int(round(t)) for t in res.x]
    ys = {(v[1], v[2]): x[j] for v, j in idx.items() if v[0] == 'y' and x[j]}
    zs = {(v[1], v[2], v[3]): x[j] for v, j in idx.items() if v[0] == 'z' and x[j]}
    for v, j in idx.items():
        if v[0] == 'x' and x[j]:
            if v[2]:
                ys[(v[1], v[2])] = x[j]
            else:
                ys[(v[1], ())] = ys.get((v[1], ()), 0) + x[j]
            for w, k in v[2]:
                zs[(v[1], w, k)] = zs.get((v[1], w, k), 0) + x[j]
    return ys, zs, -res.fun, dict(status=res.status, message=res.message, seconds=time.time() - t0, nvar=n,
                                  gap=getattr(res, 'mip_gap', None))


def expand(win, ys, zs, supply):
    """explicit owners: every hosted residual (call k, width w, helper index) -> (replica, child, instance, offset).
    Checks: one owner per helper, <= one residual per role per call, sum of widths <= rank of the child."""
    used = Counter(); slots = []; rem = {}
    for key, (r, ks) in win.items():
        pats = sorted([(pat, cnt) for (kk, pat), cnt in ys.items() if kk == key], key=lambda t: (-len(t[0]), -sum(x[0] if isinstance(x, tuple) else x for x in t[0])))
        pool = defaultdict(Counter)                       # width -> Counter(call -> count)
        for (kk, w, k), cnt in zs.items():
            if kk == key:
                pool[w][k] += cnt
        inst = 0
        if len(ks) > SHORT:                          # long window: Gale-Ryser style, call by call
            need = []
            for pat, cnt in pats:
                for _ in range(cnt):
                    need.append(Counter(pat))
            got = [[] for _ in need]
            calls_order = sorted({k for w in pool for k in pool[w] if pool[w][k]}, key=lambda k: -sum(pool[w][k] for w in pool))
            for k in calls_order:
                items = []
                for w in pool:
                    items += [w] * pool[w][k]
                    pool[w][k] = 0
                # Kuhn's augmenting paths: each instance takes at most one item of this call; instances with the
                # largest remaining need are tried first (Gale-Ryser order), scarce widths first
                order = sorted(range(len(need)), key=lambda i: (-sum(need[i].values()), i))
                wcount = Counter(items)
                items.sort(key=lambda w: (sum(1 for i in order if need[i][w] > 0) - wcount[w], -w))
                match = {}                                # instance -> item index
                def try_item(t, seen):
                    w = items[t]
                    for i in order:
                        if need[i][w] > 0 and i not in seen:
                            seen.add(i)
                            if i not in match or try_item(match[i], seen):
                                match[i] = t
                                return True
                    return False
                import sys as _s
                _s.setrecursionlimit(10000)
                for t in range(len(items)):
                    if not try_item(t, set()):
                        DROPPED[0] += 1; DROPPED[1] += items[t]       # this residual stays in a bank
                for i, t in match.items():
                    w = items[t]; need[i][w] -= 1; got[i].append((w, k))
            pats = Counter(tuple(sorted(g, reverse=True)) for g in got)
            pats = sorted(pats.items(), key=lambda t: (-len(t[0]), t[0]))
            pool = defaultdict(Counter)
            for g in got:
                for w, k in g:
                    pool[w][k] += 1
        for pat, cnt in pats:
            for _ in range(cnt):
                calls = set(); off = 0
                fixed = isinstance(pat[0], tuple) if pat else False
                for item in sorted(pat, reverse=True):
                    w = item[0] if fixed else item
                    cand = [item[1]] if fixed else [k for k, z in pool[w].items() if z > 0 and k not in calls]
                    assert not fixed or (pool[w][item[1]] > 0 and item[1] not in calls)
                    assert cand, ('expand failed', key, pat)
                    k = max(cand, key=lambda k: (sum(pool[x][k] for x in pool), pool[w][k], -k))
                    pool[w][k] -= 1; calls.add(k)
                    h = used[(k, w)]; used[(k, w)] += 1
                    assert h < supply[w]
                    slots.append((k, w, h, key, inst, off)); off += w
                assert off <= r
                rem[(key, inst)] = r - off; inst += 1
        assert inst == INST, (key, inst)
    # role-level check: a role = (replica, pair, instance); its children have disjoint windows
    pair_of = {name: pair for name, _, pair, _, _ in CH}
    seen = set()
    for k, w, h, (rp, name), inst, off in slots:
        role = (rp, pair_of[name], inst)
        assert (role, k) not in seen; seen.add((role, k))
        assert not (k // 9 == rp and PAIR[k % 9] == pair_of[name])
    return slots, rem, used


def exact_split(key, pats, pool, ks):
    """rewrite (pattern, count) as (assigned pattern, count): each width tagged with a distinct call, matching pool exactly"""
    import itertools
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    var = []
    for pat, cnt in pats:
        maps = set()
        for calls in itertools.permutations(ks, len(pat)):
            maps.add(tuple(sorted(zip(pat, calls))))
        for mp in sorted(maps):
            var.append((pat, mp))
    n = len(var); rows = []; lo = []; hi = []
    for pat, cnt in pats:
        rows.append([1 if v[0] == pat else 0 for v in var]); lo.append(cnt); hi.append(cnt)
    for w in pool:
        for k in ks:
            rows.append([sum(1 for (ww, kk) in v[1] if ww == w and kk == k) for v in var])
            lo.append(pool[w][k]); hi.append(pool[w][k])
    res = milp(np.zeros(n), constraints=LinearConstraint(np.array(rows), lo, hi), integrality=np.ones(n),
               bounds=Bounds(np.zeros(n), np.full(n, np.inf)))
    assert res.x is not None, ('exact split infeasible', key, res.message)
    out = []
    for v, x in zip(var, res.x):
        x = int(round(x))
        if x:
            out.append((v[1], x))
    return out


REPAIR = [1]


def repair(win, slots, rem, used, supply, calls=72):
    """host residuals the expansion left over: any child instance whose window holds the call, whose role is free in that
    call and whose remainder drops f (first-order gain).  Same checks as expand."""
    pair_of = {name: pair for name, _, pair, _, _ in CH}
    busy = set()
    for k, w, h, (rp, name), inst, off in slots:
        busy.add(((rp, pair_of[name], inst), k))
    inwin = defaultdict(list)                       # call -> keys whose window holds it
    for key, (r, ks) in win.items():
        for k in ks:
            inwin[k].append(key)
    n = 0
    for k in range(calls):
        for w in sorted(supply, reverse=True):
            while used[(k, w)] < supply[w]:
                best = None
                for key in inwin[k]:
                    rp, name = key
                    for inst in range(INST):
                        rr = rem[(key, inst)]
                        if rr < w or ((rp, pair_of[name], inst), k) in busy:
                            continue
                        gain = f(rr) - f(rr - w)
                        if gain > 1e-9 and (best is None or gain > best[0]):
                            best = (gain, key, inst)
                if best is None:
                    break
                gain, key, inst = best
                rp, name = key
                h = used[(k, w)]; used[(k, w)] += 1
                off = win[key][0] - rem[(key, inst)]
                slots.append((k, w, h, key, inst, off)); rem[(key, inst)] -= w
                busy.add(((rp, pair_of[name], inst), k)); n += 1
    return n


def banks(supply, used, calls=72):
    """first-fit decreasing of the unhosted residuals into banks of 73 coordinates, <= one residual per call per bank"""
    items = []
    for k in range(calls):
        for w, n in supply.items():
            items += [(w, k)] * (n - used[(k, w)])
    items.sort(key=lambda t: (-t[0], t[1]))
    bl = []                                         # [load, set(calls)]
    for w, k in items:
        for b in bl:
            if b[0] + w <= M and k not in b[1]:
                b[0] += w; b[1].add(k); break
        else:
            bl.append([w, {k}])
    return [M - b[0] for b in bl], len(items)


def run(word, copies, defer, time_limit=900):
    prof = Counter({int(k): n for k, n in word['profile'].items()})
    wid = Counter({int(k): n for k, n in word['widths'].items()})
    if copies:
        prof[1] += V; prof[8] += V; wid[9] += V
    for r, n in word.get('extra', {}).items():
        prof[int(r)] += n
    supply = dict(wid)
    win = windows(defer)
    import pickle
    cache = CACHE[0]
    if cache and os.path.exists(cache):
        ys, zs, proxy, info = pickle.load(open(cache, 'rb'))
    elif SAMERUN[0]:
        # h-L3-design: same-run windows make the programme block-diagonal (run r feeds only replica r): solve run 0, copy
        win0 = {key: v for key, v in win.items() if key[0] == 0}
        ys0, zs0, proxy0, info = solve(win0, supply, time_limit)
        ys, zs = {}, {}
        for rp in range(RUNS):
            for ((r0, name), pat), cnt in ys0.items():
                pat2 = tuple((w, k + 9 * rp) for w, k in pat) if pat and isinstance(pat[0], tuple) else pat
                ys[((rp, name), pat2)] = cnt
            for ((r0, name), w, k), cnt in zs0.items():
                zs[((rp, name), w, k + 9 * rp)] = cnt
        proxy = RUNS * proxy0
        if cache:
            pickle.dump((ys, zs, proxy, info), open(cache, 'wb'))
    else:
        ys, zs, proxy, info = solve(win, supply, time_limit)
        if cache:
            pickle.dump((ys, zs, proxy, info), open(cache, 'wb'))
    slots, rem, used = expand(win, ys, zs, supply)
    repaired = repair(win, slots, rem, used, supply) if REPAIR[0] else 0
    pads, left = banks(supply, used)
    hist = Counter({r: 72 * n for r, n in prof.items() if n})
    for (key, inst), r in rem.items():
        if r:
            hist[r] += 1
    for p in pads:
        if p:
            hist[p] += 1
    W = 64 * V + len(pads)
    mass = sum(r * n for r, n in hist.items())
    D = M * W - mass
    assert D == 2560, D
    u = dict(hist=dict(sorted(hist.items())), W=W, D=D)
    p = price(u)
    Lidle = sum(f(r) for r in rem.values())
    out = dict(copies=copies, defer=defer, dropped=list(DROPPED), repaired=repaired, W=W, banks=len(pads), unhosted=left, hosted=len(slots), pads=dict(Counter(pads)),
               hist=u['hist'], L_idle=round(Lidle, 1), proxy_saving=round(proxy, 1), solver=info,
               child_remainders={name: dict(Counter(r for (key, inst), r in rem.items() if key[1] == name)) for name, *_ in CH},
               **p)
    return out, slots


if __name__ == '__main__':
    args = sys.argv[1:]
    words = json.load(open(args[0])); name = args[1]
    opt = dict(copies=1, defer=1, out=None, tl=900, samerun=0)
    for i in range(2, len(args), 2):
        opt[args[i].lstrip('-')] = args[i + 1]
    CACHE[0] = opt['out'] + '.ilp.pkl' if opt['out'] else None
    SAMERUN[0] = int(opt['samerun'])
    res, slots = run(words[name], int(opt['copies']), int(opt['defer']), float(opt['tl']))
    res['name'] = name; res['samerun'] = SAMERUN[0]
    print(json.dumps({k: res[k] for k in ('name', 'copies', 'defer', 'dropped', 'repaired', 'W', 'banks', 'hosted', 'unhosted', 'L_idle',
                                         'root_lower', 'kernel_A', 'dft')}), flush=True)
    print('  remainders:', res['child_remainders'], ' solver:', res['solver']['status'], round(res['solver']['seconds'], 1), 's')
    if opt['out']:
        open(opt['out'], 'w').write(json.dumps(res, indent=1, default=str))
        import gzip
        with gzip.open(opt['out'] + '.slots.gz', 'wt') as g:
            json.dump(slots, g)
