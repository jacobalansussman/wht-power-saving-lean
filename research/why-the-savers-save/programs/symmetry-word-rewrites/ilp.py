"""t-egraph: exact fixed-order frame optimum (over a candidate pool) by MILP (scipy/HiGHS).

For a given gate order: P = lazy forward frames (joins), Q = eager backward bounds (meets, with end frames,
low-end cones, centre frames); every legal frame at position k lies in [P(k), Q(k)].  Candidates at k = every
P/Q/current frame seen on the chains of k's registers that lies in [P(k), Q(k)] and passes ok_at.
Variables x[k,U]; per register and consecutive pair of its positions (incl. start/end) a transport y[a,b,U,V]
with U <= V, cost FC[dim V - dim U] (0 if U == V).  Minimises sum d ln(73/d) over climbs.
"""
import sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
sys.path.insert(0, '.')
import fw, search
from search import le, join, meet, FC

def bounds_PQ(W):
    N = len(W.seq)
    cur = list(W.start); P = [None] * N
    for k in range(N):
        rr = W.regs[k]; U = ()
        for r in rr: U = join(U, cur[r])
        if k in W.frozen: U = join(U, W.f[k])
        P[k] = U
        for r in rr: cur[r] = U
    cur = list(W.end); Q = [None] * N
    for k in reversed(range(N)):
        rr = W.regs[k]; U = W.full
        for r in rr:
            U = meet(U, cur[r])
            if W.sec[k] != 'T': U = meet(U, W.ub[r])
            if k < W.nA and r in W.cframe: U = meet(U, W.cframe[r])
        if k in W.frozen: U = meet(U, W.f[k])
        Q[k] = U
        for r in rr: cur[r] = U
    return P, Q

def solve(W, time_limit=600, verbose=True, extra_pool=None):
    N = len(W.seq)
    P, Q = bounds_PQ(W)
    if any(not le(P[k], Q[k]) for k in range(N)):
        return None, 'P not <= Q (order infeasible)'
    pool = {}
    for r, ps in W.pos.items():
        s = {W.start[r], W.end[r]}
        for k in ps: s.update((P[k], Q[k], W.f[k]))
        pool[r] = s
    cand = []
    for k in range(N):
        s = set()
        for r in W.regs[k]: s |= pool[r]
        if extra_pool: s |= extra_pool.get(k, set())
        c = sorted({U for U in s if le(P[k], U) and le(U, Q[k]) and W.ok_at(k, U)}, key=lambda U: (len(U), U))
        if not c: return None, ('no candidate', k)
        cand.append(c)
    # variables
    xid = {}; nv = 0
    for k in range(N):
        for U in cand[k]: xid[(k, U)] = nv; nv += 1
    rows, cols, vals, lo, hi = [], [], [], [], []
    cost = []
    def addrow(entries, l, h):
        r = len(lo)
        for c, v in entries: rows.append(r); cols.append(c); vals.append(v)
        lo.append(l); hi.append(h)
    for k in range(N):
        addrow([(xid[(k, U)], 1) for U in cand[k]], 1, 1)
    yterms = []
    for r, ps in W.pos.items():
        seqp = [('S', None)] + [('k', k) for k in ps] + [('E', None)]
        for (ta, a), (tb, b) in zip(seqp, seqp[1:]):
            CA = [W.start[r]] if ta == 'S' else cand[a]
            CB = [W.end[r]] if tb == 'E' else cand[b]
            ys = []
            for U in CA:
                for V in CB:
                    if le(U, V): ys.append((U, V))
            base = nv + len(yterms)
            for (U, V) in ys:
                yterms.append(FC[len(V) - len(U)] if U != V else 0.0)
            # transport constraints
            idx = {uv: base + i for i, uv in enumerate(ys)}
            for U in CA:
                ent = [(idx[(U, V)], 1) for V in CB if (U, V) in idx]
                if ta == 'S': addrow(ent, 1, 1)
                else: addrow(ent + [(xid[(a, U)], -1)], 0, 0)
            for V in CB:
                ent = [(idx[(U, V)], 1) for U in CA if (U, V) in idx]
                if tb == 'E': addrow(ent, 1, 1)
                else: addrow(ent + [(xid[(b, V)], -1)], 0, 0)
    ntot = nv + len(yterms)
    c = np.zeros(ntot); c[nv:] = yterms
    A = coo_matrix((vals, (rows, cols)), shape=(len(lo), ntot)).tocsr()
    integrality = np.zeros(ntot); integrality[:nv] = 1
    if verbose: print('milp vars', ntot, 'x', nv, 'rows', len(lo), 'free positions', sum(1 for k in range(N) if len(cand[k]) > 1), flush=True)
    t = time.time()
    res = milp(c, constraints=LinearConstraint(A, lo, hi), integrality=integrality, bounds=Bounds(0, 1),
               options=dict(time_limit=time_limit, disp=False))
    if verbose: print('status', res.status, res.message, 'obj', res.fun, 'secs', round(time.time() - t, 1), flush=True)
    if res.x is None: return None, res.message
    fr = []
    for k in range(N):
        best = max(cand[k], key=lambda U: res.x[xid[(k, U)]]); fr.append(best)
    return fr, dict(obj=res.fun, status=int(res.status), mip_gap=getattr(res, 'mip_gap', None))

if __name__ == '__main__':
    W = search.Word(fw.load(sys.argv[1] if len(sys.argv) > 1 else 'K8-FULL-WORD.json.gz'))
    print('orig cost', W.total_cost())
    fr, info = solve(W)
    print(info)
    if fr:
        W.f = fr
        print('new cost', W.total_cost(), dict(sorted(W.profile().items())))
