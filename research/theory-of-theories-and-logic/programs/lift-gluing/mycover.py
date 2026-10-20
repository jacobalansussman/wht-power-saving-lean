# mycover.py m K4 TL -- independent ILP for the "Fano cover + gluing row" lower bound on b_m (all lifts of one interval).
# Own enumeration (subspaces as sorted tuples of bitmask spans), own row construction, two solvers if available.
# Rows for every 3-dim sub-interval [a,b] (edges s<t with a<=s, t<=b), with j = dim s - dim a in {0,1,2}:
#   R1: sum_{j in {0,1}} e >= 1   R2: sum_{j in {1,2}} e >= 1   R3: sum e >= 2   R4: 2*sum_{j=1} e + sum_{j in {0,2}} e >= 4
# K4 > 0: every 4-dim sub-interval has sum e >= K4.
import sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
m = int(sys.argv[1]); K4 = int(sys.argv[2]); TL = float(sys.argv[3])
def span(gs):
    S = {0}
    for g in gs:
        if g not in S: S |= {s ^ g for s in S}
    return frozenset(S)
subs = {frozenset({0})}; layer = [frozenset({0})]
while layer:
    nxt = []
    for S in layer:
        for v in range(1, 1 << m):
            if v not in S:
                T = span(list(S) + [v])
                if T not in subs: subs.add(T); nxt.append(T)
    layer = nxt
dim = {S: len(S).bit_length() - 1 for S in subs}
bydim = {}
for S in subs: bydim.setdefault(dim[S], []).append(S)
edges = [(s, t) for k in range(m) for s in bydim[k] for t in bydim[k + 1] if s <= t]
eid = {e: i for i, e in enumerate(edges)}
up = {s: [t for t in bydim.get(dim[s] + 1, []) if s <= t] for s in subs}
def sub_edges(a, b):
    out = []; seen = {a}; layer = [a]
    while layer:
        nl = []
        for s in layer:
            for t in up[s]:
                if t <= b:
                    out.append((s, t))
                    if t not in seen: seen.add(t); nl.append(t)
        layer = nl
    return out
rows = []; lb = []
for a in subs:
    for b in subs:
        if not a <= b: continue
        k = dim[b] - dim[a]
        if k == 3:
            es = sub_edges(a, b); assert len(es) == 35
            j = {e: dim[e[0]] - dim[a] for e in es}
            rows.append({eid[e]: 1 for e in es if j[e] in (0, 1)}); lb.append(1)
            rows.append({eid[e]: 1 for e in es if j[e] in (1, 2)}); lb.append(1)
            rows.append({eid[e]: 1 for e in es}); lb.append(2)
            rows.append({eid[e]: (2 if j[e] == 1 else 1) for e in es}); lb.append(4)
        elif k == 4 and K4 > 0:
            es = sub_edges(a, b); assert len(es) == 240
            rows.append({eid[e]: 1 for e in es}); lb.append(K4)
n = len(edges)
A = lil_matrix((len(rows), n))
for i, r in enumerate(rows):
    for c, v in r.items(): A[i, c] = v
A = A.tocsr(); lb = np.array(lb, float)
print(f"m={m}: {len(subs)} subspaces, {n} edges, {len(rows)} rows (K4={K4})", flush=True)
t0 = time.time()
lp = milp(np.ones(n), constraints=LinearConstraint(A, lb, np.inf), integrality=np.zeros(n), bounds=Bounds(0, np.inf))
print("LP:", lp.fun, flush=True)
res = milp(np.ones(n), constraints=LinearConstraint(A, lb, np.inf), integrality=np.ones(n), bounds=Bounds(0, np.inf),
           options={"time_limit": TL})
print("scipy/HiGHS:", res.status, res.message, "incumbent", res.fun, "dual bound", getattr(res, "mip_dual_bound", None),
      "time %.0f s" % (time.time() - t0), flush=True)
if res.x is not None:
    x = np.round(res.x).astype(int); assert (A @ x >= lb - 1e-9).all(); print("incumbent verified feasible, sum", x.sum())
# second solver: OR-tools CP-SAT if installed (exact integer arithmetic)
try:
    from ortools.sat.python import cp_model
except Exception as ex:
    print("ortools not available:", ex); sys.exit(0)
mdl = cp_model.CpModel(); X = [mdl.NewIntVar(0, max(4, K4), f"x{i}") for i in range(n)]
for r, l in zip(rows, lb): mdl.Add(sum(v * X[c] for c, v in r.items()) >= int(l))
mdl.Minimize(sum(X)); sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = TL; sol.parameters.num_workers = 8
st = sol.Solve(mdl)
print("CP-SAT:", sol.StatusName(st), "objective", sol.ObjectiveValue(), "bound", sol.BestObjectiveBound(), flush=True)
