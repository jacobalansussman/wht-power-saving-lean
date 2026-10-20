# cert3.py CUTSFILE -- turn the LP bound of lp3cut.py into an EXACT rational dual certificate (pf-duality).
# The certificate: lambda >= 0 on inequality rows (elemental polymatroid, DFZ char != 2, Ingleton), mu free on the
# equality rows (uniformity, normalisation mean delta_1 = 1), with  A_ub^T lambda + A_eq^T mu >= c  componentwise.
# Then for every feasible r:  mean delta_3 = c.r <= mu_norm.  Everything is re-verified with Fractions.
import sys, os, json
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
import lp3
from lp3 import orb, EMPTY, NV, NF, rows, lin, dist, bit
from sympy.polys.matrices import DomainMatrix
from sympy import QQ

data = json.load(open(sys.argv[1]))
UNIFORM = data["uniform"]
cuts = [(nm, {int(k): v for k, v in r.items()}) for nm, r in data["cuts"]]
singles = [bit(i) for i in range(NF)]
pairs = {d: [(bit(i) | bit(j), bit(j)) for i in range(NF) for j in range(NF) if i != j and dist(i, j) == d] for d in (1, 2, 3)}
def linF(terms):
    d = {}
    for c, mk in terms:
        v = int(orb[mk])
        if v == EMPTY: continue
        d[v] = d.get(v, F(0)) + F(c)
    return {k: v for k, v in d.items() if v != 0}
ub = [({k: F(-v) for k, v in r.items()}, 'elem') for r in rows] + [({k: F(-v) for k, v in r.items()}, nm) for nm, r in cuts]
eq = []; eqdef = []
if UNIFORM:
    for s in singles[1:]:
        e = linF([(1, s), (-1, singles[0])])
        if e: eq.append((e, F(0))); eqdef.append(["single", s, singles[0]])
    for d in (1, 2, 3):
        p0 = pairs[d][0]
        for p in pairs[d][1:]:
            e = linF([(1, p[0]), (-1, p[1]), (-1, p0[0]), (1, p0[1])])
            if e: eq.append((e, F(0))); eqdef.append(["pair", p[0], p[1], p0[0], p0[1]])
def mean(d):
    acc = {}
    for (pm, tm) in pairs[d]:
        for k, v in linF([(1, pm), (-1, tm)]).items(): acc[k] = acc.get(k, F(0)) + v / len(pairs[d])
    return acc
eq.append((mean(1), F(1))); eqdef.append(["norm"])
cobj = mean(3)
def mat(L):
    I, J, X = [], [], []
    for r_, (row, _) in enumerate(L):
        for k, v in row.items(): I.append(r_); J.append(k); X.append(float(v))
    return coo_matrix((X, (I, J)), shape=(len(L), NV)).tocsr()
c = np.zeros(NV)
for k, v in cobj.items(): c[k] = -float(v)
bounds = [(0, 0) if k == EMPTY else (0, None) for k in range(NV)]
res = linprog(c, A_ub=mat(ub), b_ub=np.zeros(len(ub)), A_eq=mat(eq), b_eq=[float(b) for _, b in eq], bounds=bounds, method='highs')
opt = -res.fun; x = res.x
lam = -res.ineqlin.marginals; mu = -res.eqlin.marginals
print("float optimum", opt, "ineq rows", len(ub), "eq rows", len(eq))
S = [i for i in range(len(ub)) if lam[i] > 1e-11]
P = [j for j in range(NV) if j != EMPTY and x[j] > 1e-9]
print("dual support", len(S), "positive primal vars", len(P))
cols = [('l', i) for i in S] + [('m', e) for e in range(len(eq))]
def entry(col, j):
    t, i = col
    return (ub[i][0].get(j, F(0))) if t == 'l' else eq[i][0].get(j, F(0))
M = [[entry(col, j) for col in cols] + [cobj.get(j, F(0))] for j in P]
dm = DomainMatrix([[QQ(v.numerator, v.denominator) for v in row] for row in M], (len(P), len(cols) + 1), QQ)
R, piv = dm.rref()
R = R.to_Matrix()
ncol = len(cols)
assert ncol not in piv, "complementary-slackness system inconsistent"
free = [k for k in range(ncol) if k not in piv]
def floatval(col):
    t, i = col
    return lam[i] if t == 'l' else mu[i]
y = [None] * ncol
for k in free: y[k] = F(floatval(cols[k])).limit_denominator(10 ** 7)
for r_, pk in enumerate(piv):
    val = F(int(R[r_, ncol].p), int(R[r_, ncol].q))
    for k in free:
        a = R[r_, k]
        if a != 0: val -= F(int(a.p), int(a.q)) * y[k]
    y[pk] = val
# verification, all exact
ok = True
for k, col in enumerate(cols):
    if col[0] == 'l' and y[k] < 0: ok = False; print("negative multiplier", col, y[k])
tot = {}
for k, col in enumerate(cols):
    row = ub[col[1]][0] if col[0] == 'l' else eq[col[1]][0]
    for j, v in row.items(): tot[j] = tot.get(j, F(0)) + y[k] * v
for j in range(NV):
    if j == EMPTY: continue
    if tot.get(j, F(0)) < cobj.get(j, F(0)): ok = False; print("dual constraint fails at var", j, tot.get(j, F(0)) - cobj.get(j, F(0)))
bound = sum(y[len(S) + e] * eq[e][1] for e in range(len(eq)))
print("EXACT CERTIFICATE VALID:", ok, " bound: mean delta_3 <= ", bound, "=", float(bound), "x mean delta_1")
used = sorted([(y[k], ub[cols[k][1]][1]) for k in range(len(S)) if y[k] > 0 and ub[cols[k][1]][1] != 'elem'], reverse=True)
print("non-Shannon / characteristic rows used (exact multipliers):")
for v, nm in used: print("  ", v, nm)
print("elemental rows used:", sum(1 for k in range(len(S)) if y[k] > 0 and ub[cols[k][1]][1] == 'elem'))

# witnesses for the elemental rows used, and a dump for the independent verifier verify3.py
used_rows = {tuple(sorted(ub[cols[k][1]][0].items())): k for k in range(len(S)) if y[k] > 0 and ub[cols[k][1]][1] == 'elem'}
wit = {}
masks = np.arange(1 << NF, dtype=np.int64)
full = (1 << NF) - 1
for i in range(NF):
    key = tuple(sorted({k: F(-v) for k, v in lin([(1, full), (-1, full ^ (1 << i))]).items()}.items()))
    if key in used_rows and key not in wit: wit[key] = ["mono", i]
for i in range(NF):
    for j in range(i + 1, NF):
        Sm = masks[((masks >> i) & 1 == 0) & ((masks >> j) & 1 == 0)]
        a = orb[Sm | (1 << i)]; b = orb[Sm | (1 << j)]; cc = orb[Sm | (1 << i) | (1 << j)]; d = orb[Sm]
        for t in range(len(Sm)):
            dd = {}
            for v, cf in ((int(a[t]), 1), (int(b[t]), 1), (int(cc[t]), -1), (int(d[t]), -1)):
                if v == EMPTY: continue
                dd[v] = dd.get(v, 0) + cf
            key = tuple(sorted((k, F(-v)) for k, v in dd.items() if v != 0))
            if key in used_rows and key not in wit: wit[key] = ["sub", int(Sm[t]), i, j]
        if len(wit) == len(used_rows): break
    if len(wit) == len(used_rows): break
assert len(wit) == len(used_rows), (len(wit), len(used_rows))
reps = [int(r) for r in lp3.reps]
dump = {"uniform": UNIFORM, "bound": str(bound), "objective": "mean delta_3", "normalisation": "mean delta_1 = 1",
        "reps": reps, "rows": [], "eqs": []}
for key, k in used_rows.items():
    dump["rows"].append({"type": "elemental", "witness": wit[key], "mult": str(y[k])})
for k in range(len(S)):
    nm = ub[cols[k][1]][1]
    if y[k] > 0 and nm != 'elem': dump["rows"].append({"type": nm.split(':')[0], "label": nm, "mult": str(y[k])})
for e in range(len(eq)):
    if y[len(S) + e] != 0: dump["eqs"].append({"def": eqdef[e], "mult": str(y[len(S) + e])})
json.dump(dump, open(sys.argv[1].replace('cuts_', 'certificate_'), 'w'))
print("certificate written:", len(dump["rows"]), "inequality rows,", len(dump["eqs"]), "equality rows")
