# cert_m.py CUTSFILE -- exact rational dual certificate for the size-restricted interval LP (pf-duality, 14 Oct 2026).
# Re-solves the LP of cutsearch.py with the cuts in CUTSFILE, takes the positive primal variables, solves the
# complementary-slackness system exactly (sympy, QQ), checks dual feasibility with Fractions, and writes a certificate in
# which every row is given by explicit FRAMES (witness sets / letters / T), so verify_m.py can rebuild it independently.
import sys, os, json, re, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from lpm import Interval, DFZ, ING, inst_row

data = json.load(open(sys.argv[1]))
m, K, UNIFORM = data["m"], data["K"], data["uniform"]
I = Interval(m, K)
NF = I.NF; allf = np.arange(NF)
# polymatroid rows WITH witnesses
rows = {}   # key -> (dict, witness)
def addrow(d, wit):
    d = {k: v for k, v in d.items() if v}
    if not d: return
    key = tuple(sorted(d.items()))
    if key not in rows: rows[key] = (d, wit)
for j in range(K - 1):
    for S in I.reps[j]:
        comp = [f for f in range(NF) if f not in S]
        for ia in range(len(comp)):
            for ib in range(ia + 1, len(comp)):
                a, b = comp[ia], comp[ib]
                d = {}
                for fr, c in ((S + (a,), 1), (S + (b,), 1), (S + (a, b), -1), (S, -1)):
                    v = I.id1(fr)
                    if v: d[v] = d.get(v, 0) + c
                addrow(d, ["sub", list(S), a, b])
for j in range(K):
    for S in I.reps[j]:
        for a in range(NF):
            if a in S: continue
            d = {}
            for fr, c in ((S + (a,), 1), (S, -1)):
                v = I.id1(fr)
                if v: d[v] = d.get(v, 0) + c
            addrow(d, ["mono", list(S), a])
print("polymatroid rows with witnesses:", len(rows)); sys.stdout.flush()
ub = [(d, wit) for d, wit in rows.values()]
for nm, r in data["cuts"]:
    ub.append(({int(k): v for k, v in r.items()}, ["cut", nm]))
# equalities with witnesses
eq = []
if UNIFORM:
    s0 = I.id1([0]); seen = {s0}
    for i in range(NF):
        v = I.id1([i])
        if v not in seen: seen.add(v); eq.append(({v: 1, s0: -1}, F(0), ["single", 0, i]))
    for d in range(1, m + 1):
        IJ = np.argwhere(I.D == d); first = None; seenp = set()
        for s, t in IJ:
            key = (I.id1([s, t]), I.id1([t]))
            if key in seenp: continue
            seenp.add(key)
            if first is None: first = (int(s), int(t), key); continue
            e = {}
            for v, c in ((key[0], 1), (key[1], -1), (first[2][0], -1), (first[2][1], 1)):
                if v: e[v] = e.get(v, 0) + c
            e = {k: v for k, v in e.items() if v}
            if e: eq.append((e, F(0), ["pair", [int(s), int(t)], [first[0], first[1]]]))
def meanF(d):
    IJ = np.argwhere(I.D == d); acc = {}
    for s, t in IJ:
        for v, c in ((I.id1([s, t]), 1), (I.id1([t]), -1)):
            if v: acc[v] = acc.get(v, F(0)) + F(c, len(IJ))
    return {k: v for k, v in acc.items() if v}
eq.append((meanF(1), F(1), ["norm"]))
cobj = meanF(m)
NV = I.NV
def mat(L):
    Ii, J, X = [], [], []
    for r_, row in enumerate(L):
        for k, v in row.items(): Ii.append(r_); J.append(k); X.append(float(v))
    return coo_matrix((X, (Ii, J)), shape=(len(L), NV)).tocsr()
A_ub = [{k: -v for k, v in d.items()} for d, _ in ub]
c = np.zeros(NV)
for k, v in cobj.items(): c[k] = -float(v)
bounds = [(0, 0) if k == 0 else (0, None) for k in range(NV)]
t0 = time.time()
res = linprog(c, A_ub=mat(A_ub), b_ub=np.zeros(len(A_ub)), A_eq=mat([e for e, _, _ in eq]), b_eq=[float(b) for _, b, _ in eq],
              bounds=bounds, method='highs')
opt = -res.fun; x = res.x; lam = -res.ineqlin.marginals; mu = -res.eqlin.marginals
print("float optimum", opt, "rows", len(ub), "eqs", len(eq), f"({time.time()-t0:.1f}s)"); sys.stdout.flush()
S = [i for i in range(len(ub)) if lam[i] > 1e-10]
P = [j for j in range(1, NV) if x[j] > 1e-9]
print("dual support", len(S), "positive primal vars", len(P)); sys.stdout.flush()
cols = [('l', i) for i in S] + [('m', e) for e in range(len(eq))]
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
def entry(col, j):
    t, i = col
    if t == 'l': return F(-ub[i][0].get(j, 0))
    return F(eq[i][0].get(j, 0))
Mrows = [[entry(col, j) for col in cols] + [cobj.get(j, F(0))] for j in P]
dm = DomainMatrix([[QQ(v.numerator, v.denominator) for v in row] for row in Mrows], (len(P), len(cols) + 1), QQ)
t0 = time.time()
R, piv = dm.to_sparse().rref() if hasattr(dm, 'to_sparse') else dm.rref()
R = R.to_Matrix(); ncol = len(cols)
print("rref done", f"({time.time()-t0:.1f}s)"); sys.stdout.flush()
assert ncol not in piv, "complementary-slackness system inconsistent"
free = [k for k in range(ncol) if k not in piv]
def fl(col): return lam[col[1]] if col[0] == 'l' else mu[col[1]]
y = [None] * ncol
for k in free: y[k] = F(fl(cols[k])).limit_denominator(10 ** 6)
for r_, pk in enumerate(piv):
    val = F(int(R[r_, ncol].p), int(R[r_, ncol].q))
    for k in free:
        a = R[r_, k]
        if a != 0: val -= F(int(a.p), int(a.q)) * y[k]
    y[pk] = val
ok = all(y[k] >= 0 for k, col in enumerate(cols) if col[0] == 'l')
tot = {}
for k, col in enumerate(cols):
    row = {kk: -vv for kk, vv in ub[col[1]][0].items()} if col[0] == 'l' else eq[col[1]][0]
    for j, v in row.items(): tot[j] = tot.get(j, F(0)) + y[k] * F(v)
bad = [j for j in range(1, NV) if tot.get(j, F(0)) < cobj.get(j, F(0))]
ok = ok and not bad
bound = sum(y[len(S) + e] * eq[e][1] for e in range(len(eq)))
print("EXACT CERTIFICATE VALID:", ok, "bound: mean delta_%d <=" % m, bound, "=", float(bound), "x mean delta_1"); sys.stdout.flush()
if not ok:
    print("negative multipliers:", sum(1 for k, col in enumerate(cols) if col[0] == 'l' and y[k] < 0), "bad vars:", len(bad)); sys.exit(1)
used = sorted([(y[k], ub[cols[k][1]][1][1]) for k in range(len(S)) if y[k] > 0 and ub[cols[k][1]][1][0] == 'cut'], reverse=True)
print("characteristic / Ingleton rows used (exact multipliers):")
for v, nm in used: print("  ", v, float(v), nm)
dump = {"m": m, "uniform": UNIFORM, "bound": str(bound), "rows": [], "eqs": []}
for k in range(len(S)):
    if y[k] > 0:
        wit = ub[cols[k][1]][1]
        if wit[0] == 'cut':
            nm = wit[1]
            mo = re.match(r"(dfz|ing):L=\[([0-9, ]*)\],T=\[([0-9, ]*)\]$", nm)
            if mo:
                L = [[int(v)] for v in mo.group(2).split(',')]; T = [int(v) for v in mo.group(3).split(',') if v.strip()]
                dump["rows"].append({"type": mo.group(1), "letters": L, "T": T, "mult": str(y[k])})
            else:   # fixed instances of cutsearch.py
                e = [1 << b for b in range(m)]
                from lpm import fano_letters
                if nm.startswith("dfz:points@bottom"):
                    L = [list(t) for t in fano_letters(I, 0, e[:3])]; T = [] if nm.endswith("none") else [0]
                else:
                    x0 = I.IDX[frozenset({0, e[0]})]
                    L = [list(t) for t in fano_letters(I, x0, e[1:4])]; T = [] if nm.endswith("none") else [x0]
                dump["rows"].append({"type": "dfz", "letters": L, "T": T, "mult": str(y[k]), "label": nm})
        else:
            dump["rows"].append({"type": wit[0], "witness": wit[1:], "mult": str(y[k])})
for e in range(len(eq)):
    if y[len(S) + e] != 0: dump["eqs"].append({"def": eq[e][2], "mult": str(y[len(S) + e])})
out = sys.argv[1].replace("cuts_", "certificate_")
json.dump(dump, open(out, "w"))
print("certificate written:", out, len(dump["rows"]), "rows,", len(dump["eqs"]), "equalities")
