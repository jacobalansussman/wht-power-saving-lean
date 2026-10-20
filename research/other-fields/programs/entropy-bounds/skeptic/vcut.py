# vcut.py -- v-nc-gain's own cut LPs for a fixed layout of sessions on an arbitrary graph (independent code).
#  cut_flow(n, edges, pairs, directed): the cut LP in compact flow form, one flow per terminal pattern under shared
#    capacities, solved by HiGHS (scipy), with a RIGOROUS Lagrangian lower bound computed here:
#      every variable boxed in [0, U], U = W (some optimal solution respects it: capacities can be capped at the
#      largest requirement <= W, and each pattern's flow taken acyclic of value = requirement <= W);
#      LB = y_ub.b_ub + y_eq.b_eq + sum_j min(0, red_j) * U  for any y_ub <= 0, y_eq.
#    undirected: c_e per edge, f_uv + f_vu <= c_e, requirement sep(Q) = #pairs split;
#    directed:   c_a per arc,  f_a <= c_a,         requirement #{(s,t): t in Q, s not in Q}.
#  cut_all(n, edges, pairs): undirected cut LP with EVERY cut listed (small graphs only), dual rationalised and checked
#    in exact fractions (a cut packing with every edge load <= 1).
import numpy as np, itertools
import scipy.sparse as sp
from scipy.optimize import linprog
from fractions import Fraction

def patterns(pairs, directed):
    X = sorted({v for p in pairs for v in p})
    out = []
    if directed:
        for mask in range(1, (1 << len(X)) - 1):
            P = {X[i] for i in range(len(X)) if mask >> i & 1}
            need = sum(1 for (s, t) in pairs if t in P and s not in P)
            if need > 0: out.append(([x for x in X if x not in P], sorted(P), need))
    else:
        for mask in range(0, 1 << (len(X) - 1)):
            A = {X[-1]} | {X[i] for i in range(len(X) - 1) if mask >> i & 1}
            if len(A) == len(X): continue
            need = sum(1 for (s, t) in pairs if (s in A) != (t in A))
            if need > 0: out.append((sorted(A), [x for x in X if x not in A], need))
    return out

def cut_flow(n, edges, pairs, directed=False, method='highs-ds', U=None, model_only=False):
    W = len(pairs); U = float(W if U is None else U)
    E = len(edges); eu = np.array([u for u, v in edges]); ev = np.array([v for u, v in edges])
    tail = np.concatenate([eu, ev]); head = np.concatenate([ev, eu]); A2 = 2 * E
    ncap = A2 if directed else E
    pats = patterns(pairs, directed)
    Iu, Ju, Vu, bu = [], [], [], []      # A_ub rows
    Ie, Je, Ve, be = [], [], [], []      # A_eq rows
    nv = ncap; ru = 0; re_ = 0
    for (src, snk, need) in pats:
        f0 = nv; g0 = f0 + A2; h0 = g0 + len(src); nv = h0 + len(snk)
        ar = np.arange(A2)
        if directed:          # f_a - c_a <= 0
            Iu += [ru + ar, ru + ar]; Ju += [f0 + ar, ar]; Vu += [np.ones(A2), -np.ones(A2)]; bu.append(np.zeros(A2)); ru += A2
        else:                 # f_uv + f_vu - c_e <= 0
            er = np.arange(E)
            Iu += [ru + er, ru + er, ru + er]; Ju += [f0 + er, f0 + E + er, er]; Vu += [np.ones(E), np.ones(E), -np.ones(E)]
            bu.append(np.zeros(E)); ru += E
        # conservation: inflow - outflow + g - h = 0 at every node
        Ie += [re_ + head, re_ + tail, re_ + np.array(src), re_ + np.array(snk)]
        Je += [f0 + ar, f0 + ar, g0 + np.arange(len(src)), h0 + np.arange(len(snk))]
        Ve += [np.ones(A2), -np.ones(A2), np.ones(len(src)), -np.ones(len(snk))]
        be.append(np.zeros(n)); re_ += n
        # value: -sum g <= -need
        Iu.append(np.full(len(src), ru)); Ju.append(g0 + np.arange(len(src))); Vu.append(-np.ones(len(src)))
        bu.append(np.array([-float(need)])); ru += 1
    Aub = sp.csr_matrix((np.concatenate(Vu), (np.concatenate(Iu), np.concatenate(Ju))), shape=(ru, nv))
    Aeq = sp.csr_matrix((np.concatenate(Ve), (np.concatenate(Ie), np.concatenate(Je))), shape=(re_, nv))
    bub = np.concatenate(bu); beq = np.concatenate(be)
    c = np.zeros(nv); c[:ncap] = 1.0
    if model_only: return c, Aub, bub, Aeq, beq, U, len(pats)
    res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=(0, U), method=method)
    if res.status != 0: return None, None, res.message
    yu = np.minimum(np.array(res.ineqlin.marginals), 0.0); ye = np.array(res.eqlin.marginals)
    red = c - Aub.T @ yu - Aeq.T @ ye
    LB = float(yu @ bub + ye @ beq + np.minimum(red, 0).sum() * U)
    return float(res.fun), LB, len(pats)

def all_cuts_rows(n, edges, pairs):
    eu = np.array([u for u, v in edges]); ev = np.array([v for u, v in edges])
    rows = []; rhs = []
    for mask in range(1, 1 << (n - 1)):
        q = np.array([(mask >> i) & 1 for i in range(n)], bool)
        sep = sum(1 for s, t in pairs if q[s] != q[t])
        if sep == 0: continue
        rows.append(q[eu] != q[ev]); rhs.append(sep)
    return np.array(rows, float), np.array(rhs, float)

def cut_all(n, edges, pairs, rows_rhs=None):
    R, rhs = rows_rhs if rows_rhs is not None else all_cuts_rows(n, edges, pairs)
    if rows_rhs is not None:   # recompute rhs for these pairs on the stored cut list
        pass
    E = len(edges)
    res = linprog(np.ones(E), A_ub=-sp.csr_matrix(R), b_ub=-rhs, bounds=(0, None), method='highs')
    lam = -np.array(res.ineqlin.marginals)
    lamF = [Fraction(max(0.0, float(x))).limit_denominator(10 ** 4) for x in lam]
    load = [Fraction(0)] * E
    for i, l in enumerate(lamF):
        if l:
            for e in np.nonzero(R[i] > 0.5)[0]: load[e] += l
    mx = max(load); sc = 1 if mx <= 1 else 1 / mx
    val = sum(l * int(rhs[i]) for i, l in enumerate(lamF)) * sc
    return float(res.fun), val
