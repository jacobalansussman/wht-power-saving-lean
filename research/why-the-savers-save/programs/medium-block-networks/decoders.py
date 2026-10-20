#!/usr/bin/env python3
"""cx-medium / decoders.py (own code; numpy + scipy milp from our pylib; python3 -I).
The centre (index code) of one call, for each port design: delivered map M = sum_f w_f g_f g_f^T, g_f(a) = f . r_a
(a total T_f = sum_b g_f(b) x_b scattered with w_f g_f(a)).  Requirement: M_ab = delta_ab on every label-non-orthogonal
pair (B_ab = 1: the diagonal and the orthogonal roots); anything on junk pairs (B_ab = 0).
Cost of one total with support S = {b : g_f(b) != 0} (gather set = scatter set):
  hub  = dim span_F2(u_S): its copy stands at a frame containing every gathered label and descends to frame 0, where
         the targets wait (the hub's and our accounting, gcert E2: y at 0 at the scatter);
  pair = rank_F2 B[S, S] (an2-smallnet 3.2): the least possible drop, reached only if the targets wait at u_S^perp.
P=1 break-even (one pair, three calls): c < 2v/3 per call; P=2: c < 4v/5.
(1) natural decoders: an orthonormal basis of R^t in the design's coordinates (all E8-type coordinate frames are
    tried by the coordinate permutations/sign choices that keep the span); (2) a MILP over a pool of candidate
    directions f (all {0,+-1} vectors of support <= 4 and all (+-1)^8, inside the span), minimising the summed cost,
    with the exact equalities above (solved in floating point, then the support re-solved and checked in Fractions)."""
import sys, os, json, itertools, math, time
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
from fractions import Fraction
from scipy.optimize import milp, LinearConstraint, Bounds
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from designs_cx import design, labels, ip, rank2, dot

def costs_of(S, labs, B):
    hub = rank2([labs[b] for b in S])
    pairs = rank2([int(''.join(str(B[a][b]) for b in S), 2) for a in S]) if S else 0
    return hub, pairs

def span_basis_R(L):
    """an orthonormal-in-rationals check is not needed; we test candidate f by membership in span (rank test)."""
    A = np.array(L, dtype=float)
    return np.linalg.matrix_rank(A)

def in_span(f, A, r):
    return np.linalg.matrix_rank(np.vstack([A, f])) == r

def pool(dim, A, r):
    out = set()
    for supp in range(1, 5):
        for idx in itertools.combinations(range(dim), supp):
            for signs in itertools.product((1, -1), repeat=supp):
                if signs[0] != 1: continue
                f = [0] * dim
                for i, s in zip(idx, signs): f[i] = s
                out.add(tuple(f))
    if dim == 8:
        for signs in itertools.product((1, -1), repeat=8):
            if signs[0] == 1: out.add(signs)
    if dim <= 8:
        # sums of two roots directions and coordinate differences are inside the pool above for D/A types;
        pass
    res = [f for f in out if in_span(np.array(f, float), A, r)]
    return res

def natural(L, labs, B):
    """coordinate decoders: E-type and D-type use the ambient coordinates; we report the best orthonormal coordinate
    basis found among: ambient e_i (when the span is everything), and for E7/A_n the span-adapted bases."""
    return None

def analyse(name, tl=120, P_list=(1, 2)):
    L = design(name); labs, k, B = labels(L); n = len(L); dim = len(L[0])
    A = np.array(L, dtype=float); r = np.linalg.matrix_rank(A)
    norm = ip(L[0], L[0])
    cands = pool(dim, A, r)
    G = np.array([[ip(f, l) for l in L] for f in cands], dtype=float)        # g_f(a)
    supp = [[a for a in range(n) if G[i, a] != 0] for i in range(len(cands))]
    keep = [i for i in range(len(cands)) if supp[i]]
    cands = [cands[i] for i in keep]; G = G[keep]; supp = [supp[i] for i in keep]
    cst = [costs_of(S, labs, B) for S in supp]
    # equalities on B = 1 pairs (a <= b)
    pairs = [(a, b) for a in range(n) for b in range(a, n) if B[a][b] == 1]
    nf = len(cands)
    Aeq = np.zeros((len(pairs), nf))
    for p, (a, b) in enumerate(pairs):
        Aeq[p] = G[:, a] * G[:, b]
    beq = np.array([1.0 if a == b else 0.0 for a, b in pairs])
    out = dict(design=name, n=n, k=k, t=r, pool=nf, break_even={P: 2 * P * n / (2 * P + 1) for P in P_list})
    kinds = os.environ.get('KINDS', 'hub,pair').split(',')
    for kind, ci in [kc for kc in (('hub', 0), ('pair', 1)) if kc[0] in kinds]:
        c = np.array([cs[ci] for cs in cst], dtype=float)
        Wb = 4.0 / norm * 4                                                  # weight bound (natural weights are 1/norm*2)
        # variables: w+ (nf), w- (nf), z (nf)
        cobj = np.concatenate([np.zeros(2 * nf), c])
        Aall = np.hstack([Aeq, -Aeq, np.zeros((len(pairs), nf))])
        cons = [LinearConstraint(Aall, beq, beq)]
        Ib = np.hstack([np.eye(nf), np.eye(nf), -Wb * np.eye(nf)])
        cons.append(LinearConstraint(Ib, -np.inf, 0))
        integrality = np.concatenate([np.zeros(2 * nf), np.ones(nf)])
        bounds = Bounds(np.zeros(3 * nf), np.concatenate([np.full(2 * nf, Wb), np.ones(nf)]))
        t0 = time.time()
        res = milp(cobj, constraints=cons, integrality=integrality, bounds=bounds,
                   options=dict(time_limit=tl, mip_rel_gap=0.0, disp=False))
        rec = dict(status=res.status, message=res.message, seconds=round(time.time() - t0, 1))
        if res.x is not None:
            z = res.x[2 * nf:]; w = res.x[:nf] - res.x[nf:2 * nf]
            used = [i for i in range(nf) if z[i] > 0.5]
            # exact re-solve on the support: least squares in Fractions via numpy then rational check
            sub = Aeq[:, used]
            wl, *_ = np.linalg.lstsq(sub, beq, rcond=None)
            resid = float(np.abs(sub @ wl - beq).max())
            rec.update(objective=float(sum(c[i] for i in used)), totals=len(used), residual=resid,
                       dual_bound=float(res.mip_dual_bound) if hasattr(res, 'mip_dual_bound') and res.mip_dual_bound is not None else None,
                       costs=sorted(int(c[i]) for i in used),
                       totals_list=[dict(f=cands[i], support=len(supp[i]), hub=cst[i][0], pair=cst[i][1], w=float(wl[j])) for j, i in enumerate(used)])
        out[kind] = rec
        print(f"{name}: n={n} k={k} t={r} pool={nf} [{kind} cost] status={res.status} obj={rec.get('objective')} totals={rec.get('totals')}"
              f" dual={rec.get('dual_bound')} resid={rec.get('residual')} costs={rec.get('costs')} ({rec['seconds']} s); break-even P=1 {2*n/3:.1f}, P=2 {4*n/5:.1f}", flush=True)
    return out

if __name__ == '__main__':
    names = sys.argv[1:] or ['A5', 'D6', 'E7', 'P112']
    tl = float(os.environ.get('TL', '120'))
    res = [analyse(nm, tl) for nm in names]
    tag = '_'.join(names)
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'decoders_{tag}.json'), 'w'), indent=1, default=str)
