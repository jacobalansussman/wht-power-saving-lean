# capacity.py IRR[,IRR..] -- the module game of a2-gct's relaxation (relax2.load) read as a subspace game, per P_0-type tau
# with multiplicity m: (i) depth: is k in S + S h S (the depth-1 space)? (ii) capacity: need_tau = how much the target's
# tau-rows exceed the inputs' tau-rows (per array), cap_tau = the most one paid move can add to the tau row space (the rank
# of h's rows into tau, from the other coordinates, in the tau-view); bound M >= W * need / cap for every tau.
import sys, numpy as np, relax2
np.set_printoptions(linewidth=200)
js = [int(a) for a in sys.argv[1].split(',')]
n, h, k, tl = relax2.load(js, np.random.default_rng(12345))
def rk(X, tol=1e-8):
    s = np.linalg.svd(X, compute_uv=False); return int((s > tol * max(1, s[0])).sum()), s
print('irreps %s: n = %d, types (dim, multiplicity):' % (js, n), [(d, len(sl)) for d, sl in tl])
if all(len(sl) == 1 for d, sl in tl):
    # multiplicity-free: S h S contains Hom(tau', tau) iff the block h[tau, tau'] is nonzero, so S + ShS is the sum of
    # End(tau) and those Hom blocks (no span computation needed)
    nz = [[np.linalg.norm(h[np.ix_(s1[0], s2[0])]) for (_, s2) in tl] for (_, s1) in tl]
    print('multiplicity-free; norms of the blocks h[tau, tau\']:', np.round(nz, 3).tolist())
    full = all(nz[a][b] > 1e-8 for a in range(len(tl)) for b in range(len(tl)) if a != b)
    print('S + ShS = End(V): %s  (so k has depth 1: %s)' % (full, full))
elif len(sys.argv) > 2 and sys.argv[2] == 'span':
    # S basis: for each type, matrix units repeated diagonally on its slots
    Sb = []
    for d, slots in tl:
        for a in range(d):
            for b in range(d):
                M = np.zeros((n, n))
                for s_ in slots: M[s_[a], s_[b]] = 1
                Sb.append(M)
    dimS = len(Sb)
    F1 = [M.ravel() for M in Sb]
    # S h S spanned by Ea h Eb over matrix units (it is a bimodule: products of basis elements span it)
    for A_ in Sb:
        Ah = A_ @ h
        if not np.any(Ah): continue
        for B_ in Sb:
            P_ = Ah @ B_
            if np.any(P_): F1.append(P_.ravel())
    F1 = np.array(F1); r1, _ = rk(F1)
    r1k, _ = rk(np.vstack([F1, k.ravel()[None]]))
    print('dim S = %d, dim (S + ShS) = %d of %d; k in S + ShS: %s' % (dimS, r1, n * n, r1k == r1))

else: print('not multiplicity-free: depth test skipped (pass span to run it)')
# capacity per type (W = 1 view; the tau-view of a source x is the d x (m n) matrix of its tau-rows over all slots)
for t, (d, slots) in enumerate(tl):
    idx = np.concatenate(slots); other = np.setdiff1d(np.arange(n), idx)
    # inputs' tau-rows: e restricted to tau rows = coordinates idx (rows); rowspace = coordinate space on idx
    # target tau-rows: k[idx, :]; need = rank of k[idx, other]  (the part outside the inputs' rows)
    need, sk = rk(k[np.ix_(idx, other)])
    # a move adds rows (h c)[idx] = h[idx, idx] c[idx] + h[idx, other] c[other]; the first part lies in the old row space
    # only if the mixing is free on the tau rows -- with multiplicity m the mixing acts as X (x) 1, so c[idx] rows are
    # constrained; we report rank h[idx, other] as the (upper bound on the) gain from other types, per move
    cap, sh = rk(h[np.ix_(idx, other)])
    gap = lambda s, r: '%.1e | %.1e' % (s[r - 1] if r else 0, s[r] if r < len(s) else 0)
    print('  type %d: dim %d x mult %d: need %d per array (rows of k outside the inputs), cap <= %d per move -> M/W >= %.3f   [sv gaps: need %s; cap %s; generic need would be %d]'
          % (t, d, len(slots), need, cap, need / cap if cap else float('inf'), gap(sk, need), gap(sh, cap), min(d * len(slots), len(other))), flush=True)
