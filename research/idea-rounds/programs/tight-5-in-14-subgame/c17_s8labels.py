# c17_s8labels.py -- labelled certificates in the singular subgame across PAIRS of S8 irreducibles.
# For every S8 irreducible, its AGL(3,2) isotypic components with the AGL(3,2) irreducible identified by its character.
# For a pair (lam, mu) (lam = mu allowed: two copies), a certificate is Y = sum_tau U_tau (x) L_tau with L_tau a subspace
# of the multiplicity space of tau in V_lam + V_mu. Search: every dimension pattern (d_tau), several random subspaces
# each (generic position), plus coordinate (cut-like) choices. Ratio f(K)/f(h). Calibration: the hook pattern on
# (2,1^6)+(3,1^5)+(2,2,2,1,1) is a triple; check pairs reproduce <= 2.8 and report the maximum.
import numpy as np, itertools, time
from lib514 import *
from libs8 import *
import c8_s4 as C8
t0 = time.time(); rng = np.random.default_rng(17)
PBg = [to_perm(g) for g in C8.Bp_gens]
ph = to_perm(h1); pK = to_perm(K)
Bperms = None
def bprime_elements():
    el = {tuple(range(8)): None}; fr = [tuple(range(8))]
    while fr:
        nx = []
        for k in fr:
            for p in PBg:
                q = tuple(int(x) for x in compose(p, np.array(k)))
                if q not in el: el[q] = None; nx.append(q)
        fr = nx
    return [np.array(k) for k in el]
Bperms = bprime_elements(); assert len(Bperms) == 1344
tests = [Bperms[j] for j in rng.choice(1344, 30, replace=False)]
info = {}
for lam in partitions(8):
    Yr = YoungRep(lam)
    mats = [Yr(p) for p in Bperms]
    S = rng.standard_normal((Yr.dim, Yr.dim)); S = S + S.T
    Rm = sum(M @ S @ M.T for M in mats) / len(mats)
    w, V = np.linalg.eigh((Rm + Rm.T) / 2)
    pieces = []; st = 0
    for i in range(1, Yr.dim + 1):
        if i == Yr.dim or w[i] - w[i - 1] > 1e-7 * max(1, np.abs(w).max()): pieces.append(V[:, st:i]); st = i
    # group pieces by character; within a type, align bases by intertwiners so that labels make sense
    types = {}
    for Q in pieces:
        chi = tuple(np.round([np.trace(Q.T @ Yr(t) @ Q) for t in tests], 4))
        types.setdefault(chi, []).append(Q)
    info[lam] = (Yr, types)
# global type list
alltypes = {}
for lam, (Yr, types) in info.items():
    for chi, Qs in types.items(): alltypes.setdefault(chi, []).append((lam, len(Qs)))
print('AGL(3,2) types (dim: list of (S8 irreducible, multiplicity)):')
for chi, lst in alltypes.items(): print('  dim', int(round(chi[0] if False else 0)) or '', [(l, m) for l, m in lst])
def align(Q1, Q2, Yr1, Yr2):
    """B'-isometry: columns of Q2 X represent the same basis as Q1"""
    d = Q1.shape[1]
    eqs = [np.kron(np.eye(d), (Q1.T @ Yr1(p) @ Q1).T) - np.kron(Q2.T @ Yr2(p) @ Q2, np.eye(d)) for p in PBg]
    _, sv, vt = np.linalg.svd(np.vstack(eqs)); X = vt[-1].reshape(d, d); X = X / np.sqrt(np.trace(X @ X.T) / d)
    assert sv[-1] < 1e-7, sv[-1]
    return Q2 @ X
def pair_search(lam, mu, ntrials=6):
    Yl, tl = info[lam]; Ym, tm = info[mu]
    dl, dm = Yl.dim, Ym.dim; n = dl + dm
    Hh = np.zeros((n, n)); Hh[:dl, :dl] = Yl(ph); Hh[dl:, dl:] = Ym(ph)
    KK = np.zeros((n, n)); KK[:dl, :dl] = Yl(pK); KK[dl:, dl:] = Ym(pK)
    # per type: list of aligned embeddings (n x d) for each copy
    copies = {}
    for chi, Qs in tl.items():
        for Q in Qs:
            E = np.zeros((n, Q.shape[1])); E[:dl] = Q; copies.setdefault(chi, []).append((E, Yl, 'l', Q))
    for chi, Qs in tm.items():
        for Q in Qs:
            E = np.zeros((n, Q.shape[1])); E[dl:] = Q; copies.setdefault(chi, []).append((E, Ym, 'm', Q))
    # align all copies of a type to the first
    aligned = {}
    for chi, cl in copies.items():
        ref = cl[0]; out = [ref[0]]
        for (E, Yr, side, Q) in cl[1:]:
            Qa = align(ref[3], Q, ref[1], Yr)
            Ea = np.zeros_like(E)
            if side == 'l': Ea[:dl] = Qa
            else: Ea[dl:] = Qa
            out.append(Ea)
        aligned[chi] = out
    shared = [chi for chi, cl in aligned.items() if len(cl) >= 2]
    best = (0, None)
    chis = list(aligned)
    mults = [len(aligned[c]) for c in chis]
    for dims in itertools.product(*[range(m + 1) for m in mults]):
        for trial in range(ntrials if any(0 < d < m for d, m in zip(dims, mults)) else 1):
            cols = []
            for c, d, m in zip(chis, dims, mults):
                if d == 0: continue
                Lsub = np.linalg.qr(rng.standard_normal((m, d)))[0]   # d-dim subspace of C^m (real generic)
                for j in range(d):
                    cols.append(sum(Lsub[i, j] * aligned[c][i] for i in range(m)))
            if not cols: continue
            Yc = np.hstack(cols)
            if Yc.shape[1] in (0, n): continue
            cap = rank(np.hstack([Yc, Hh @ Yc])) - Yc.shape[1]
            if cap == 0: continue
            dem = rank(np.hstack([Yc, KK @ Yc])) - Yc.shape[1]
            r = dem / cap
            if r > best[0] + 1e-9: best = (r, (dims, cap, dem))
    return best, len(shared)
if __name__ == '__main__':
    lams = [l for l in partitions(8) if info[l][0].dim > 1]
    results = []
    for i, lam in enumerate(lams):
        for mu in lams[i:]:
            # only pairs that share at least one type (or lam = mu)
            sh = set(info[lam][1]) & set(info[mu][1])
            if not sh: continue
            best, nsh = pair_search(lam, mu)
            results.append((best[0], lam, mu, best[1]))
    results.sort(key=lambda t: -t[0])
    print('best generic-label pair certificates (ratio, lam, mu, (dims, cap, demand)):')
    for r in results[:12]: print('   %.4f' % r[0], r[1], r[2], r[3])
    print('pairs examined:', len(results), '[%.0fs]' % (time.time() - t0))
