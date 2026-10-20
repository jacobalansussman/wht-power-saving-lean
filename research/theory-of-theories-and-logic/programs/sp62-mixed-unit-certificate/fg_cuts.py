# fg_cuts.py -- for each irreducible of Sp(6,2) (from fg_build.py): sanity checks, P_0-pieces (via A = Sym(3,2) and the
# Levi GL(3,2), never a sum over all 10,752 elements), Theorem S's spectral ceiling, and every cut (sum of P_0-isotypic
# components): profile (psi1, psi2, psi3) = ranks of (1-Pi_S) g Pi_S for g = h, h2, k; Pareto front; mixed criteria.
# Usage: fg_cuts.py irreps.npz out.json [only_j,...]
import sys, json, itertools, time
import numpy as np
from fractions import Fraction as Fr
D = np.load(sys.argv[1]); TABLE = [tuple(int(v) for v in r) for r in D['table']]
ct = json.load(open('gap/ctbl.json'))
cls = ct['classes']; iX = {c: cls.index(c) for c in ['2a', '2c', '2d']}
rng = np.random.default_rng(12345)
PRUNE, PRUNE_DIM = 1.95, 280     # prune (only for irreducibles of dimension >= 280) cuts whose Theorem-S bound is <= 88/45 - eps
only = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else None
def rank(M, tol=1e-7):
    if M.size == 0: return 0, None, None
    s = np.linalg.svd(M, compute_uv=False)
    r = int((s > tol).sum()); lo = float(s[r - 1]) if r else None; hi = float(s[r]) if r < len(s) else None
    return r, lo, hi
# GL(3,2) as words in the two Levi generators (F_2 level), for the averaging
def levi_f2():
    g1 = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]]); g2 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    seen = {np.eye(3, dtype=int).tobytes(): ()}; fr = [(np.eye(3, dtype=int), ())]
    while fr:
        nx = []
        for M, w in fr:
            for gi, g in enumerate((g1, g2)):
                N = (g @ M) % 2; kk = N.tobytes()
                if kk not in seen: seen[kk] = w + (gi,); nx.append((N, w + (gi,)))
        fr = nx
    return list(seen.values())
WORDS = levi_f2(); assert len(WORDS) == 168
results = {}
for j in (only if only else [int(x) for x in D['found']]):
    t0 = time.time()
    d, va, vb = TABLE[j]
    A, L, H, G = D[f'irr{j}_A'], D[f'irr{j}_L'], D[f'irr{j}_H'], D[f'irr{j}_G']
    # sanity: orthogonal, characters of h, h2, k
    err = max(np.abs(M @ M.T - np.eye(d)).max() for M in list(A) + list(L) + list(H) + list(G))
    chis = [float(np.trace(H[i])) for i in range(3)]
    expect = [ct['irr'][j][iX[c]] for c in ['2a', '2c', '2d']]
    assert err < 1e-8 and np.allclose(chis, expect, atol=1e-6), (j, err, chis, expect)
    # elements of A (64: products of 6 commuting generators) and L (168 words)
    Ael = []
    for bits in range(64):
        M = np.eye(d)
        for i in range(6):
            if bits >> i & 1: M = M @ A[i]
        Ael.append(M)
    Lel = {(): np.eye(d)}
    for w in sorted(WORDS, key=len):
        if w: Lel[w] = L[w[-1]] @ Lel[w[:-1]]
    Lel = list(Lel.values())
    S = rng.standard_normal((d, d)); S = S + S.T
    SL = sum(M @ S @ M.T for M in Lel) / 168
    R = sum(M @ SL @ M.T for M in Ael) / 64
    # check R commutes with P_0 generators
    cerr = max(np.abs(M @ R - R @ M).max() for M in list(A) + list(L))
    w, V = np.linalg.eigh(R)
    cl = []; st = 0
    for i in range(1, d + 1):
        if i == d or w[i] - w[i - 1] > 1e-7 * max(1, abs(w).max()): cl.append((st, i)); st = i
    # signature of each copy: dim and traces on test elements of P_0
    tests = [A[0], A[3], A[0] @ A[1], A[0] @ A[1] @ A[2], A[3] @ A[5], L[0], L[1], L[0] @ L[1], L[1] @ L[1] @ L[0],
             A[0] @ L[1], A[3] @ L[0], A[2] @ A[4] @ L[1] @ L[0]]
    copies = []
    for (a, b) in cl:
        Q = V[:, a:b]
        sig = (b - a,) + tuple(round(float(np.trace(Q.T @ T @ Q)), 5) + 0.0 for T in tests)
        copies.append((sig, Q))
    types = {}
    for sig, Q in copies: types.setdefault(sig, []).append(Q)
    pieces = [(sig, np.concatenate(Qs, axis=1), len(Qs)) for sig, Qs in types.items()]
    pieces.sort(key=lambda p: (p[0][0], p[0]))
    dims = sorted(p[0][0] for p in pieces for _ in range(p[2]))
    expdims = sorted(ct['P0_degrees'][t] for t, mm in enumerate(ct['restriction'][j]) for _ in range(mm))
    ok_dims = dims == expdims
    mults = [p[2] for p in pieces]; mfree = max(mults) == 1
    r = len(pieces)
    # Theorem S on isotypic pieces (for multiplicity-free modules this is pf-duality's statement)
    Hm = H[0]
    Pi = [p[1] @ p[1].T for p in pieces]; dd = np.array([p[1].shape[1] for p in pieces], float)
    Mw = np.array([[np.trace(Pi[a] @ Hm @ Pi[b] @ Hm) / dd[a] for b in range(r)] for a in range(r)])
    Ms = np.diag(np.sqrt(dd)) @ Mw @ np.diag(1 / np.sqrt(dd)); ev = np.sort(np.linalg.eigvalsh((Ms + Ms.T) / 2))[::-1]
    lam2 = float(ev[1]) if r > 1 else None
    ceil_mod = 2 / (1 - lam2) if r > 1 else 0.0
    # cuts
    Bk = [[[p1[1].T @ H[g] @ p2[1] for p2 in pieces] for p1 in pieces] for g in range(3)]
    cuts = []; gaps = [np.inf, 0.0]; skipped = 0
    for bits in range(1, 1 << r):
        if bits >> (r - 1) & 1: continue            # complement symmetry: f_N(g) = f_{N^perp}(g) (g orthogonal)
        Sset = [i for i in range(r) if bits >> i & 1]; Sc = [i for i in range(r) if not bits >> i & 1]
        dimN = int(sum(dd[i] for i in Sset)); p = dimN / d
        tS = 1 / ((1 - lam2) * max(p, 1 - p)) if r > 1 else 0
        if d >= PRUNE_DIM and tS <= PRUNE:           # Theorem S: psi_j/psi_1 <= tS for j = 2, 3, so this cut cannot reach
            skipped += 1; continue                    # 2.8, nor psi3+4psi2 > 10psi1 (needs tS > 2), nor 7psi3+38psi2 > 88psi1 (tS > 88/45)
        psi = []
        for g in range(3):
            M = np.block([[Bk[g][a][b] for b in Sset] for a in Sc])
            rk, lo, hi = rank(M); psi.append(rk)
            if lo is not None: gaps[0] = min(gaps[0], lo)
            if hi is not None: gaps[1] = max(gaps[1], hi)
        cuts.append(dict(S=Sset, dimN=dimN, psi=psi, thmS=tS))
    best = max((c for c in cuts if c['psi'][0] > 0), key=lambda c: c['psi'][2] / c['psi'][0], default=None)
    mix1 = max((c['psi'][2] + 4 * c['psi'][1] - 10 * c['psi'][0] for c in cuts), default=None)
    mix2 = max((7 * c['psi'][2] + 38 * c['psi'][1] - 88 * c['psi'][0] for c in cuts), default=None)
    pareto = sorted({(Fr(c['psi'][1], c['psi'][0]), Fr(c['psi'][2], c['psi'][0])) for c in cuts if c['psi'][0] > 0})
    front = [pt for pt in pareto if not any(q[0] >= pt[0] and q[1] >= pt[1] and q != pt for q in pareto)]
    res = dict(deg=d, chi2a=va, mult_free=mfree, piece_dims=[int(p[0][0]) for p in pieces], piece_mults=mults,
               dims_match_table=ok_dims, comm_err=float(cerr), orth_err=float(err), lam2=lam2, thmS_module=ceil_mod,
               ncuts=len(cuts), pruned_by_thmS=skipped, best_ratio=(str(Fr(best['psi'][2], best['psi'][0])) if best else None),
               best_cut=best, max_psi3_plus_4psi2_minus_10psi1=mix1, max_7psi3_38psi2_minus_88psi1=mix2,
               pareto_front=[(str(a), str(b), float(b)) for a, b in front], sv_gap=[float(gaps[0]), float(gaps[1])],
               cuts=cuts, secs=round(time.time() - t0, 1))
    results[j] = res
    print(f"[{j:2d}] {d:3d}({va:+d}) {'MF ' if mfree else 'REP'} pieces {res['piece_dims']} mult {mults} ok={ok_dims} lam2={lam2 if lam2 is None else round(lam2,4)} "
          f"thmS={ceil_mod:.3f} cuts={len(cuts)}+{skipped}pruned best={res['best_ratio']} {best['psi'] if best else ''} mix1={mix1} mix2={mix2} "
          f"front={[(a, b) for a, b, _ in res['pareto_front']]} svgap={gaps[0]:.2e}/{gaps[1]:.1e} ({res['secs']}s)", flush=True)
    json.dump(results, open(sys.argv[2], 'w'), default=lambda o: o if not isinstance(o, np.integer) else int(o))
