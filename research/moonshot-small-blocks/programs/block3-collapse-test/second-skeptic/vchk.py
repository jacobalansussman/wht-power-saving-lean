"""v-xw3-extract (skeptic): independent re-test of xw3-extract's classification of one batch-B/C flag.

Uses only the claimed collapse vectors (fu_<tag>.npz: keep, X, labels, red) and my own letters/paths (vlet.py).
 1. labels: my enumeration of top words must equal the npz labels in order.
 2. f = 3,4 stacked, FRESH inputs (seed 7000+): every claimed collapse vector must be in the complete lower span
    (residual ~1e-15), otherwise the extraction (or my reading of the conventions) is wrong.
    Also: column dependencies null(A34) of ALL top columns, and the trivial columns (each in the lower span).
 3. f = 5, fresh inputs, my own row subsample: candidate space C = span(X) + null(A34) + e_trivial (every
    coefficient vector that collapses at f = 3,4 with my inputs); report dim{c in C : B5 c in L5, B5 c != 0},
    the X-subspace spectrum, the trivial columns at f = 5, and the IRREDUCIBLE part: rank of B5[:, irr] c_irr over
    the f = 5 collapse space, with irr from MY reducibility test (dense pulled-back moves at f = 3, confirmed at f = 4)
    and, separately, from the npz 'red' array.  A class (c) collapse would show as a positive irreducible rank.
 4. f = 6 for the f = 5 collapse space (rows subsampled).
usage: vchk.py <tag> <npz> <moves> <R> <W1> ... <W_{n-1}> <B> <S> <rows5> <rows6>
"""
import sys, os, time, json, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlet import P, F, parse_moves, top_labels, path_vec, lower_list, A_, B_

T0 = time.time()
def log(*a): print(f"[{time.time()-T0:6.0f}s]", *a, flush=True)

def orth(A, tol=1e-9):
    if A.shape[1] == 0: return A[:, :0]
    U, s, _ = np.linalg.svd(A, full_matrices=False); r = int((s > s[0] * tol).sum()); return U[:, :r]

def proj(Q, Y): return Y - Q @ (Q.conj().T @ Y)

def null(A, tol=1e-9):
    _, s, Vh = np.linalg.svd(A, full_matrices=False)
    r = int((s > s[0] * tol).sum()); return Vh[r:].conj().T, s

MENU = {"R9": "1,X1,X2,X4,Z1,Z2,Z4,P01,P12", "K9": "1,T2,CS01,A,Ai,A.T2.Ai,Pr0,P01,X1", "K5": "1,T2,A,A.T2.Ai,Pr0",
        "H15": "1,P01,P12,C01,S1,X4,Z1,A,Ai,B,lX1,fS1,Pr0,T2,CS01",
        "H30": "1,P01,P12,C01,C12,C20,S1,S2,S4,X4,Z1,CZ01,A,Ai,a0,a0i,E01,B,Bi,lS1,lX1,fS1,Pr0,Pr2,T2,T0,CS01,CS12,A.T2.Ai,Ai.CS01.A"}
def M(s): return MENU.get(s, s).split(",")


def dfs_cols(Pp, mv, sub, mixl, fronts, col):
    """columns front . M_sub[-1] . mix ... M_sub[0] . mix . v in nested-product order (mixes outer, front inner), DFS"""
    out = []
    def rec(t, V):
        if t == len(sub):
            for g in fronts: out.append(col(Pp.word(g, V)))
            return
        i = sub[t]
        for k in mixl[t]:
            rec(t + 1, Pp.move(mv[i][0], mv[i][1], Pp.word(k, V)))
    rec(0, Pp.v)
    return out


def mats(fs, S, seed, mv, R, W, B, tl, low, rows=None):
    tops, lows = [], []
    n = len(mv)
    for k, f in enumerate(fs):
        Pp = P(f, S, seed + 17 * k)
        if rows is None: col = lambda V: V.reshape(-1).copy()
        else:
            rr = rows[k]; col = lambda V, rr=rr: V[rr, :].reshape(-1).copy()
        tops.append(np.stack(dfs_cols(Pp, mv, list(range(n)), [["1"]] + W, R, col), 1))
        lw = []
        for kk in range(n):
            for sub in itertools.combinations(range(n), kk):
                front = R if (kk > 0 and sub[-1] == n - 1) else B
                mixl = [["1"] if i == 0 else W[i - 1] for i in sub]
                lw += dfs_cols(Pp, mv, list(sub), mixl, front, col)
        lows.append(np.stack(lw, 1))
    return np.concatenate(tops, 0), np.concatenate(lows, 0)


# ---------------- my reducibility test (pulled-back moves applied to all basis vectors at f = 3, no matmul)
_PP = {}
def _pp(f):
    if f not in _PP: _PP[f] = P(f, 1, 0)
    return _PP[f]

def app_pulled(Pp, mv, lab, t, V):
    """M'_t V = P_t^-1 M_t P_t V, P_t = the mixes before move t in order of action"""
    pref = [w for w in lab[1:t + 1] if w not in ("", "1")]
    for w in pref: V = Pp.word(w, V)
    V = Pp.move(mv[t][0], mv[t][1], V)
    for w in reversed(pref): V = Pp.word(w, V, inv=True)
    return V

def is_mono(Mx, tol=1e-8): return bool(((np.abs(Mx) > tol).sum(0) == 1).all())

_RC = {}
def reducible_dense(f, mv, lab, n):
    key = (f, tuple(lab[:n]))
    if key in _RC: return _RC[key]
    Pp = _pp(f); I = np.eye(Pp.n, dtype=complex)
    rng = np.random.default_rng(1); Vr = rng.standard_normal((Pp.n, 3)) + 1j * rng.standard_normal((Pp.n, 3))
    res = None
    for j in range(n):
        for k in range(j + 1, n):
            if not is_mono(app_pulled(Pp, mv, lab, k, app_pulled(Pp, mv, lab, j, I))): continue
            mid = range(j + 1, k)
            def comm(p, l):
                return np.abs(app_pulled(Pp, mv, lab, p, app_pulled(Pp, mv, lab, l, Vr)) -
                              app_pulled(Pp, mv, lab, l, app_pulled(Pp, mv, lab, p, Vr))).max() < 1e-9
            if all(comm(k, l) for l in mid) or all(comm(j, l) for l in mid):
                res = (j, k); break
        if res: break
    _RC[key] = res
    return res


def main():
    a = sys.argv[1:]
    tag, npz = a[0], a[1]; mv = parse_moves(a[2]); n = len(mv)
    R = M(a[3]); W = [M(a[4 + j]) for j in range(n - 1)]; B = M(a[3 + n]); S = int(a[4 + n])
    rows5, rows6 = int(a[5 + n]), int(a[6 + n])
    d = np.load(npz); keep = d["keep"]; X = d["X"]; labs = [tuple(s.split("|")) for s in d["lab"]]; red_np = d["red"]
    rec = dict(tag=tag)
    tl = top_labels(n, R, W); low = lower_list(n, R, W, B)
    assert tl == labs, "label order differs"
    log(f"{tag}: labels match ({len(tl)} top words); lower words {len(low)}; claimed collapse vectors {X.shape[1]}")
    nt = X.shape[1]; Xf = np.zeros((len(tl), nt), complex); Xf[keep] = X
    # ---- 2. f = 3,4 fresh inputs
    S34 = int(os.environ.get('V_S34', S)); A34, L34 = mats([3, 4], S34, 7001, mv, R, W, B, tl, low)
    Q34 = orth(L34); del L34
    log(f"f=3,4 fresh: rows {A34.shape[0]}, lower rank {Q34.shape[1]}")
    Y = A34 @ Xf
    res34 = [float(np.linalg.norm(proj(Q34, Y[:, j])) / np.linalg.norm(Y[:, j])) for j in range(nt)]
    log(f"claimed collapse vectors at f=3,4 (fresh inputs): residuals {[f'{x:.1e}' for x in res34]}")
    rec["res34"] = res34
    cn = np.linalg.norm(proj(Q34, A34), axis=0) / np.linalg.norm(A34, axis=0)
    triv = np.flatnonzero(cn < 1e-8)
    log(f"trivial columns at f=3,4 (fresh): {triv.size}; npz kept {keep.size} -> implied trivial {len(tl)-keep.size}")
    Nd, sA = null(A34)
    log(f"column dependencies of all top columns at f=3,4: {Nd.shape[1]}")
    # full f=3,4 collapse count with my inputs (projected SVD) -- completeness of the extraction
    if os.environ.get("V_NT34", "1") == "1":
      Ap = proj(Q34, A34[:, np.setdiff1d(np.arange(len(tl)), triv)])
      rk = int((np.linalg.svd(A34[:, np.setdiff1d(np.arange(len(tl)), triv)], compute_uv=False) > sA[0] * 1e-9).sum())
      sp = np.linalg.svd(Ap, compute_uv=False); rp = int((sp > sp[0] * 1e-9).sum())
      log(f"non-trivial collapse count at f=3,4 with my inputs: {rk - rp} (claimed {nt}); gap {sp[rp-1]/sp[0]:.1e} / "
        f"{(sp[rp]/sp[0]) if rp < sp.size else 0:.1e}")
      rec["nt34_mine"] = rk - rp; del Ap
    del A34, Q34
    E = np.zeros((len(tl), triv.size), complex); E[triv, np.arange(triv.size)] = 1
    Cb = orth(np.concatenate([Xf, Nd, E], 1), 1e-10)
    log(f"candidate space dim {Cb.shape[1]} (= X {nt} + dependencies {Nd.shape[1]} + trivial {triv.size}, overlaps removed)")
    # ---- my reducibility
    red = np.zeros(len(tl), bool); pairs = {}
    for i, lab in enumerate(tl):
        p = reducible_dense(3, mv, lab, n)
        if p is not None: red[i] = True; pairs[i] = p
    red4 = red.copy(); conf = {}
    for i in np.flatnonzero(red):       # confirm at f = 4 on 32 random basis vectors
        key = tuple(tl[i][:n])
        if key not in conf:
            Pp = _pp(4); rng = np.random.default_rng(len(conf))
            cols = rng.choice(Pp.n, 32, replace=False); V = np.zeros((Pp.n, 32), complex); V[cols, np.arange(32)] = 1
            j, k = pairs[i]
            Yk = app_pulled(Pp, mv, tl[i], k, app_pulled(Pp, mv, tl[i], j, V))
            conf[key] = bool(((np.abs(Yk) > 1e-8).sum(0) == 1).all())
        red4[i] = conf[key]
    log(f"my reducibility: {int(red.sum())} top words reducible at f=3, {int(red4.sum())} confirmed at f=4; "
        f"npz red: {int(red_np.sum())}; agreement on kept columns: {int((red4[keep] == red_np[keep]).sum())}/{keep.size}")
    rec["red_mine"] = int(red4.sum()); rec["red_npz"] = int(red_np.sum())
    rec["red_disagree"] = [tl[i] for i in keep if red4[i] != red_np[i]][:10]
    # ---- 3. f = 5
    def stage(f, nrows, seed):
        rng = np.random.default_rng(seed + 99)
        rows = np.sort(rng.choice(8 ** f, min(nrows, 8 ** f), replace=False))
        Bf, Lf = mats([f], 1, seed, mv, R, W, B, tl, low, rows=[rows])
        Qf = orth(Lf); del Lf
        log(f"f={f} seed {seed}: rows {Bf.shape[0]} of {8**f}, lower rank {Qf.shape[1]}")
        out = {}
        Yx = Bf @ Xf
        per = [float(np.linalg.norm(proj(Qf, Yx[:, j])) / np.linalg.norm(Yx[:, j])) for j in range(nt)]
        Uy = orth(Yx); spec = np.sort(np.linalg.svd(proj(Qf, Uy), compute_uv=False))
        log(f"f={f}: claimed X per-vector residuals {[f'{x:.2e}' for x in per]}; X-subspace spectrum {[f'{x:.2e}' for x in spec]}")
        out["per"] = per; out["spec"] = [float(x) for x in spec]
        ct = np.linalg.norm(proj(Qf, Bf[:, triv]), axis=0) / np.linalg.norm(Bf[:, triv], axis=0) if triv.size else np.array([0.0])
        log(f"f={f}: f=3,4-trivial columns, max residual at f={f}: {ct.max():.1e}")
        out["triv_max"] = float(ct.max())
        Yc = Bf @ Cb; sY = np.linalg.svd(Yc, compute_uv=False); rY = int((sY > sY[0] * 1e-9).sum())
        Nc, sP = null(proj(Qf, Yc), 1e-9)
        Z = Yc @ Nc                      # images of the coefficient vectors that land in L_f (incl. kernel)
        sZ = np.linalg.svd(Z, compute_uv=False) if Z.shape[1] else np.array([0.0])
        coll = int((sZ > sY[0] * 1e-9).sum())
        rP = Cb.shape[1] - Nc.shape[1]
        log(f"f={f}: candidate space image rank {rY}, mod lower {rP}: collapses surviving {rY - rP} (from images: {coll}); "
            f"gap {sP[rP-1]/sP[0]:.1e} / {(sP[rP]/sP[0]) if rP < sP.size else 0:.1e}")
        out["coll"] = rY - rP
        Cn = Cb @ Nc
        # collapses without the trivial directions (compare with the claimed survivor count)
        Cb2 = orth(np.concatenate([Xf, Nd], 1), 1e-10); Y2 = Bf @ Cb2
        s2 = np.linalg.svd(Y2, compute_uv=False); r2 = int((s2 > s2[0] * 1e-9).sum())
        p2 = np.linalg.svd(proj(Qf, Y2), compute_uv=False); q2 = int((p2 > s2[0] * 1e-9).sum())
        log(f"f={f}: span(X) + f=3,4 dependencies (no trivial cols): image rank {r2}, mod lower {q2}: collapses {r2 - q2}")
        out["coll_noE"] = r2 - q2
        # B+ (independent): irreducible coordinates of the f=3,4 candidate space against L+ = lower + reducible + trivial
        for nm, rd in (("mine", red4), ("npz", red_np)):
            irr = np.flatnonzero(~rd & ~np.isin(np.arange(len(tl)), triv))
            redc = np.flatnonzero(rd | np.isin(np.arange(len(tl)), triv))
            Qp = np.concatenate([Qf, orth(proj(Qf, Bf[:, redc]))], 1) if redc.size else Qf
            Ci = orth(Cb[irr], 1e-10)
            Yi = Bf[:, irr] @ Ci
            si = np.linalg.svd(Yi, compute_uv=False); ri = int((si > si[0] * 1e-9).sum()) if si.size else 0
            pi_ = np.linalg.svd(proj(Qp, Yi), compute_uv=False) if Yi.shape[1] else np.array([])
            qi = int((pi_ > si[0] * 1e-9).sum()) if si.size else 0
            Ui = orth(Yi) if Yi.shape[1] else Yi
            spi = np.sort(np.linalg.svd(proj(Qp, Ui), compute_uv=False)) if Ui.shape[1] else np.array([])
            log(f"f={f}: [B+ {nm} reducibility] irreducible cols {irr.size}, candidate dim {Ci.shape[1]}, image rank {ri}, "
                f"mod lower+reducible {qi}: IRREDUCIBLE collapses {ri - qi}; smallest residuals {[f'{x:.1e}' for x in spi[:4]]}")
            out["irr_coll_" + nm] = ri - qi
        return out, Cn
    rec["f5"], Cn5 = stage(5, rows5, 5101)
    if rec["f5"]["coll"] and rows6:
        Pp = None
        rng = np.random.default_rng(66)
        rows = np.sort(rng.choice(8 ** 6, rows6, replace=False))
        Bf, Lf = mats([6], 1, 6101, mv, R, W, B, tl, low, rows=[rows])
        Qf = orth(Lf); del Lf
        Yc = Bf @ Cn5; sY = np.linalg.svd(Yc, compute_uv=False); rY = int((sY > sY[0] * 1e-9).sum())
        Uy = orth(Yc); spec = np.sort(np.linalg.svd(proj(Qf, Uy), compute_uv=False))
        log(f"f=6: f=5 collapse space image rank {rY}; spectrum vs lower span {[f'{x:.1e}' for x in spec]}")
        rec["f6_spec"] = [float(x) for x in spec]
    print("RESULT", json.dumps(rec, default=str), flush=True)
    log("done")


if __name__ == "__main__":
    main()
