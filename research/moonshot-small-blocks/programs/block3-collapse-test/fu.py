"""xw3-extract (v-ms-xw3 follow-up, 14 Oct 2026): extract the non-trivial collapses that path_ext.py flagged at
f = 3, 4 against the COMPLETE lower span, re-test the WHOLE collapse SUBSPACE at f = 5 (fresh inputs, same
coefficients), and classify.

Stage A (f = 3, 4, seed 31, exactly as path_ext.py): collapse space X34 (nt vectors over the non-trivial top columns,
orthogonal to the column dependencies), plus the candidate space N34 (= X34 + column dependencies at f = 3, 4).
Stage B (f = 5, fresh seed): for every vector of X34 its relative distance from the complete lower span at f = 5;
the relative-residual SPECTRUM of the subspace (min over combinations, i.e. singular values of P_perp U_Y where
Y = B5 X34 = U_Y S V^H); number of survivors (< 1e-6); and the count dim(span(B5 N34) cap L5) (kernel included).
Stage R (reducibility, dense operators at f = 3, confirmed on basis vectors at f = 4): a top column is REDUCIBLE when
two of its paid moves, pulled back through the mixes (M'_i = P_i^-1 M_i P_i), satisfy M'_k M'_j^-1 = monomial and
one of them commutes with every pulled move in between: then the word equals a free product times a word with two
fewer paid moves (D_z D_z = X^z), i.e. it is formally a lower term.
Stage A+/B+: the same extraction for the IRREDUCIBLE columns against L+ = lower span + reducible columns.  A collapse
of class (c) must survive here; survivors of B but not B+ are explained by reducible columns (class (b)).
Stage F6 (only for survivors): same coefficients at f = 6, fresh input, rows subsampled (membership is preserved
exactly by any row subset; a non-member keeps a positive residual).
usage: fu.py <tag> <moves> <R> <W1> ... <W_{n-1}> <B> <f1,f2> <S> [seed5]
"""
import sys, os, time, json, itertools
import numpy as np
np.seterr(all="ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xtalk_ext import L, fl, zc
from path_ext import build
from xw3lib import Game

T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.0f}s]", *a, flush=True)


def orth(A, tol=1e-9):
    if A.shape[1] == 0: return np.zeros((A.shape[0], 0), A.dtype), np.array([])
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    r = int((s > s[0] * tol).sum()); return U[:, :r], s


def proj_off(Q, Y):
    return Y - Q @ (Q.conj().T @ Y)


def norm2(A, it=30, seed=0):
    rng = np.random.default_rng(seed); x = rng.standard_normal(A.shape[1]) + 0j
    for _ in range(it):
        y = A @ x; x = A.conj().T @ y; x /= np.linalg.norm(x)
    return float(np.linalg.norm(A @ x))


def extract(A3, Q, cols, tag):
    """collapse space of the columns `cols` of A3 against span(Q): returns (keep, N, X, info)"""
    if cols.size == 0:
        log(f"[{tag}] no columns"); z = np.zeros((0, 0), complex)
        return cols, z, z, dict(tag=tag, ncols=0, nt=0)
    A = A3[:, cols]; Ap = proj_off(Q, A)
    cn = np.linalg.norm(Ap, axis=0) / np.maximum(np.linalg.norm(A, axis=0), 1e-300)
    triv = cols[cn < 1e-8]; keep = cols[cn >= 1e-8]
    if keep.size == 0:
        log(f"[{tag}] all {cols.size} columns trivial"); z = np.zeros((0, 0), complex)
        return keep, z, z, dict(tag=tag, ncols=int(cols.size), trivial=int(triv.size), nt=0)
    Ak = A3[:, keep]; Akp = proj_off(Q, Ak)
    _, sP, Vh = np.linalg.svd(Akp, full_matrices=False)
    rP = int((sP > sP[0] * 1e-9).sum()); N = Vh[rP:].conj().T
    Wm = Ak @ N
    nA = norm2(A3)
    if N.shape[1]:
        _, sw, vh = np.linalg.svd(Wm, full_matrices=False)
        nt = int((sw > nA * 1e-8).sum()); X = N @ vh[:nt].conj().T
    else:
        sw = np.array([]); nt = 0; X = np.zeros((keep.size, 0), complex)
    gap = (sP[rP - 1] / sP[0], sP[rP] / sP[0] if rP < sP.size else 0.0)
    res34 = [float(np.linalg.norm(proj_off(Q, Ak @ X[:, j])) / np.linalg.norm(Ak @ X[:, j])) for j in range(nt)]
    info = dict(tag=tag, ncols=int(cols.size), trivial=int(triv.size), keep=int(keep.size), dimN=int(N.shape[1]),
                nt=nt, kernel=int(N.shape[1] - nt), gap=[float(gap[0]), float(gap[1])],
                sw=[float(x / nA) for x in sw[:nt + 3]], res34_max=max(res34) if res34 else 0.0)
    log(f"[{tag}] cols {cols.size}: trivial {triv.size}, kept {keep.size}; null space {N.shape[1]} = collapses {nt} + "
        f"column dependencies {N.shape[1] - nt}; sv gap {gap[0]:.1e} / {gap[1]:.1e}; collapse sv/|A| "
        f"{[f'{x / nA:.1e}' for x in sw[max(0, nt - 2):nt + 2]]}; f=3,4 residual of the extracted vectors <= "
        f"{info['res34_max']:.1e}")
    return keep, N, X, info


def test_subspace(B, Q, keep, X, N, tag, lab):
    """B: top columns at the new f (all columns); X: collapse vectors over `keep`; returns survivors"""
    out = {}
    if X.shape[1] == 0:
        log(f"[{tag}] nothing to test"); return out, np.zeros((keep.size, 0), complex)
    Y = B[:, keep] @ X; Yp = proj_off(Q, Y)
    per = np.linalg.norm(Yp, axis=0) / np.linalg.norm(Y, axis=0)
    Uy, sy, Vyh = np.linalg.svd(Y, full_matrices=False)
    ry = int((sy > sy[0] * 1e-9).sum())
    Uy, sy, Vyh = Uy[:, :ry], sy[:ry], Vyh[:ry]
    Up = proj_off(Q, Uy)
    _, sp, Vph = np.linalg.svd(Up, full_matrices=False)
    order = np.argsort(sp); sp = sp[order]; Vph = Vph[order]
    surv = int((sp < 1e-6).sum())
    log(f"[{tag}] {X.shape[1]} collapse vectors (image rank {ry}); per-vector relative distance from the complete "
        f"lower span: {[f'{x:.2e}' for x in per]}")
    log(f"[{tag}] subspace spectrum (min relative distance over combinations, ascending): "
        f"{[f'{x:.2e}' for x in sp]}; SURVIVORS {surv}")
    # N-based count (column dependencies of f = 3, 4 included)
    YN = B[:, keep] @ N
    sN = np.linalg.svd(YN, compute_uv=False); rN = int((sN > sN[0] * 1e-9).sum())
    sNp = np.linalg.svd(proj_off(Q, YN), compute_uv=False); rNp = int((sNp > sN[0] * 1e-9).sum())
    log(f"[{tag}] candidate space incl. f=3,4 column dependencies: dim {N.shape[1]}, image rank {rN}, "
        f"mod lower {rNp}: collapses surviving {rN - rNp}")
    out = dict(per=[float(x) for x in per], spectrum=[float(x) for x in sp], survivors=surv, N_count=rN - rNp)
    # survivor coefficient vectors: c = V_y S_y^-1 d, d = right singular vector of Up
    Xs = np.zeros((keep.size, surv), complex)
    for j in range(surv):
        d = Vph[j].conj(); c = Vyh.conj().T @ (d / sy); x = X @ c; x /= np.linalg.norm(x); Xs[:, j] = x
        big = np.argsort(-np.abs(x))[:14]
        log(f"[{tag}]  survivor {j}: residual {sp[j]:.1e}; terms {[(lab[keep[i]], np.round(x[i], 3)) for i in big if abs(x[i]) > 1e-3]}")
    return out, Xs


# ------------------------------------------------------------------ reducibility (dense pulled-back moves)
class Dense:
    def __init__(self, mv, f=3):
        self.G = Game([f], W=1, S=1, seed=5); self.f = f; self.n = 8 ** f; self.mv = mv; self.cache = {}

    def I(self):
        return {self.f: np.eye(self.n, dtype=complex)}

    def apply_word(self, w, X, inv=False):
        if w in ("", "1"): return X
        return self.G.mono(self.G.winv(w) if inv else self.G.word(w), X)

    def apply_move(self, i, X, inv=False):
        z, C = self.mv[i]
        if C in ("", "1"): return self.G.D(z, X, inverse=inv)
        return self.G.mono(self.G.word(C), self.G.D(z, self.G.mono(self.G.winv(C), X), inverse=inv))

    def pulled(self, i, prefix, inv=False):
        key = (i, prefix, inv)
        if key not in self.cache:
            X = self.I()
            for w in prefix: X = self.apply_word(w, X)
            X = self.apply_move(i, X, inv)
            for w in reversed(prefix): X = self.apply_word(w, X, inv=True)
            self.cache[key] = X[self.f]
        return self.cache[key]


def is_mono(M, tol=1e-8):
    nz = (np.abs(M) > tol).sum(axis=0)
    return bool((nz == 1).all())


RCACHE = {}


def reducible(Dn, lab, n):
    """lab: (mix before move 0 = '1', mix before move 1, ..., mix before move n-1, readout)"""
    key = tuple(lab[:n])
    if key not in RCACHE: RCACHE[key] = _reducible(Dn, lab, n)
    return RCACHE[key]


def _reducible(Dn, lab, n):
    pref = [tuple(x for x in lab[1:i + 1] if x not in ("", "1")) for i in range(n)]
    Mp = [Dn.pulled(i, pref[i]) for i in range(n)]
    for j in range(n):
        for k in range(j + 1, n):
            Op = Mp[k] @ Dn.pulled(j, pref[j], inv=True)
            if not is_mono(Op): continue
            mid = range(j + 1, k)
            ok_k = all(np.abs(Mp[k] @ Mp[l] - Mp[l] @ Mp[k]).max() < 1e-9 for l in mid)
            ok_j = all(np.abs(Mp[j] @ Mp[l] - Mp[l] @ Mp[j]).max() < 1e-9 for l in mid)
            if ok_k or ok_j: return (j, k)
    return None


def confirm_f4(mv, lab, n, pair, rng):
    """check M'_k M'_j^-1 maps random basis vectors to multiples of basis vectors at f = 4"""
    Dn = Dense(mv, f=4) if not hasattr(confirm_f4, "D4") else confirm_f4.D4
    confirm_f4.D4 = Dn
    j, k = pair
    pref = [tuple(x for x in lab[1:i + 1] if x not in ("", "1")) for i in range(n)]
    cols = rng.choice(Dn.n, 24, replace=False)
    X = {4: np.zeros((Dn.n, cols.size), complex)}; X[4][cols, np.arange(cols.size)] = 1
    def app(i, X, inv):
        for w in pref[i]: X = Dn.apply_word(w, X)
        X = Dn.apply_move(i, X, inv)
        for w in reversed(pref[i]): X = Dn.apply_word(w, X, inv=True)
        return X
    Y = app(k, app(j, X, True), False)[4]
    return bool(((np.abs(Y) > 1e-8).sum(axis=0) == 1).all())


# ------------------------------------------------------------------ f = 6 with subsampled rows
def build_rows(Gf, mv, R, W, B, rows):
    """as path_ext.build, but every column is cut to `rows` (single f)"""
    import path_ext
    f = Gf.FS[0]
    orig = path_ext.fl
    path_ext.fl = lambda G, X: X[f].reshape(-1)[rows]
    try:
        top, low = build(Gf, mv, R, W, B)
    finally:
        path_ext.fl = orig
    return np.stack([t[0] for t in top], 1), np.stack([t[0] for t in low], 1)


def main():
    tag = sys.argv[1]
    a = sys.argv[1:]
    mv = [zc(m) for m in a[1].split(";")]; n = len(mv)
    R = L(a[2]); W = [L(a[3 + j]) for j in range(n - 1)]; B = L(a[2 + n])
    FS = [int(x) for x in a[3 + n].split(",")]; S = int(a[4 + n])
    seed5 = int(a[5 + n]) if len(a) > 5 + n else 501
    rec = dict(tag=tag, moves=a[1], R=a[2], W=a[3:2 + n], B=a[2 + n], FS=FS, S=S, seed5=seed5)
    log(f"# {tag}: moves={a[1]} |R|={len(R)} |W|={[len(w) for w in W]} |B|={len(B)} FS={FS} S={S}")
    # ---------------- stage A
    G = Game(FS, W=1, S=S, seed=31)
    top, low = build(G, mv, R, W, B)
    lab = [t[1] for t in top]
    A3 = np.stack([t[0] for t in top], 1); A2 = np.stack([t[0] for t in low], 1); del top, low
    log(f"rows {A3.shape[0]}, top {A3.shape[1]}, lower {A2.shape[1]}")
    Q, s2 = orth(A2); del A2
    log(f"lower rank {Q.shape[1]}")
    allc = np.arange(A3.shape[1])
    keep, N, X, infoA = extract(A3, Q, allc, "A f=3,4")
    rec["A"] = infoA
    # ---------------- stage R (reducibility of the kept columns)
    Dn = Dense(mv, f=3); rng = np.random.default_rng(3)
    red = np.zeros(A3.shape[1], bool); pairs = {}
    for i in keep:
        p = reducible(Dn, lab[i], n)
        if p is not None and confirm_f4(mv, lab[i], n, p, rng):
            red[i] = True; pairs[i] = p
    log(f"reducible kept columns: {int(red[keep].sum())} of {keep.size}; examples "
        f"{[(lab[i], pairs[i]) for i in list(pairs)[:8]]}")
    rec["reducible"] = int(red[keep].sum())
    # irreducible against L+ = lower + reducible
    redc = keep[red[keep]]; irr = keep[~red[keep]]
    Qp = np.concatenate([Q, orth(proj_off(Q, A3[:, redc]))[0]], 1) if redc.size else Q
    keepP, NP, XP, infoP = extract(A3, Qp, irr, "A+ f=3,4 irreducible vs lower+reducible")
    rec["Aplus"] = infoP
    del A3
    # ---------------- stage B (f = 5)
    if X.shape[1] == 0 and XP.shape[1] == 0:
        log("no collapse to test"); rec["verdict"] = "nothing at f=3,4"; print("RESULT", json.dumps(rec)); return
    F5 = int(os.environ.get("FU_F5", "5")); F6 = int(os.environ.get("FU_F6", "6"))
    G5 = Game([F5], W=1, S=1, seed=seed5)
    R5 = int(os.environ.get("FU_ROWS5", "0"))
    if R5 and R5 < 8 ** F5:      # row subsample: membership is kept exactly; a distance is a LOWER bound (absolute)
        rows5 = np.sort(np.random.default_rng(seed5 + 3).choice(8 ** F5, R5, replace=False))
        B5, L5 = build_rows(G5, mv, R, W, B, rows5)
        log(f"f={F5}: rows subsampled {R5} of {8 ** F5}")
    else:
        top5, low5 = build(G5, mv, R, W, B)
        B5 = np.stack([t[0] for t in top5], 1); L5 = np.stack([t[0] for t in low5], 1); del top5, low5
    Q5, _ = orth(L5); del L5
    log(f"f={F5} (seed {seed5}): rows {B5.shape[0]}, lower rank {Q5.shape[1]}")
    rec["B"], Xs = test_subspace(B5, Q5, keep, X, N, "B f=5", lab)
    red5 = keep[red[keep]]
    Q5p = np.concatenate([Q5, orth(proj_off(Q5, B5[:, red5]))[0]], 1) if red5.size else Q5
    rec["Bplus"], XsP = test_subspace(B5, Q5p, keepP, XP, NP, "B+ f=5 irreducible vs lower+reducible", lab)
    del B5, Q5, Q5p
    # ---------------- stage F6 for survivors
    if Xs.shape[1] or XsP.shape[1]:
        G6 = Game([F6], W=1, S=1, seed=seed5 + 1000)
        rr = np.random.default_rng(seed5 + 7)
        nrows = min(int(os.environ.get("FU_ROWS6", "40000")), 8 ** F6)
        rows = np.sort(rr.choice(8 ** F6, nrows, replace=False))
        B6, L6 = build_rows(G6, mv, R, W, B, rows)
        log(f"f={F6}: rows sampled {nrows} of {8 ** F6}, top {B6.shape[1]}, lower {L6.shape[1]}")
        Q6, _ = orth(L6); del L6
        res6 = []
        for j in range(Xs.shape[1]):
            y = B6[:, keep] @ Xs[:, j]; res6.append(float(np.linalg.norm(proj_off(Q6, y)) / np.linalg.norm(y)))
        log(f"[F6] survivors of B: relative distance at f=6 {[f'{x:.1e}' for x in res6]}")
        rec["F6"] = res6
        if XsP.shape[1]:
            Q6p = np.concatenate([Q6, orth(proj_off(Q6, B6[:, red5]))[0]], 1) if red5.size else Q6
            res6p = []
            for j in range(XsP.shape[1]):
                y = B6[:, keepP] @ XsP[:, j]
                res6p.append(float(np.linalg.norm(proj_off(Q6p, y)) / np.linalg.norm(y)))
            log(f"[F6] survivors of B+: relative distance at f=6 {[f'{x:.1e}' for x in res6p]}")
            rec["F6plus"] = res6p
    np.savez_compressed(f"fu_{tag}.npz", keep=keep, X=X, Xs=Xs, keepP=keepP, XP=XP, XsP=XsP, red=red,
                        lab=np.array(["|".join(l) for l in lab]))
    print("RESULT", json.dumps(rec), flush=True)
    log("done")


if __name__ == "__main__":
    main()
