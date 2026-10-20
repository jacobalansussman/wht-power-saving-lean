"""tri_als.py -- w7-xwborder (15 Oct 2026).  EXACT search for a TRIANGULAR 2 in 3 at block 2 in XW.

Why (PROOF.md s4, Lemma Lambda1 and the triangular reduction): an order-1 border 2 in 3 whose base is driven by x only
and whose exploding readouts sit on array 2 only is the same thing as an EXACT 3-move program with
      out1(x, y) = K x + T y   (T arbitrary),        out2(x, y) = K y,
read with inputs (x, eps y) and outputs (out1, out2 / eps).  Theorem 2W's rank argument allows it only if the layer-1
part of T has full rank.  This script searches for such programs with exact block least squares (every block of
coefficients enters the outputs affinely), fitting f in FS jointly, then re-scores on FRESH probes at f in FCHECK.
Residual = sqrt( |out1(x,0) - Kx|^2 + |out2(x,0)|^2 + |out2(0,y) - Ky|^2 ) / sqrt( |Kx|^2 + |Ky|^2 ).
Controls: 'plain4' (pattern of 4 moves, target as above: solvable), 'one' (array 2's target D_1 y instead of K y:
solvable in 3 moves), 'diag' (out1's y-part constrained to 0 too: the ordinary 2 in 3, excluded by Theorem T').
usage: python3 -I tri_als.py pattern menu FS sweeps seed target [S] [FCHECK]
"""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import xw2lib

MENUS = {
    "small": ["1", "X1", "Z1", "X2", "Z2", "lZ1", "lX1", "P", "C01", "CS"],
    "chain": ["1", "a1", "a1i", "a0", "a0i", "X1", "X2", "Z1", "Z2", "lZ1", "lX1", "lZ2", "lX2", "P", "C01", "CS"],
    "braid": ["1", "S1", "S2", "S3", "a1", "a1i", "a0", "a0i", "X1", "X2", "Z1", "Z2", "lZ1", "lX1", "P", "C01", "CS",
              "Pr0", "Pr1"],
    "wide": ["1", "X1", "X2", "Z1", "Z2", "S1", "S2", "S3", "P", "C01", "C10", "CZ", "CS", "CSi", "A", "Ai", "B", "Bi",
             "a0", "a1", "a0i", "a1i", "b1", "E01", "E10", "Pr0", "Pr1", "lCS", "fCS", "lZ1", "lX1", "fX2"],
    "lamp": ["1", "a1", "a1i", "a0", "a0i", "fX2", "a1.fX2.a1i", "fZ2", "a1i.fZ2.a1", "fX1", "a0.fX1.a0i", "fZ1",
             "a0i.fZ1.a0", "X1", "X2", "Z1", "Z2", "S1", "S2", "P", "C01", "CS", "b1", "Pr0", "Pr1"],
}


class Prog:
    def __init__(self, FS, zs, menu, S, seed, target):
        self.FS, self.zs, self.menu, self.M = FS, zs, menu, len(zs)
        self.G = xw2lib.Game(FS)
        self.target = target
        rng = np.random.default_rng(seed)
        self.mono = {b: self.G.mono(b) for b in menu}
        self.nb = len(menu)
        self.In, self.T, self.mask = {}, {}, {}
        for f in FS:
            n = 4 ** f
            X = (rng.standard_normal((n, S)) + 1j * rng.standard_normal((n, S))) / np.sqrt(2)
            Y = (rng.standard_normal((n, S)) + 1j * rng.standard_normal((n, S))) / np.sqrt(2)
            Z = np.zeros((n, S), complex)
            In0 = np.concatenate([X, Z], 1); In1 = np.concatenate([Z, Y], 1)
            self.In[f] = [In0, In1]
            KX, KY = self.G.K(X, f), self.G.K(Y, f)
            t2 = KY if target != "one" else self.G.apply_D(1, Y, f)
            # targets for out1 (only x-columns constrained, unless target 'diag') and out2 (all columns)
            self.T[f] = [np.concatenate([KX, Z], 1), np.concatenate([Z, t2], 1)]
            mu = float(os.environ.get("MU", "0"))      # weight on out1's y-part (keeps the cross-term T bounded)
            m1 = np.ones(2 * S) if target == "diag" else np.r_[np.ones(S), np.sqrt(mu) * np.ones(S)]
            self.mask[f] = [m1, np.ones(2 * S)]
        self.WIn = {f: [[self.G.apply_mono(self.mono[b][f], self.In[f][l]) for l in range(2)] for b in menu] for f in FS}
        nb = self.nb
        self.shapes = {}
        for k in range(self.M):
            self.shapes[('g', k)] = (nb, 2)
            for j in range(k):
                self.shapes[('h', k, j)] = (nb,)
        for l in range(2):
            self.shapes[('R', l)] = (nb, 2)
            for j in range(self.M):
                self.shapes[('r', l, j)] = (nb,)
        self.keys = list(self.shapes)

    def wapply(self, V, f):
        return [self.G.apply_mono(self.mono[b][f], V) for b in self.menu]

    def forward(self, P, f):
        u = []; Wu = []
        for k in range(self.M):
            c = np.zeros_like(self.In[f][0])
            g = P[('g', k)]
            for bi in range(self.nb):
                for l in range(2):
                    if g[bi, l] != 0: c = c + g[bi, l] * self.WIn[f][bi][l]
            for j in range(k):
                h = P[('h', k, j)]
                for bi in range(self.nb):
                    if h[bi] != 0: c = c + h[bi] * Wu[j][bi]
            uk = self.G.apply_D(self.zs[k], c, f)
            u.append(uk); Wu.append(self.wapply(uk, f))
        outs = []
        for l in range(2):
            o = np.zeros_like(self.In[f][0])
            R = P[('R', l)]
            for bi in range(self.nb):
                for l2 in range(2):
                    if R[bi, l2] != 0: o = o + R[bi, l2] * self.WIn[f][bi][l2]
            for j in range(self.M):
                r = P[('r', l, j)]
                for bi in range(self.nb):
                    if r[bi] != 0: o = o + r[bi] * Wu[j][bi]
            outs.append(o)
        return outs, u, Wu

    def resid_vec(self, outs, f):
        return [((outs[l] - self.T[f][l]) * self.mask[f][l]) for l in range(2)]

    def score(self, P, FS=None):
        FS = FS or self.FS
        num = den = 0.0
        for f in FS:
            outs, _, _ = self.forward(P, f)
            for l in range(2):
                e = (outs[l] - self.T[f][l]) * self.mask[f][l]
                num += np.vdot(e, e).real; den += np.vdot(self.T[f][l] * self.mask[f][l], self.T[f][l] * self.mask[f][l]).real
        return np.sqrt(num / den)

    def prop(self, P, f, k, dc, Wu):
        """linear propagation of a content perturbation dc at move k to the two (masked) outputs."""
        du = {k: self.G.apply_D(self.zs[k], dc, f)}
        Wdu = {k: self.wapply(du[k], f)}
        for m in range(k + 1, self.M):
            c = None
            for j in range(k, m):
                h = P[('h', m, j)]
                for bi in range(self.nb):
                    if h[bi] != 0:
                        t = h[bi] * Wdu[j][bi]; c = t if c is None else c + t
            du[m] = self.G.apply_D(self.zs[m], c if c is not None else np.zeros_like(dc), f)
            Wdu[m] = self.wapply(du[m], f)
        outs = []
        for l in range(2):
            o = np.zeros_like(dc)
            for j in range(k, self.M):
                r = P[('r', l, j)]
                for bi in range(self.nb):
                    if r[bi] != 0: o = o + r[bi] * Wdu[j][bi]
            outs.append(o * self.mask[f][l])
        return outs

    def tri_report(self, P, FS=None):
        """unpenalised triangular residual and |T y| / |K y| (T = out1's y-part)"""
        FS = FS or self.FS
        num = den = tn = td = 0.0
        for f in FS:
            outs, _, _ = self.forward(P, f)
            S = outs[0].shape[1] // 2
            e1 = outs[0][:, :S] - self.T[f][0][:, :S]; e2 = outs[1] - self.T[f][1]
            if self.target == "diag":
                e1 = outs[0] - self.T[f][0]
            num += np.vdot(e1, e1).real + np.vdot(e2, e2).real
            den += np.vdot(self.T[f][0], self.T[f][0]).real + np.vdot(self.T[f][1], self.T[f][1]).real
            tn += np.vdot(outs[0][:, S:], outs[0][:, S:]).real; td += np.vdot(self.T[f][1][:, S:], self.T[f][1][:, S:]).real
        return np.sqrt(num / den), np.sqrt(tn / td)

    def block_update(self, P, key, ridge):
        """exact least squares in one block (outputs are affine in it)."""
        shp = self.shapes[key]; nvar = int(np.prod(shp))
        cur = {f: self.forward(P, f) for f in self.FS}
        Jcols = [[] for _ in range(nvar)]
        for f in self.FS:
            outs, u, Wu = cur[f]
            if key[0] in ('r', 'R'):
                l = key[1]
                for i in range(nvar):
                    if key[0] == 'r':
                        v = Wu[key[2]][i]
                    else:
                        bi, l2 = divmod(i, 2); v = self.WIn[f][bi][l2]
                    o = [np.zeros_like(v), np.zeros_like(v)]; o[l] = v * self.mask[f][l]
                    Jcols[i].append(np.concatenate([o[0].ravel(), o[1].ravel()]))
            else:
                k = key[1]
                for i in range(nvar):
                    if key[0] == 'g':
                        bi, l2 = divmod(i, 2); dc = self.WIn[f][bi][l2]
                    else:
                        dc = Wu[key[2]][i]
                    o = self.prop(P, f, k, dc, Wu)
                    Jcols[i].append(np.concatenate([o[0].ravel(), o[1].ravel()]))
        J = np.stack([np.concatenate(c) for c in Jcols], 1)
        Ocur = np.concatenate([np.concatenate([(cur[f][0][l] * self.mask[f][l]).ravel() for l in range(2)]) for f in self.FS])
        Tm = np.concatenate([np.concatenate([(self.T[f][l] * self.mask[f][l]).ravel() for l in range(2)]) for f in self.FS])
        th0 = P[key].ravel()
        b = Tm - Ocur + J @ th0
        A = J.conj().T @ J; y = J.conj().T @ b
        lam = ridge * (np.trace(A).real / nvar + 1e-30)
        th = np.linalg.solve(A + lam * np.eye(nvar), y)
        P[key] = th.reshape(shp)
        return P


def main():
    pat = [int(t) for t in sys.argv[1].split(",")]
    menu = MENUS[sys.argv[2]]
    FS = [int(t) for t in sys.argv[3].split(",")]
    sweeps = int(sys.argv[4]); seed = int(sys.argv[5]); target = sys.argv[6]
    S = int(sys.argv[7]) if len(sys.argv) > 7 else 8
    FCHECK = [int(t) for t in sys.argv[8].split(",")] if len(sys.argv) > 8 else [2, 3, 4, 5]
    pr = Prog(FS, pat, menu, S, seed, target)
    rng = np.random.default_rng(seed + 1000)
    P = {k: 0.3 * (rng.standard_normal(s) + 1j * rng.standard_normal(s)) / np.sqrt(2) for k, s in pr.shapes.items()}
    t0 = time.time(); hist = []
    for sw in range(sweeps):
        order = list(pr.keys); rng.shuffle(order)
        for key in order:
            P = pr.block_update(P, key, 1e-7)
        r = pr.score(P)
        hist.append(r)
        mx = max(np.abs(v).max() for v in P.values())
        tr, tt = pr.tri_report(P)
        print("sweep %3d  fit %.6f  tri %.6f  |T| %.3f  max|coef| %.2f  t %.0fs" % (sw, r, tr, tt, mx, time.time() - t0), flush=True)
        if r < 1e-9: break
    # fresh-probe re-score (new probes; same coefficients) at FCHECK
    chk = Prog(FCHECK, pat, menu, S, seed + 7777, target)
    per_f = {f: chk.score(P, [f]) for f in FCHECK}
    print("FRESH re-score:", "  ".join("f=%d %.6f" % (f, per_f[f]) for f in FCHECK), flush=True)
    print("FRESH tri/|T|:", "  ".join("f=%d %.6f/%.3f" % ((f,) + chk.tri_report(P, [f])) for f in FCHECK), flush=True)
    tag = "%s_%s_%s_s%d_%s_mu%s" % (sys.argv[1], sys.argv[2], sys.argv[3].replace(",", ""), seed, target, os.environ.get("MU", "0"))
    np.savez("tri_%s.npz" % tag, **{"|".join(map(str, k)): v for k, v in P.items()}, menu=np.array(menu))
    with open("tri_%s.json" % tag, "w") as fh:
        json.dump({"pattern": pat, "menu": sys.argv[2], "FS": FS, "fit": hist[-1], "fresh": per_f,
                   "max_coef": float(max(np.abs(v).max() for v in P.values()))}, fh)


if __name__ == "__main__":
    main()
