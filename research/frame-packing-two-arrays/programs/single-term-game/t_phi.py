#!/usr/bin/env python3
"""t_phi.py: test the potential Phi(V) = sum_{k>=1} dim eps(G_k), G_k = {v in V : every term at level >= k},
eps = augmentation per input (sum of coefficients), along random plays in toy frame games.
If every move raised Phi by at most 1, then M >= L1 + L2 (Conjecture T) would follow, since an available single of
input l at level L_l puts e_l into eps(G_k) for all k <= L_l.  Plays: single-term (stg-style moves, as in
t_stgcal.py) and, for comparison, MULTI-term (moved vectors random elements of V)."""
import sys, random
import numpy as np
import scipy.linalg as _sl
_np_svd = np.linalg.svd
def _safe_svd(A, full_matrices=True, **kw):
    try:
        return _np_svd(A, full_matrices=full_matrices, **kw)
    except np.linalg.LinAlgError:
        return _sl.svd(A, full_matrices=full_matrices, lapack_driver='gesvd')
np.linalg.svd = _safe_svd
from toylab import *
from sweep import GAMES
from stg import *
from t_stg_toys import paid_classes

def run(name, ci, plays=30, steps=8, seed=1, multi=False):
    if name.startswith("H"):
        g = hamming_game(int(name[1])); G = g.fg.G; P = g.fg.P
    else:
        G, P, _ = GAMES[name]()
        g = from_toylab(G, P, paid_classes(G, P)[ci], name=name)
    n = G.n; H = g.fg.H
    lvl = np.array([g.lv[g.elem_pos[e]] for e in range(n)] * 2)
    rng = random.Random(seed); nrng = np.random.default_rng(seed)
    coset_elems = {}
    for e in range(n): coset_elems.setdefault(g.elem_pos[e], []).append(e)
    def tr(v, h):
        w = np.zeros(2 * n, complex)
        for l in range(2):
            idx = np.nonzero(v[l * n:(l + 1) * n])[0]
            for e in idx: w[l * n + G.m(h, int(e))] += v[l * n + e]
        return w
    def basis(gens):
        rows = []
        for l in range(2):
            for p in P:
                v = np.zeros(2 * n, complex); v[l * n + p] = 1; rows.append(v)
        for v in gens:
            for p in P: rows.append(tr(v, p))
        A = np.array(rows); u, s, vh = np.linalg.svd(A, full_matrices=False)
        return vh[: int((s > 1e-9 * s[0]).sum())]
    def phi(B):
        tot = 0; per = []
        for k in range(1, g.m + 1):
            out = np.nonzero(lvl < k)[0]
            if len(out) == 0: Gk = B
            else:
                Bo = B[:, out]
                # null space of c -> c Bo
                u, s, vh = np.linalg.svd(Bo.T)
                r = int((s > 1e-9 * max(1, s[0] if len(s) else 1)).sum())
                Nsp = vh[r:].conj()            # rows c with c Bo = 0  (c in C^rank(B))
                Gk = Nsp @ B
            if Gk.shape[0] == 0: d = 0
            else:
                E = np.stack([Gk[:, :n].sum(1), Gk[:, n:].sum(1)], 1)
                d = int(np.linalg.matrix_rank(E, tol=1e-7)) if E.size else 0
            per.append(d); tot += d
        return tot, per
    worst = -9; bad = 0; nmoves = 0; ex = None
    for play in range(plays):
        gens = []; pairs = []
        A = [{0}, {0}]           # available positions (true, recomputed)
        B = basis(gens); ph0, _ = phi(B)
        for st in range(steps):
            # true availability
            av = [set(), set()]
            for l in range(2):
                for x, els in coset_elems.items():
                    v = np.zeros(2 * n, complex); v[l * n + els[0]] = 1
                    if np.linalg.norm(v - (v @ B.conj().T) @ B) < 1e-8: av[l].add(x)
            h = rng.choice(H)
            if multi == "G":
                # a random element of G_k(V) (all terms at level >= k) for a random k: packed high-level vectors
                k = rng.randint(1, g.m)
                out = np.nonzero(lvl < k)[0]
                u_, s_, vh_ = np.linalg.svd(B[:, out].T)
                r_ = int((s_ > 1e-9 * max(1, s_[0])).sum())
                Nsp = vh_[r_:].conj()
                if Nsp.shape[0] == 0 or rng.random() < 0.3:
                    c = nrng.normal(size=B.shape[0]) + 1j * nrng.normal(size=B.shape[0])
                    v = (c * (nrng.random(B.shape[0]) < 0.3)) @ B
                else:
                    c = nrng.normal(size=Nsp.shape[0]) + 1j * nrng.normal(size=Nsp.shape[0])
                    v = (c @ Nsp) @ B
                v[np.abs(v) < 1e-10] = 0
                if np.linalg.norm(v) < 1e-9: continue
                w = tr(v, h)
            elif multi:
                # random element of V (combination of a few basis vectors), moved
                c = nrng.normal(size=B.shape[0]) + 1j * nrng.normal(size=B.shape[0])
                mask = nrng.random(B.shape[0]) < 0.3
                v = (c * mask) @ B
                if np.linalg.norm(v) < 1e-9: continue
                w = tr(v, h)
            else:
                kind = rng.choice(["s1", "s2", "p", "p", "c", "c"])
                if kind == "c" and not pairs: kind = "p"
                if kind in ("s1", "s2"):
                    l = int(kind[1]) - 1; x = rng.choice(sorted(av[l])); e = rng.choice(coset_elems[x])
                    v = np.zeros(2 * n, complex); v[l * n + e] = 1
                elif kind == "p":
                    x0 = rng.choice(sorted(av[0])); y0 = rng.choice(sorted(av[1]))
                    v = np.zeros(2 * n, complex)
                    v[rng.choice(coset_elems[x0])] = nrng.normal() + 1j * nrng.normal()
                    v[n + rng.choice(coset_elems[y0])] = nrng.normal() + 1j * nrng.normal()
                else:
                    v = rng.choice(pairs)
                w = tr(v, h)
                if np.count_nonzero(np.abs(w[:n]) > 1e-12) == 1 and np.count_nonzero(np.abs(w[n:]) > 1e-12) == 1: pairs.append(w)
            gens.append(w); nmoves += 1
            B = basis(gens); ph1, per = phi(B)
            if ph1 - ph0 > worst: worst = ph1 - ph0
            if ph1 - ph0 > 1:
                bad += 1
                if ex is None: ex = (play, st, ph0, ph1, per)
            ph0 = ph1
    print("%s class %d %s: %d moves, max increase of Phi per move = %d, moves with increase > 1: %d %s" % (
        name, ci, "MULTI-term" if multi else "single-term", nmoves, worst, bad, ("example " + str(ex)) if ex else ""), flush=True)

if __name__ == "__main__":
    for spec in sys.argv[1:]:
        nm, ci, mode, pl, stp = spec.split(":"); run(nm, int(ci), plays=int(pl), steps=int(stp), multi=("G" if mode == "G" else (mode == "multi")))
