"""symsolve.py -- mw-sym4: the H-invariant unit equations as a MULTILINEAR system, solved numerically (LM with variable
projection), structure fixed in advance.
Skeleton: a list of orbit-moves (rank 1); move t has a menu of summands (source generator, free letter b) with unknown
complex coefficients alpha_{t,k}; its content is sum_k alpha_k b_k u_{s_k} ('free' move: the tau-orbit of u is added)
or sum_k alpha_k sum_{h in H} tau_h(b_k u_{s_k}) ('sym' move: an orbit-sum content, u usually tau-fixed).  The
orbit images carry the same coefficients (tau is linear and exact).  End: for each array l, [K]x_l must lie in the
C[B]-span of all generators: closure over frames (bm.Game.closure, structure from generic alpha), then
min_beta ||A(alpha) beta - e_K|| (variable projection); LM over alpha.  A residual -> 0 is then reconstructed and
checked EXACTLY (rational/Gaussian-rational coefficients, bm span test) -- not done automatically here.
Usage: symsolve.py m W Hname shape menu starts seed      shape e.g. 'F,F,F' (F = free orbit move, S = sym move)"""
import sys, time
from symsearch import *

def rand_free_b(G, rng, m):
    n = 1 << m
    while True:
        A = [rng.randrange(1, n) for _ in range(m)]
        span = {0}; ok = True
        for a in A:
            if a in span: ok = False; break
            span |= {x ^ a for x in span}
        if ok: break
    c = rng.randrange(n)
    sig = []
    for x in range(n):
        y = c
        for j in range(m):
            if (x >> j) & 1: y ^= A[j]
        sig.append(y)
    a = [rng.randrange(4) for _ in range(m)]
    B = {(j, k): rng.randrange(2) for j in range(m) for k in range(j + 1, m)}
    q = []
    for x in range(n):
        b = [(x >> j) & 1 for j in range(m)]
        v = sum(a[j] * b[j] for j in range(m)) + 2 * sum(B[(j, k)] * b[j] * b[k] for j in range(m) for k in range(j + 1, m))
        q.append(v % 4)
    return G.mono(sig, q)

def fixed_dirs(R, rng, samples=4000):
    """rank-1 move elements delta = b D1 b^-1 with c delta c^-1 delta^-1 free for every c in H (tau-fixed directions),
    one per frame of delta"""
    G = R.G
    out = {}
    for _ in range(samples):
        b = rand_free_b(G, rng, R.m)
        d = G.mul(b, G.D1, np.conj(b).T)
        if all(G.is_monomial(G.mul(c, d, np.conj(c).T, np.conj(d).T)) for (c, pi, cp) in R.Hc[1:]):
            fk = G.frame_key(np.conj(d).T.copy())
            if fk not in out:
                out[fk] = d
    return list(out.values())

class Skel:
    def __init__(self, m, W, Hname, seed, R=None):
        self.R = Run(m, W, Hname, seed) if R is None else R
        self.G = self.R.G; self.W = W; self.m = m
        self.rng = random.Random(seed)
        G = self.G
        # generators: list of (terms: list of (l,gid)), coef source: ('inp',) or ('move', t, perm/None)
        self.gens = []          # term lists
        self.coef_of = []       # (kind, data)
        for l in range(W):
            self.gens.append([(l, G.ID)]); self.coef_of.append(('inp',))
        self.moves = []         # per move: dict(summands=[(src, b, hlist)], D, terms, maps)
        self.nalpha = 0
        self.fixed = []

    def dict_of(self, i, vals):
        return {t: v for t, v in zip(self.gens[i], vals)}

    def add_move(self, kind, summ, Dm_override=None):
        """summ: list of (src_index, b matrix).  kind 'F' or 'S'.  Builds the structure: the representative's term list,
        gather maps (for each summand and each h in the sum: src term index -> rep term index), then the tau-orbit."""
        G, R = self.G, self.R
        Hs = R.Hc if kind in ('S', 'A', 'B', 'C', 'N', 'L', 'T') else R.Hc[:1]
        nH = len(R.Hc)
        rpow = {'S': 0, 'A': 1, 'B': 2, 'C': 3}.get(kind, 0)
        chi = [np.exp(2j * np.pi * rpow * j / nH) for j in range(nH)]
        Dm = G.D1 if Dm_override is None else Dm_override
        terms = {}; order = []
        maps = []                               # list of (alpha_index, src, src_term_idx array, rep_term_idx array)
        for k, (src, b) in enumerate(summ):
            a_idx = self.nalpha + k
            for hj, h in enumerate(Hs):
                # tau_h(b u_src): term (l,g) -> tau_h( b g ) ; then D on the left
                si, ri = [], []
                for j, (l, g) in enumerate(self.gens[src]):
                    g1 = G.intern(G.mul(b, G.MATS[g]))
                    l2, g2 = R.tau(h, {(l, g1): 1}).popitem()[0] if h is not R.Hc[0] else (l, g1)
                    if kind == 'T':          # w3-blueprint: c-twisted orbit sum, tau(u) = c u for ANY direction
                        g2 = G.intern(G.mul(np.conj(h[0]).T, G.MATS[g2]))
                    g3 = G.intern(G.mul(Dm, G.MATS[g2]))
                    key = (l2, g3)
                    if key not in terms:
                        terms[key] = len(order); order.append(key)
                    si.append(j); ri.append(terms[key])
                maps.append((a_idx, src, np.array(si), np.array(ri), chi[hj] if kind in ('S', 'A', 'B', 'C', 'N', 'L', 'T') else 1.0))
        self.fixed.append(self.nalpha)
        self.nalpha += len(summ)
        rep = len(self.gens)
        self.gens.append(order); self.coef_of.append(('move', len(self.moves)))
        mv = {"kind": kind, "maps": maps, "rep": rep, "n": len(order), "images": []}
        self.moves.append(mv)
        # orbit images (structure check with generic coefficients carried through the actual structure)
        if not hasattr(self, 'testcv'):
            self.testcv = [np.ones(1, dtype=complex) for _ in range(self.W)]
            self.testal = {}
        v = np.zeros(len(order), dtype=complex)
        for (a, src, si, ri, ch) in maps:
            if a not in self.testal:
                self.testal[a] = complex(self.rng.gauss(0, 1), self.rng.gauss(0, 1))
            np.add.at(v, ri, ch * self.testal[a] * self.testcv[src][si])
        self.testcv.append(v)
        test = {t: v[j] for j, t in enumerate(order)}
        reps = [test]
        for h in R.Hc[1:]:
            img = R.tau(h, test)
            if any(self.R.S.same_up_to_free(img, w) for w in reps):
                continue
            reps.append(img)
            # term permutation: image term list in the same order as the representative's
            imgterms = [R.tau_term_cached(h, l, g) for (l, g) in order]
            self.gens.append(imgterms); self.coef_of.append(('img', len(self.moves) - 1)); self.testcv.append(v)
            mv["images"].append(len(self.gens) - 1)
        return 1 + len(mv["images"])

    def coefs(self, alpha):
        """coefficient vectors of every generator for complex alpha"""
        out = []
        for i, (kind, *d) in enumerate(self.coef_of):
            if kind == 'inp':
                out.append(np.ones(1, dtype=complex))
            elif kind == 'move':
                mv = self.moves[d[0]]
                v = np.zeros(mv["n"], dtype=complex)
                for (a, src, si, ri, ch) in mv["maps"]:
                    np.add.at(v, ri, ch * alpha[a] * out[src][si])
                out.append(v)
            else:
                out.append(out[self.moves[d[0]]["rep"]])
        return out

    def setup_final(self, cap=3000):
        """closure per array with index-coefficients; stores (nrows, ncols, entries as (r, c, gen, termidx))"""
        G = self.G
        self.final = []
        idx_gens = [{t: (i, j) for j, t in enumerate(T)} for i, T in enumerate(self.gens)]
        # bm closure wants dict gens {(l,gid): co}; we pass co = a token we can decode
        tok = []
        gdicts = []
        for i, T in enumerate(self.gens):
            d = {}
            for j, t in enumerate(T):
                d[t] = len(tok) + 1; tok.append((i, j))
            gdicts.append(d)
        for l in range(self.W):
            cl = G.closure(gdicts, (l, G.KID), cap)
            if cl is None:
                return False
            rows, cols, ent = cl
            # ent: (r, c, co) with co the token of the generator term (gen i, term j); the column's free multiplier only
            # moves terms, coefficient = coefficient of that term
            E = [(r, c, tok[co - 1][0], tok[co - 1][1]) for (r, c, co) in ent]
            self.final.append((len(rows), len(cols), E))
        return True

    def residual(self, alpha, want_beta=False):
        cv = self.coefs(alpha)
        res = []; betas = []
        for (nr, nc, E) in self.final:
            A = np.zeros((nr, nc), dtype=complex)
            for (r, c, i, j) in E:
                A[r, c] += cv[i][j]
            t = np.zeros(nr, dtype=complex); t[0] = 1
            if not np.isfinite(A).all():
                return (np.full(nr, 1e3, dtype=complex), None) if want_beta else np.full(nr, 1e3, dtype=complex)
            try:
                beta, *_ = np.linalg.lstsq(A, t, rcond=None)
            except np.linalg.LinAlgError:
                beta = np.zeros(nc, dtype=complex)
            res.append(A @ beta - t); betas.append(beta)
        r = np.concatenate(res)
        return (r, betas) if want_beta else r

def lm(f, x0, iters=60, lam=1e-2, eps=1e-6, tol=1e-12):
    """Levenberg-Marquardt on real vector x (f returns complex residual)"""
    x = x0.copy()
    def F(x):
        r = f(x); return np.concatenate([r.real, r.imag])
    r = F(x); cost = r @ r
    for it in range(iters):
        J = np.zeros((len(r), len(x)))
        for k in range(len(x)):
            xp = x.copy(); xp[k] += eps
            J[:, k] = (F(xp) - r) / eps
        g = J.T @ r; Hm = J.T @ J
        while True:
            try:
                dx = -np.linalg.solve(Hm + lam * np.diag(np.diag(Hm) + 1e-9), g)
            except np.linalg.LinAlgError:
                lam *= 10; continue
            xn = x + dx; rn = F(xn); cn = rn @ rn
            if cn < cost:
                x, r, cost = xn, rn, cn; lam = max(lam / 3, 1e-9); break
            lam *= 4
            if lam > 1e8:
                return x, cost
        if cost < tol:
            break
    return x, cost

def build_skeleton(m, W, Hname, shape, menu, seed, steps_hint=None):
    """shape: string of F/S per orbit move; menu: number of summands per move (candidates: progress letters of lead
    terms of recent generators, preferring other arrays, plus exact alignments onto the first summand's new term)"""
    sk = Skel(m, W, Hname, seed)
    G, R, rng = sk.G, sk.R, sk.rng
    cost = 0
    prev_ids = list(range(W))
    for tok in shape:
        kind, mode = (tok.split(':') + ['any'])[:2]
        if mode == 'inp':
            pool = list(range(W))
        elif mode == 'prev':
            pool = list(prev_ids)
        elif mode == 'prev2':
            pool = list(prev_ids) + [i for i in range(W)]
        else:
            pool = list(range(len(sk.gens)))[::-1]
        cands = []
        for i in pool:
            gd = {t: 1 for t in sk.gens[i]}
            ld = R.lead(gd)
            for (d0, l, g) in ld[:3]:
                if G.dist(g) < 1:
                    continue
                for b in R.plets(g, 1)[:3]:
                    cands.append((i, b, g))
            if mode == 'any' and len(cands) > 4 * menu:
                break
        if not cands:
            return None, None
        rng.shuffle(cands)
        cands.sort(key=lambda c: G.dist(c[2]))      # deepest (closest to the kernel frame) first
        Dm = G.D1
        if kind == 'P':
            src = pool[0]
            gd = {t: 1 for t in sk.gens[src]}
            ld = R.lead(gd)
            per_l = {}
            for (d0, l, g) in ld:
                if G.dist(g) < 1 or len(per_l.get(l, [])) >= 1:
                    continue
                per_l.setdefault(l, []).append(g)
            summ = []
            for l, gs in sorted(per_l.items()):
                for g in gs:
                    pls = R.plets(g, 1)
                    for b in rng.sample(pls, min(len(pls), menu)):
                        summ.append((src, b))
            if not summ:
                return None, None
            n0 = len(sk.gens)
            cost += sk.add_move('F', summ, None)
            prev_ids = list(range(n0, len(sk.gens)))
            continue
        if kind in ('S', 'A', 'B', 'C', 'L'):
            if not hasattr(sk, 'fdirs'):
                sk.fdirs = fixed_dirs(R, rng)
            if not sk.fdirs:
                return None, None
            # a fixed direction delta and a letter b with delta b g one step nearer the kernel frame, g a deep lead term
            found = None
            for (i, b0, g) in cands[:6]:
                for d in rng.sample(sk.fdirs, min(len(sk.fdirs), 6)):
                    pl = G.progress_letters(g, rng, want=(24 if kind == 'L' else 3), Dmat=d, steps=1)
                    if kind == 'L':
                        Ki = np.conj(G.K).T
                        pl = [b for b in pl if not any(G.is_monomial(G.mul(Ki, np.conj(b).T, c, b, G.K)) for (c, pi, cp) in R.Hc[1:])]
                    if pl:
                        found = (i, rng.choice(pl), g, d); break
                if found:
                    break
            if found is None:
                return None, None
        elif kind == 'N':
            if not hasattr(sk, 'fdirs'):
                sk.fdirs = fixed_dirs(R, rng)
            if not sk.fdirs:
                return None, None
            i, b0, g = cands[0]
            found = (i, b0, g, rng.choice(sk.fdirs))
        if kind in ('S', 'A', 'B', 'C', 'N', 'L'):
            first = (found[0], found[1], found[2]); Dm = found[3]
        else:
            first = cands[0]
        summ = [(first[0], first[1])]
        T = G.intern(G.mul(Dm, first[1], G.MATS[first[2]]))
        # aligned summands onto T from other generators (the pool first)
        for i in list(pool) + [j for j in range(len(sk.gens)) if j not in pool]:
            if len(summ) >= menu:
                break
            for (l, g) in sk.gens[i]:
                if i == first[0] and g == first[2]:
                    continue
                b = G.mul(np.conj(Dm).T, G.MATS[T], np.conj(G.MATS[g]).T)
                if G.is_monomial(b) and G.dist(g) >= 1:
                    summ.append((i, b)); break
        for c in cands[1:]:
            if len(summ) >= menu:
                break
            summ.append((c[0], c[1]))
        n0 = len(sk.gens)
        cost += sk.add_move(kind, summ, Dm if kind in ('S', 'A', 'B', 'C', 'N', 'L') else None)
        prev_ids = list(range(n0, len(sk.gens)))
    return sk, cost

def main():
    m, W, Hname, shape, menu, starts, seed = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4].split(','), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    nskel = int(sys.argv[8]) if len(sys.argv) > 8 else 1
    t0 = time.time()
    best_all = []
    for s in range(nskel):
        sk, cost = build_skeleton(m, W, Hname, shape, menu, seed * 100 + s)
        if sk is None:
            print("skeleton failed"); continue
        ok = sk.setup_final()
        if not ok:
            print(f"skel {s}: cost {cost} closure cap"); continue
        dims = [(nr, nc) for (nr, nc, E) in sk.final]
        N = sk.nalpha
        free_idx = [a for a in range(N) if a not in sk.fixed]
        nf = len(free_idx)
        def full(x):
            al = np.ones(N, dtype=complex)
            al[free_idx] = x[:nf] + 1j * x[nf:]
            return al
        print(f"skel {s}: cost {cost}, unknowns {nf} complex, final dims {dims}", flush=True)
        best = None
        for st in range(starts):
            x0 = np.array([sk.rng.gauss(0, 1) for _ in range(2 * nf)])
            f = lambda x: sk.residual(full(x))
            if nf == 0:
                c = float(np.sum(np.abs(f(x0)) ** 2)); best = (c, x0); break
            x, c = lm(f, x0)
            if best is None or c < best[0]:
                best = (c, x)
        print(f"skel {s}: shape={shape} cost={cost} unknowns={nf} final dims={dims} best residual^2={best[0]:.3e} ({time.time()-t0:.0f}s)", flush=True)
        best_all.append((best[0], cost))
    print("summary:", sorted(best_all)[:10])

if __name__ == "__main__" and not (len(sys.argv) > 1 and sys.argv[1] == 'fast'):
    main()

# ---------------------------------------------------------------- fast variable projection (Kaufman Jacobian)
def prep_fast(sk):
    sizes = [len(T) for T in sk.gens]
    off = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
    sk.off = off
    sk.fin = []
    for (nr, nc, E) in sk.final:
        Ea = np.array(E, dtype=np.int64).reshape(-1, 4)
        sk.fin.append((nr, nc, Ea[:, 0], Ea[:, 1], off[Ea[:, 2]] + Ea[:, 3]))

def flat(sk, alpha):
    return np.concatenate(sk.coefs(alpha))

def vp_eval(sk, alpha, free_idx, want_J=True, eps=1e-7):
    fl = flat(sk, alpha)
    if not np.isfinite(fl).all():
        return None, None
    D = []
    if want_J:
        for a in free_idx:
            al = alpha.copy(); al[a] += eps
            D.append((flat(sk, al) - fl) / eps)
    rs, Js = [], []
    for (nr, nc, R, C, F) in sk.fin:
        A = np.zeros((nr, nc), dtype=complex)
        np.add.at(A, (R, C), fl[F])
        t = np.zeros(nr, dtype=complex); t[0] = 1
        try:
            beta = np.linalg.lstsq(A, t, rcond=1e-10)[0]
        except np.linalg.LinAlgError:
            return None, None
        r = A @ beta - t
        rs.append(r)
        if want_J:
            V = np.zeros((nr, len(D)), dtype=complex)
            for k, d in enumerate(D):
                V[:, k] = np.bincount(R, weights=(d[F] * beta[C]).real, minlength=nr) + 1j * np.bincount(R, weights=(d[F] * beta[C]).imag, minlength=nr)
            X = np.linalg.lstsq(A, V, rcond=1e-10)[0]
            P = V - A @ X                      # (I - A A^+) dA_k beta
            Js.append(P)
    r = np.concatenate(rs)
    if not want_J:
        return r, None
    Jc = np.concatenate(Js, axis=0)            # complex derivative wrt alpha_k (holomorphic direction)
    # real parameters (Re a_k, Im a_k): d r/d Re = Jc, d r/d Im = i Jc ; stack real/imag parts of r
    Jr = np.concatenate([np.concatenate([Jc.real, (1j * Jc).real], axis=1),
                         np.concatenate([Jc.imag, (1j * Jc).imag], axis=1)], axis=0)
    return r, Jr

def vp_lm(sk, N, free_idx, x0, iters=80, tol=1e-20, log=None):
    nf = len(free_idx)
    def full(x):
        al = np.ones(N, dtype=complex); al[free_idx] = x[:nf] + 1j * x[nf:]; return al
    x = x0.copy(); lam = 1e-3
    r, J = vp_eval(sk, full(x), free_idx)
    if r is None:
        return x, np.inf
    rr = np.concatenate([r.real, r.imag]); cost = rr @ rr
    for it in range(iters):
        g = J.T @ rr; Hm = J.T @ J
        improved = False
        for _ in range(12):
            try:
                dx = -np.linalg.solve(Hm + lam * (np.diag(np.diag(Hm)) + 1e-12 * np.eye(len(x))), g)
            except np.linalg.LinAlgError:
                lam *= 10; continue
            xn = x + dx
            rn, _ = vp_eval(sk, full(xn), free_idx, want_J=False)
            if rn is None:
                lam *= 10; continue
            rrn = np.concatenate([rn.real, rn.imag]); cn = rrn @ rrn
            if cn < cost:
                x, cost = xn, cn; lam = max(lam / 5, 1e-12); improved = True; break
            lam *= 5
        if not improved or cost < tol:
            break
        r, J = vp_eval(sk, full(x), free_idx)
        rr = np.concatenate([r.real, r.imag])
    return x, cost

def main_fast():
    m, W, Hname, shape, menu, starts, seed = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5].split(','), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8])
    nskel = int(sys.argv[9]) if len(sys.argv) > 9 else 1
    cap = int(sys.argv[10]) if len(sys.argv) > 10 else 3000
    t0 = time.time()
    summary = []
    for s in range(nskel):
        sk, cost = build_skeleton(m, W, Hname, shape, menu, seed * 100 + s)
        if sk is None:
            print(f"skel {s}: build failed", flush=True); continue
        if not sk.setup_final(cap):
            print(f"skel {s}: cost {cost} closure cap", flush=True); continue
        prep_fast(sk)
        N = sk.nalpha
        free_idx = [a for a in range(N) if a not in sk.fixed]
        nf = len(free_idx)
        dims = [(nr, nc) for (nr, nc, E) in sk.final]
        best = (np.inf, None)
        for st in range(starts):
            x0 = np.array([sk.rng.gauss(0, 1) for _ in range(2 * nf)])
            if nf == 0:
                r, _ = vp_eval(sk, np.ones(N, dtype=complex), [], want_J=False)
                best = (float(np.sum(np.abs(r) ** 2)), x0); break
            x, c = vp_lm(sk, N, free_idx, x0)
            if c < best[0]:
                best = (c, x)
        print(f"skel {s}: cost={cost} moves/orbits={[1 + len(mv['images']) for mv in sk.moves]} unknowns={nf} dims={dims} best res^2={best[0]:.3e} ({time.time()-t0:.0f}s)", flush=True)
        summary.append((best[0], cost, s))
        if best[0] < 1e-14:
            np.save(f"sol_{m}_{W}_{Hname}_{''.join(shape)}_{seed}_{s}.npy", best[1])
    print("summary (res^2, cost, skel):", sorted(summary)[:12], flush=True)

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == 'fast':
    main_fast()
