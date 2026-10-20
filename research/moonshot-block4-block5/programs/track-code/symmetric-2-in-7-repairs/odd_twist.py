"""odd_twist.py -- mw2-absorb4 (resumed run): STRUCTURED odd-order twists for the absorbed shared move.

v-ms-fano's sigma (an order-3 address relabelling) separates two pieces with (1 + sigma)^-1 = (1 - sigma + sigma^2)/2.
In a Z2-symmetric program the place for it is the absorbed shared move: tau(u) = beta_tot u with beta_tot of ODD order,
so that the two arrays' parts are related by an odd-order letter instead of an involution.  The random twist sampler of
twist.beta_candidates (Paulis + random engine letters) almost never proposes such letters.  Here the candidates are ALL
address maps addr(P), P in GL(m,2) (and the same times Pauli-Z phases), filtered by the twist rules of twist.py:
norm g_W scalar (tau inverts beta up to a scalar when W = 2) and D beta D^-1 free; kept when ord(beta_tot) is odd.

Mode 'recipe': odd_twist.py recipe recipe.pkl move_index maxbeta starts [cap]
   twist the shared move of a recipe (e.g. the seed-66 limit) by every admissible odd twist; full exact test each.
Mode 'search': odd_twist.py search m W H shape menu seed nstarts steps cap
   fullsearch.py with the structured odd candidates (patched beta_candidates)."""
import sys, itertools, time, pickle
import twist
from twist import *

_orig = twist.beta_candidates

def gl_all(m):
    n = 1 << m
    out = []
    def rec(cols, span):
        if len(cols) == m:
            out.append(list(cols)); return
        for v in range(1, n):
            if v in span: continue
            rec(cols + [v], span | {x ^ v for x in span})
    rec([], {0})
    return out

def addr_letters(G, m, with_z=True):
    n = 1 << m
    S = Sym(m)
    out = []
    for P in gl_all(m):
        sig = [gl_apply(P, x, m) for x in range(n)]
        if sig == list(range(n)):
            continue
        A = G.mono(sig, [0] * n)
        out.append(A)
        if with_z:
            for b in range(1, n):
                Z = G.mono(list(range(n)), [2 * (bin(x & b).count('1') & 1) for x in range(n)])
                out.append(G.mul(Z, A))
    return out

_cache = {}
def beta_candidates_odd(G, R, Dm, rng, nrand=0):
    m, n = G.m, G.n
    key = (id(G), G.key(Dm))
    if key in _cache:
        return _cache[key]
    c = R.Hc[1][0]
    Dmi = np.conj(Dm).T
    beta1 = G.mul(c, Dm, np.conj(c).T, Dmi)
    S = TSkel.__new__(TSkel); S.G = G; S.R = R
    out = []; seen = set()
    if 'letters' not in _cache:
        _cache['letters'] = addr_letters(G, m, with_z=False)
    for bt in _cache['letters']:
        dbd = G.mul(Dm, bt, Dmi)
        if not G.is_monomial(dbd):
            continue
        if S.gammas(bt) is None:
            continue
        btot = G.mul(beta1, dbd)
        if is_scalar(btot):
            continue
        o = order_of(G, btot)
        if o % 2 == 0:
            continue
        k = G.key(bt)
        if k in seen: continue
        seen.add(k)
        out.append((bt, o, btot))
    _cache[key] = out
    return out

def main():
    mode = sys.argv[1]
    if mode == 'search':
        twist.beta_candidates = beta_candidates_odd
        import fullsearch
        fullsearch.beta_candidates = beta_candidates_odd
        sys.argv = [sys.argv[0]] + sys.argv[2:] + ['odd']
        fullsearch.main()
        return
    d = pickle.load(open(sys.argv[2], 'rb'))
    t_sh = int(sys.argv[3]); maxb = int(sys.argv[4]); starts = int(sys.argv[5]); cap = int(sys.argv[6]) if len(sys.argv) > 6 else 4000
    m, W, H = d.get('m', 4), d.get('W', 2), d.get('H', 'Z2swap')
    rec = d['recipe']
    base = build_T(m, W, H, rec, d['skseed'], {})
    G, R = base.G, base.R
    kind, summ, Dm = rec[t_sh]
    Dmat = G.D1 if Dm is None else Dm
    t0 = time.time()
    cands = beta_candidates_odd(G, R, Dmat, None)
    print(f"odd admissible twists for the shared move: {len(cands)}; orders of beta_tot {sorted(set(o for _, o, _ in cands))} ({time.time()-t0:.0f}s)", flush=True)
    rng = random.Random(5)
    res = []
    for k, (bt, o, btot) in enumerate(cands[:maxb]):
        rec2 = list(rec); rec2[t_sh] = ('T', summ, Dm)
        try:
            sk = build_T(m, W, H, rec2, d['skseed'], {t_sh: bt}, R)
        except AssertionError:
            continue
        if not sk.moves[t_sh]['absorbed']:
            print(f"twist {k}: orbit check failed"); continue
        if not sk.setup_final(cap):
            print(f"twist {k} order {o}: closure cap"); continue
        prep_fast(sk)
        N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
        best = np.inf
        for st in range(starts):
            if nf == 0:
                r, _ = vp_eval(sk, np.ones(N, dtype=complex), [], want_J=False); best = float(np.sum(np.abs(r) ** 2)); break
            x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
            x, cst = vp_lm(sk, N, free_idx, x0, iters=100)
            best = min(best, cst)
        dims = [(nr, nc) for (nr, nc, E) in sk.final]
        res.append((best, o, k))
        print(f"twist {k}: order(beta_tot) {o}: closure {dims}: full res^2 {best:.4f} ({time.time()-t0:.0f}s)", flush=True)
        if best < 1e-12:
            print("*** FULL HIT ***", flush=True)
            pickle.dump({"recipe": rec2, "twist": {t_sh: bt}, "skseed": d['skseed'], "m": m, "W": W, "H": H}, open(f"HIT_oddtw_{k}.pkl", "wb"))
    res.sort(); print("best", res[:10])

if __name__ == "__main__":
    main()
