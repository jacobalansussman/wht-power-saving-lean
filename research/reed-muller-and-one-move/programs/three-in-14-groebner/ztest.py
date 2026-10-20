"""ztest.py: my own implementation of the clique test (C) and rm-theory's zero-block test 8.7 (Z), with
certificates, for a multi-target configuration (targets T_l, moves z_i; any time order).

Facts used (rm-theory 8.5, 8.7; re-derived here only as far as needed to state the certificate):
- a form Q is tight when |D_Q| = r_Q >= 1; then E[D_Q, D_Q] = 0;
- two moves are *linked* when one tight support contains both (a move is linked to itself when it lies in a tight
  support); a set Y of pairwise linked moves (self-links included) has E[Y, Y] = 0;
- C: |Y| <= W;   Z: every form F with Y inside D_F has |D_F| + r_F >= 2|Y|.
A certificate is (Y; for each pair in Y a tight form containing both; a form F with |D_F| + r_F < 2|Y|).
The certificate's forms are the subsystem that gb.py hands to Groebner.
"""
import itertools
from f2q import popc


def links(cfg):
    tight = cfg.tight_supports()            # D -> (Q, p)
    M = cfg.M
    wit = {}                                # (i, j) with i <= j -> (size of D, Q, D)
    for D, (Q, p) in tight.items():
        mv = [i for i in range(M) if D >> i & 1]
        s = popc(D)
        for a in range(len(mv)):
            for b in range(a, len(mv)):
                key = (mv[a], mv[b])
                if key not in wit or s < wit[key][0]:
                    wit[key] = (s, Q, D)
    return tight, wit


def cliques(cfg, wit, maxsize=None):
    """all cliques Y (pairwise linked incl. self) of size >= 1, as sorted tuples."""
    M = cfg.M
    live = [i for i in range(M) if (i, i) in wit]
    adj = {i: {j for j in live if j != i and (min(i, j), max(i, j)) in wit} for i in live}
    out = []
    def rec(cur, cand):
        out.append(tuple(cur))
        if maxsize and len(cur) >= maxsize: return
        for v in sorted(cand):
            if cur and v < cur[-1]: continue
            rec(cur + [v], cand & adj[v])
    for v in live:
        rec([v], adj[v])
    return sorted(set(out))


def zero_block(cfg, want_all=False):
    """returns dict(clique=max clique size, C=bool, Z=bool, cert=..., nfail=number of (Y, F) failures)."""
    tight, wit = links(cfg)
    cl = cliques(cfg, wit)
    cmax = max((len(y) for y in cl), default=0)
    res = dict(ntight_supports=len(tight), clique=cmax, C=cmax <= cfg.W)
    if not res["C"]:
        Y = max(cl, key=len)
        res.update(Z=None, cert=dict(kind="C", Y=list(Y)))
        return res
    # maximal cliques only matter for Z; test each against every form whose support contains it
    big = [y for y in cl if len(y) >= 1]
    best = None; nfail = 0
    ymask = [(sum(1 << i for i in y), y) for y in big]
    for Q in range(len(cfg.Dmask)):
        d = cfg.Dmask[Q]; r = popc(cfg.pmask[Q]); sd = popc(d)
        if sd + r >= 2 * cmax: continue          # cannot fail for any clique
        for ym, y in ymask:
            if d & ym == ym and sd + r < 2 * len(y):
                nfail += 1
                # prefer the smallest support, then the fewest distinct moves in the certificate
                pairs = [(y[a], y[b]) for a in range(len(y)) for b in range(a, len(y))]
                U = d
                for pr in pairs: U |= wit[pr][2]
                key = (sd, popc(U), r)
                if best is None or key < best[0]:
                    best = (key, dict(kind="Z", Y=list(y), F=Q, DF=d, rF=r,
                                      pairs=[(i, j, wit[(i, j)][1], wit[(i, j)][2]) for (i, j) in pairs], U=U))
    res["Z"] = (nfail == 0); res["nfail"] = nfail
    if best: res["cert"] = best[1]
    return res


def check_cert(cfg, cert):
    """re-check a Z certificate from scratch (forms recomputed from their bitmasks)."""
    F = cfg.F
    zm = [F.dir_mask(z) for z in cfg.moves]; tm = [F.sym_mask(T) for T in cfg.Ts]
    def supp(Q): return sum(1 << i for i, x in enumerate(zm) if popc(Q & x) & 1)
    def r(Q): return sum(popc(Q & x) & 1 for x in tm)
    Y = cert["Y"]
    for (i, j, Q, D) in cert["pairs"]:
        assert supp(Q) == D and (D >> i & 1) and (D >> j & 1)
        assert popc(D) == r(Q) >= 1, "witness form not tight"
    covered = {(i, j) for (i, j, _, _) in cert["pairs"]}
    assert all((min(a, b), max(a, b)) in covered for a in Y for b in Y)
    DF = supp(cert["F"]); assert DF == cert["DF"]
    assert all(DF >> y & 1 for y in Y)
    assert popc(DF) + r(cert["F"]) < 2 * len(Y)
    return True
