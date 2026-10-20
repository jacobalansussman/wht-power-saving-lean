# fc-reopt: an alternative per-tick placement for #3815's build_fill.build(), applied as a text patch to an in-memory
# copy of the hub file (the hub file is not changed). Everything before and after the placement loop is #3815's:
# the class-level pattern lists from the Hall LP/MIP, the bins, the helper widths, the FILL/PROFILE writer.
#
# #3815's loop places roles tick by tick in reverse time, width by width, earliest-deadline class first, with the pool
# classes as the residual. That greedy can run short near the end of the clock when the full-cut mix concentrates a
# width in few bins. This loop places ALL roles of a tick at once by a min-cost transportation problem (scipy milp;
# the matrix is totally unimodular, so it is solved at the root):
#   source -> width w (capacity: roles of width w in the call) -> bin group g (bins of one class with the same free-slot
#   vector; capacity = number of bins, one role per bin per call) -> sink,
# cost of a bin = its slack = (ticks left to the class's window start, t - lo) - (its free slots); a bin with slack 0
# must take a role now. Any infeasibility stops the build (it never places a role outside its window or twice in a bin
# at one call); the result is then checked by the hub's own auditors as before.
START = "    # À rebours : une classe devient disponible à t = hi et doit être servie avant t = lo.\n"
END = "    assert demand.sum() == 0\n    rows = [None] * (H.N * 408)\n"
BLOCK = '''    from scipy.optimize import milp as _milp, LinearConstraint as _LC, Bounds as _Bd
    from scipy.sparse import coo_matrix as _coo
    _need = np.array([len(helpers_of[wi]) for wi in range(6)], dtype=np.int64)
    _E = int(__import__('os').environ.get('FC_ENDGAME', '0'))
    if _E:
        # exact end game: the last _E ticks as one integer program at group level (bins of one class with the same
        # pattern are interchangeable), then a balanced realisation per group (each bin at most one role per tick and
        # at most its slots of each width; König's theorem guarantees one exists).
        N_ = H.N; TT = list(range(N_ - _E, N_))
        gk, gv, gb = [], [], []
        for k in np.where((hi_k >= N_ - _E) & (lo_k < N_ - 1))[0]:
            vecs, invb = np.unique(free[k].astype(np.int64), axis=0, return_inverse=True)
            for gi_, v in enumerate(vecs):
                if v.sum() == 0: continue
                gk.append(int(k)); gv.append(v); gb.append(np.where(invb.ravel() == gi_)[0])
        var = []                                    # (g, t, wi)
        for g, k in enumerate(gk):
            for t in TT:
                if lo_k[k] < t <= hi_k[k]:
                    for wi in range(6):
                        if gv[g][wi] > 0: var.append((g, t, wi))
        nvar = len(var); vidx = {x: i for i, x in enumerate(var)}
        rows, lor, hir = [], [], []
        def addrow(ix, lo_, hi_):
            rows.append(ix); lor.append(lo_); hir.append(hi_)
        for t in TT:
            for wi in range(6):
                addrow([vidx[x] for x in var if x[1] == t and x[2] == wi], float(_need[wi]), float(_need[wi]))
        byg = {}
        for i, (g, t, wi) in enumerate(var): byg.setdefault(g, []).append(i)
        for g, ix in byg.items():
            n_ = len(gb[g]); k = gk[g]
            for t in TT:
                sub = [i for i in ix if var[i][1] == t]
                if sub: addrow(sub, 0.0, float(n_))
            for wi in range(6):
                sub = [i for i in ix if var[i][2] == wi]
                if not sub: continue
                full = float(n_ * gv[g][wi])
                if lo_k[k] >= N_ - _E: addrow(sub, full, full)          # the class closes inside the end game
                else: addrow(sub, 0.0, full)
            if lo_k[k] < N_ - _E:                                       # what is left must fit in the earlier ticks
                left_ticks = (N_ - _E - 1) - int(lo_k[k])
                need_lo = n_ * (int(gv[g].sum()) - left_ticks)
                if need_lo > 0: addrow(ix, float(need_lo), float(n_ * int(gv[g].sum())))
        r_, c_ = [], []
        for ri, ix in enumerate(rows):
            r_ += [ri] * len(ix); c_ += ix
        Aeg = _coo((np.ones(len(r_)), (r_, c_)), shape=(len(rows), nvar)).tocsr()
        res = _milp(np.zeros(nvar), constraints=_LC(Aeg, np.array(lor), np.array(hir)), integrality=np.ones(nvar),
                    bounds=_Bd(np.zeros(nvar), np.full(nvar, np.inf)), options={'time_limit': 1800})
        print('END GAME: %d ticks, %d groups, %d variables, status %s' % (_E, len(gk), nvar, res.message), flush=True)
        if res.status != 0 or res.x is None:
            raise SystemExit(f'fin de partie infaisable : {res.message}')
        Z = np.rint(res.x).astype(np.int64)
        def _colour(edges, n):
            """bipartite edge colouring with n colours (max degree <= n): alternating-path recolouring."""
            L, R = {}, {}
            col = [None] * len(edges)
            for e, (u, v) in enumerate(edges):
                Lu = L.setdefault(u, {}); Rv = R.setdefault(v, {})
                a_ = next(c for c in range(n) if c not in Lu)
                b_ = next(c for c in range(n) if c not in Rv)
                if a_ not in Rv:
                    c = a_
                else:
                    # flip the a/b path starting at v (it alternates a, b, a, ... and cannot reach u)
                    path = []; x = ('R', v); cc = a_
                    while True:
                        side, node = x
                        tab = R if side == 'R' else L
                        if cc not in tab.get(node, {}): break
                        f = tab[node][cc]; path.append((f, cc))
                        uu, vv = edges[f]
                        x = ('L', uu) if side == 'R' else ('R', vv)
                        cc = b_ if cc == a_ else a_
                    for f, cc0 in path:
                        uu, vv = edges[f]
                        del L[uu][cc0]; del R[vv][cc0]
                    for f, cc0 in path:
                        uu, vv = edges[f]; nc = b_ if cc0 == a_ else a_
                        L[uu][nc] = f; R[vv][nc] = f; col[f] = nc
                    c = a_
                assert c not in L[u] and c not in R[v]
                L[u][c] = e; R[v][c] = e; col[e] = c
            return col
        place = {}                                  # (t, wi) -> list of (k, bin)
        for g in range(len(gk)):
            k = gk[g]; n_ = len(gb[g]); q = gv[g]
            edges = []; meta = []
            copyfill = {wi: 0 for wi in range(6)}
            for t in sorted(TT, reverse=True):
                for wi in range(6):
                    if (g, t, wi) not in vidx: continue
                    for _ in range(int(Z[vidx[(g, t, wi)]])):
                        cp = copyfill[wi] % int(q[wi]); copyfill[wi] += 1     # round robin over the width's slot copies
                        edges.append((t, (wi, cp))); meta.append((t, wi))
            if not edges: continue
            col = _colour(edges, n_)
            for (t, wi), c in zip(meta, col):
                place.setdefault((t, wi), []).append((k, int(gb[g][c])))
        for t in sorted(TT, reverse=True):
            for wi in range(6):
                lst = place.get((t, wi), [])
                assert len(lst) == _need[wi], ('end game count', t, wi, len(lst))
                hs = list(helpers_of[wi]); pos = 0
                byk = {}
                for k, b in lst: byk.setdefault(k, []).append(b)
                for k, bl in byk.items():
                    take(k, np.array(bl), t, wi, hs[pos:pos + len(bl)]); pos += len(bl)
        print('END GAME placed', flush=True)
    for t in range(H.N - 1 - _E, -1, -1):
        ks = np.where((lo_k < t) & (t <= hi_k))[0]
        fr = free[ks].astype(np.int64)
        tot = fr.sum(2)
        rem = (t - lo_k[ks]).astype(np.int64)
        if (tot > rem[:, None]).any():
            raise SystemExit(f'case en retard (flot) : t={t}')
        ki, bi = np.nonzero(tot > 0)
        keys = np.concatenate([ks[ki][:, None], fr[ki, bi]], axis=1)
        uniq, inv, cnt = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
        inv = inv.ravel()
        ng = len(uniq)
        slack = (t - lo_k[uniq[:, 0]]) - uniq[:, 1:].sum(1)
        gw = np.nonzero(uniq[:, 1:] > 0)            # arcs (group, width) with a free slot of that width
        na = len(gw[0])
        # transportation: x[a] roles of width gw[1][a] into bins of group gw[0][a]; rows: widths (=), groups (<=)
        r = np.concatenate([gw[1], 6 + gw[0]]); c = np.concatenate([np.arange(na), np.arange(na)])
        A = _coo((np.ones(2 * na), (r, c)), shape=(6 + ng, na)).tocsr()
        forced = (slack <= 0).astype(float) * cnt          # a bin with no slack must take a role at this call
        lo_r = np.concatenate([_need.astype(float), forced]); hi_r = np.concatenate([_need.astype(float), cnt.astype(float)])
        _ca, _cb = [float(x) for x in __import__('os').environ.get('FC_COST', '1,0').split(',')]
        _cost = _ca * slack[gw[0]].astype(float) - _cb * uniq[gw[0], 1 + gw[1]].astype(float)
        res = _milp(_cost, constraints=_LC(A, lo_r, hi_r), integrality=np.ones(na),
                    bounds=_Bd(np.zeros(na), cnt[gw[0]].astype(float)))
        if res.status != 0:
            capw = [int(cnt[uniq[:, 1 + wi] > 0].sum()) for wi in range(6)]
            print('PLACER DIAG t', t, 'need', _need.tolist(), 'bins with a slot of each width', capw, 'bins with any slot', int(cnt.sum()),
                  'classes open', len(ks), flush=True)
            import itertools as _it
            for r_ in range(2, 7):
                for S_ in _it.combinations(range(6), r_):
                    m_ = (uniq[:, [1 + x for x in S_]] > 0).any(1)
                    if int(cnt[m_].sum()) < int(_need[list(S_)].sum()):
                        print('   Hall violation at this tick: widths', [int(H.W[x]) for x in S_], 'bins', int(cnt[m_].sum()), 'need', int(_need[list(S_)].sum()), flush=True)
            raise SystemExit(f'flot infaisable : t={t} statut {res.status} {res.message}')
        flows = np.rint(res.x).astype(np.int64)
        order_g = np.argsort(inv, kind='stable')
        bounds = np.searchsorted(inv[order_g], np.arange(ng + 1))
        hs_left = [list(helpers_of[wi]) for wi in range(6)]
        usedcnt = np.zeros(ng, dtype=np.int64)
        for a_, f in enumerate(flows):
            if not f: continue
            g = int(gw[0][a_]); wi = int(gw[1][a_]); k = int(uniq[g, 0])
            mem = order_g[bounds[g]:bounds[g + 1]]
            pick = bi[mem[usedcnt[g]:usedcnt[g] + f]]
            usedcnt[g] += f
            hs = hs_left[wi][:f]; hs_left[wi] = hs_left[wi][f:]
            take(k, pick, t, wi, hs)
        assert all(not x for x in hs_left), ('roles left', t)
        for k in np.where(lo_k == t)[0]:
            if demand[k].any():
                raise SystemExit(f'classe {k} fermée incomplète à t={t} : {demand[k]}')
'''


def patch(src):
    i = src.index(START); j = src.index(END)
    assert src.count(START) == 1 and src.count(END) == 1
    return src[:i] + BLOCK + src[j:]
