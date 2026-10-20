# fc-reopt round 2: #3815's reverse-time EDF placement (construction/build_fill.py, its loop copied here) made
# TOLERANT, applied as a text patch to an in-memory copy of the hub file (the hub file is not changed).
#
# Round 1 found that the class-level optimum of the full-cut Hall LP (1.7823049e-3) cannot be placed EXACTLY: at a few
# calls the one-role-per-bin-per-call rule leaves some roles of one width without a planned slot. Here such a role is
# not a stop: it goes to a SPARE place, i.e. any bin open at that call that is not a head bin (pair class 4,
# children 8/10, stays empty), has not been touched at this call, and whose final capacity (final widths, every planned
# slot counted as filled, plus earlier spare roles) still has room; at most one widened (width 4 -> 5) role per bin
# (#3945's verify.py 'changed bin increment'). Among those, the bin whose full-cut cost phi_24 rises least.
# A planned slot whose class closes before it is used is simply left empty (more idle in that bin).
# Every placement still respects the window, the per-call distinctness and the capacity; the hub's auditors check
# the result as before. The price is whatever the realised fill gives (evaluated by our tools afterwards).
START = "    # À rebours : une classe devient disponible à t = hi et doit être servie avant t = lo.\n"
END = "    assert demand.sum() == 0\n    rows = [None] * (H.N * 408)\n"
BLOCK = '''    import os as _os, json as _json
    _tau = 1.0 - float(_os.environ.get('FC_A', '0.0017807'))
    _C = 24
    def _phi(l):
        if l <= 0: return 0.0
        q, m = divmod(int(l), _C)
        return q * _C ** _tau + (m ** _tau if m else 0.0)
    _PHI = np.array([_phi(l) for l in range(0, 65)])
    _capk = np.array([int(lp.cap[k]) for k in range(K)], dtype=np.int64)
    _W4 = _os.environ.get('FC_WIDEN4', '1') == '1'                              # line 3945: yes; line 2979: no
    _wnew = np.array([w + (1 if (w == 4 and _W4) else 0) for w in W], dtype=np.int64)     # #3945's offset step: 4 -> 5
    _plan_new = (free.astype(np.int64) * _wnew[None, None, :]).sum(2)            # final load if every slot is used
    _has4 = free[:, :, 0] > 0
    _extra = np.zeros((K, H.BINS), dtype=np.int64)
    _extra4 = np.zeros((K, H.BINS), dtype=bool)
    _nothead = (np.arange(K) % 9) != 4
    _log = dict(dropped={}, spare={}, spare_ticks=[], subst={})
    def _drop(k, wi, t, why):
        n = int(demand[k, wi])
        if n:
            _log['dropped'][str(W[wi])] = _log['dropped'].get(str(W[wi]), 0) + n
            demand[k, wi] = 0
    def _spare(t, wi, hs):
        """place roles hs (width W[wi]) at call t into spare places (see the header)."""
        if not hs: return
        wn = int(_wnew[wi])
        ks = np.where((lo_k < t) & (t <= hi_k) & _nothead)[0]
        sp = _capk[ks, None] - _plan_new[ks] - _extra[ks]
        ok = (sp >= wn) & (used_by[ks] != t)
        if W[wi] == 4 and _W4:
            ok &= ~_has4[ks] & ~_extra4[ks]
        spc = np.clip(sp, 0, 64)
        cost = np.where(ok, _PHI[np.clip(spc - wn, 0, 64)] - _PHI[spc], np.inf).ravel()
        n = min(len(hs), int(np.isfinite(cost).sum()))
        rest = hs[n:]; hs = hs[:n]
        if rest:
            _subst(t, wi, rest, ks)
        if not n: return
        top = np.argpartition(cost, n - 1)[:n]
        for f, h in zip(top, hs):
            k = int(ks[f // H.BINS]); b = int(f % H.BINS)
            used_by[k, b] = t
            _extra[k, b] += wn
            if W[wi] == 4: _extra4[k, b] = True
            placed[k, b].append((t, h, W[wi]))
        _log['spare'][str(W[wi])] = _log['spare'].get(str(W[wi]), 0) + n
        _log['spare_ticks'].append([int(t), int(W[wi]), int(n)])
    def _subst(t, wi, hs, ks):
        """no spare place left at this call: the role takes a planned slot of a LARGER width w' in an open bin not touched at
        this call (the bin's plan changes from w' to w, its final idle grows by w' - w); the class then has one w' slot
        fewer, so one w' role finds no planned slot later and goes the same way (spare place or substitution)."""
        wn = int(_wnew[wi])
        for h in hs:
            done = False
            for wj in range(wi + 1, 6):
                ok = (free[ks, :, wj] > 0) & (used_by[ks] != t)
                if W[wi] == 4 and _W4:
                    ok &= ~_has4[ks] & ~_extra4[ks]
                if not ok.any(): continue
                sp = _capk[ks, None] - _plan_new[ks] - _extra[ks]
                d = int(_wnew[wj]) - wn
                spc = np.clip(sp, 0, 64)
                cost = np.where(ok, _PHI[np.clip(spc + d, 0, 64)] - _PHI[spc] + 1e-3 * jitter[ks], np.inf).ravel()
                f = int(np.argmin(cost)); k = int(ks[f // H.BINS]); b = int(f % H.BINS)
                free[k, b, wj] -= 1; ftot[k, b] -= 1; demand[k, wj] -= 1
                _plan_new[k, b] -= d
                if W[wi] == 4: _has4[k, b] = True
                used_by[k, b] = t
                placed[k, b].append((t, h, W[wi]))
                key = '%d>%d' % (W[wj], W[wi]); _log['subst'][key] = _log['subst'].get(key, 0) + 1
                done = True
                break
            if not done:
                raise SystemExit(f'pas de place de secours ni de substitution : t={t} largeur {W[wi]}')
    for t in range(H.N - 1, -1, -1):
        for k in closing.get(t, ()):
            if lo_k[k] < 0:
                in_pool[pool_pos[k]] = True
            else:
                for wi in range(6):
                    if demand[k, wi]:
                        heapq.heappush(heaps[wi], (-int(lo_k[k]), int(hi_k[k]), k))
        for wi in range(6):
            hs = list(helpers_of[wi])
            hp = heaps[wi]
            aside = []
            while hs and hp:
                neglo, _, k = hp[0]
                if demand[k, wi] == 0:
                    heapq.heappop(hp)
                    continue
                if -neglo >= t:
                    _drop(k, wi, t, 'deadline')
                    heapq.heappop(hp)
                    continue
                elig = np.where((free[k, :, wi] > 0) & (used_by[k] != t))[0]
                n = min(len(hs), int(demand[k, wi]), len(elig))
                if n == 0:
                    aside.append(heapq.heappop(hp))
                    continue
                pick = elig[np.argsort(-(ftot[k, elig] + jitter[k, elig]))[:n]]
                take(k, pick, t, wi, hs[:n])
                hs = hs[n:]
                if demand[k, wi] and hs:
                    aside.append(heapq.heappop(hp))
            for item in aside:
                heapq.heappush(hp, item)
            if hs:
                ks = pool[in_pool]
                fr = free[ks, :, wi]
                ok = (fr > 0) & (used_by[ks] != t)
                score = np.where(ok, ftot[ks] + jitter[ks], -1.0)
                flat = score.ravel()
                n = min(len(hs), int((flat >= 0).sum()))
                if n:
                    top = np.argpartition(-flat, n - 1)[:n]
                    by_class = defaultdict(list)
                    for f in top:
                        by_class[int(ks[f // H.BINS])].append(int(f % H.BINS))
                    pos = 0
                    for k, bins in by_class.items():
                        take(k, bins, t, wi, hs[pos:pos + len(bins)])
                        pos += len(bins)
                    hs = hs[n:]
            if hs:
                _spare(t, wi, hs)
        for k in np.where(lo_k == t)[0]:
            for wi in range(6):
                _drop(k, wi, t, 'closed')
    for k in range(K):
        for wi in range(6):
            _drop(k, wi, -1, 'end')
    _log['spare_total'] = int(sum(_log['spare'].values())); _log['dropped_total'] = int(sum(_log['dropped'].values()))
    print('TOLERANT PLACEMENT: spare roles %s, substitutions %s, dropped planned slots %s' % (_log['spare'], _log['subst'], _log['dropped']), flush=True)
    _p = _os.environ.get('FC_FALLBACK_LOG')
    if _p:
        open(_p, 'w').write(_json.dumps(_log))
'''


def patch(src):
    assert src.count(START) == 1 and src.count(END) == 1
    i = src.index(START); j = src.index(END)
    return src[:i] + BLOCK + src[j:]
