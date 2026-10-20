"""p-own (wht45/lever, 19 Oct 2026): independent audit of a p-own host table (second code path, numpy, vectorised).
Reads ONLY the table, plus the hub's PROFILE.json / ORDER.json / DISPATCH.json as REFERENCE values to compare with.
It does not reuse own_to_host.py (different code, different method: whole-array arithmetic, sort-based prefixes,
the window predicate in BOTH forms: interval [lo, hi] from the table AND p-a3's progT form recomputed from time).

Checks (all must hold): shape; global order (time a permutation, order_k its inverse, increasing along each replica);
windows recomputed from time == table windows; every entry inside its window (interval form and progT form:
KZ[c] <= #{l : time[9 rho_h + l] < time[k]} <= KZ[c] + KK[c] - 1); one host per slot per call (rows injective);
one family per call (owner, hub pair, side, port distinct in every row); windows of kinds sharing a hub pair are
disjoint in every replica; no banks (every entry < NH, nb = 0); capacities and remainders (rem + load == rank, rem >= 0);
compact prefixes (in every bin, offsets in time order are 0, w1, w1 + w2, ..., and end at most the rank);
remainder histogram == {3: 23328, 5: 4320, 48: 137952}; full profile (sum of the providers' local profiles over the
calls + remainders) == PROFILE.json; mass = 65088 T; cross-replica roles == 1,161,423 (bB-data copy for design B; p-own's original: 1,152,628); inventory per provider
{5: 28, 6: 18, 7: 53, 8: 96, 9: 213}; dispatch == DISPATCH.json (1,672 extra, 1,568 gauge44) and order == ORDER.json.
Negatives (each must be REJECTED, and by the check it targets): remainder off by one, window violation, provider
flipped, family repeated, offset hole.
usage: python3 -I check_own.py HOSTTAB.json.gz PROFILE.json ORDER.json DISPATCH.json OUT.json"""
import sys, json, gzip, hashlib, collections, time as _time
import numpy as np
t0 = _time.time()
tabf, proff, ordf, dispf, outf = sys.argv[1:6]
d = json.loads(gzip.decompress(open(tabf, 'rb').read()))
PROF = {int(a): b for a, b in json.load(open(proff)).items()}
ORD = json.load(open(ordf)); DSP = json.load(open(dispf))
T, R, NH = d['T'], d['R'], d['NH']; NC = 9 * T
RANK = np.array(d['ranks']); KZ = np.array(d['kz']); KK = np.array(d['kk']); HP = np.array(d['hpair'])
PROV = {'extra': 0, 'gauge44': 1}
assert d['providers'] == {'0': 'extra', '1': 'gauge44'}
LP = {int(p): {int(a): b for a, b in v.items()} for p, v in d['local_profile'].items()}
H = np.array(d['hrows'], dtype=np.int64); O = np.array(d['offs'], dtype=np.int64)
REM = np.array(d['hrem'], dtype=np.int64); TM = np.array(d['time'], dtype=np.int64)
DISP = np.array(d['disp'], dtype=np.int64); WIN = np.array(d['windows'], dtype=np.int64)
RES = np.array([d['res']['0'], d['res']['1']], dtype=np.int64); OK_ = np.array(d['order_k'], dtype=np.int64)
ref_disp = np.array([PROV[tag] for s, r, tag in sorted(DSP, key=lambda x: 9 * x[1] + x[0])])
ref_k = np.array([9 * r + s for s, r in ORD])
def check(H, O, REM, TM, DISP, WIN, RES):
    ok = {}
    ok['shape'] = (H.shape == (NC, R) and O.shape == (NC, R) and REM.shape == (NH,) and TM.shape == (NC,)
                   and DISP.shape == (NC,) and WIN.shape == (9 * T, 2) and RES.shape == (2, R) and d['nb'] == 0
                   and NH == 2016 * T)
    ok['no_banks'] = bool(((H >= 0) & (H < NH)).all())
    ok['order'] = (bool((np.sort(TM) == np.arange(NC)).all()) and bool((TM[OK_] == np.arange(NC)).all())
                   and bool((np.diff(TM.reshape(T, 9), axis=1) > 0).all()))
    ok['order_eq_ORDER_json'] = bool((OK_ == ref_k).all())
    tr = TM.reshape(T, 9)
    lo = np.where(KZ[None, :] == 0, 0, tr[:, np.maximum(KZ - 1, 0)] + 1); hi = tr[:, KZ + KK - 1]
    ok['windows_recomputed'] = bool((WIN[:, 0] == lo.ravel()).all() and (WIN[:, 1] == hi.ravel()).all())
    rho = H // 2016; kind = H % 2016 // 224; side = H % 224 // 112; port = H % 112
    tk = TM[:, None]
    w = WIN[9 * rho + kind]
    ok['window_interval'] = bool(((w[..., 0] <= tk) & (tk <= w[..., 1])).all())
    prog = (tr[rho] < tk[..., None]).sum(axis=2)
    ok['window_progT'] = bool(((KZ[kind] <= prog) & (prog <= KZ[kind] + KK[kind] - 1)).all())
    Hs = np.sort(H, axis=1); ok['one_host_per_slot'] = bool((Hs[:, 1:] != Hs[:, :-1]).all())
    fam = ((rho * 4 + HP[kind]) * 2 + side) * 112 + port
    Fs = np.sort(fam, axis=1); ok['one_family_per_call'] = bool((Fs[:, 1:] != Fs[:, :-1]).all())
    dis = True
    for pr in range(4):
        ks = [c for c in range(9) if HP[c] == pr]
        for a in ks:
            for b in ks:
                if a < b:
                    dis &= bool((np.maximum(WIN[a::9, 0], WIN[b::9, 0]) > np.minimum(WIN[a::9, 1], WIN[b::9, 1])).all())
    ok['family_windows_disjoint'] = dis
    W = RES[DISP][:, :]          # (NC, R): width of slot q of call k
    load = np.bincount(H.ravel(), weights=W.ravel(), minlength=NH).astype(np.int64)
    rk = RANK[np.arange(NH) % 2016 // 224]
    ok['rem_plus_load_eq_rank'] = bool((REM + load == rk).all()) and bool((REM >= 0).all())
    hf, wf, of_, tf = H.ravel(), W.ravel(), O.ravel(), np.repeat(TM, R)
    srt = np.lexsort((tf, hf)); hs, ws, os_ = hf[srt], wf[srt], of_[srt]
    cs = np.cumsum(ws) - ws                      # exclusive prefix over the whole sorted list
    start = np.ones(len(hs), dtype=bool); start[1:] = hs[1:] != hs[:-1]
    base = np.maximum.accumulate(np.where(start, np.arange(len(hs)), 0))
    pref = cs - cs[base]
    ok['compact_prefix'] = bool((os_ == pref).all()) and bool((os_ + ws <= RANK[hs % 2016 // 224]).all())
    rh = collections.Counter(int(x) for x in REM[REM > 0])
    ok['remainders_exact'] = rh == {3: 23328, 5: 4320, 48: 137952}
    prof = collections.Counter()
    for p in (0, 1):
        n = int((DISP == p).sum())
        for a, b in LP[p].items(): prof[a] += n * b
    prof.update(rh)
    ok['profile_eq_PROFILE_json'] = dict(prof) == PROF
    ok['mass'] = sum(a * b for a, b in prof.items()) == 65088 * T
    cross = int((rho != (np.arange(NC) // 9)[:, None]).sum())
    ok['cross_replica_1161423'] = cross == 1161423   # design B (#3945's own fill: 1,152,628)
    ok['inventory'] = all(collections.Counter(RES[p].tolist()) == {5: 28, 6: 18, 7: 53, 8: 96, 9: 213} for p in (0, 1))
    ok['dispatch_eq_DISPATCH_json'] = bool((DISP == ref_disp).all()) and (int((DISP == 0).sum()), int((DISP == 1).sum())) == (1672, 1568)
    info = dict(cross=cross, remainders=dict(sorted(rh.items())), profile=dict(sorted(prof.items())))
    return ok, info
ok, info = check(H, O, REM, TM, DISP, WIN, RES)
out = dict(table=tabf.split('/')[-1], table_sha256=hashlib.sha256(open(tabf, 'rb').read()).hexdigest(),
           checks=ok, all_checks_pass=all(ok.values()), **info,
           cross_share=round(info['cross'] / (NC * R), 6), provider_slot_diffs=[int(q) for q in np.nonzero(RES[0] != RES[1])[0]],
           negatives={})
def neg(name, target, fn, undo):
    fn(); okn, _ = check(H, O, REM, TM, DISP, WIN, RES); undo()
    failed = sorted(k for k, v in okn.items() if not v)
    out['negatives'][name] = dict(rejected=bool(failed), failed_checks=failed, target=target, target_failed=target in failed)
# 1. one remainder off by one
i0 = int(np.nonzero(REM > 0)[0][0])
neg('rem_off_by_one', 'rem_plus_load_eq_rank', lambda: REM.__setitem__(i0, REM[i0] - 1), lambda: REM.__setitem__(i0, REM[i0] + 1))
# 2. window violation: entry (k, q) moved to the same bin of another replica whose window misses the call, no clash
def find_move():
    for k in range(NC):
        for q in range(R):
            e = int(H[k, q]); c = e % 2016 // 224
            for o in range(T):
                lo, hi = WIN[9 * o + c]
                e2 = 2016 * o + e % 2016
                if not (lo <= TM[k] <= hi) and e2 not in H[k]: return k, q, e, e2
k2, q2, e2a, e2b = find_move()
neg('window_violation', 'window_interval', lambda: H.__setitem__((k2, q2), e2b), lambda: H.__setitem__((k2, q2), e2a))
# 3. provider of one call flipped
k3 = 0
neg('provider_flipped', 'dispatch_eq_DISPATCH_json', lambda: DISP.__setitem__(k3, 1 - DISP[k3]), lambda: DISP.__setitem__(k3, 1 - DISP[k3]))
# 4. family repeated: slot 1 of call k4 put on the host of slot 0 (the hub's own control: row[1] host := row[0] host)
k4 = int(OK_[0]); s4 = (int(H[k4, 1]), int(O[k4, 1]))
def f4(): H[k4, 1] = H[k4, 0]; O[k4, 1] = O[k4, 0]
def u4(): H[k4, 1], O[k4, 1] = s4
neg('family_repeated', 'one_family_per_call', f4, u4)
# 5. offset hole
neg('offset_hole', 'compact_prefix', lambda: O.__setitem__((k4, 0), O[k4, 0] + 1), lambda: O.__setitem__((k4, 0), O[k4, 0] - 1))
ok2, _ = check(H, O, REM, TM, DISP, WIN, RES)
out['restored_table_passes_again'] = all(ok2.values())
out['all_ok'] = (all(ok.values()) and out['restored_table_passes_again']
                 and all(v['rejected'] and v['target_failed'] for v in out['negatives'].values()))
out['seconds'] = round(_time.time() - t0, 1)
print(json.dumps(out)); json.dump(out, open(outf, 'w'), indent=1)
print('ALL OK' if out['all_ok'] else 'SOME CHECK FAILED')
