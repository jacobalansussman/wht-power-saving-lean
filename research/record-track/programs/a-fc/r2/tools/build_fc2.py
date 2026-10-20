# fc-reopt round 2: copy of round 1's tools/build_fc.py with one more placer, 'tolerant' (r2/tools/fallback_patch.py).
# fc-reopt (round wht45/fullcut): an integral literal fill built with the hub's construction, priced for a separated
# cut of paid-rank cap C (C = 24: the full cut). Follows o-model's P3b (#3945 line) and o-param's build_any.py
# (#2979 line); only the objective changes (hall_fc.py):
#   1. #3815's per-tick Hall LP with the cut objective (hall_fc.load), linearised at a = --a (Dinkelbach point);
#   2. a MIP on the generated columns (y integral), re-separating Hall rows on the MIP point (o-model's hall_mip.py);
#   3. #3815's reverse-time placement construction/build_fill.py, UNCHANGED (for #2979's words through o-param's stub
#      package whose NATIVE frames encode #2979's helper widths);
#   4. line 3945: PR3734's owner_audit.py (copied by build_fill into the output) and #3945's offset step
#      (o-model's offsets3945.py re-implementation); line 2979: the PR2839 owner engine vendored by #2979 with its
#      rejection controls (o-param's code).
# Usage: build_fc.py LINE(3945|2979) HUBROOT_3815 HUBROOT_LINE ORDER.json OUT_DIR --C 24 [--ext] --a A [--mipsec S]
import argparse, gzip, json, os, shutil, subprocess, sys, time, types
from collections import Counter
from pathlib import Path

R2 = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(os.path.dirname(os.path.dirname(R2)), 'tools')   # round-1 tools (hall_fc, fc_cost, offsets3945)
sys.path.insert(0, TOOLS); sys.path.insert(0, R2)
import hall_fc as HF  # noqa: E402

LOCAL3945 = {1: 2482, 2: 712, 3: 218, 4: 109, 5: 20, 8: 8}
LOCAL3815 = {1: 2440, 2: 725, 3: 224, 4: 102, 5: 20, 8: 8}
SUP3815 = [26, 2, 18, 53, 96, 213]


def rd(p):
    b = Path(p).read_bytes()
    return json.loads(gzip.decompress(b) if b[:2] == b'\x1f\x8b' else b)


def stub_package(root, template_path, widths):
    """o-param's stand-in for PR3734's package directory as build_fill.build() reads it (widths in NATIVE frames)."""
    inp = Path(root) / 'inputs'
    inp.mkdir(parents=True, exist_ok=True)
    frames = [[1 << i for i in range(k)] for k in range(10)]
    nat = {'frames': frames, 'start': [0] * 632, 'final': [0] * 224 + [int(w) for w in widths]}
    (inp / 'NATIVE-base.json.gz').write_bytes(gzip.compress(json.dumps(nat).encode()))
    shutil.copy2(template_path, inp / 'TEMPLATE.json.gz')
    for f in ('DISPATCH.json', 'NATIVE-extra.json.gz', 'PLAN.json.gz'):
        (inp / f).write_bytes(b'{}')
    (Path(root) / 'owner_audit.py').write_text('# stub: build_fill copies this file; not run\n')
    return Path(root)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('line', choices=['3945', '2979']); ap.add_argument('hub3815'); ap.add_argument('hubline')
    ap.add_argument('order'); ap.add_argument('out')
    ap.add_argument('--C', type=int, default=24); ap.add_argument('--ext', action='store_true')
    ap.add_argument('--a', type=float, required=True); ap.add_argument('--mipsec', type=float, default=900.0)
    ap.add_argument('--distinct', type=int, default=40); ap.add_argument('--lp-only', action='store_true')
    ap.add_argument('--suffix', type=int, default=0, help='suffix-interval rows for A in the last SUFFIX ticks')
    ap.add_argument('--suffix-step', type=int, default=1)
    ap.add_argument('--placer', choices=['hub', 'flow', 'tolerant'], default='hub')
    ap.add_argument('--suffix-sep', type=int, default=0, help='after each MIP, add the violated suffix rows over the last N ticks and re-solve')
    ap.add_argument('--fix-late', help='SOL json of a placeable design (e.g. B-like): classes whose window ends in the last --fix-ticks calls keep exactly its patterns')
    ap.add_argument('--fix-ticks', type=int, default=200)
    ap.add_argument('--from-sol', help='skip the LP/MIP: use this SOL-*.json class solution (integral) as the LP result')
    a = ap.parse_args()
    import numpy as np, highspy
    construction = Path(a.hub3815) / 'research/staged544-clock-ceiling/construction'
    order = rd(a.order); T = 1 + max(r for _, r in order)
    if a.line == '3945':
        sup, local, final_local, widen4 = SUP3815, LOCAL3815, LOCAL3945, True
    else:
        P = Path(a.hubline) / 'research/rank7-joint-endframes'
        providers = rd(P / 'allocation/PROVIDERS.json')
        widths = providers['gauge44']['widths']
        assert widths == providers['extra']['widths'] and len(widths) == 408
        inv = Counter(widths)
        sup = [inv.get(w, 0) for w in (4, 5, 6, 7, 8, 9)]
        local = {int(k): v for k, v in providers['gauge44']['actual_paid_profile'].items()}
        final_local, widen4 = local, False
    H = HF.load(str(construction / 'hall_lp.py'), T, sup, local, a.C, a.ext, widen4)
    BF = types.ModuleType('build_fill_hub'); BF.__file__ = str(construction / 'build_fill.py')
    bsrc = (construction / 'build_fill.py').read_text()
    if a.placer == 'flow':
        import place_patch
        bsrc = place_patch.patch(bsrc)
    elif a.placer == 'tolerant':
        import fallback_patch
        bsrc = fallback_patch.patch(bsrc)
        os.environ.setdefault('FC_A', str(a.a))
    exec(compile(bsrc, 'build_fill_hub', 'exec'), BF.__dict__)
    assert BF.H is H
    out = Path(a.out); out.mkdir(parents=True, exist_ok=False)
    if a.line == '2979':
        BF.PKG = stub_package(out / 'stubpkg', P / 'allocation/vendor/TEMPLATE.json.gz', widths)
    log = dict(line=a.line, C=a.C, ext=a.ext, a=a.a, T=T, order=a.order)
    orig_init, orig_run = H.HallLP.__init__, H.HallLP.run

    fixsol = None
    if a.fix_late:
        fixsol = {}
        for k, j, v in json.load(open(a.fix_late))['cols']:
            if round(v): fixsol.setdefault(int(k), {})[int(j)] = int(round(v))
    orig_add = H.HallLP.add_columns

    def add_columns(self, new):
        if fixsol is None or not hasattr(self, 'fixed'):
            return orig_add(self, new)
        new = [(k, j) for k, j in new if k not in self.fixed or j in fixsol.get(k, {})]
        n0 = len(self.cols); r = orig_add(self, new)
        for idx in range(n0, len(self.cols)):
            k, j = self.cols[idx]
            if k in self.fixed:
                v = float(fixsol[k][j]); col = self.first_y + idx
                self.h.changeColsBounds(1, np.array([col], dtype=np.int32), np.array([v]), np.array([v]))
        return r
    H.HallLP.add_columns = add_columns

    def init(self, order, a_=None, seeds=None, verbose=True, **kw):
        if fixsol is not None:
            cap_, lo_, hi_ = H.classes(order)
            self.fixed = set(int(k) for k in np.where(hi_ >= H.N - a.fix_ticks)[0])
            print('fixing %d late classes to the given placeable design' % len(self.fixed), flush=True)
            # seeds: for fixed classes only the given patterns (others are never created)
        orig_init(self, order, a=a.a, seeds=seeds, verbose=False)
        if fixsol is not None:
            have = set(self.cols)
            missing = [(k, j) for k in self.fixed for j in fixsol.get(k, {}) if (k, j) not in have]
            self.add_columns(missing)
            for idx, (k, j) in enumerate(self.cols):
                if k in self.fixed:
                    v = float(fixsol.get(k, {}).get(j, 0)); col = self.first_y + idx
                    self.h.changeColsBounds(1, np.array([col], dtype=np.int32), np.array([v]), np.array([v]))

    def report(self, tag):
        s = HF.summarize(H, self, T, a.C, a.ext, widen4, final_local)
        y = self.x[self.first_y:]
        sol = [[int(k), int(j), float(v)] for (k, j), v in zip(self.cols, y) if v > 1e-9]
        gm = {}
        for k, j, v in sol:
            c = int(self.cap[k]); q = [int(x) for x in H.PATS[c][0][j]]
            key = '%d:%s:%s' % (c, ','.join(map(str, q)), 'H' if (k % 9 in ({4} | ({0, 2} if a.ext else set())) and int(H.PATS[c][1][j]) == c) else '')
            gm[key] = gm.get(key, 0.0) + v
        (out / ('SOL-%s%s.json' % (tag.split()[0], ('%d' % len(log.get('mip', [])) if tag.startswith('MIP') else '')))).write_text(json.dumps(dict(cols=sol, global_mix=gm, widths=[int(w) for w in H.W])))
        print('%s: cut price %.13e heads %s fractional classes %d slack %.3g final remainders %s' % (
            tag, s['cut_price'], {k: round(v, 3) for k, v in s['heads'].items()}, s['fractional_classes'], s['slack'],
            {k: round(v, 2) for k, v in s['final_remainders'].items()}), flush=True)
        return s

    def run(self, max_rounds=200):
        if a.from_sol:
            sol = json.load(open(a.from_sol))
            cols = [(int(k), int(j)) for k, j, v in sol['cols']]
            self.add_columns([c for c in cols if c not in self.col_set])
            self.x = np.zeros(self.ncols)
            pos = {c: n for n, c in enumerate(self.cols)}
            for k, j, v in sol['cols']:
                self.x[self.first_y + pos[(int(k), int(j))]] = round(v)
            log['from_sol'] = a.from_sol
            report(self, 'FROMSOL')
            return 0.0
        if a.suffix:
            self.add_suffix(list(range(H.N - a.suffix, H.N, a.suffix_step)))
            print('suffix rows added: %d' % len(self.suffix), flush=True)
        t0 = time.time(); r = orig_run(self, max_rounds)
        log['lp'] = report(self, 'LP'); log['lp']['seconds'] = time.time() - t0
        if a.lp_only:
            (out / 'BUILD.json').write_text(json.dumps(log, indent=1)); raise SystemExit('LP only')
        if not log['lp']['fractional_classes']:
            return r
        ycols = np.arange(self.first_y, self.first_y + len(self.cols), dtype=np.int32)
        for it in range(8):
            self.h.changeColsIntegrality(len(ycols), ycols, np.array([highspy.HighsVarType.kInteger] * len(ycols)))
            self.h.setOptionValue('time_limit', self.h.getRunTime() + a.mipsec); self.h.setOptionValue('mip_rel_gap', 1e-9)
            t1 = time.time(); self.h.run(); st = self.h.modelStatusToString(self.h.getModelStatus())
            self.x = np.array(self.h.getSolution().col_value)
            rec = report(self, 'MIP %d (%s, %.0fs)' % (it, st, time.time() - t1)); rec['status'] = st
            added = self.separate(); rec['new_hall_rows'] = added
            sadd = suffix_separate(self) if a.suffix_sep else 0
            rec['new_suffix_starts'] = sadd
            log.setdefault('mip', []).append(rec)
            print('   new Hall rows on the MIP point: %d, new suffix starts: %d' % (added, sadd), flush=True)
            if added == 0 and sadd == 0:
                break
            # the new rows go straight into the next MIP round (HiGHS's timer is cumulative, see time_limit above)
        return H.root(self.profile())
    def suffix_separate(self):
        """add suffix-interval rows for every start A in the last a.suffix_sep ticks whose condition the current point
        violates (vectorised over the columns)."""
        y = self.x[self.first_y:]
        Nn = H.N; have = {A for A, i, row in getattr(self, 'suffix', [])}
        nz = np.where(y > 1e-9)[0]
        CK = np.array([self.cols[i][0] for i in nz]); CJ = np.array([self.cols[i][1] for i in nz]); CN = y[nz]
        CS = np.array([HF_SLOTS[int(self.cap[k])][j] for k, j in zip(CK, CJ)])
        clo = self.lo[CK] + 1; chi = np.minimum(self.hi[CK], Nn - 1)
        bad = []
        for A in range(Nn - a.suffix_sep, Nn):
            if A in have: continue
            ov = np.maximum(chi - np.maximum(clo, A) + 1, 0)
            capS = (CN[:, None] * np.minimum(CS, ov[:, None])).sum(0)
            if (capS < H.SUBSET_SUP * (Nn - A) - 1e-6).any(): bad.append(A)
        if bad:
            self.add_suffix(bad)
            print('   suffix rows added for starts', bad[:30], '...' if len(bad) > 30 else '', flush=True)
        return len(bad)
    import fc_cost as _fcc
    HF_SLOTS = _fcc.SLOTS
    H.HallLP.__init__ = init; H.HallLP.run = run
    t0 = time.time()
    try:
        built, r_int = BF.build(a.order, out / 'built', distinct=a.distinct)
    except SystemExit as e:
        log['build_error'] = str(e); (out / 'BUILD.json').write_text(json.dumps(log, indent=1, default=str))
        print('BUILD STOPPED:', e, flush=True); raise
    log['build_seconds'] = time.time() - t0
    if a.line == '3945':
        res = subprocess.run([sys.executable, '-B', str(built / 'owner_audit.py')], capture_output=True, text=True, cwd=str(built))
        log['owner_audit'] = dict(rc=res.returncode, out=res.stdout[-3000:], err=res.stderr[-3000:])
        print('owner_audit (PR3734):', res.returncode, res.stdout[:800], res.stderr[-1500:], flush=True)
        res = subprocess.run([sys.executable, '-I', '-B', os.path.join(TOOLS, 'offsets3945.py'),
                              str(Path(a.hubline) / 'research/staged544-terminal-gap-fill/owner'),
                              str(built / 'inputs/FILL.json.gz'), str(out / 'final3945')], capture_output=True, text=True)
        log['offsets3945'] = dict(rc=res.returncode, out=res.stdout[-3000:], err=res.stderr[-3000:])
        print('offsets3945:', res.returncode, res.stdout[-800:], res.stderr[-1500:], flush=True)
    else:
        fill = rd(built / 'inputs/FILL.json.gz')
        dispatch_raw = rd(P / 'allocation/DISPATCH.json')
        dispatch = {(s, r): p for s, r, p in dispatch_raw}
        plan, template = rd(P / 'allocation/vendor/PLAN.json.gz'), rd(P / 'allocation/vendor/TEMPLATE.json.gz')
        prov = {k: {'widths': v['widths'], 'actual_paid_profile': {int(r): n for r, n in v['actual_paid_profile'].items()}}
                for k, v in providers.items()}
        engine = P / 'allocation/vendor/independent_owners.py'
        E = types.ModuleType('pr2839_engine'); E.__file__ = str(engine)
        exec(compile(engine.read_bytes(), 'pr2839_engine', 'exec'), E.__dict__)
        result, win, lastrow = E.check(T, order, fill, dispatch, None, plan, template, prov)
        rejected = []

        def reject(name, fn, reason):
            try:
                fn()
            except ValueError as e:
                assert str(e) == reason, (name, str(e)); rejected.append(name)
            else:
                raise AssertionError('accepted corruption ' + name)
        rerun = lambda T_=T, f=fill, d=dispatch, c=None: E.check(T_, order, f, d, c, plan, template, prov)
        reject('false_replica_count', lambda: rerun(T_=T + 1), 'complete T call coverage')
        reject('missing_literal_helper_role', lambda: rerun(f={**fill, 'rows': fill['rows'][:-1]}), 'complete role count')
        bad = dict(dispatch); bad[tuple(order[0])] = 'unknown-provider'
        reject('unbound_provider_at_actual_call', lambda: rerun(d=bad), 'provider dispatch coverage')

        def mutate(i, new):
            rows = list(fill['rows']); rows[i] = new; return {**fill, 'rows': rows}
        t, h, o, c, p, z = fill['rows'][0]
        owner, child = next(k for k, (cap, lo, hi) in win.items() if not lo < t <= hi)
        reject('unavailable_source_window', lambda: rerun(f=mutate(0, [t, h, owner, child, p, z])), 'literal source window')
        j = next(i for i, r in enumerate(fill['rows'][1:], 1) if r[0] == t); r = fill['rows'][j]
        reject('same_call_physical_family_collision', lambda: rerun(f=mutate(j, [r[0], r[1], o, c, p, r[5]])), 'same-call family collision')
        key, (zz, i, w) = next((k, v) for k, v in lastrow.items() if v[0] + v[2] < win[k[:2]][0]); r = list(fill['rows'][i]); r[5] += 1
        reject('unpaid_hole_in_compact_pending_prefix', lambda: rerun(f=mutate(i, r)), 'compact bin prefix')
        wrong = dict(result['profile']); wrong[1] = wrong.get(1, 0) + 1
        reject('changed_total_paid_profile', lambda: rerun(c=wrong), 'complete claimed profile')
        final = out / 'design'; final.mkdir()
        (final / 'ORDER.json').write_text(json.dumps(order))
        (final / 'DISPATCH.json').write_text(json.dumps(dispatch_raw))
        (final / 'FILL.json.gz').write_bytes(gzip.compress(json.dumps(fill).encode()))
        (final / 'PROFILE.json').write_text(json.dumps({str(k): v for k, v in sorted(result['profile'].items(), key=lambda x: int(x[0]))}, indent=2))
        log.update(engine='PASS', engine_controls_rejected=rejected, idle=result['idle'], profile=result['profile'])
        print('ENGINE PASS (PR2839 vendored by #2979), controls rejected:', rejected, flush=True)
        print('idle', result['idle'], flush=True)
    (out / 'BUILD.json').write_text(json.dumps(log, indent=1, default=str))
    print('DONE', flush=True)


if __name__ == '__main__':
    main()
