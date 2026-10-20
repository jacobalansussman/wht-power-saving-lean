"""afc-cut (round wht45/afc, 20 Oct 2026): Work/L2979Fc/Prov.lean, the data of the #2979 warp network for A-fc.

aA-final's tools/mk_prov2979H.py made Work/L2979H/Prov.lean from x-final's Work/L2979/Prov.lean by renames.  This
generator makes Work/L2979Fc/Prov.lean from aA-final's Work/L2979H/Prov.lean (sha-pinned) by:
  * the renames of afc-data's INTERFACE (schedule `sched2979Fc`, host table `K.DH2979Fc`, facts `sched2979Fc_*`,
    `resK2979Fc_*`, namespace `MX.L2979Fc`, imports S2979FcFacts / Hist2979Fc) and the unit profile `MX.unit2979Fc`
    (afc-cut's Hist2979Fc) in the place of `PB.unit2979H`;
  * the remainder histogram of `hU2979FcW`: A-fc's {1: 2269, 6: 5, 7: 3, 17: 8, 19: 1, 23: 25905, 48: 16576}
    (= PROFILE.json 5df79411 minus 666 calls of #2979's local profile), in the EXACT text of afc-data's
    `sched2979Fc_hcost` right-hand side when its S2979FcFacts.lean is given (the generator parses it and refuses if
    the numbers differ from PROFILE.json's), else in the planned v0 form (placeholder);
  * texts.
It refuses (exit 2) if any design-A fill-dependent name or number is left over (2979H, 764395fc, 17310, 4956, 1728,
28190, 16576-free design A hcost), if a substitution does not hit exactly the expected number of times, or if the
hcost of the given facts file differs from PROFILE.json minus the calls.

usage: python3 -I -B mk_prov2979Fc.py OUT_DIR [S2979FcFacts.lean]
"""
import sys, os, re, json, hashlib

SRC = os.environ['GEN_ROOT'] + '/wht45/allheads/aA-final/lean/Work/L2979H/Prov.lean'
SRC_SHA = 'b7894a73b24e1e3820e3b896f6546f913c8372abd1210b271451450543cb7f1e'
PROFILE = os.environ['GEN_ROOT'] + '/wht45/fullcut/fc-reopt/r2/designs/u1-2979sep/design/PROFILE.json'
PROFILE_SHA = '5df7941120150492c64e62dc3362095eff3dfac9facff041ebb5c25b2e78f75e'
LP = {1: 2445, 2: 706, 3: 232, 4: 108, 5: 18, 8: 8}     # #2979's local word (x-tables lp2979), 666 calls
OLD_HCOST = '17310 * φ 1 + 4956 * φ 3 + 1728 * φ 5 + 28190 * φ 48'


def die(m):
    print('REFUSED:', m)
    sys.exit(2)


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def hcost_from_profile():
    if sha(PROFILE) != PROFILE_SHA:
        die('PROFILE.json hash')
    prof = {int(k): v for k, v in json.load(open(PROFILE)).items()}
    for r, n in LP.items():
        prof[r] -= 666 * n
    assert all(v >= 0 for v in prof.values())
    return {r: v for r, v in sorted(prof.items()) if v}


def parse_hcost(text):
    """the right-hand side of `theorem sched2979Fc_hcost (φ : ℕ → ℝ) : ... = RHS := by` and its numbers"""
    mo = re.search(r'theorem sched2979Fc_hcost \(φ : ℕ → ℝ\) :\s*∑ ρ : Fin 74, ∑ c : Fin 9, hcostT sched2979Fc ρ c φ\s*'
                   r'=\s*(.*?)\s*:= by', text, re.S)
    if not mo:
        die('sched2979Fc_hcost statement not found in the facts file')
    rhs = ' '.join(mo.group(1).split())
    terms = rhs.split(' + ')
    got = {}
    for t in terms:
        m2 = re.fullmatch(r'(?:(\d+) \* )?φ (\d+)', t)
        if not m2:
            die('unexpected hcost term %r' % t)
        got[int(m2.group(2))] = int(m2.group(1) or 1)
    return rhs, got


def main():
    out = sys.argv[1]
    facts = sys.argv[2] if len(sys.argv) > 2 else None
    if sha(SRC) != SRC_SHA:
        die('aA-final Prov.lean hash')
    s = open(SRC).read()
    want = hcost_from_profile()
    if facts:
        rhs, got = parse_hcost(open(facts).read())
        if got != want:
            die('hcost of %s = %s differs from PROFILE.json minus calls %s' % (facts, got, want))
        mode = 'afc-data facts %s (sha %s)' % (os.path.basename(facts), sha(facts)[:8])
    else:
        rhs = ' + '.join('%d * φ %d' % (n, r) for r, n in want.items())
        mode = 'PLACEHOLDER (afc-data INTERFACE v0 form)'

    def rep(old, new, cnt):
        nonlocal s
        c = s.count(old)
        if c != cnt:
            die('%r occurs %d times, expected %d' % (old, c, cnt))
        s = s.replace(old, new)

    # header doc
    i, j = s.index('/-!'), s.index('-/') + 2
    doc = '''/-!
# afc-cut: the data of the #2979 warp network for A-fc (fill u1-2979sep, 7fc766c6; full cut, suffix L2979Fc)

* `I0W`, `I1W`: x-words' two #2979 provider packages (`PWG2979.ginvE` = extra, provider 0; `PWG2979.ginvB` = gauge44,
  provider 1; windowed `KLW.GInvEW .. GXD.S2979B.okW`) as p-a3's warp providers `XHostTW.GInvW` (the same seven
  fields);
* `Ic2979W`: the per-call providers by x-tables' dispatch table `D2979.disp` (x-inst's `D2979W.IcW`), with
  `hIc2979W` (one target-vector family, x-words' `PWG2979.tv_eq`);
* `hOk2979`: every window of the network (Shea's `win401N`, conjugation windows `gwinC`) has a warp meaning
  (x-words' `GXD.S2979B.okC`, by p-warpsem's `WLift.ofConj`);
* `resIcW2979Fc`: the entrance residual of every slot of every call is the host table's (x-words' `hent` against
  afc-data's residual table `DH2979Fc.resT` (= x-tables'), two kernel checks `resE_ok`, `resB_ok`, and afc-data's
  `dp_disp`); hence `hfit2979FcW` (afc-data's `sched2979Fc_fits`) and `hpos2979FcW`;
* `hU2979FcW`: the unit price of the network is `listCost φ MX.unit2979Fc` (666 calls of the local profile
  {1: 2445, 2: 706, 3: 232, 4: 108, 5: 18, 8: 8} plus the remaining climbs {%s}).
aA-final's `Work/L2979H/Prov.lean` (x-final's `Work/L2979/Prov.lean` for design A) with the schedule, host table and
unit profile of A-fc (afc-data's `sched2979Fc`, `K.DH2979Fc`; afc-cut's `MX.unit2979Fc`: 16,576 heads and no other
rank-48 block; the providers, dispatch and residuals do not depend on the fill).  afc-cut (wht45/afc, 20 Oct 2026;
made by tools/mk_prov2979Fc.py); prepared with assistance from Claude.  No `sorry`.
-/''' % ', '.join('%d: %d' % (r, n) for r, n in want.items())
    s = s[:i] + doc + s[j:]
    rep('import Work.XHostT.Inst.S2979HFacts', 'import Work.XHostT.Inst.S2979FcFacts', 1)
    rep('import Work.Mixed.Bridge.Inst2979H', 'import Work.Mixed.Hist2979Fc', 1)
    rep('namespace L2979H', 'namespace L2979Fc', 1)
    rep('end L2979H', 'end L2979Fc', 1)
    rep('MX.L2979H.', 'MX.L2979Fc.', 6)
    rep('K.DH2979H.', 'K.DH2979Fc.', s.count('K.DH2979H.'))
    rep('Inst.resK2979H_disp', 'Inst.resK2979Fc_disp', 1)
    rep('Inst.resK2979H_pos', 'Inst.resK2979Fc_pos', 1)
    rep('Inst.sched2979H_fits', 'Inst.sched2979Fc_fits', 1)
    rep('Inst.sched2979H_hcost', 'Inst.sched2979Fc_hcost', 1)
    rep('Inst.sched2979H', 'Inst.sched2979Fc', s.count('Inst.sched2979H'))
    rep('resIcW2979H', 'resIcW2979Fc', s.count('resIcW2979H'))
    rep('hfit2979HW', 'hfit2979FcW', s.count('hfit2979HW'))
    rep('hpos2979HW', 'hpos2979FcW', s.count('hpos2979HW'))
    rep('hU2979HW', 'hU2979FcW', s.count('hU2979HW'))
    rep('PB.unit2979H', 'MX.unit2979Fc', 3)
    rep('(' + OLD_HCOST + ')', '(' + rhs + ')', 1)
    rep("(aA-data's `sched2979H_fits`)", "(afc-data's `sched2979Fc_fits`)", 1)
    rep('its profile** `unit2979H`', 'its profile** `MX.unit2979Fc`', 1)
    for bad in ('2979H', '764395fc', '17310', '4956', '1728 ', '28190', 'design A', 'aA-data'):
        for ln in s.split('\n'):
            if bad in ln and 'aA-final' not in ln:
                die('left over %r in: %s' % (bad, ln))
    path = os.path.join(out, 'Work', 'L2979Fc', 'Prov.lean')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, 'w').write(s)
    print('written %s (%s); hcost %s; sha256 %s' % (path, mode, rhs, sha(path)))


if __name__ == '__main__':
    main()
