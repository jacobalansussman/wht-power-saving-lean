"""e-r45m1 (round wht45/lever, 19 Oct 2026): R5, the MIXED rate fact of a cut word, as a Lean data module.

Rate fact (Theorem 1G, rl-prove; Work.Mixed.Phi `RateQ z 1` in normalised form, e-r12's `twoTarget_rateQ_hist`):

    sum_r n_r (r/m)^z1  <=  (W - H) + H (kc/m)^z1          (children left, arrays right)

at the rational exponent z1 = 1 - eps1, eps1 = (N*K + d1) / (10^10 * K), where N is the DFT grid numerator
(final saving N/10^10); the seam exponent is alphaZ = 1 - (N*K + d2)/(10^10*K) with 0 < d2 < d1 < K, so that
z1 < alphaZ < 1 - N/10^10 = target.

Enclosures (all checked here EXACTLY with Fractions, then re-checked by Lean's norm_num):
  * each child term: Work.BlockAccounting.Moment4.term_le4 (upper bound, quartic) with a rational L >= log(m/r)
    certified by `log_le_of_le_expsum` (x-pack tools4 gen_rates.log_upper);
  * the cut-target term: Work.BlockAccounting.Sharp.term_ge (lower bound, cubic) with a rational L <= log(m/kc)
    certified by `le_log_of_expbound` (x-pack tools4 gen_sharp.log_lower).
The full-target term is (m/m)^z1 = 1 exactly.

The generator REFUSES (exit 2, nothing written) when the enclosed inequality fails.  `--scan` finds the largest N
that the enclosure proves and checks that N + 1 is refused.

usage:
  python3 -I -B mk_rate_mixed.py emit  NAME N OUTDIR      (writes OUTDIR/Work/Mixed/Rate<NAME>.lean + OUTDIR/../out/rate-<NAME>.json)
  python3 -I -B mk_rate_mixed.py check NAME N             (exit 0 iff the enclosure proves N)
  python3 -I -B mk_rate_mixed.py scan  NAME               (largest N, and the refusal of N + 1)
NAME in PROFILES below (half-cut profiles of #3945 and #3815, lever-reprice/tools/lever_price.py), or NAME=profile.json
for any other separated two-target profile (e.g. the port's exact histogram): {"children": {...}, "W", "H", "kc", "m"}.
"""
import sys, os, json
from fractions import Fraction
from decimal import Decimal, getcontext
sys.path.insert(0, os.environ['GEN_ROOT'] + '/wht43/bound-check/x-pack/tools4')
import gen_rates, gen_sharp
from gen_rates import frac_lean
getcontext().prec = 80

K, D1, D2 = 100, 2, 1        # z1 at N + 0.02 grid steps, alphaZ at N + 0.01, target at N
DD = 12                      # decimals of the log bounds

# half cut (rl-plan / lever-reprice): heads (rank 48) cut by 24 axes -> rank-24 remnants; targets {73: W-H, 49: H}
def half(prof, W, heads, non):
    ch = dict(prof)
    ch[48] -= heads + non
    assert ch[48] == 0
    del ch[48]
    ch[24] = ch.get(24, 0) + heads
    ch[48] = non
    return dict(children=dict(sorted(ch.items())), W=W, H=heads, kc=49, m=73)

PROFILES = {
    # #3945 (Dugongue, staged544 T360 terminal-gap fill): lever-reprice P3945, heads 57,312 (heads_fill.py), others 80,640
    '3945': half({1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920, 48: 137952}, 322560, 57312, 80640),
    # #3815 (ChaloupeMusic, staged544 T360 clock): same idle-48 bins as #3945 (lever-reprice bounds_rows, lower end)
    '3815': half({1: 7989840, 2: 2349000, 3: 749088, 4: 330480, 5: 69120, 8: 25920, 48: 137952}, 322560, 57312, 80640),
}


def ceil_frac(x, G):
    y = x * 10 ** G
    return Fraction(-((-y.numerator) // y.denominator), 10 ** G)


def floor_frac(x, G):
    y = x * 10 ** G
    return Fraction(y.numerator // y.denominator, 10 ** G)


def root_dec(p):
    """exact-ish saving root (Decimal bisection), same formula as lever-reprice mixed.py"""
    m = p['m']
    ch = [(Decimal(n), (Decimal(r) / m).ln()) for r, n in p['children'].items()]
    tg = [(Decimal(p['W'] - p['H']), Decimal(0)), (Decimal(p['H']), (Decimal(p['kc']) / m).ln())]
    lo, hi = Decimal('0.97'), Decimal(1)
    for _ in range(200):
        z = (lo + hi) / 2
        lhs = sum(n * (l * z).exp() for n, l in ch)
        rhs = sum(n * (l * z).exp() for n, l in tg)
        if lhs > rhs:
            lo = z
        else:
            hi = z
    return 1 - (lo + hi) / 2


def plan(p, N):
    """the rational certificate at eps1 = (N K + D1)/(10^10 K); None if the enclosure fails"""
    eps = Fraction(N * K + D1, 10 ** 10 * K)
    m = p['m']
    rows = []
    for r, n in p['children'].items():
        assert 1 <= r < p['kc'] <= m
        L, Nt = gen_rates.log_upper(m, r, DD)
        x = eps * L
        assert x <= 1 and L >= 0
        val = n * Fraction(r, m) * (1 + x + x * x / 2 + x ** 3 / 6 + x ** 4 * Fraction(5, 96))
        rows.append(dict(r=r, n=n, L=L, N=Nt, val=val))
    Ll, s, nn = gen_sharp.log_lower(m, p['kc'], DD)
    x = eps * Ll
    low = p['H'] * Fraction(p['kc'], m) * (1 + x + x * x / 2 + x ** 3 / 6)
    margin = (p['W'] - p['H']) + low - sum(rw['val'] for rw in rows)
    if margin <= 0:
        return None
    G = 0
    while Fraction(len(rows) + 1, 10 ** G) > margin / 50:
        G += 1
    for rw in rows:
        rw['q'] = ceil_frac(rw['val'], G)
        assert rw['q'] >= rw['val']
    qlo = floor_frac(low, G)
    assert qlo <= low
    tot = sum(rw['q'] for rw in rows)
    assert tot <= (p['W'] - p['H']) + qlo          # the inequality Lean's linarith closes
    return dict(eps=eps, rows=rows, lower=dict(L=Ll, s=s, n=nn, val=low, q=qlo), G=G,
                margin=margin, slack_rel=float(margin / p['W']))


def plan_fail(p, N1):
    """certificate that the rate fact FAILS at the grid exponent 1 - N1/10^10 (lower bounds of the children with
    term_ge, upper bound of the cut target with term_le4); None if the enclosure cannot show it"""
    eps = Fraction(N1, 10 ** 10)
    m = p['m']
    rows = []
    for r, n in p['children'].items():
        L, s, nn = gen_sharp.log_lower(m, r, DD)
        x = eps * L
        rows.append(dict(r=r, n=n, L=L, s=s, nn=nn, val=n * Fraction(r, m) * (1 + x + x * x / 2 + x ** 3 / 6)))
    Lu, Nt = gen_rates.log_upper(m, p['kc'], DD)
    x = eps * Lu
    assert x <= 1
    up = p['H'] * Fraction(p['kc'], m) * (1 + x + x * x / 2 + x ** 3 / 6 + x ** 4 * Fraction(5, 96))
    margin = sum(rw['val'] for rw in rows) - (p['W'] - p['H']) - up
    if margin <= 0:
        return None
    G = 0
    while Fraction(len(rows) + 1, 10 ** G) > margin / 50:
        G += 1
    for rw in rows:
        rw['q'] = floor_frac(rw['val'], G)
    qu = ceil_frac(up, G)
    assert sum(rw['q'] for rw in rows) > (p['W'] - p['H']) + qu
    return dict(eps=eps, rows=rows, upper=dict(L=Lu, N=Nt, q=qu), G=G, margin=margin)


SHARP_HEAD = """import Work.Mixed.Rate{name}
import Work.Mixed.RateMono

/-!
# R5 negative control for #{name}: the mixed rate fact is FALSE one grid step above the figure
(generated by e-r45m1's tools/mk_rate_mixed.py; do not edit)

* `rate{name}_sharp`: at `1 - {N1}/10^10` the children's moment exceeds `{Wf} + T {m} z {kc} {H}`
  (lower bounds `term_ge` for the children, upper bound `term_le4` for the cut target; margin {mg:.3e});
* `rate{name}_fails`: the same at EVERY exponent `z ≤ 1 - {N1}/10^10` (`twoTarget_fail_mono`), so the rate fact of
  `Rate{name}` cannot be moved to the saving {N1}/10^10 or beyond.
-/

set_option maxRecDepth 4000

namespace OAI.PowerSaving.MX
open OAI.PowerSaving.BlockAccounting Finset

theorem rate{name}_sharp :
    (({Wf}:ℕ):ℝ) + T {m} (1 - {epsL}) {kc} {H} < moment {m} (1 - {epsL}) cut{name} := by
  have hsum : moment {m} (1 - {epsL}) cut{name}
      = {chain} := rfl
"""

SHARP_ROW = """  have b{i} : {q} ≤ T {m} (1 - {epsL}) {r} {n} :=
    term_ge {m} {r} {n} {epsL} {L} _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (le_log_of_expbound _ _ {s} {nn} (by norm_num) (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num)
"""

SHARP_TAIL = """  have c0 : T {m} (1 - {epsL}) {kc} {H} ≤ {q} :=
    term_le4 {m} {kc} {H} {epsL} {L} _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ {Nt} (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  rw [hsum]
  have e0 : (({Wf}:ℕ):ℝ) = {Wf} := by norm_num
  rw [e0]
  linarith

theorem rate{name}_fails {{z : ℝ}} (hz : z ≤ 1 - {epsL}) :
    (({Wf}:ℕ):ℝ) + T {m} z {kc} {H} < moment {m} z cut{name} :=
  twoTarget_fail_mono {m} {kc} {Wf} {H} cut{name} (by norm_num) (by norm_num) (by norm_num)
    (fun b hb => by have := cut{name}_proper b hb; constructor <;> omega) hz rate{name}_sharp

end OAI.PowerSaving.MX

#print axioms OAI.PowerSaving.MX.rate{name}_sharp
#print axioms OAI.PowerSaving.MX.rate{name}_fails
"""


def lean_sharp(name, p, N, pf):
    m, kc, W, H = p['m'], p['kc'], p['W'], p['H']
    N1 = N + 1
    epsL = '(%d/(10:ℝ)^10)' % N1
    out = [SHARP_HEAD.format(name=name, N1=N1, Wf=W - H, m=m, kc=kc, H=H, mg=float(pf['margin']), epsL=epsL,
                             chain=nested(m, epsL, p['children']))]
    for i, rw in enumerate(pf['rows']):
        out.append(SHARP_ROW.format(i=i, q=frac_lean(rw['q']), m=m, epsL=epsL, r=rw['r'], n=rw['n'],
                                    L=frac_lean(rw['L']), s=rw['s'], nn=rw['nn']))
    u = pf['upper']
    out.append(SHARP_TAIL.format(m=m, epsL=epsL, kc=kc, H=H, q=frac_lean(u['q']), L=frac_lean(u['L']), Nt=u['N'],
                                 Wf=W - H, name=name))
    return ''.join(out)


def scan(p):
    s = root_dec(p)
    N0 = int(s * 10 ** 10)
    N = N0 + 2
    while plan(p, N) is None:
        N -= 1
        assert N > N0 - 20
    assert plan(p, N + 1) is None
    return N, s


def lean_module(name, p, N, pl):
    m, kc, W, H = p['m'], p['kc'], p['W'], p['H']
    epsL = '(%d/(10:ℝ)^%d)' % (N * K + D1, 12)
    assert K == 100
    lst = ', '.join('(%d, %d)' % (r, n) for r, n in p['children'].items())
    rmax = max(p['children'])
    out = []
    out.append('''import Work.BlockAccounting.Moment4
import Work.BlockAccounting.Sharp

/-!
# R5: the mixed rate fact of the half cut of #{name} (generated by e-r45m1's tools/mk_rate_mixed.py; do not edit)

Round wht45/lever (19 Oct 2026).  Cut word: the hub's word of PR #{name} with each of its {H} literal rank-48 heads cut
by 24 of its 48 axes (rank-24 remnant kept in the body, target rank {kc} on the head array); the other arrays keep the
full target {m}.  Children (paid blocks of one unit, by rank): `cut{name}`.  Arrays: {Wf} full, {H} cut.

* `zRate{name} = 1 - {e1}/10^12`: the exponent of the rate fact (DFT grid numerator {N}, plus {d1}/{K} of a grid step);
* `rate{name}_moment`: `moment {m} zRate{name} cut{name} ≤ {Wf} + T {m} zRate{name} {kc} {H}`;
* `rate{name}`: the same in the exact shape of e-r12's `Rule.twoTarget_rateQ_hist` hypothesis `h`
  (`m = {m}`, `kc = {kc}`, `W = {W}`, `H = {H}`, `rmax = {rmax}`), for any paid-round enumeration `rk` whose
  histogram is `histFn cut{name}`.
Enclosures: `term_le4` (Moment4, upper, quartic) for every child; `term_ge` (Sharp, lower, cubic) for the cut target.
Relative slack of the rational certificate: {slack:.3e}.  The generator refuses N = {N1} (one grid step higher).
No `sorry`, no `native_decide`.
-/

set_option maxRecDepth 4000

namespace OAI.PowerSaving.MX
open OAI.PowerSaving.BlockAccounting Finset

/-- paid blocks of one unit of the half-cut word of #{name}, `(rank, count)` -/
def cut{name} : List (ℕ × ℕ) := [{lst}]

/-- exponent of the rate fact: `1 - {e1}/10^12` -/
noncomputable def zRate{name} : ℝ := 1 - {epsL}

theorem cut{name}_count : (cut{name}.map Prod.snd).sum = {P} := by norm_num [cut{name}]

theorem cut{name}_proper : ∀ b ∈ cut{name}, 1 ≤ b.1 ∧ b.1 ≤ {rmax} := by decide

theorem rate{name}_moment :
    moment {m} zRate{name} cut{name} ≤ (({Wf}:ℕ):ℝ) + T {m} zRate{name} {kc} {H} := by
  unfold zRate{name}
  have hsum : moment {m} (1 - {epsL}) cut{name}
      = {chain} := rfl
'''.format(name=name, H=H, kc=kc, m=m, Wf=W - H, W=W, e1=N * K + D1, N=N, d1=D1, K=K, rmax=rmax,
           slack=pl['slack_rel'], N1=N + 1, epsL=epsL, lst=lst, P=sum(p['children'].values()),
           chain=nested(m, epsL, p['children'])))
    for i, rw in enumerate(pl['rows']):
        out.append('''  have b{i} : T {m} (1 - {epsL}) {r} {n} ≤ {q} :=
    term_le4 {m} {r} {n} {epsL} {L} _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ {Nt} (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
'''.format(i=i, m=m, epsL=epsL, r=rw['r'], n=rw['n'], q=frac_lean(rw['q']), L=frac_lean(rw['L']), Nt=rw['N']))
    lo = pl['lower']
    out.append('''  have c0 : {q} ≤ T {m} (1 - {epsL}) {kc} {H} :=
    term_ge {m} {kc} {H} {epsL} {L} _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (le_log_of_expbound _ _ {s} {n} (by norm_num) (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num)
  rw [hsum]
  have e0 : (({Wf}:ℕ):ℝ) = {Wf} := by norm_num
  rw [e0]
  linarith

/-- **R5 for #{name}, half cut**, in the shape of the hypothesis `h` of `Rule.twoTarget_rateQ_hist`
(e-r12, Work/Mixed/TwoTarget.lean) with `m = {m}`, `kc = {kc}`, `W = {W}`, `H = {H}`, `rmax = {rmax}`, `z = zRate{name}`:
for every enumeration `rk` of the paid rounds whose histogram is that of `cut{name}`. -/
theorem rate{name} {{P : ℕ}} (rk : Fin P → ℕ)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut{name} r) :
    ∑ r ∈ Finset.range ({rmax} + 1), ((Finset.univ.filter fun j => rk j = r).card : ℝ) * ((r : ℝ) / (({m}:ℕ):ℝ)) ^ zRate{name}
      ≤ ((({W}:ℕ) - ({H}:ℕ) : ℕ) : ℝ) + (({H}:ℕ) : ℝ) * ((({kc}:ℕ) : ℝ) / (({m}:ℕ):ℝ)) ^ zRate{name} := by
  simp_rw [hcnt]
  rw [sum_histFn cut{name} (fun r => ((r:ℝ) / (({m}:ℕ):ℝ)) ^ zRate{name}) (Finset.range ({rmax} + 1))
    (fun b hb => by
      have := cut{name}_proper b hb
      exact Finset.mem_range.mpr (by omega))]
  have h := rate{name}_moment
  have e1 : ((({W}:ℕ) - ({H}:ℕ) : ℕ) : ℝ) = (({Wf}:ℕ):ℝ) := by norm_num
  rw [e1]
  exact h

end OAI.PowerSaving.MX

#print axioms OAI.PowerSaving.MX.rate{name}_moment
#print axioms OAI.PowerSaving.MX.rate{name}
'''.format(q=frac_lean(lo['q']), m=m, epsL=epsL, kc=kc, H=H, L=frac_lean(lo['L']), s=lo['s'], n=lo['n'],
           Wf=W - H, name=name, W=W, rmax=rmax))
    return ''.join(out)


def nested(m, epsL, ch):
    items = list(ch.items())
    s = '0'
    for r, n in reversed(items):
        s = 'T %d (1 - %s) %d %d + %s' % (m, epsL, r, n, s if s == '0' else '(' + s + ')')
    return s


def load_profile(name):
    """NAME from PROFILES, or NAME=path/to/profile.json with {"children": {rank: count}, "W", "H", "kc", "m"}
    (e.g. the port's exact paid histogram; NAME must be a Lean identifier suffix)"""
    if '=' in name:
        name, path = name.split('=', 1)
        j = json.load(open(path))
        p = dict(children={int(k): int(v) for k, v in sorted(j['children'].items(), key=lambda kv: int(kv[0]))},
                 W=int(j['W']), H=int(j['H']), kc=int(j['kc']), m=int(j['m']))
        assert all(1 <= r < p['kc'] <= p['m'] for r in p['children']) and 1 <= p['H'] < p['W']
        assert name.isalnum()
        return name, p
    return name, PROFILES[name]


def main():
    cmd = sys.argv[1]
    name, p = load_profile(sys.argv[2])
    if cmd == 'scan':
        N, s = scan(p)
        print('profile %s: root saving %s' % (name, s))
        print('largest DFT numerator the enclosure proves: %d (saving %.7e); N+1 = %d refused' % (N, N / 1e10, N + 1))
        pl = plan(p, N)
        print('  slack_rel %.3e, root - N = %.4f grid steps, rate point N + %d/%d' % (pl['slack_rel'], float(s * 10 ** 10 - N), D1, K))
        return
    N = int(sys.argv[3])
    pl = plan(p, N)
    if pl is None:
        print('REFUSED: the enclosure does not prove the rate fact at N = %d (z1 = 1 - %d/10^12)' % (N, N * K + D1))
        sys.exit(2)
    if cmd == 'check':
        print('OK N = %d slack_rel %.3e' % (N, pl['slack_rel']))
        return
    assert cmd == 'emit'
    outdir = sys.argv[4]
    path = os.path.join(outdir, 'Work/Mixed/Rate%s.lean' % name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, 'w').write(lean_module(name, p, N, pl))
    pf = plan_fail(p, N + 1)
    assert pf is not None, 'the enclosure cannot show the failure one step higher'
    open(os.path.join(outdir, 'Work/Mixed/Rate%sSharp.lean' % name), 'w').write(lean_sharp(name, p, N, pf))
    js = dict(name=name, dft_numerator=N, rate_eps=[N * K + D1, 10 ** 12], alpha_eps=[N * K + D2, 10 ** 12],
              children={str(k): v for k, v in p['children'].items()}, W=p['W'], H=p['H'], kc=p['kc'], m=p['m'],
              slack_rel=pl['slack_rel'], G=pl['G'], root=str(root_dec(p)))
    jo = os.path.join(outdir, '..', 'out')
    os.makedirs(jo, exist_ok=True)
    json.dump(js, open(os.path.join(jo, 'rate-%s.json' % name), 'w'), indent=1)
    print('written', path, 'slack_rel %.3e' % pl['slack_rel'])


if __name__ == '__main__':
    main()
