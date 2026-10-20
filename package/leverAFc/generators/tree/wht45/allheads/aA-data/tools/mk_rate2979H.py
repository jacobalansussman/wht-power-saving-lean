"""aA-data (wht45/allheads, 19 Oct 2026): R5 for design A (#2979, every head kept) with e-r45m1's generator, unchanged
(imported by path, sha-pinned), plus ONE difference in the negative control, forced by where the root lies:

  the mixed root of the half cut is b* = 1.75647080068e-3, i.e. 0.0068 grid steps ABOVE N + 1 = 17564708.  The
  chain's convention puts the rate point at N + 0.02 grid steps, so N + 1 is refused (its rate point 17564708.02 lies
  above the root) and the largest numerator the enclosure proves is N = 17564707.  But the rate fact ITSELF still
  holds at the saving 1 - 17564708/10^10 (it is 0.0068 steps below the root), so "the rate fact is false one grid
  step higher" is a false statement and e-r45m1's emit (which asks for exactly that) stops.  This wrapper proves the
  failure at the first grid point where it is true and certifiable, N + 2 = 17564709 (0.993 steps above the root):
  `rate2979H_fails` for every z <= 1 - 17564709/10^10.  The refusal of N + 1 itself is the generator's (exit 2).

usage: python3 -I -B mk_rate2979H.py PROFILE.json N LEANROOT OUTJSON
  (writes LEANROOT/Work/Mixed/Rate2979H.lean, Rate2979HSharp.lean and OUTJSON; refuses N + 1)"""
import sys, os, json, hashlib, importlib.util
from fractions import Fraction
from decimal import Decimal
GEN = os.environ['GEN_ROOT'] + '/wht45/lever/e-r45m1/tools/mk_rate_mixed.py'
GEN_SHA = 'bb56d9998536b506e49df349edbb40801d812a34d442d489dac361fc79c55757'
prof, N, root, outj = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
def die(m):
    print('REFUSED:', m); sys.exit(2)
if hashlib.sha256(open(GEN, 'rb').read()).hexdigest() != GEN_SHA:
    die('e-r45m1 generator hash')
spec = importlib.util.spec_from_file_location('mk_rate_mixed', GEN); g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)
name, p = g.load_profile('2979H=' + prof)
pl = g.plan(p, N)
if pl is None:
    die('the enclosure does not prove the rate fact at N = %d' % N)
if g.plan(p, N + 1) is not None:
    die('N + 1 = %d is also proved: N is not the largest' % (N + 1))
s = g.root_dec(p)
steps = lambda n: float(s * 10 ** 10 - n)
# where the failure can first be shown (lower enclosure of the children, upper of the cut target)
N1 = None
for d in (1, 2, 3):
    if g.plan_fail(p, N + d) is not None:
        N1 = N + d; break
if N1 is None:
    die('no failure certificate within three grid steps')
pf = g.plan_fail(p, N1)
if N1 == N + 1:
    print('note: the failure one step higher is certifiable; this wrapper is not needed')
os.makedirs(os.path.join(root, 'Work/Mixed'), exist_ok=True)
open(os.path.join(root, 'Work/Mixed/Rate%s.lean' % name), 'w').write(g.lean_module(name, p, N, pl))
sharp = g.lean_sharp(name, p, N1 - 1, pf)        # lean_sharp states the failure at (its N) + 1 = N1
old = '# R5 negative control for #%s: the mixed rate fact is FALSE one grid step above the figure' % name
assert sharp.count(old) == 1
new = ('# R5 negative control for design A (#2979, every head kept; suffix %s): the mixed rate fact is FALSE at %d/10^10\n'
       '(%d grid steps above the figure N = %d; at N + 1 the rate fact still holds, %.4f steps below the mixed root,\n'
       'but the chain\'s rate point N + 1 + %d/%d lies above the root, so the generator refuses N + 1)'
       % (name, N1, N1 - N, N, steps(N + 1), g.D1, g.K))
sharp = sharp.replace(old, new)
open(os.path.join(root, 'Work/Mixed/Rate%sSharp.lean' % name), 'w').write(sharp)
js = dict(name=name, dft_numerator=N, rate_eps=[N * g.K + g.D1, 10 ** 12], alpha_eps=[N * g.K + g.D2, 10 ** 12],
          children={str(k): v for k, v in p['children'].items()}, W=p['W'], H=p['H'], kc=p['kc'], m=p['m'],
          slack_rel=pl['slack_rel'], G=pl['G'], root=str(s), root_minus_N_steps=steps(N),
          root_minus_N1_steps=steps(N + 1), refused=N + 1, sharp_numerator=N1, sharp_margin=float(pf['margin']),
          generator=GEN, generator_sha256=GEN_SHA, wrapper_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(),
          profile_sha256=hashlib.sha256(open(prof, 'rb').read()).hexdigest())
json.dump(js, open(outj, 'w'), indent=1)
print('OK N = %d (saving %.7e), N + 1 = %d refused (root - (N+1) = %.4f steps < %d/%d), failure proved at %d; slack_rel %.3e'
      % (N, N / 1e10, N + 1, steps(N + 1), g.D1, g.K, N1, pl['slack_rel']))
