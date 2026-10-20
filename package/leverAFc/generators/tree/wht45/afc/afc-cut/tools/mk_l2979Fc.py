"""afc-cut (round wht45/afc, 20 Oct 2026): aA-data's tools/mk_l2979H.py (sha256 f4d3f2021bd7ba2a) for A-fc, the FULL cut on
#2979's words with the re-solved fill u1-2979sep (suffix L2979Fc, N 17887045, profile 2979Fc = tools/profile2979Fc.json,
rate module Rate2979Fc.lean of afc-cut, cut target kc = 25 instead of 49).  Only defaults, the profile, paths, the kc
in the profile check and texts change; the chain text and every check are x-inst's.
aA-data's description: x-inst's tools/mk_l2979.py for design A (#2979 with every head kept,
FILL 764395fc; suffix L2979H, N 17564707, profile 2979H = tools/profile2979H.json, rate module Rate2979H.lean of
aA-data).  Only defaults, the profile, paths and texts change; the chain text and every check are x-inst's.
x-inst's description: generator of the integration skeleton of the #2979 lever instance
L2979.  It is i-skel's tools/mk_l3945.py (wht45/lever, sha256 c5ba66d8...) with the #2979 defaults and profile;
the chain text and the checks are unchanged.  Original description, with the #2979 values:

It writes, for the suffix SUF (default L2979Fc) and the DFT grid numerator N (default 17887045, afc-cut's R5 figure;
final saving N/10^10):

  stage/Work/<SUF>/Stmt.lean      `MX.<SUF>.EngineStmt`: OpenAI's `hills_program` statement (TensorProgram.lean:95-99)
                                  with the envelope at the seam exponent alpha = 1 - (100N+1)/10^12, as a proposition
                                  (index types in universe 0), and `MX.<SUF>.EngineStmt.lift` (every universe, by
                                  OpenAI's `Able.reindex`);
  stage/Work/Fourier<SUF>/*, stage/Work/Fourier<SUF>Neg/*, stage/comparator/fourier-<suf>{,-neg}.json
                                  the Fourier chain, the statement files (Goal: decimal exponent 1 - N/10^10) and the
                                  challenge files at N (positive) and N+1 (negative), by e-r45m1's mk_fourier_mixed.py
                                  (copied byte for byte into tools/, sha256 4df517b4...), with the seam's engine call
                                  set by MODE;
  stage/Work/Mixed/SeamConst.lean e-r45m1's R4 lemmas (bundle copy, sha256 9c623573...), imported by every mixed seam;
  meta/                           provenance, engine call, gate report, the extra axioms file, sha256 list.

MODE
  hyp    the seam's engine fact enters as an explicit hypothesis: every chain file downstream of the seam gets
         `variable [hE_<SUF> : Fact OAI.PowerSaving.MX.<SUF>.EngineStmt.{0}]` and `include hE_<SUF>` right after its
         `namespace OAI`, and the seam's engine theorem is `MX.<SUF>.EngineStmt.lift hE_<SUF>.out cl m k v`.
         No sorry, no axiom: the chain compiles now; its final theorems take the instance as an argument.
  final  the seam's engine theorem is `MX.<SUF>.EngineStmt.lift MX.<SUF>.engine cl m k v`, imported from the module
         --engine-module (default Work.<SUF>.Engine): THE ONE THEOREM x-final supplies,
             theorem OAI.PowerSaving.MX.<SUF>.engine : OAI.PowerSaving.MX.<SUF>.EngineStmt.{0}
         The chain text is then exactly e-r45m1's mk_fourier_mixed.py output (MXT control: Compare and comparator).
  stub   final, plus stage/Work/<SUF>/Engine.lean proving `MX.<SUF>.engine` from our BUILT U8X0HQ engine (envelope
         1 - 16256630/10^10) through `MX.can_envelope_of_real`; only for N <= 16256629 (dry run of the final path).

GATE (the generator refuses, exit 2, nothing written, when it fails):
  hyp/final: (1) the mixed rate fact of the profile holds at the rate point z1 = 1 - (100N+2)/10^12, by my own
             Decimal bisection of the root (formula of lever-reprice mixed.py), and (2) e-r45m1's enclosure
             (`mk_rate_mixed.py check PROFILE N`, the certificate the Lean rate fact is generated from) proves it.
             For A-fc (2979Fc) this accepts 17887045 and refuses 17887046 (one grid step higher).
  stub:      the U8X0HQ engine's exponent 1 - 16256630/10^10 is at most alpha: accepts 16256629, refuses 16256630.
Also checked: alpha < target = 1 - N/10^10 exactly; both challenge files byte-identical to my own derivation from
OpenAI's statement (UniformFourier-fd4aeeb2..., sha256 4c571f27...); the negative differs only in the numerator.

usage: python3 -I -B mk_l2979Fc.py MODE OUT [--suffix L2979Fc] [--num 17887045] [--profile 2979Fc|NAME=profile.json]
                                           [--engine-module Work.L2979Fc.Engine]
"""
import sys, os, re, json, hashlib, shutil, subprocess, argparse
from fractions import Fraction
from decimal import Decimal, getcontext

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.environ['GEN_ROOT'] + '/wht43/bound-check/fsp-build/stage-u8x0hq'
MKF = os.path.join(HERE, 'mk_fourier_mixed.py')
MKF_SHA = '4df517b4d0011d142d08cbea485e98aca15ec024b1758b9873375330064d1d0b'
RATE = os.path.join(HERE, 'mk_rate_mixed.py')      # fc-rate's copy of e-r45m1's generator (doc keys only)
RATE_SHA = '58caa520891c3a3628f73bcb74037e38009c1dac035082bb56d72ec7e2f301ce'
RATELEAN = os.path.join(HERE, '..', 'lean', 'Work', 'Mixed', 'Rate%s.lean')
# e-r45m1's rate script knows only 3945/3815 by name: the #2979 profile goes to it as NAME=profile.json
PROFILE_FILES = {'2979Fc': os.path.join(HERE, 'profile2979Fc.json')}
SEAMCONST = os.environ['GEN_ROOT'] + '/wht45/lever/e-r45m1/out/bundle/Work/Mixed/SeamConst.lean'
SEAMCONST_SHA = '9c623573443be64eaab1da26ca0ca75d2a5cd11b8798e5024f21675ff1f135b0'
REF = os.environ['GEN_ROOT'] + '/wht43/bound-check/k-build/ref/UniformFourier-fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb.lean'
REF_SHA = '4c571f275bd7506c3ccd1e6d78f2fca77e07112c969c1c62ad8359d61544ec74'
U8X0HQ_KERNEL = 16256630          # the built U8X0HQ engine's envelope exponent 1 - 16256630/10^10 (stub mode)
# the files of the chain downstream of the seam (all of them import Seam.lean transitively); Goal.lean, the
# challenge and Axioms.lean are not touched by the hypothesis
DOWNSTREAM = ['Seam', 'SectorAlgorithm', 'SynchronizedAlgorithm', 'WorkingTransform', 'WorkingCompiler',
              'WorkingPreparation', 'ArbitraryLength', 'TransformProgram', 'ConvolutionProgram', 'UniformBounds',
              'Asymptotics', 'Main']
# A-fc, FULL cut (afc-cut): PROFILE.json of u1-2979sep minus the 16,576 literal heads (cut by all 48 axes, target 25);
# no other rank-48 block, nothing split; every paid rank <= 23 < 25 (= afc-cut's cut2979Fc, Rate2979Fc.lean)
PROFILES = {
    '2979Fc': dict(children={1: 1630639, 2: 470196, 3: 154512, 4: 71928, 5: 11988, 6: 5, 7: 3, 8: 5328, 17: 8, 19: 1,
                             23: 25905},
                   W=66304, H=16576, kc=25, m=73),
}


def die(msg):
    print('REFUSED:', msg)
    sys.exit(2)


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def load_profile(name):
    if '=' in name:
        nm, path = name.split('=', 1)
        j = json.load(open(path))
        p = dict(children={int(k): int(v) for k, v in j['children'].items()}, W=int(j['W']), H=int(j['H']),
                 kc=int(j['kc']), m=int(j['m']))
        return nm, p
    return name, PROFILES[name]


def root_saving(p):
    """saving b* = 1 - z*, z* the root of  sum_r n_r (r/m)^z = (W - H) + H (kc/m)^z  (Decimal, 60 digits)."""
    getcontext().prec = 60
    m = Decimal(p['m'])
    ch = [(Decimal(n), (Decimal(r) / m).ln()) for r, n in p['children'].items()]
    full, cut, lc = Decimal(p['W'] - p['H']), Decimal(p['H']), (Decimal(p['kc']) / m).ln()
    f = lambda z: sum(n * (l * z).exp() for n, l in ch) - full - cut * (lc * z).exp()
    lo, hi = Decimal('0.9'), Decimal(1)          # f(lo) > 0 (fails), f(1) < 0 (holds)
    assert f(lo) > 0 and f(hi) < 0, 'rate fact must fail at 0.9 and hold at 1'
    for _ in range(190):
        z = (lo + hi) / 2
        if f(z) > 0:
            lo = z
        else:
            hi = z
    return 1 - hi


def gate(mode, N, prof):
    rep = []
    alpha = 1 - Fraction(100 * N + 1, 10 ** 12)
    target = 1 - Fraction(N, 10 ** 10)
    if not (Fraction(1, 2) < alpha < target < 1):
        die('alpha < target fails')
    rep.append('alpha = 1 - %d/10^12 < target = 1 - %d/10^10 (exact)' % (100 * N + 1, N))
    if mode == 'stub':
        if not Fraction(100 * N + 1, 10 ** 12) <= Fraction(U8X0HQ_KERNEL, 10 ** 10):
            die('stub: the U8X0HQ engine (exponent 1 - %d/10^10) does not reach alpha = 1 - %d/10^12' % (U8X0HQ_KERNEL, 100 * N + 1))
        rep.append('stub: U8X0HQ engine exponent 1 - %d/10^10 <= alpha' % U8X0HQ_KERNEL)
        return rep
    name, p = load_profile(prof)
    if name in PROFILES:
        # my copy of the profile must be the histogram of e-r45m1's Lean rate module (cut<NAME>, W - H, H, kc, 73)
        lean = open(RATELEAN % name).read()
        mo = re.search(r'def cut%s : List \(ℕ × ℕ\) := \[([^\]]*)\]' % name, lean)
        cut = dict((int(x), int(y)) for x, y in re.findall(r'\((\d+), (\d+)\)', mo.group(1)))
        if cut != p['children'] or ('((%d:ℕ):ℝ)' % (p['W'] - p['H'])) not in lean or ('%d %d' % (p['kc'], p['H'])) not in lean:
            die('profile %s differs from Rate%s.lean' % (name, name))
        rep.append('profile %s = the histogram cut%s of the Lean rate module Rate%s.lean (afc-cut, by fc-rate\'s copy of e-r45m1\'s generator; sha %s)' % (name, name, name, sha(RATELEAN % name)[:8]))
    b = root_saving(p)
    pt = Decimal(100 * N + 2) / Decimal(10 ** 12)
    rep.append('profile %s: mixed root saving %s (Decimal bisection), %.4f grid steps above N' % (name, str(b)[:16], float((b - Decimal(N) / Decimal(10 ** 10)) * 10 ** 10)))
    if not pt < b:
        die('the rate point 1 - %d/10^12 lies above the mixed root (saving %s): the rate fact is false at N = %d' % (100 * N + 2, str(b)[:14], N))
    rep.append('rate point 1 - %d/10^12 below the root: the rate fact holds (independent check)' % (100 * N + 2))
    rp = ('%s=%s' % (prof, PROFILE_FILES[prof])) if prof in PROFILE_FILES else prof
    if prof in PROFILE_FILES:
        _, pj = load_profile(rp)
        if pj != p:
            die('profile file %s differs from the built-in profile %s' % (PROFILE_FILES[prof], prof))
    if sha(RATE) != RATE_SHA:
        die('tools/mk_rate_mixed.py is not fc-rate\'s 1224229c')
    r = subprocess.run([sys.executable, '-I', '-B', RATE, 'check', rp, str(N)], capture_output=True, text=True)
    rep.append('e-r45m1 enclosure: ' + (r.stdout.strip() or r.stderr.strip()))
    if r.returncode != 0:
        die('e-r45m1 enclosure does not prove the rate fact at N = %d: %s' % (N, r.stdout.strip()))
    return rep


STMT = '''import Work.Scratch.Engine
import OAI.Computability.FourierTransform.TensorProgram

/-!
# The engine statement of the instance `{S}` (generated by afc-cut's tools/mk_l2979Fc.py, round wht45/afc)

New file, not a copy of an OpenAI file.  `MX.{S}.EngineStmt` is OpenAI's `hills_program` statement (lines 95-99 of
lean/OAI/Computability/FourierTransform/TensorProgram.lean of github.com/openai/math, commit fd4aeeb, Apache-2.0) as a
proposition, with the envelope `envelope z k = ⌈(k+1)^z⌉` at the seam exponent `z = 1 - {A}/10^12` in the place of
`hills`, for index types `ι` in one universe.  DFT figure of the chain: `1 - {N}/10^10` (e-r45m1's exponent convention:
rate point `1 - {A2}/10^12` < seam exponent < statement exponent).

The integration supplies ONE theorem, in the module `{EM}`:

    theorem OAI.PowerSaving.MX.{S}.engine : OAI.PowerSaving.MX.{S}.EngineStmt.{{0}}

and the seam of `Work.Fourier{S}` uses it through `EngineStmt.lift`, which gives every universe from universe 0 by
OpenAI's `Able.reindex` along `i ↦ ⟨k i, v i⟩ : Σ n, (Bits n → ℂ)`.
-/

namespace OAI
namespace PowerSaving
namespace MX
namespace {S}
open Binary Matrix Finset RAM RAM.Ty
universe u

/-- The seam's engine fact at `1 - {A}/10^12`, index types in universe `w`. -/
def EngineStmt.{{w}} : Prop :=
  ∀ (cl : Paint) (m : Bool) {{ι : Type w}} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ),
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope (1 - {A}/(10:ℝ)^12) (k i))

/-- The engine fact for index types in universe 0 gives it in every universe (`Able.reindex`). -/
theorem EngineStmt.lift (h : EngineStmt.{{0}}) (cl : Paint) (m : Bool) {{ι : Type u}} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope (1 - {A}/(10:ℝ)^12) (k i)) :=
  Able.reindex (h cl m (fun j : (Σ n : ℕ, (Bits n → ℂ)) => j.1) (fun j => j.2))
    (fun i => (⟨k i, v i⟩ : Σ n : ℕ, (Bits n → ℂ)))

end {S}
end MX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.MX.{S}.EngineStmt.lift
'''

STUB = '''import Work.{S}.Stmt
import Work.FourierU8X0HQ.Seam
import Work.Mixed.SeamConst

/-!
# STUB engine of the dry run `{S}` (afc-cut, round wht45/afc; NOT a result)

Generated by afc-cut's tools/mk_l2979Fc.py (mode stub).  It proves `MX.{S}.EngineStmt` (seam exponent `1 - {A}/10^12`)
from our BUILT U8X0HQ engine `SeamU8X0HQ.engine_U8X0HQx` (envelope exponent `1 - 16256630/10^10`, below it) through
e-r45m1's `MX.can_envelope_of_real` with the constant `A = 2`, exactly as e-r45m1's MXT control did.  Its only purpose is
to run the FINAL path of the integration (this module in the place of the real `Work.L2979Fc.Engine`) end to end.
-/

namespace OAI
namespace PowerSaving
namespace MX
namespace {S}
open Binary Matrix Finset RAM RAM.Ty

theorem engine : EngineStmt.{{0}} := by
  intro cl m ι k v
  exact MX.can_envelope_of_real k (fun i => 2^k i*envelope (1 - 16256630/(10:ℝ)^10) (k i)) (A := 2)
    (fun i => MX.envelope_nat_real_le (by norm_num) (by norm_num) (k i))
    (SeamU8X0HQ.engine_U8X0HQx cl m k v)

end {S}
end MX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.MX.{S}.engine
'''

DOC = {
    'hyp': ("SKELETON of the A-fc full-cut lever instance (afc-cut, round wht45/afc; NOT a result).  The kernel program "
            "`Seam{S}.engine_{S}` is `MX.{S}.EngineStmt.lift` (Work/{S}/Stmt.lean) of the instance HYPOTHESIS "
            "`hE_{S} : Fact MX.{S}.EngineStmt` (OpenAI's `hills_program` statement at `1 - {A}/10^12`), a section "
            "variable included in every theorem of this chain downstream of the seam; no `sorry`, no axiom.  The "
            "integration replaces the hypothesis by the theorem `MX.{S}.engine` (mode final of the generator)."),
    'final': ("The kernel program of the full-cut #2979 word (mdumitrean, rank7 joint end frames, PR #2979 of "
              "CrocSwap/integer-mult-bounds, on Jeremy Kothe's #2957 and Mason Hale's #2917 T74 words, with the re-solved fill u1-2979sep (A-fc): every child-8/10 bin a literal head, no other idle rank-48 block; the 16,576 "
              "literal rank-48 heads cut by all 48 of their axes, target 25, every paid rank at most 23) by the "
              "mixed-target recursion (Theorem 1G): `MX.{S}.engine` (module `{EM}`) proves `MX.{S}.EngineStmt`, "
              "OpenAI's `hills_program` statement at the seam exponent `1 - {A}/10^12` for index types in universe 0, "
              "and `MX.{S}.EngineStmt.lift` (Work/{S}/Stmt.lean) gives every universe by `Able.reindex`."),
    'stub': ("DRY RUN of the final integration path (afc-cut, round wht45/afc; NOT a result): `MX.{S}.engine` "
             "(Work/{S}/Engine.lean, a stub) proves `MX.{S}.EngineStmt` (seam exponent `1 - {A}/10^12`) from our "
             "built engine of fsp-build's 14 Oct stack (envelope exponent `1 - 16256630/10^10`) through "
             "`MX.can_envelope_of_real`; `MX.{S}.EngineStmt.lift` (Work/{S}/Stmt.lean) gives every universe by "
             "`Able.reindex`.  The DFT figure is that stack's own, `1 - {N}/10^10`."),
}

HYP_LINES = ('-- x-inst HYPOTHESIS MODE (skeleton, not a result): the seam\'s engine fact enters as this instance hypothesis,\n'
             '-- included in every theorem below (no `sorry`, no axiom); the final chain has neither line.\n'
             'variable [hE_{S} : Fact OAI.PowerSaving.MX.{S}.EngineStmt.{{0}}]\n'
             'include hE_{S}\n')


def derive(ref_text, num, suf):
    """k-build's derive_challenges.sh (perl) in Python: exponent 1 - num/10^10, names with suffix suf."""
    t = ref_text.replace('1 - 1/(10^(13:ℕ))', '1 - %d/(10^(10:ℕ))' % num)
    return re.sub(r'\b(decimalExponent|decimalTime|TimeBounds|DFTGoal|ConvGoal|transform_main|convolution_main)\b',
                  lambda mo: mo.group(1) + suf, t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['hyp', 'final', 'stub'])
    ap.add_argument('out')
    ap.add_argument('--suffix', default=None)
    ap.add_argument('--num', type=int, default=None)
    ap.add_argument('--profile', default='2979Fc')
    ap.add_argument('--engine-module', default=None)
    a = ap.parse_args()
    mode = a.mode
    S = a.suffix or ('L2979FcS' if mode == 'stub' else 'L2979Fc')
    N = a.num if a.num is not None else (16256629 if mode == 'stub' else 17887045)
    EM = a.engine_module or 'Work.%s.Engine' % S
    if not re.fullmatch(r'L[A-Za-z0-9]+', S):
        die('suffix must look like L2979Fc')
    A, A2 = 100 * N + 1, 100 * N + 2
    rep = gate(mode, N, a.profile)
    if sha(MKF) != MKF_SHA:
        die('tools/mk_fourier_mixed.py is not e-r45m1\'s 4df517b4')
    if sha(SEAMCONST) != SEAMCONST_SHA:
        die('SeamConst.lean is not e-r45m1\'s 9c623573')
    if sha(REF) != REF_SHA:
        die('OpenAI statement file hash')
    if os.path.exists(a.out):
        die('refusing to overwrite %s' % a.out)
    st, meta = os.path.join(a.out, 'stage'), os.path.join(a.out, 'meta')
    os.makedirs(meta)
    fmt = dict(S=S, N=N, A=A, A2=A2, EM=EM)
    # 1. the statement module (and the stub engine)
    os.makedirs(os.path.join(st, 'Work', S))
    open(os.path.join(st, 'Work', S, 'Stmt.lean'), 'w').write(STMT.format(**fmt))
    if mode == 'stub':
        if EM != 'Work.%s.Engine' % S:
            die('stub engine module name')
        open(os.path.join(st, 'Work', S, 'Engine.lean'), 'w').write(STUB.format(**fmt))
    # 2. the chain through e-r45m1's generator
    if mode == 'hyp':
        engine = dict(imports=['Work.%s.Stmt' % S], call='  MX.%s.EngineStmt.lift hE_%s.out cl m k v\n' % (S, S))
    else:
        engine = dict(imports=['Work.%s.Stmt' % S, EM], call='  MX.%s.EngineStmt.lift MX.%s.engine cl m k v\n' % (S, S))
    engine['doc'] = DOC[mode].format(**fmt)
    ej = os.path.join(meta, 'engine-%s.json' % S)
    json.dump(engine, open(ej, 'w'), indent=1)
    r = subprocess.run([sys.executable, '-I', '-B', MKF, SRC, st, S, str(N), ej], capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        shutil.rmtree(a.out)
        die('mk_fourier_mixed.py: ' + r.stdout.strip() + r.stderr.strip())
    os.makedirs(os.path.join(st, 'Work', 'Mixed'))
    shutil.copy2(SEAMCONST, os.path.join(st, 'Work', 'Mixed', 'SeamConst.lean'))
    shutil.move(os.path.join(st, 'fourier-%s-provenance.json' % S), os.path.join(meta, 'fourier-%s-provenance.json' % S))
    chain = os.path.join(st, 'Work', 'Fourier' + S)
    # 3. hypothesis mode: the instance variable in every downstream file
    if mode == 'hyp':
        for fn in DOWNSTREAM:
            p = os.path.join(chain, fn + '.lean')
            t = open(p).read()
            if len(re.findall(r'(?m)^namespace OAI$', t)) != 1:
                die('namespace OAI in %s' % fn)
            t = re.sub(r'(?m)^namespace OAI\n', 'namespace OAI\n' + HYP_LINES.format(**fmt), t, count=1)
            open(p, 'w').write(t)
    # 4. checks of the written text
    seam = open(os.path.join(chain, 'Seam.lean')).read()
    if seam.count(engine['call'].strip()) != 1 or ('def alpha%s : ℝ := 1 - %d/(10:ℝ)^12' % (S, A)) not in seam:
        die('seam text')
    for fn in os.listdir(chain):
        t = open(os.path.join(chain, fn)).read()
        if any('U8X0HQ' in ln and 'from Work/FourierU8X0HQ/Seam.lean (fsp-build)' not in ln for ln in t.split('\n')):
            die('U8X0HQ left in %s' % fn)
        if fn != 'UniformFourierChallenge.lean' and any(re.search(r'\b(sorry|admit)\b', ln) and '`sorry`' not in ln for ln in t.split('\n')):
            die('sorry in %s' % fn)
    ref = open(REF).read()
    for num, sub in ((N, 'Fourier' + S), (N + 1, 'Fourier' + S + 'Neg')):
        got = open(os.path.join(st, 'Work', sub, 'UniformFourierChallenge.lean')).read()
        if got != derive(ref, num, S):
            die('challenge %s differs from my derivation from OpenAI\'s statement' % sub)
        rep.append('challenge %s (1 - %d/10^10): byte-identical to my derivation from OpenAI\'s statement' % (sub, num))
    goal = open(os.path.join(chain, 'Goal.lean')).read()
    if ('def decimalExponent%s : ℝ := 1 - %d/(10^(10:ℕ))' % (S, N)) not in goal:
        die('Goal.lean exponent')
    for nm, cm in (('fourier-%s.json' % S.lower(), 'Work.Fourier%s.UniformFourierChallenge' % S),
                   ('fourier-%s-neg.json' % S.lower(), 'Work.Fourier%sNeg.UniformFourierChallenge' % S)):
        c = json.load(open(os.path.join(st, 'comparator', nm)))
        if not (c['challenge_module'] == cm and c['solution_module'] == 'Work.Fourier%s.Main' % S and
                c['theorem_names'] == ['OAI.PowerSaving.transform_main%s' % S, 'OAI.PowerSaving.convolution_main%s' % S]
                and c['permitted_axioms'] == ['propext', 'Quot.sound', 'Classical.choice'] and c['enable_nanoda']):
            die('comparator config %s' % nm)
    # 5. meta: extra axioms file, gate report, hashes
    ax = ['import Work.Fourier%s.Main' % S] + ['#print axioms %s' % x for x in (
        'OAI.PowerSaving.MX.%s.EngineStmt.lift' % S, 'OAI.PowerSaving.Seam%s.engine_%s' % (S, S),
        'OAI.PowerSaving.RAM.hills%s_program' % S, 'OAI.PowerSaving.RAM.envelope_program_ge%s' % S,
        'OAI.PowerSaving.Seam%s.alpha%s_lt_target' % (S, S), 'OAI.PowerSaving.Skies.alpha_theta%s' % S,
        'OAI.PowerSaving.transform_main%s' % S, 'OAI.PowerSaving.convolution_main%s' % S)]
    if mode != 'hyp':
        ax.insert(1, '#print axioms OAI.PowerSaving.MX.%s.engine' % S)
    if mode == 'hyp':
        ax += ['#check @OAI.PowerSaving.transform_main%s' % S, '#check @OAI.PowerSaving.convolution_main%s' % S]
    open(os.path.join(meta, 'AxExtra%s.lean' % S), 'w').write('\n'.join(ax) + '\n')
    files = sorted(os.path.relpath(os.path.join(d, f), st) for d, _, fs in os.walk(st) for f in fs)
    with open(os.path.join(meta, 'stage.sha256'), 'w') as fh:
        for f in files:
            fh.write('%s  %s\n' % (sha(os.path.join(st, f)), f))
    info = dict(mode=mode, suffix=S, dft_numerator=N, statement_exponent='1 - %d/10^10' % N,
                seam_exponent='1 - %d/10^12' % A, rate_point='1 - %d/10^12' % A2, negative='1 - %d/10^10' % (N + 1),
                profile=a.profile, engine_module=EM if mode != 'hyp' else None,
                engine_theorem=('OAI.PowerSaving.MX.%s.engine : OAI.PowerSaving.MX.%s.EngineStmt.{0}' % (S, S)),
                solution_module='Work.Fourier%s.Main' % S,
                theorems=['OAI.PowerSaving.transform_main%s' % S, 'OAI.PowerSaving.convolution_main%s' % S],
                build_targets=['Work.%s.Stmt' % S] + ([EM] if mode != 'hyp' else []) +
                              ['Work.Fourier%s.Main' % S, 'Work.Fourier%s.UniformFourierChallenge' % S,
                               'Work.Fourier%sNeg.UniformFourierChallenge' % S, 'Work.Fourier%s.Axioms' % S],
                gate=rep, files=len(files), generator_sha256=sha(os.path.abspath(__file__)))
    json.dump(info, open(os.path.join(meta, 'info.json'), 'w'), indent=1)
    open(os.path.join(meta, 'targets.txt'), 'w').write(' '.join(info['build_targets']) + '\n')
    for line in rep:
        print('  ' + line)
    print('ok %s mode %s N %d: %d files in %s' % (S, mode, N, len(files), st))


if __name__ == '__main__':
    main()
