"""e-r45m1 (round wht45/lever, 19 Oct 2026): R4, the Fourier chain for a MIXED-engine instance.

Copies the U8X0HQ chain (fsp-build stage-u8x0hq: Work/FourierU8X0HQ/* (15 files), Work/FourierU8X0HQNeg/*, comparator
configs) to the suffix NEW, mechanically, as fsp-cc's mk_fourier_hx.py did for every earlier instance, with:
  * the seam exponent OFF the 1e-10 grid:  alphaNEW = 1 - (100 N + 1)/10^12   (was 1 - (N+1)/10^10);
    the same in `alpha_thetaNEW`'s local constant `let d : ℝ := (100 N + 1)/(10:ℝ)^12` (d must equal 1 - alpha);
  * the target / decimal exponent 1 - N/10^10 (Goal, Seam's targetExponent, the challenge, docstrings);
  * the negative challenge at 1 - (N+1)/10^10 (one grid step above);
  * THE SEAM's engine theorem `SeamNEW.engine_NEW` (section 3) proved by the Lean term ENGINE.json["call"], which
    must have type  Can m (fun i => 2^k i) (simIn cl) (bundle cl) (fun i => givenBits cl (k i) (v i))
                      (fun i => (bits (k i)).tape (ripple (k i) (v i))) (fun i => 2^k i*envelope alphaNEW (k i))
    with `cl m k v` in scope (e.g. the mixed engine through MX.can_envelope_of_real).
Everything else (hillsNEW_pos / _mono / _le / _real, the Fourier chain, the goal statement) is the U8X0HQ text with
names changed: the seam with a constant needs no other change (Can carries the constant; lg-g2 DESIGN finding 1).

The generator REFUSES (exit 2) if a number it must replace is missing or occurs where it should not, and if
alpha < target fails in exact rational arithmetic.

usage: python3 -I -B mk_fourier_mixed.py SRC OUT NEW N ENGINE.json
  SRC = folder with Work/FourierU8X0HQ, Work/FourierU8X0HQNeg, comparator/ (fsp-build/stage-u8x0hq)
  ENGINE.json = {"imports": [...], "call": "<Lean term, indented 2>", "doc": "<one paragraph>"}
"""
import sys, os, re, json, hashlib
from fractions import Fraction

OLD = 'U8X0HQ'
OLD_K, OLD_N = '16256630', '16256629'


def die(msg):
    print('REFUSED:', msg)
    sys.exit(2)


def main():
    src, out, NEW, N, spec = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), json.load(open(sys.argv[5]))
    if not re.fullmatch(r'[A-Z][A-Za-z0-9]*', NEW) or OLD in NEW:
        die('bad suffix')
    A = 100 * N + 1
    alpha = 1 - Fraction(A, 10 ** 12)
    target = 1 - Fraction(N, 10 ** 10)
    if not (alpha < target and Fraction(1, 2) < alpha):
        die('alpha < target fails')
    prov, written = {}, []
    for sub in ('Work/Fourier' + OLD, 'Work/Fourier' + OLD + 'Neg'):
        d = os.path.join(src, sub)
        neg = sub.endswith('Neg')
        dsub = sub.replace('Fourier' + OLD, 'Fourier' + NEW)
        os.makedirs(os.path.join(out, dsub), exist_ok=True)
        for fn in sorted(os.listdir(d)):
            if not fn.endswith('.lean'):
                continue
            t = open(os.path.join(d, fn)).read()
            prov[sub + '/' + fn] = hashlib.sha256(t.encode()).hexdigest()
            if fn == 'Seam.lean':
                t = seam(t, NEW, N, A, spec)
            else:
                t = t.replace(OLD, NEW)
            if neg:
                if t.count(OLD_K) != 1 or OLD_N in t:
                    die('negative challenge numbers in %s' % fn)
                t = t.replace(OLD_K, str(N + 1))
            else:
                for a, b in (('%s/(10:ℝ)^10' % OLD_K, '%d/(10:ℝ)^12' % A), ('%s/10^10' % OLD_K, '%d/10^12' % A)):
                    t = t.replace(a, b)
                if OLD_K in t:
                    die('unreplaced kernel numerator in %s' % fn)
                t = t.replace(OLD_N, str(N))
            if OLD in t:
                die('old suffix left in %s' % fn)
            if fn == 'Seam.lean' and not neg:
                if t.count('__IMPORTS__') != 1 or t.count('__CALL__') != 1 or t.count('__DOC__') != 1:
                    die('seam placeholders')
                imps = list(spec['imports']) + [i for i in ('Work.Block.WHT', 'Work.Mixed.SeamConst') if i not in spec['imports']]
                t = t.replace('__IMPORTS__\n', ''.join('import %s\n' % i for i in imps))
                t = t.replace('__DOC__', spec['doc'].strip()).replace('__OLD__', OLD)
                t = t.replace('__CALL__', spec['call'].rstrip('\n'))
            p = os.path.join(out, dsub, fn)
            if os.path.exists(p):
                die('refusing to overwrite %s' % p)
            open(p, 'w').write(t)
            written.append(os.path.join(dsub, fn))
    os.makedirs(os.path.join(out, 'comparator'), exist_ok=True)
    for fn in ('fourier-%s.json' % OLD.lower(), 'fourier-%s-neg.json' % OLD.lower()):
        t = open(os.path.join(src, 'comparator', fn)).read().replace(OLD, NEW)
        q = os.path.join(out, 'comparator', fn.replace(OLD.lower(), NEW.lower()))
        open(q, 'w').write(t)
        written.append(os.path.relpath(q, out))
    meta = dict(new=NEW, dft_numerator=N, alpha=[A, 10 ** 12], target=[N, 10 ** 10], negative=[N + 1, 10 ** 10],
                engine=spec, source_sha256=prov, written=written)
    json.dump(meta, open(os.path.join(out, 'fourier-%s-provenance.json' % NEW), 'w'), indent=1)
    print('ok', NEW, 'files', len(written), 'alpha 1 - %d/10^12, target 1 - %d/10^10' % (A, N))


def seam(t, NEW, N, A, spec):
    # imports
    m = re.match(r'import Work\.XHost\.%s\n' % OLD, t)
    if not m:
        die('seam import')
    t = '__IMPORTS__\n' + t[m.end():]
    # header docstring: replace the engine description paragraph (from "Twin of" to the line before "OpenAI's Fourier proof")
    i0 = t.index('/-!\n# THE SEAM')
    i1 = t.index("OpenAI's Fourier proof uses")
    head = ('/-!\n# THE SEAM for the Fourier transfer: our kernel program in OpenAI\'s `hills_program` form\n\n'
            'New file, not a copy of an OpenAI file: generated by e-r45m1\'s tools/mk_fourier_mixed.py (round wht45/lever)\n'
            'from Work/Fourier__OLD__/Seam.lean (fsp-build), suffix `%s`.  It contains lines 95-99 of\n'
            'lean/OAI/Computability/FourierTransform/TensorProgram.lean (the statement of `hills_program`) and line 114 of\n'
            'lean/OAI/Computability/FourierTransform/TensorSaving.lean (the definition of `hills`) of github.com/openai/math\n'
            '(commit fd4aeeb, Apache-2.0), with `hills`, `alpha` replaced by `hills%s`, `alpha%s`.\n\n'
            '__DOC__\n\n') % (NEW, NEW, NEW)
    t = t[:i0] + head + t[i1:]
    t = t.replace(OLD, NEW)
    # the bullet list of the header
    j0 = t.index('* `GX.engine_of_ginv%s`' % NEW)
    j1 = t.index('* `RAM.hills%s_program`' % NEW)
    t = t[:j0] + '* `Seam%s.engine_%s`   the kernel program at `alpha%s` (section 3; see the paragraph above).\n' % (NEW, NEW, NEW) + t[j1:]
    t = t.replace('* `RAM.hills%s_program`   THE SEAM, closed with the kernel-checked scalar half `GXD.U8X0.scalOK`.' % NEW,
                  '* `RAM.hills%s_program`   THE SEAM, closed with `Seam%s.engine_%s`.' % (NEW, NEW, NEW))
    # section 3
    k0 = t.index('/-! ## 3.')
    k1 = t.index('end Seam%s\n' % NEW) + len('end Seam%s\n' % NEW)
    sec3 = ('/-! ## 3. The kernel program at `alpha%s` -/\n\n'
            'namespace Seam%s\n'
            'open Binary Matrix Finset RAM RAM.Ty\n'
            'universe U\n\n'
            '/-- **Kernel program at `alpha%s = 1 - %d/10^12`**: every colour, both modes, every batch of inputs. -/\n'
            'theorem engine_%s (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :\n'
            '    Can m (fun i => 2^k i) (simIn cl) (bundle cl)\n'
            '      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))\n'
            '      (fun i => 2^k i*envelope alpha%s (k i)) :=\n'
            '__CALL__\n\n'
            'end Seam%s\n') % (NEW, NEW, NEW, A, NEW, NEW, NEW)
    t = t[:k0] + sec3 + t[k1:]
    # section 4: the seam's docstring and proof
    t = re.sub(r'the kernel program of\nthe certificate `GXD\.%s` \(the one behind `WHT\.wht_main_block_%sx`\), for every colour,'
               % (NEW, NEW), 'the kernel program `Seam%s.engine_%s`, for every colour,' % (NEW, NEW), t)
    old_call = '  Seam%s.engine_%sx cl m k v' % (NEW, NEW)
    if t.count(old_call) != 1:
        die('seam call')
    t = t.replace(old_call, '  Seam%s.engine_%s cl m k v' % (NEW, NEW))
    return t


if __name__ == '__main__':
    main()
