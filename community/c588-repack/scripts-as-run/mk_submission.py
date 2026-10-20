# c-pack: build a bigospeedrun-style submission folder for x-pack's B-quartic overlay, from
#   (a) nmamano's published packaging submissions/dft/nmamano-shea-c588 (verifier c92b9b6) and
#   (b) Shea's lean/complex588 tree (the files that submission copies) and the overlay files.
# Rule of the packaging (checked on the control first): submission file = Shea's file with every line
# "import X" turned into "import Submission.X".  Usage:
#   python3 -I mk_submission.py <nmamano submission dir> <shea complex588 dir> <overlay dir (has Work/)> <out dir> <num> <den>
import sys, os, re, shutil, hashlib
sub, shea, ov, out, num, den = sys.argv[1:7]
def tx(text):
    return re.sub(r'(?m)^import (?!Mathlib)(\S+)', r'import Submission.\1', text)
def rd(p): return open(p, encoding='utf-8').read()
files = []
for root, _, fs in os.walk(os.path.join(sub, 'Submission')):
    for f in fs:
        if f.endswith('.lean'):
            files.append(os.path.relpath(os.path.join(root, f), os.path.join(sub, 'Submission')))
files.sort()
# control: every submission file equals tx(Shea's file)
bad = [f for f in files if not os.path.exists(os.path.join(shea, f)) or tx(rd(os.path.join(shea, f))) != rd(os.path.join(sub, 'Submission', f))]
print(f"control: {len(files)} submission modules; {len(files)-len(bad)} equal tx(Shea's file); mismatches: {bad[:10]}")
assert not bad, "packaging rule does not hold"
# Mathlib imports untouched? check rule on one file with Mathlib import
ovfiles = []
for root, _, fs in os.walk(os.path.join(ov, 'Work')):
    for f in fs:
        if f.endswith('.lean'):
            ovfiles.append(os.path.relpath(os.path.join(root, f), ov))
ovfiles.sort()
if os.path.exists(out): shutil.rmtree(out)
os.makedirs(os.path.join(out, 'Submission'))
changed = new = 0
extra = [f for f in ovfiles if f not in files]
print('overlay files outside the submission closure:', extra)
keep = set(files) | {f for f in extra if f == 'Work/BlockAccounting/Moment4.lean'}
for f in sorted(keep):
    src = os.path.join(ov, f) if f in ovfiles else os.path.join(shea, f)
    dst = os.path.join(out, 'Submission', f)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    open(dst, 'w', encoding='utf-8').write(tx(rd(src)))
    if f in ovfiles:
        if f in files: changed += 1
        else: new += 1
print(f"overlay: {changed} files replaced, {new} new (Moment4)")
# The new module must be reached from the closure: RateC588x imports it in the overlay.
a = rd(os.path.join(sub, 'Submission.lean'))
old_n, old_d = '2466337', '2000000000'
print('bound occurrences in adapter:', a.count(f'({old_n} : ℝ) / {old_d}')); assert a.count(f'({old_n} : ℝ) / {old_d}') == 3
a = a.replace(f'({old_n} : ℝ) / {old_d}', f'({num} : ℝ) / {den}')
a = a.replace(f'Saving: {old_n}/{old_d}.', f'Saving: {num}/{den}. LOCAL CHECK by c-pack (13 Oct 2026), not a submission: x-pack 5,009-bank table + term_le4 on this packaging.')
open(os.path.join(out, 'Submission.lean'), 'w', encoding='utf-8').write(a)
open(os.path.join(out, 'manifest.json'), 'w').write('{"problem":"dft","bound":"%s/%s"}' % (num, den))
n = sum(1 for r, _, fs in os.walk(os.path.join(out, 'Submission')) for f in fs if f.endswith('.lean'))
print(f"out: {n} modules under Submission/ + Submission.lean; manifest bound {num}/{den}")
