#!/usr/bin/python3
"""fs-build (13 Oct 2026): k-build's kb_submission.py (m-build's mb_submission.py), unchanged except the adapter's credit comment.
k-build (13 Oct 2026): m-build's mb_submission.py, unchanged except the adapter's credit comment.
m-build (13 Oct 2026): package the C588M proof closure the way nmamano packaged Shea's C588 for Big O Speedrun
(bigospeedrun/verifier submissions/dft/nmamano-shea-c588): every local module of the import closure of the
solution module is copied under Submission/, with ONLY its local import lines prefixed by `Submission.`; an adapter
Submission.lean derives `Race.bound` at the given bound from transform_main<SUF>.
usage: mb_submission.py <tree lean/complex588> <solution module> <suffix, e.g. C588M> <bound n/d> <stage dir>"""
import re, sys
from fractions import Fraction
from pathlib import Path

T, ROOT, SUF, BOUND, S = Path(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], Path(sys.argv[5])
NAME = re.compile(r'[A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)*\Z')

def strip_comments(text):  # bigospeedrun/verifier scripts/verifier.py, verbatim logic
    out, depth, i = [], 0, 0
    while i < len(text):
        if text.startswith('/-', i):
            depth += 1; out.append('  '); i += 2
        elif depth and text.startswith('-/', i):
            depth -= 1; out.append('  '); i += 2
        elif not depth and text.startswith('--', i):
            end = text.find('\n', i)
            i = len(text) if end < 0 else end
        else:
            out.append(text[i] if not depth or text[i] == '\n' else ' '); i += 1
    if depth:
        raise ValueError('unclosed block comment')
    return ''.join(out)

def path_of(m):
    return T / (m.replace('.', '/') + '.lean')

def header_imports(text):
    """(line index, [modules]) of the import header, as the verifier's plan_sources reads it."""
    res = []
    for k, line in enumerate(strip_comments(text).splitlines()):
        s = line.strip()
        if not s or s in ['prelude', 'module']:
            continue
        h = re.fullmatch(r'(?:(?:public|private|meta)\s+)*import\s+(.+)', s)
        if not h:
            break
        names = h[1].split()
        assert names and all(NAME.fullmatch(n) for n in names), (k, s)
        res.append((k, names))
    return res

closure, order = set(), []
def visit(m):
    if m in closure or not path_of(m).is_file():
        return
    closure.add(m)
    for _, names in header_imports(path_of(m).read_text()):
        for n in names:
            visit(n)
    order.append(m)
visit(ROOT)
assert ROOT in closure, 'solution module missing'
S.mkdir(parents=True, exist_ok=True)
changed = 0
for m in sorted(closure):
    text = path_of(m).read_text()
    lines = text.split('\n')
    for k, names in header_imports(text):
        orig = lines[k]
        for n in names:
            if n in closure:
                lines[k] = re.sub(r'(?<![\w.])' + re.escape(n) + r'(?![\w.])', 'Submission.' + n, lines[k], count=1)
        if lines[k] != orig:
            changed += 1
    out = S / 'Submission' / (m.replace('.', '/') + '.lean')
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text('\n'.join(lines))
n, d = Fraction(BOUND).as_integer_ratio()
adapter = f'''/-
Local packaging by fs-build (13 Oct 2026) for an offline re-creation of the Big O Speedrun verifier check;
not a submission. Same packaging as nmamano's submissions/dft/nmamano-shea-c588 (m-build's script): the proof
closure of {ROOT} is copied under Submission/, changing local import paths only.
Same-run hosting: DaysSky's #442 hosting argument (x-hosting; h-L3's hosted network Work/XHost); mutable sources
with 5 centre-source copies kept (x-mutable); exit frames (x-exit-lean, k-lower, k-combine).
Credits: the hub's U8 word (Jeremy Kothe #869/#870, Chafik Boukhalfa #890, Capalbo #942, on Kothe's #799/#588
lineage), FanPUYen/Hao Yan (owners), the hub authors (exit frames and tail #591/#593/#599); Ryan Shea's Lean
instantiation (generators, bank and seam pipeline); Nil Mamano's gauged Lean extension and record files;
Dugongue (#401) and jacklightChen (#377), Zhihao Chen (#425); OpenAI's model and Fourier reduction; Jacob Sussman's
engine, checkers and generators, with the exit-frame (XExit, KLow, KMX) and mutable-source (XMut) engine sides.
-/
import Submission.{ROOT}

open OAI.PowerSaving Filter Asymptotics

theorem Race.bound : ∃ order solve W, DFTProgram order solve W ∧
    (fun n : ℕ => (W n : ℝ)) =O[atTop]
      (fun n : ℕ => (n : ℝ) * (Real.log (n : ℝ)) ^
        (1 - ({n} : ℝ) / {d})) := by
  obtain ⟨order, solve, W, hprogram, htime⟩ := transform_main{SUF}
  refine ⟨order, solve, W, hprogram, ?_⟩
  have hdecimal : decimalTime{SUF} =
      (fun n : ℕ => (n : ℝ) * (Real.log (n : ℝ)) ^
        (1 - ({n} : ℝ) / {d})) := by
    funext n
    have hexponent : decimalExponent{SUF} = 1 - ({n} : ℝ) / {d} := by
      norm_num [decimalExponent{SUF}]
    simp only [decimalTime{SUF}, hexponent]
  exact hdecimal ▸ htime.2.1
'''
(S / 'Submission.lean').write_text(adapter)
print(f'closure modules: {len(closure)}; import lines rewritten: {changed}; bound {n}/{d}')
