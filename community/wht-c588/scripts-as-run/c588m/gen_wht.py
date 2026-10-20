#!/usr/bin/env python3
"""gen_wht.py: generate the Walsh-Hadamard corollary files for one kernel certificate (swappable).

r-wht, 13 Oct 2026.  For a kernel program `ENGINE` (a theorem of the `hills_program` shape at
exponent 1 - NUM/10^10, as `SeamC588M.engine_C588Mx` of Work/FourierC588M/Seam.lean), write

  Work/WHT<TAG>/Main.lean       the one-term proof: BlockWHT.wht_main_of_engine applied to ENGINE
  Work/WHT<TAG>/Solution.lean   = published Work/GCert/Data/SolutionB2Ge8x.lean, names/exponent/import changed
  Work/WHT<TAG>/Challenge.lean  = published Work/GCert/Data/ChallengeB2Ge8x.lean, exponent and theorem name changed
  Work/WHT<TAG>/Axioms.lean     #print axioms of both theorems
  Work/WHT<TAG>Neg/Challenge.lean   negative control, exponent one unit stronger (NUM+1)
  comparator/WHT<TAG>.json, WHT<TAG>-nanoda.json, WHT<TAG>Neg.json   official comparator configurations
  MANIFEST.sha256

The challenge files are produced from OUR published challenge by textual substitution, and the script
CHECKS that, after undoing the substitution of the theorem name and the exponent, every line from
`namespace OAI` to the end is byte-identical to the published file (so the statement differs from the
published one only in the exponent and the theorem name).

usage: gen_wht.py --tag C588M --seam Work.FourierC588M.Seam --engine SeamC588M.engine_C588Mx \
         --num 13552607 --published <clone of wht-power-saving-lean> --note <kernel-note.txt> --out <dir>
"""
import argparse, hashlib, json, sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument('--tag', required=True)
ap.add_argument('--seam', required=True, help='Lean module that proves the engine theorem')
ap.add_argument('--engine', required=True, help='engine theorem, relative to namespace OAI.PowerSaving')
ap.add_argument('--num', required=True, type=int, help='kernel exponent is 1 - NUM/10^10')
ap.add_argument('--published', required=True)
ap.add_argument('--note', required=True, help='text file: the paragraph describing the kernel program and its credit')
ap.add_argument('--out', required=True)
a = ap.parse_args()

PUB = Path(a.published)
pubC = (PUB / 'Work/GCert/Data/ChallengeB2Ge8x.lean').read_text()
pubS = (PUB / 'Work/GCert/Data/SolutionB2Ge8x.lean').read_text()
note = Path(a.note).read_text().rstrip('\n')
T, N = a.tag, a.num
D, DN = f'WHT{T}', f'WHT{T}Neg'
INNER = f'wht_main_{T}x'          # OAI.PowerSaving.WHT<TAG>.wht_main_<TAG>x
OUTER = f'wht_main_block_{T}x'    # OAI.PowerSaving.WHT.wht_main_block_<TAG>x
pubC_sha = hashlib.sha256(pubC.encode()).hexdigest()
pubS_sha = hashlib.sha256(pubS.encode()).hexdigest()

def must(cond, msg):
    if not cond:
        sys.exit('gen_wht: ' + msg)

# the published files have exactly these occurrences (guards against a changed base)
must(pubC.count('8762479') == 4 and pubC.count('B2Ge8x') == 4, 'published challenge is not the expected file')
must(pubS.count('8762479') == 3 and pubS.count('B2Ge8') == 7, 'published solution is not the expected file')

def challenge(num, mod):
    s = pubC
    s = s.replace('`Work.GCert.Data.ChallengeB2Ge8x`', f'`Work.{mod}.Challenge`')
    s = s.replace('`Work.GCert.Data.SolutionB2Ge8x`.', f'`Work.{D}.Solution`.')
    s = s.replace('B2Ge8x', f'{T}x').replace('8762479', str(num))
    head, sep, rest = s.partition('\n-/\n')
    must(sep != '', 'header comment end not found')
    head += (f'\n\nDerived 13 Oct 2026 from the published `Work/GCert/Data/ChallengeB2Ge8x.lean` of\n'
             f'github.com/jacobalansussman/wht-power-saving-lean (commit 9c94857, sha256 {pubC_sha[:8]}...{pubC_sha[-2:]}) by four textual\n'
             f'substitutions: the two module names in this comment, `B2Ge8x` -> `{T}x`, `8762479` -> `{num}`.')
    s = head + sep + rest
    # check: undo the name and exponent; the body must be the published body byte for byte
    body = s[s.index('\nnamespace OAI\n'):]
    pbody = pubC[pubC.index('\nnamespace OAI\n'):]
    must(body.replace(f'{T}x', 'B2Ge8x').replace(str(num), '8762479') == pbody, f'challenge body for {num} differs beyond name/exponent')
    nd = sum(1 for x, y in zip(body.split('\n'), pbody.split('\n')) if x != y)
    must(len(body.split('\n')) == len(pbody.split('\n')) and nd == 3, f'expected 3 changed body lines, got {nd}')
    return s

def solution():
    s = pubS
    s = s.replace('import Work.GCert.Data.B2Ge8Main', f'import Work.{D}.Main')
    s = s.replace('`Work.GCert.Data.SolutionB2Ge8x`', f'`Work.{D}.Solution`')
    s = s.replace('`Work.GCert.Data.B2Ge8Main`', f'`Work.{D}.Main`')
    s = s.replace('    theorem wht_main_B2Ge8x :', f'    theorem {D}.{INNER} :')
    s = s.replace('\n  wht_main_B2Ge8x\n', f'\n  {D}.{INNER}\n')
    s = s.replace('B2Ge8x', f'{T}x').replace('8762479', str(N))
    must('B2Ge8' not in s, 'solution substitution incomplete (B2Ge8 left)')
    s = s.replace('come from `WHTCheck.Solution`.\n-/',
                  'come from `WHTCheck.Solution`.\n'
                  'Same text as `Work/GCert/Data/SolutionB2Ge8x.lean` of wht-power-saving-lean (9c94857), with the\n'
                  'exponent, the theorem name and the import changed.\n-/')
    must('Same text as' in s, 'solution substitution incomplete (header note)')
    return s

def main():
    z = f'1 - {N}/(10:ℝ)^10'
    return f'''import {a.seam}
import Work.Block.WHT

/-!
# Walsh-Hadamard transform from the {T} kernel program

New file (13 Oct 2026); not part of github.com/openai/math and not part of Ryan Shea's tree
(github.com/shea256/fourier-transform-below-nlogn, lean/complex588 at 5239484), which it is laid over
without changing any of his files.

{note}

The Walsh-Hadamard corollary is `BlockWHT.wht_main_of_engine` (Work/Block/WHT.lean, Jacob Sussman's engine
file, byte-identical in Shea's tree and in wht-power-saving-lean 9c94857): pre-scale, kernel engine,
post-scale.

This file adds one term and nothing else.  No `sorry`, no `native_decide`, no axiom.
-/

namespace OAI
namespace PowerSaving
namespace {D}
open RAM

/-- **Walsh-Hadamard transform at `1 - {N}/10^10`**: `BlockWHT.wht_main_of_engine` applied to the
kernel program `{a.engine}` (colour `Paint.left`, transform mode `false`). -/
theorem {INNER} :
    ∃ solve W, WHT.WHTProgram solve W ∧ BlockWHT.WHTTimeBoundsAt ({z}) W :=
  BlockWHT.wht_main_of_engine ({z}) (by norm_num) (by norm_num)
    (fun k v => {a.engine} Paint.left false k v)

end {D}
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.{D}.{INNER}
'''

def cfg(chal, nanoda):
    return json.dumps({"challenge_module": chal, "solution_module": f"Work.{D}.Solution",
                       "theorem_names": [f"OAI.PowerSaving.WHT.{OUTER}"], "definition_names": [],
                       "permitted_axioms": ["propext", "Quot.sound", "Classical.choice"],
                       "enable_nanoda": nanoda}, indent=2) + '\n'

out = Path(a.out)
files = {
    f'Work/{D}/Main.lean': main(),
    f'Work/{D}/Solution.lean': solution(),
    f'Work/{D}/Challenge.lean': challenge(N, D),
    f'Work/{DN}/Challenge.lean': challenge(N + 1, DN),
    f'Work/{D}/Axioms.lean': f'import Work.{D}.Solution\n#print axioms OAI.PowerSaving.{D}.{INNER}\n#print axioms OAI.PowerSaving.WHT.{OUTER}\n',
    f'comparator/{D}.json': cfg(f'Work.{D}.Challenge', False),
    f'comparator/{D}-nanoda.json': cfg(f'Work.{D}.Challenge', True),
    f'comparator/{DN}.json': cfg(f'Work.{DN}.Challenge', False),
}
man = []
for p, txt in files.items():
    f = out / p
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(txt)
    man.append(f'{hashlib.sha256(txt.encode()).hexdigest()}  {p}')
(out / 'MANIFEST.sha256').write_text('\n'.join(man) + '\n')
(out / 'PARAMS.json').write_text(json.dumps({'tag': T, 'seam': a.seam, 'engine': a.engine, 'num': N,
    'theorem': f'OAI.PowerSaving.WHT.{OUTER}', 'inner': f'OAI.PowerSaving.{D}.{INNER}',
    'dir': D, 'negdir': DN, 'published_challenge_sha256': pubC_sha, 'published_solution_sha256': pubS_sha}, indent=2) + '\n')
print(f'wrote {len(files)} files to {out} (tag {T}, exponent 1 - {N}/10^10, negative {N+1})')
