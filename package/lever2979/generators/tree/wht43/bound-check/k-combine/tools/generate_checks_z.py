#!/usr/bin/env python3
"""m-gen (13 Oct 2026): the full Fourier challenge at the new exponent, the one-unit-stronger negative challenge,
the axiom module, the Lean negative controls, the comparator configuration, for the certificate `C588M`.

Ryan Shea's research/complex-transfer-588-596/experiments/generate_checks.py (shea256/fourier-transform-below-nlogn
@ 5239484), parametrized by the name, with these changes:
* the controls: Shea's corrupted first A gate is NOT rejected by the copy-free word at batch (0, 28) (its first A
  gates move entrance dirt between slots); the corrupted datum is the sign of the first coefficient of the first
  scatter row (`bad_scalar_rejected`, GS.scalCheckZ); new `classic_rejected` (Sussman's classic GS.scalCheck
  refuses this certificate: it writes sources); `duplicate_owner_rejected` as Shea, with this bank count.
  All three values are predicted FALSE by the exact mirrors (control_facts.py) before any build;
* a comparator configuration comparator/fourier-<name>.json next to Shea's comparator/fourier.json.
usage: generate_checks_z.py --out <dir> --original-challenge <UniformFourier.lean> --comparator <Compare.lean>
       [--name C588M]"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--original-challenge", type=Path, required=True)
    ap.add_argument("--comparator", type=Path, required=True)
    ap.add_argument("--name", default="C588M")
    a = ap.parse_args()
    NAME = a.name
    info = json.loads((a.out / ("logs/global-%s.json" % NAME)).read_text())
    numerator = info["dft_numerator"]
    nb = info["banks"]
    source = a.original_challenge.read_text()
    needle = "1 - 1/(10^(13:ℕ))"
    assert source.count(needle) == 1
    names = ["decimalExponent", "decimalTime", "TimeBounds", "DFTGoal", "ConvGoal", "transform_main", "convolution_main"]
    for negative in (False, True):
        txt = source.replace(needle, f"1 - {numerator + int(negative)}/(10^(10:ℕ))")
        for name in names:
            txt = re.sub(r"\b" + name + r"\b", name + NAME, txt)
        directory = a.out / "Work" / (("Fourier%sNeg" if negative else "Fourier%s") % NAME)
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "UniformFourierChallenge.lean").write_text(txt)
    p = a.out / ("Work/Fourier%s/Axioms.lean" % NAME)
    p.write_text("import Work.Fourier%s.Main\n" % NAME + "\n".join("#print axioms " + n for n in [
        "OAI.PowerSaving.transform_main" + NAME, "OAI.PowerSaving.convolution_main" + NAME,
        "OAI.PowerSaving.RAM.hills%s_program" % NAME]) + "\n")
    controls = f'''import Work.GCert.Data.Gen.{NAME}Scal0
import Work.GCert.Data.Gen.{NAME}Sched
set_option maxRecDepth 100000
namespace GXD.{NAME}Controls
open GXD
/-- the sign of the first coefficient of the first scatter row flipped -/
def flipRows : List (List (Nat × Co)) → List (List (Nat × Co))
  | ((k, c) :: r) :: rows => ((k, {{c with neg := !c.neg}}) :: r) :: rows
  | rows => rows
def flipScat : Scat → Scat
  | .table rows => .table (flipRows rows)
  | s => s
noncomputable def badRaw : Raw := {{{NAME}.raw with scat := flipScat {NAME}.raw.scat}}
/-- the mutable-source scalar check rejects a corrupted scatter coefficient (y_0 no longer receives exactly x_0) -/
theorem bad_scalar_rejected : GS.scalCheckZ badRaw ⟨10, 21, 4, 8, 8⟩ 0 28 = false := by decide +kernel
/-- Sussman's classic scalar check refuses the copy-free certificate itself: it writes sources -/
theorem classic_rejected : GS.scalCheck {NAME}.raw ⟨10, 21, 4, 8, 8⟩ 0 28 = false := by decide +kernel
noncomputable def badRows : List (List Nat) := match {NAME}.schedRows with
  | (x::_y::xs)::rs => (x::x::xs)::rs
  | other => other
theorem duplicate_owner_rejected : OAI.PowerSaving.GX.schedK {nb} {NAME}.raw.R badRows {NAME}.schedPad {NAME}.schedRes = false := by decide +kernel
end GXD.{NAME}Controls
#print axioms GXD.{NAME}Controls.bad_scalar_rejected
#print axioms GXD.{NAME}Controls.classic_rejected
#print axioms GXD.{NAME}Controls.duplicate_owner_rejected
'''
    (a.out / ("Work/GCert/Data/%sControls.lean" % NAME)).write_text(controls)
    (a.out / "tools").mkdir(exist_ok=True)
    shutil.copy2(a.comparator, a.out / "tools/Compare.lean")
    (a.out / "comparator").mkdir(exist_ok=True)
    (a.out / ("comparator/fourier-%s.json" % NAME.lower())).write_text(json.dumps(dict(
        challenge_module="Work.Fourier%s.UniformFourierChallenge" % NAME, solution_module="Work.Fourier%s.Main" % NAME,
        theorem_names=["OAI.PowerSaving.transform_main" + NAME, "OAI.PowerSaving.convolution_main" + NAME],
        definition_names=[], permitted_axioms=["propext", "Quot.sound", "Classical.choice"], enable_nanoda=True),
        indent=2) + "\n")
    (a.out / ("logs/challenge-provenance-%s.json" % NAME)).write_text(json.dumps(dict(
        source=a.original_challenge.name, source_sha256=hashlib.sha256(a.original_challenge.read_bytes()).hexdigest(),
        changes={"exponent": f"1 - {numerator}/10^10", "seven_names_suffix": NAME},
        negative_control=f"one unit stronger: {numerator + 1}/10^10"), indent=2) + "\n")


if __name__ == "__main__":
    main()
