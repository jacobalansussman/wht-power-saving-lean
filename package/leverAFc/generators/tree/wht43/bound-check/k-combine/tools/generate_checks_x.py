#!/usr/bin/env python3
"""k-lower: the Fourier challenge at the new exponent, the one-unit-stronger negative challenge, the axiom module and
the comparator configuration of an exit-frame certificate (m-gen's generate_checks_z.py / Shea's generate_checks.py,
without the scalar controls).  usage: generate_checks_x.py --out <global dir> --original-challenge <UniformFourier.lean>
--comparator <Compare.lean> [--name X4X]"""
import argparse, json, re, shutil
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--original-challenge", type=Path, required=True)
    ap.add_argument("--comparator", type=Path, required=True)
    ap.add_argument("--name", default="X4X")
    a = ap.parse_args()
    NAME = a.name
    info = json.loads((a.out / ("logs/global-%s.json" % NAME)).read_text())
    numerator = info["dft_numerator"]
    source = a.original_challenge.read_text()
    needle = "1 - 1/(10^(13:ℕ))"
    assert source.count(needle) == 1
    names = ["decimalExponent", "decimalTime", "TimeBounds", "DFTGoal", "ConvGoal", "transform_main", "convolution_main"]
    for negative in (False, True):
        txt = source.replace(needle, f"1 - {numerator + int(negative)}/(10^(10:ℕ))")
        for name in names:
            txt = re.sub(r"\b" + name + r"\b", name + NAME, txt)
        d = a.out / "Work" / (("Fourier%sNeg" if negative else "Fourier%s") % NAME)
        d.mkdir(parents=True, exist_ok=True)
        (d / "UniformFourierChallenge.lean").write_text(txt)
    (a.out / ("Work/Fourier%s/Axioms.lean" % NAME)).write_text("import Work.Fourier%s.Main\n" % NAME + "\n".join(
        "#print axioms " + n for n in ["OAI.PowerSaving.transform_main" + NAME, "OAI.PowerSaving.convolution_main" + NAME,
                                       "OAI.PowerSaving.RAM.hills%s_program" % NAME]) + "\n")
    (a.out / "tools").mkdir(exist_ok=True)
    shutil.copy2(a.comparator, a.out / "tools/Compare.lean")
    (a.out / "comparator").mkdir(exist_ok=True)
    (a.out / ("comparator/fourier-%s.json" % NAME.lower())).write_text(json.dumps(dict(
        challenge_module="Work.Fourier%s.UniformFourierChallenge" % NAME, solution_module="Work.Fourier%s.Main" % NAME,
        theorem_names=["OAI.PowerSaving.transform_main" + NAME, "OAI.PowerSaving.convolution_main" + NAME],
        definition_names=[], permitted_axioms=["propext", "Quot.sound", "Classical.choice"], enable_nanoda=True), indent=2) + "\n")
    print("challenge exponent 1 - %d/10^10 (negative %d)" % (numerator, numerator + 1))

if __name__ == "__main__":
    main()
