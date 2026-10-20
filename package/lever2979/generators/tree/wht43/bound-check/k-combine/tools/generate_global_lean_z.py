#!/usr/bin/env python3
"""m-gen (13 Oct 2026): the bank witness, rate, glue and all-length Fourier instance of `C588M` (copy-free #588).

Ryan Shea's research/complex-transfer-588-596/experiments/generate_global_lean.py (shea256/fourier-transform-below-nlogn
@ 5239484), parametrized by the certificate name (C588M) and the price file, with these changes only:
* R = 504 helpers (no source copies), the bank count, W, the pad ranks and the kernel numerator come from the price
  file (make_sched.py: annealed or cyclic schedule);
* the glue (from the G5 template of Nil Mamano's record submission, as Shea) builds the gauged invocation package
  with x-mutable's mutable-source glue: `GXD.ginvOfZ raw scalOKZ es es_ok aux validE`, `GXD.ginvOfZ_cost`,
  `GXD.ginvOfZ_ent_dim` (Work/XMut/Glue.lean), E7 batches `GS.eCheckZ` at SPar (10, 8, 4, 8, 8);
* the Fourier files are the record's FourierC401G5 files renamed to C588M (Shea's substitution map, target C588M).
usage: generate_global_lean_z.py <price.json> --toolkit <toolkit/tools/gen> --record <record-templates> --out <dir>
       [--name C588M]"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

EPAR = "⟨10, 8, 4, 8, 8⟩"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("price", type=Path)
    ap.add_argument("--toolkit", type=Path, required=True)
    ap.add_argument("--record", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--name", default="C588M")
    a = ap.parse_args()
    NAME = a.name
    p = json.loads(a.price.read_text())
    sched = json.loads(gzip.decompress(a.price.with_suffix(".schedule.json.gz").read_bytes()))
    sys.path.insert(0, str(a.toolkit))
    import foldlib as fl
    import foldgen as fg
    u = dict(m=p["m"], W=p["W"], D=p["deficit"], H={int(k): v for k, v in p["profile"].items()})
    R = p["local_residual_count"]
    assert (u["m"], p["v"], len(sched["widths"]), R) == (73, 112, 504, 504)
    from fractions import Fraction as Q
    A = int(Q(p["exact_rate"]["root_lower"]) * 10 ** 10)
    while fl.plan_block(u, A, 10, 40) is None:
        A -= 1
    assert A == p["kernel_numerator"]
    delta = A - 1
    base = a.out / "Work/GCert/Data"
    (base / "Gen").mkdir(parents=True, exist_ok=True)
    body, plan = fg.block_text(NAME + "x", u, A, 10, 40)
    rate = "import Work.BlockAccounting.Moment\nset_option maxRecDepth 4000\nnamespace OAI.PowerSaving.BlockAccounting\nopen Finset\n"
    rate += fg.list_text(NAME + "x", u, "72 copies of PR588 atm3 WITHOUT source copies (mutable sources), #401 dedicated banks, "
                                        "and an explicit %s bank schedule." % p["schedule"])
    rate += "\n" + body + "\nend OAI.PowerSaving.BlockAccounting\n#print axioms OAI.PowerSaving.BlockAccounting.fold_%sx\n" % NAME
    (base / ("Rate%sx.lean" % NAME)).write_text(rate)
    rows, pads, widths = sched["rows"], sched["pads"], sched["widths"]
    nb = len(pads)
    top = max(pads)
    counts = Counter(pads)
    assert nb == p["banks"] and p["W"] == 64 * 112 + nb
    imports = []
    parts = []
    for i, lo in enumerate(range(0, 72, 4)):
        mod = "Work.GCert.Data.Gen.%sRows%d" % (NAME, i)
        imports.append("import " + mod)
        parts.append(f"rows{i}")
        txt = "import Work.GCert.Data.Gen.%s\nnamespace GXD.%s\n" % (NAME, NAME)
        txt += f"noncomputable def rows{i} : List (List Nat) := {rows[lo:lo + 4]}\nend GXD.{NAME}\n"
        (base / f"Gen/{NAME}Rows{i}.lean").write_text(txt)
    txt = "\n".join(imports) + "\nimport Work.Bridge401.SchedK\nimport Work.GCert.Data.Count\nset_option maxRecDepth 100000\nnamespace GXD.%s\nopen GXD OAI.PowerSaving\n" % NAME
    txt += "noncomputable def schedRows : List (List Nat) := " + " ++ ".join(parts) + "\n"
    txt += f"noncomputable def schedPad : List Nat := {pads}\nnoncomputable def schedRes : List Nat := {widths}\n"
    txt += f"theorem sched_ok : GX.schedK {nb} raw.R schedRows schedPad schedRes = true := by decide +kernel\n"
    txt += f"theorem rankBoundP : wsumK (SSC.wGt {top}) schedPad 0 = 0 := by decide +kernel\n"
    txt += "\n".join(f"theorem rankP{r}_ok : wsumK (SSC.wEq {r}) schedPad 0 = {counts[r]} := by decide +kernel" for r in range(top + 1))
    txt += "\nend GXD.%s\n#print axioms GXD.%s.sched_ok\n" % (NAME, NAME)
    (base / f"Gen/{NAME}Sched.lean").write_text(txt)
    inv = " + ".join(f"{n} * φ {r}" for r, n in sorted((int(r), n) for r, n in p["local_profile"].items()))
    pad = " + ".join(f"{n} * bcost φ {r}" for r, n in sorted(counts.items()) if r and n)
    txt = f'''import Work.Bridge377.Bank
import Work.FoldRate.Group
import Work.GCert.Data.Rate{NAME}x
namespace OAI.PowerSaving.BridgeRate.{NAME}
open Finset Binary Matrix RAM CB FoldRate
noncomputable section
def invPrice (φ : ℕ → ℝ) : ℝ := {inv}
def padPrice (φ : ℕ → ℝ) : ℝ := {pad}
def unitPrice (φ : ℕ → ℝ) : ℝ := 72 * invPrice φ + ((112:ℕ):ℝ) * GX.bankPrice377 φ + padPrice φ
theorem unit_cost_x (φ : ℕ → ℝ) : unitPrice φ = listCost φ BlockAccounting.blocks{NAME}x := by
  unfold unitPrice invPrice padPrice
  simp only [GX.bankPrice377, listCost, BlockAccounting.blocks{NAME}x, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  norm_num [bcost]
  try ring
theorem unit_tally : unitPrice (fun r => (r:ℝ)) = (({p['W']}:ℕ):ℝ) * 73 - 2560 := by
  rw [unit_cost_x]
  exact BlockAccounting.blocks{NAME}x_tally
end
end OAI.PowerSaving.BridgeRate.{NAME}
#print axioms OAI.PowerSaving.BridgeRate.{NAME}.unit_cost_x
'''
    (base / f"Inst{NAME}.lean").write_text(txt)
    # ---- the glue: Shea's substitutions on the G5 template, then the mutable-source package
    src = a.record / "Work/GCert/Data/G5Glue401.lean"
    txt = src.read_text().replace("import Submission.", "import ").replace("C401Gg5", NAME).replace("G5", NAME)
    txt = txt.replace("6673", str(nb)).replace("14353", str(p["W"])).replace("Fin 772", "Fin %d" % R).replace("2 * 120", "2 * 112")
    txt = txt.replace("bcost φ 9 :=", "bcost φ 8 :=")
    txt = re.sub(r"def es : List \(GS.SPar × List Nat\) := .*",
                 "def es : List (GS.SPar × List Nat) := [(%s, glb0), (%s, glb1), (%s, glb2)]" % (EPAR, EPAR, EPAR), txt)
    txt = txt.replace("  exact absurd he List.not_mem_nil",
                      "  rcases List.mem_cons.mp he with h | he\n  · subst h; exact echk2_ok\n  exact absurd he List.not_mem_nil")
    txt = txt.replace("schedPad 46 rankBoundP", f"schedPad {top} rankBoundP")
    txt = re.sub(r"rankP0_ok, .*rankP46_ok", ", ".join(f"rankP{r}_ok" for r in range(top + 1)), txt)
    txt = txt.replace("(SPar ⟨10, 6, 1, 4, 4⟩, ⟨10, 5, 1, 4, 4⟩)", "(three batches, SPar %s, sources mutable)" % EPAR)
    txt = txt.replace("Generated by tools/gx/gxgen401.py (DFT Lean).",
                      "Adapted from the G5 glue by DFT Lean in Nil Mamano's record submission as in Ryan Shea's "
                      "generate_global_lean.py; mutable sources: m-gen's generate_global_lean_z.py.")
    for old, new in [
            ("import Work.Gauge.Glue\n", "import Work.Gauge.Glue\nimport Work.XMut.Glue\n"),
            ("`ginv` is DFT Bridge3's `GXD.ginvOf` (Work/Gauge/Glue.lean) on the\nkernel checks of this certificate: scalar half `scalOK`",
             "`ginv` is x-mutable's `GXD.ginvOfZ` (Work/XMut/Glue.lean, on DFT Bridge3's `GXD.ginvOf`) on the\nkernel checks of this "
             "certificate: scalar half `scalOKZ` (sources mutable)"),
            ("GS.eCheck raw e.1 e.2", "GS.eCheckZ raw e.1 e.2"),
            ("GXD.ginvOf raw scalOK.hv scalOK.hAi scalOK.hBi scalOK.hx scalOK.hy scalOK.hid es es_ok aux (es_gl.symm ▸ validE)",
             "GXD.ginvOfZ raw scalOKZ es es_ok aux (es_gl.symm ▸ validE)"),
            ("GXD.ginvOf_cost", "GXD.ginvOfZ_cost"),
            ("GXD.ginvOf_ent_dim", "GXD.ginvOfZ_ent_dim")]:
        assert txt.count(old) == 1, old
        txt = txt.replace(old, new)
    assert "ginvOf " not in txt and "GS.eCheck " not in txt and "scalOK." not in txt
    (base / f"{NAME}Glue401.lean").write_text(txt)
    # ---- the Fourier seam and transfer (Shea's map, target NAME)
    mapping = {"C401Gg5": NAME, "C401G5": NAME, "GXD.G5": "GXD." + NAME, ".G5": "." + NAME, "/G5": "/" + NAME,
               "10221077": str(A), "10221076": str(delta), "14353": str(p["W"])}
    dst = a.out / ("Work/Fourier%s" % NAME)
    dst.mkdir(parents=True, exist_ok=True)
    provenance = {}
    for source in sorted((a.record / "Work/FourierC401G5").glob("*.lean")):
        txt = source.read_text().replace("import Submission.", "import ")
        for k, v in mapping.items():
            txt = txt.replace(k, v)
        (dst / source.name).write_text(txt)
        provenance[source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
    (a.out / "logs").mkdir(exist_ok=True)
    (a.out / ("logs/global-%s.json" % NAME)).write_text(json.dumps(dict(
        kernel_numerator=A, dft_numerator=delta, denominator=10 ** 10, schedule=p["schedule"], banks=nb, W=p["W"],
        pad_hist={str(k): v for k, v in sorted(counts.items())}, rational_slack=plan["slack_rel"],
        fourier_source_sha256=provenance, fourier_substitutions=mapping), indent=2) + "\n")
    print(f"{NAME}: kernel saving {A}/10^10; DFT saving {delta}/10^10; banks {nb}; W {p['W']}; rational slack {plan['slack_rel']}")


if __name__ == "__main__":
    main()
