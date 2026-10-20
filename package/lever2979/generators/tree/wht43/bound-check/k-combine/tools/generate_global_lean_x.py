#!/usr/bin/env python3
"""k-lower (13 Oct 2026): the bank witness, rate, glue and all-length Fourier instance of an EXIT-FRAME certificate.

m-gen's generate_global_lean_z.py (itself Ryan Shea's research/complex-transfer-588-596/experiments/
generate_global_lean.py, shea256/fourier-transform-below-nlogn), with these changes only:
* R = 616 slots (504 helpers + 112 source copies, as Shea's route), counts and the kernel numerator from the price
  file (make_sched_x.py: x-exit-words' pad-minimal table remapped);
* the glue builds the EXIT-FRAME package: `ginvX := GXD.ginvXOf raw scalOK.hv .. es es_ok aux el validX exitRC.1
  exitRC.2` (k-lower's Work/KLow/GlueX.lean on x-exit-lean's GLabX / GInvX), and Shea's gauged package
  `ginv := ginvX.toE'` (x-exit-lean's reduction XExit/Pack: same price, stand-in entrance labels of dimension
  9 - (dim exit - dim start)); the bank residuals are the EXIT residuals (`res_ok` on the start and end states,
  `GXD.ginvXOf_res`, `GInvX.entS_hd`).  Everything downstream (rate, unit price, Fourier seam) is Shea's/the
  record's, unchanged: it consumes `ginv`, `sched`, `fitsS`, `price`, `width`;
* the Fourier files are the record's FourierC401G5 files renamed to the certificate name (Shea's map).
usage: generate_global_lean_x.py <price.json> --toolkit <toolkit/tools/gen> --record <record-templates> --out <dir>
       [--name X4X]"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

EPAR = None   # set from the local preflight (E7 SPar)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("price", type=Path)
    ap.add_argument("--toolkit", type=Path, required=True)
    ap.add_argument("--record", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--name", default="X4X")
    ap.add_argument("--local", type=Path, required=True, help="the out dir of generate_local_lean_x.py")
    ap.add_argument("--cert", type=Path, required=True)
    a = ap.parse_args()
    NAME = a.name
    p = json.loads(a.price.read_text())
    pre = json.loads((a.local / "logs" / ("preflight-%s.json" % NAME)).read_text())
    d_, sw_, ux_, us_, uy_ = pre["SPar"]
    EPAR = "⟨%d, %d, %d, %d, %d⟩" % (d_, pre["e7_sw"], ux_, us_, uy_)
    NEB = -(-len(json.loads(gzip.decompress(a.cert.read_bytes()))["ext"]) // 80)
    sched = json.loads(gzip.decompress(a.price.with_suffix(".schedule.json.gz").read_bytes()))
    sys.path.insert(0, str(a.toolkit))
    import foldlib as fl
    import foldgen as fg
    u = dict(m=p["m"], W=p["W"], D=p["deficit"], H={int(k): v for k, v in p["profile"].items()})
    R = p["local_residual_count"]
    assert (u["m"], p["v"], len(sched["widths"]), R) == (73, 112, 616, 616)
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
    rate += fg.list_text(NAME + "x", u, "72 copies of the exit-frame certificate %s (hub tail in phase B), #401 dedicated banks, "
                                        "and an explicit %s bank schedule." % (NAME, p["schedule"]))
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
    # ---- the glue: the exit-frame package and its reduction to Shea's gauged package
    rp = ", ".join(f"rankP{r}_ok" for r in range(top + 1))
    e7 = json.loads((a.out.parent / "local" / "logs" / ("preflight-%s.json" % NAME)).read_text()) if False else None
    glue = f"""import Work.KLow.GlueX
import Work.Bridge401.GGeom
import Work.GCert.Data.{NAME}Cert
import Work.GCert.Data.{NAME}ScalOK
import Work.GCert.Data.Gen.{NAME}EChk
import Work.GCert.Data.Gen.{NAME}ZOK
import Work.GCert.Data.Inst{NAME}
import Work.GCert.Data.Gen.{NAME}Sched

/-!
# (key: k-lower) the EXIT-FRAME invocation package and the bank schedule of `GXD.{NAME}` (PR #401 layout)
Generated by k-lower's generate_global_lean_x.py (on Ryan Shea's C588Glue401 / Nil Mamano's G5 glue).
`ginvX` is k-lower's `GXD.ginvXOf` (Work/KLow/GlueX.lean: x-exit-lean's `GLabX` over DFT Bridge3's gauged scalar
half) on the kernel checks of this certificate: scalar half `scalOK`, E7 `echk<j>_ok`, label half `validX`
(exit-aware end check `endKX`), the two exit conditions `exitRC` (zero checks on the y-free twin, `Gen.{NAME}ZOK`).
`ginv := ginvX.toE'` is x-exit-lean's reduction (XExit/Pack) to Shea's gauged package: same price; the bank loads are
the exit residuals.
-/

set_option maxRecDepth 100000
namespace GXD.{NAME}
open GXD OAI.PowerSaving OAI.PowerSaving.CB
noncomputable section

/-- the E7 batches: parameters and gauge registers -/
def es : List (GS.SPar × List Nat) := [{', '.join('(%s, glb%d)' % (EPAR, i) for i in range(NEB))}]

theorem es_gl : es.flatMap Prod.snd = gl := by first | rfl | decide

theorem es_ok : ∀ e ∈ es, GS.eCheck raw e.1 e.2 = true := by
  intro e he
  unfold es at he
{''.join('  rcases List.mem_cons.mp he with h | he' + chr(10) + '  · subst h; exact echk%d_ok' % i + chr(10) for i in range(NEB))}  exact absurd he List.not_mem_nil

/-- **the exit-frame invocation package** of the certificate -/
def ginvX : OAI.PowerSaving.XEX.GInvX (Fin aux.p.h) (Fin raw.v) (Fin raw.R) (Fin raw.ret.length) :=
  GXD.ginvXOf raw scalOK.hv scalOK.hAi scalOK.hBi scalOK.hx scalOK.hy scalOK.hid es es_ok aux el
    (es_gl.symm ▸ validX) exitRC.1 exitRC.2

/-- **the gauged invocation package** (x-exit-lean's reduction `GInvX.toE'`) -/
def ginv : GX.GInvE (Fin aux.p.h) (Fin raw.v) (Fin raw.R) (Fin raw.ret.length) := ginvX.toE'

theorem h9 : Fintype.card (Fin aux.p.h) = 9 := by rw [Fintype.card_fin]; rfl

/-- its price is the invocation price `invPrice` (kernel counts) -/
theorem ginv_cost (φ : ℕ → ℝ) : ginv.cost φ = invPrice φ := by
  have h5 : ∀ k : Fin raw.ret.length,
      bcost φ (GLab.cdim (aux.p.tab.get (aux.SA.get (GXD.reg raw k)))) = bcost φ 8 :=
    fun k => congrArg (bcost φ) (cen_dim _ (GXD.reg_mem raw k))
  show ginvX.cost φ = _
  rw [ginvX, GXD.ginvXOf_cost, Finset.sum_congr rfl (fun k _ => h5 k), Finset.sum_const, Finset.card_univ,
    nsmul_eq_mul, Fintype.card_fin, ← inv_cost φ, show raw.ret.length = 8 from by decide]
  norm_num [bcost]
  try rfl

/-- the bank schedule of the certificate (tables `Gen.{NAME}Sched`, checker `GX.schedK` of DFT New Sums) -/
def sched : GX.BankSched (Fin raw.R) := GX.BankSched.ofTab {nb} raw.R schedRows schedPad schedRes sched_ok

/-- the EXIT residuals of the table (start and end states of the label check) -/
theorem res_ok : ∀ q : Fin raw.R, Bool.rec (motive := fun _ => Nat)
    (aux.p.h - GLab.cdim (aux.p.tab.get (aux.S0.get (2 * raw.v + q.val))))
    (GLab.cdim (aux.p.tab.get (aux.SF.get (2 * raw.v + q.val)))
      - GLab.cdim (aux.p.tab.get (aux.S0.get (2 * raw.v + q.val))))
    (GLab.memK (2 * raw.v + q.val) el) = schedRes.getD q.val 0 := by
  decide +kernel

theorem hent : ∀ q : Fin raw.R, Fintype.card (Fin aux.p.h) - (ginv.ent q).dim = schedRes.getD q.val 0 := fun q => by
  show Fintype.card (Fin aux.p.h) - (ginvX.entS q).dim = _
  rw [ginvX.entS_hd q]
  exact (GXD.ginvXOf_res raw scalOK.hv scalOK.hAi scalOK.hBi scalOK.hx scalOK.hy scalOK.hid es es_ok aux el
    (es_gl.symm ▸ validX) exitRC.1 exitRC.2 q).trans (res_ok q)

theorem fitsS : sched.fits ginv.ent := GX.BankSched.ofTab_fits _ _ _ _ _ sched_ok ginv.ent hent

theorem width : 64 * Fintype.card (Fin raw.v) + sched.nb = {p['W']} := by
  rw [Fintype.card_fin, sched, GX.BankSched.ofTab_nb, raw_v]

theorem pad_price (φ : ℕ → ℝ) : ∑ b, bcost φ (sched.pad b) = BridgeRate.{NAME}.padPrice φ := by
  rw [sched, GX.BankSched.ofTab_padPrice, GXD.rcost_wsumK φ schedPad {top} rankBoundP]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, {rp}]
  unfold BridgeRate.{NAME}.padPrice
  norm_num [bcost]

/-- the unit price identity of the PR #401 unit (`GX.unitPrice401` shape) -/
theorem price (φ : ℕ → ℝ) : 72 * ginv.cost φ + (Fintype.card (Fin raw.v) : ℝ) * GX.bankPrice377 φ
    + ∑ b, bcost φ (sched.pad b) = BridgeRate.{NAME}.unitPrice φ := by
  rw [ginv_cost, pad_price, Fintype.card_fin, raw_v]
  unfold BridgeRate.{NAME}.unitPrice BridgeRate.{NAME}.invPrice invPrice
  push_cast
  ring

end
end GXD.{NAME}
#print axioms GXD.{NAME}.ginv_cost
#print axioms GXD.{NAME}.price
#print axioms GXD.{NAME}.fitsS
"""
    (base / f"{NAME}Glue401.lean").write_text(glue)
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
