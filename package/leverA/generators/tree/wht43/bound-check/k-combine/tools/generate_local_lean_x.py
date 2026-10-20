#!/usr/bin/env python3
"""k-lower (13 Oct 2026): the LOCAL Lean modules of an EXIT-FRAME certificate (exit helpers, hub tail in phase B).

Ryan Shea's research/complex-transfer-588-596/experiments/generate_local_lean.py (shea256/fourier-transform-below-nlogn
@ 0e61af6), which drives Jacob Sussman's gxgen.py / gxcore.py / gxbind.py UNCHANGED (copies in ../toolkit/tools/gx),
with these changes:
* Python mirror: E1 lets an exit slot (final frame != 1) end in a frame containing its start frame, E6 uses the
  invoice of the exit route (N + entrance ranks = sum over slots of dim(final) + 2 v (h-1) + cst) -- source
  substitution of gxcore.chk_replay, nothing else changes;
* the end check is k-lower's exit-aware `GLab.endKX ... gl el ...` (Work/KLow/LabelsX.lean), the label half
  `validX : GLab.ValidX raw aux gl el`;
* the movement invoice of the Price module: N (as Shea; the exit ranks are the final dims);
* the gauge list module also defines the exit list `el`;
* NEW: the zero checks of the two exit conditions (`GS.zCheck` on the y-free twin `GXD.nyRaw raw`, Work/KLow/ZCheck,
  Work/KLow/NoY) from zdesign.py's batches, in modules <Name>ZChk<i>, and their coverage facts in <Name>ZOK.
usage: generate_local_lean_x.py <cert.gcert.json.gz> --toolkit <toolkit/tools/gx> --verification <toolkit/verification>
       --zdesign <zdesign.json> --out <dir> [--name X4X]
"""
import argparse, gzip, json, os, sys, inspect
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificate", type=Path)
    ap.add_argument("--toolkit", type=Path, required=True)
    ap.add_argument("--verification", type=Path, required=True)
    ap.add_argument("--zdesign", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--name", default="X4X")
    a = ap.parse_args()
    NAME = a.name
    c = json.loads(gzip.decompress(a.certificate.read_bytes()))
    zd = json.loads(a.zdesign.read_text())
    sys.path.insert(0, str(a.verification.resolve()))
    from scalar_audit import audit
    result = audit(c)
    os.environ["GX_OUT"] = str(a.out.resolve())
    sys.path.insert(0, str(a.toolkit.resolve()))
    sys.path.insert(0, str(HERE))
    import gxcore as gc
    import gxgen as gen
    import gxbind as gb
    import gsmirror as sm
    import zmirror as zm
    gc.SC = 8
    gc.chk_shape = lambda D: result["checks"][1]["updates"]
    d, sw, ux, us, uy = result["SPar"]
    sm.replay(c, (ux, us, uy), "x")
    sm.replay(c, (ux, us, uy), "y")
    sm.bigreplay(c, (d, sw, ux, us, uy), 0, c["v"], [c["A"]], [c["B"]])
    gb.spar = lambda D, path: dict(d=d, sw=sw, ux=ux, us=us, uy=uy, source="exact rational and packed preflight",
                                   digits=gc.digits(D))
    # ---- exit-aware mirror (E1, E6)
    src = inspect.getsource(gc.chk_replay)
    a1 = '''        need(start[r] == ext.get(r, 0) and final[r] == 1, "E1: slot start / final")'''
    b1 = '''        need(start[r] == ext.get(r, 0) and (final[r] == 1 or (r in EXITS and all(reduce(u, fr[final[r]]) == 0 for u in fr[start[r]]))), "E1x: slot start / exit")'''
    a2 = '''    need(N == D["N"] and N + sum(dim[f] for f in ext.values()) == R * h + 2 * v * (h - 1) + D["cst"], "E6: N")'''
    b2 = '''    need(N == D["N"] and N + sum(dim[f] for f in ext.values()) == sum(dim[final[r]] for r in range(2 * v, n)) + 2 * v * (h - 1) + D["cst"], "E6x: N")'''
    assert src.count(a1) == 1 and src.count(a2) == 1
    src = src.replace(a1, b1).replace(a2, b2)
    ns = vars(gc)
    ns["EXITS"] = set(c["exits"])
    exec(src, ns)
    assert gc.chk_replay.__code__.co_filename == "<string>"
    # ---- the exit-aware end check
    gb.LAB["end_imports"] = ["Work.KLow.LabelsX", "Work.GCert.Data.Gen.%sGl" % NAME]
    gb.LAB["end"] = "GLab.endKX par %(v)d ports (raw.ret.map Prod.fst) gl el st0 st%(a)d st%(f)d"
    sys.argv = ["gxgen.py", str(a.certificate.resolve()), NAME, "gseg=2000", "batches=4", "frames=1"]
    gen.main()
    base = a.out / "Work/GCert/Data"
    cert = base / ("%sCert.lean" % NAME)
    s = cert.read_text()
    s = s.replace("import Work.GCert.Labels.Main", "import Work.KLow.LabelsX\nimport Work.GCert.Data.Gen.%sGl" % NAME)
    assert "theorem valid : GLab.Valid raw aux where" in s
    s = s.replace("theorem valid : GLab.Valid raw aux where", "theorem validX : GLab.ValidX raw aux gl el where")
    s = s.split("/-- price of the label half over ANY scalar half -/")[0]
    cert.write_text(s + "end GXD.%s\n#print axioms GXD.%s.validX\n" % (NAME, NAME))
    price = base / ("%sPrice.lean" % NAME)
    s = price.read_text()
    old = f'{c["R"]} * {c["h"]} + 2 * {c["v"]} * ({c["h"]} - 1) + {c["cst"]}'
    assert s.count(old) == 1
    s = s.replace(old, str(c["N"]))
    s = s.replace("R h + 2 v (h - 1) + cst =", "sum_q dim(final q) + 2 v (h - 1) + cst - sum_q rank(ent q) =")
    price.write_text(s)
    # ---- gauge list, exit list, E7 checks (Shea's), with an exact mirror of each E7 batch
    gl = [r for r, f in c["ext"]]
    el = sorted(c["exits"])
    batches = [gl[i:i + 80] for i in range(0, len(gl), 80)]
    lc = zm.lean_data(c)
    esw = 8
    while not all(zm.evaluate(zm.eCheck, lc, (d, esw, ux, us, uy), b, False)["value"] for b in batches):
        esw += 1
        assert esw < 64
    head = "namespace GXD.%s\nopen GXD\n" % NAME
    defs = [f"noncomputable def gl : List Nat := {gl}"]
    defs += [f"noncomputable def glb{i} : List Nat := {b}" for i, b in enumerate(batches)]
    defs += [f"/-- the exit registers (exit helpers of the hub word) -/\nnoncomputable def el : List Nat := {el}"]
    (base / ("Gen/%sGl.lean" % NAME)).write_text("import Work.GCert.Data.Gen.%s\n" % NAME + head + "\n".join(defs)
                                                 + "\nend GXD.%s\n" % NAME)
    checks = [f"theorem echk{i}_ok : GS.eCheck raw ⟨{d}, {esw}, {ux}, {us}, {uy}⟩ glb{i} = true := by decide +kernel"
              for i in range(len(batches))]
    (base / ("Gen/%sEChk.lean" % NAME)).write_text(
        "import Work.Gauge.ECheck\nimport Work.GCert.Data.Gen.%sGl\nset_option maxRecDepth 100000\n" % NAME + head
        + "\n".join(checks) + "\nend GXD.%s\n" % NAME
        + "\n".join(f"#print axioms GXD.{NAME}.echk{i}_ok" for i in range(len(batches))) + "\n")
    # ---- the zero checks (exit conditions) on the y-free twin
    zp = zd["spar"]
    P = "⟨%d, %d, %d, %d, %d⟩" % tuple(zp)
    v = c["v"]
    lel, hel = zd["light"], zd["heavy"]
    assert sorted(lel + hel) == el and not set(lel) & set(hel)
    items = []   # (theorem name, statement)
    for b, val, k in zd["bits"]:
        pre = list(range(v, v + k))
        items.append(("zb%d%s" % (b, "t" if val else "f"),
                      f"GS.zCheck (GXD.nyRaw raw) {P} ({pre} ++ GXD.bitB lel {b} {'true' if val else 'false'}) (GXD.slotsL (GXD.nyRaw raw)) = true"))
    for j, (e, k) in enumerate(zd["heavy_batches"]):
        pre = list(range(v, v + k))
        items.append(("zh%d" % j, f"GS.zCheck (GXD.nyRaw raw) {P} ({pre} ++ [{e}]) (GXD.slotsL (GXD.nyRaw raw)) = true"))
    for j, ch in enumerate(zd["chunks"]):
        items.append(("zc%d" % j, f"GS.zCheck (GXD.nyRaw raw) {P} chunk{j} el = true"))
    per = 8
    mods = []
    zdefs = [f"noncomputable def lel : List Nat := {lel}", f"noncomputable def hel : List Nat := {hel}"]
    zdefs += [f"noncomputable def chunk{j} : List Nat := {ch}" for j, ch in enumerate(zd["chunks"])]
    zdefs += ["noncomputable def chunks : List (List Nat) := [%s]" % ", ".join("chunk%d" % j for j in range(len(zd["chunks"])))]
    (base / ("Gen/%sZL.lean" % NAME)).write_text(
        "import Work.GCert.Data.Gen.%sGl\n" % NAME + head + "\n".join(zdefs) + "\nend GXD.%s\n" % NAME)
    for i in range(0, len(items), per):
        m = "%sZChk%d" % (NAME, i // per)
        mods.append(m)
        body = "\n".join(f"theorem {nm} : {st} := by decide +kernel" for nm, st in items[i:i + per])
        (base / ("Gen/%s.lean" % m)).write_text(
            "import Work.KLow.NoY\nimport Work.KLow.GlueX\nimport Work.GCert.Data.Gen.%sZL\nset_option maxRecDepth 100000\n" % NAME
            + head + body + "\nend GXD.%s\n" % NAME
            + "\n".join(f"#print axioms GXD.{NAME}.{nm}" for nm, _ in items[i:i + per]) + "\n")
    nb = zd["nb"]
    bitcases = []
    for b in range(nb):
        for val in ("f", "t"):
            bitcases.append(f"⟨_, _, by rw [ny_v]; decide, zb{b}{val}⟩")
    nbits = len(bitcases)
    ors = " ∨ ".join(f"b = {b}" for b in range(nb))
    rfls = " | ".join(["rfl"] * nb)
    nh = len(zd["heavy_batches"])
    hcases = ", ".join(f"⟨_, _, by rw [ny_v]; decide, zh{j}⟩" for j in range(nh))
    nc = len(zd["chunks"])
    ccases = ", ".join(f"⟨_, zc{j}⟩" for j in range(nc))
    N = 2 * v + c["R"]
    imps = ''.join('import Work.GCert.Data.Gen.%s\n' % m for m in mods)
    zok = f"""{imps}import Work.GCert.Data.{NAME}ScalOK
set_option maxRecDepth 100000
/-! # (key: k-lower) the coverage facts of the zero checks of `GXD.{NAME}` and the two exit conditions
Generated by k-lower's generate_local_lean_x.py.  Bit batches over the light exit list `lel` ({len(lel)} exits, {nb} bits),
one batch per heavy exit ({nh}), {nc} chunks of the other non-y registers; all on the y-free twin `GXD.nyRaw raw`. -/
namespace GXD.{NAME}
open GXD OAI.PowerSaving

theorem ny_v : (GXD.nyRaw raw).v = {v} := rfl

theorem hbit : ∀ b, b < {nb} → ∀ val : Bool, ∃ p pre, (∀ x ∈ pre, x < 2 * (GXD.nyRaw raw).v) ∧
    GS.zCheck (GXD.nyRaw raw) p (pre ++ GXD.bitB lel b val) (GXD.slotsL (GXD.nyRaw raw)) = true := by
  intro b hb val
  have h : {ors} := by omega
  rcases h with {rfls} <;> cases val
  exacts [{', '.join(bitcases)}]

theorem hhv : ∀ e ∈ hel, ∃ p pre, (∀ x ∈ pre, x < 2 * (GXD.nyRaw raw).v) ∧
    GS.zCheck (GXD.nyRaw raw) p (pre ++ [e]) (GXD.slotsL (GXD.nyRaw raw)) = true := by
  intro e he
  simp only [hel, List.mem_cons, List.not_mem_nil, or_false] at he
  rcases he with {' | '.join(['rfl'] * nh)}
  exacts [{hcases}]

theorem hch : ∀ l ∈ chunks, ∃ p, GS.zCheck (GXD.nyRaw raw) p l el = true := by
  intro l hl
  simp only [chunks, List.mem_cons, List.not_mem_nil, or_false] at hl
  rcases hl with {' | '.join(['rfl'] * nc)}
  exacts [{ccases}]

theorem hcov' : ∀ s, s < {N} → (s < {v} ∨ {2 * v} ≤ s) → GLab.memK s el = false → ∃ l ∈ chunks, s ∈ l := by
  decide +kernel
theorem hdisj : ∀ l ∈ chunks, ∀ s ∈ l, GLab.memK s el = false := by decide +kernel
theorem hlen : lel.length ≤ 2 ^ {nb} := by decide +kernel
theorem hel2' : ∀ r ∈ el, {2 * v} ≤ r := by decide +kernel
theorem hlel : ∀ r ∈ lel, GLab.memK r el = true := by decide +kernel
theorem hhel : ∀ r ∈ hel, r ∉ lel := by decide +kernel
theorem hsplit : ∀ r ∈ el, r ∈ lel ∨ r ∈ hel := by decide +kernel

/-- **the two exit conditions** of the certificate (on its own total matrix) -/
theorem exitRC : OAI.PowerSaving.XEX.ExitRows (GLab.exitSet raw el) (GS.Mt raw scalOK.hv) ∧
    OAI.PowerSaving.XEX.ExitCols (GLab.exitSet raw el) (GS.Mt raw scalOK.hv) :=
  GXD.exit_of_twin raw scalOK.hv (fun x hx => (scalOK.addB x hx).2.2.2.2) scalOK.hy el
    (GXD.exit_of_zchecks (GXD.nyRaw raw) scalOK.hv el
      (GXD.twin_hy raw scalOK.hv (fun x hx => (scalOK.addB x hx).2.2.2.2) scalOK.hy)
      lel hel {nb} chunks hlen (by decide) (fun r hr => hel2' r hr) hlel hhel hsplit hbit hhv hch
      (fun s hs h1 h2 => hcov' s hs h1 h2) hdisj)

end GXD.{NAME}
#print axioms GXD.{NAME}.exitRC
"""
    (base / ("Gen/%sZOK.lean" % NAME)).write_text(zok)
    result["zero_checks"] = dict(spar=zp, modules=mods, items=[nm for nm, _ in items])
    result["e7_sw"] = esw
    (a.out / ("logs/preflight-%s.json" % NAME)).write_text(json.dumps(result, indent=2, default=str) + "\n")
    print("ok", NAME, "E7 sw", esw, "zero checks", len(items), "modules", mods)


if __name__ == "__main__":
    main()
