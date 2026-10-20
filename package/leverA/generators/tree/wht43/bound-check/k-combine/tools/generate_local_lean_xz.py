#!/usr/bin/env python3
"""k-combine (13 Oct 2026): the LOCAL Lean modules of an EXIT-FRAME certificate with MUTABLE SOURCES (exit helpers,
hub tail in phase B, source copies folded back).  k-lower's generate_local_lean_x.py with m-gen's mutable-source scalar
side (generate_local_lean_z.py: kinds okAZK/okBZK in the mirror, scalCheckZ/yCheckZ/eCheckZ, scalOKZ), SPar from the
exact rational audit scalar_audit_z.py and the exact packed mirror zmirror.suite (asserted), and the zero checks
`GS.zCheckZ` (Work/KMX/ZCheckZ) with rows x + slots for the exit batches (ExitColsX), exitRC from
`GXD.exit_of_twinZ` / `GXD.exit_of_zchecksZ` (Work/KMX/GlueXZ).  k-lower's text follows:


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
    sys.path.insert(0, str(HERE))
    from scalar_audit_z import audit
    import zmirror as zm
    result = audit(c)
    SPAR = tuple(result["SPar_by_Shea_rule"])
    lc = zm.lean_data(c)
    lc["ext"] = c["ext"]
    gl0 = [r for r, f in c["ext"]]
    eb0 = [gl0[i:i + 80] for i in range(0, len(gl0), 80)]
    esw = 8
    while not all(zm.evaluate(zm.eCheck, lc, (SPAR[0], esw) + SPAR[2:], b, True)["value"] for b in eb0):
        esw += 1
        assert esw < 64
    EPAR = (SPAR[0], esw) + SPAR[2:]
    exact = zm.suite(lc, SPAR, 4, EPAR, 80, True)
    assert exact["all_true"], exact
    result["SPar"] = list(SPAR)
    os.environ["GX_OUT"] = str(a.out.resolve())
    sys.path.insert(0, str(a.toolkit.resolve()))
    sys.path.insert(0, str(HERE))
    import gxcore as gc
    import gxgen as gen
    import gxbind as gb
    gc.SC = 8
    def chk_shape_z(D):
        """E3 with the mutable-source kinds (okAZK / okBZK; scatter okSK is checked in scatter_cut / replay)"""
        n = 0
        for ph in "AB":
            for i, g in enumerate(D[ph]):
                for t, s, _ in gc.adds(g):
                    ct, cs = gc.cls(D, t), gc.cls(D, s)
                    gc.need(t != s, "E3: t = s (gate %s%d)" % (ph, i))
                    if ph == "A":
                        gc.need(ct != "y" and cs != "y", "E3z: phase A names a target (gate A%d)" % i)
                    else:
                        gc.need(cs != "y" or ct == "y", "E3z: a target is read by a non-target (gate B%d)" % i)
                    n += 1
        return n

    def chk_scalar_z(D, batches=1):
        """gxcore.chk_scalar without the X3 assertion (x rows need not come back: mutable sources)"""
        v = D["v"]
        from collections import Counter
        count = Counter()
        reg, afterA = gc.replay(D, {t: {t: gc.SC} for t in range(v)}, count)
        gc.need(all(reg[v + t] == {t: gc.SC} for t in range(v)), "E5: a target does not receive exactly its source (hid)")
        if D["ext"]:
            rg, _ = gc.replay(D, {r: {r: gc.SC} for r, _ in D["ext"]})
            gc.need(all(not rg[v + t] for t in range(v)), "E7: the dirt of an exterior-gauged slot reaches a y role")
        step = -(-v // batches)
        per = [sum(count[k] for k in range(lo, min(lo + step, v))) for lo in range(0, v, step)]
        return dict(updates=sum(count.values()), per_batch=per, afterA=afterA,
                    batches=[(lo, min(lo + step, v)) for lo in range(0, v, step)],
                    x_rows_restored=all(reg[t] == {t: gc.SC} for t in range(v)))

    gc.chk_shape = chk_shape_z
    gc.chk_scalar = chk_scalar_z
    d, sw, ux, us, uy = SPAR
    gb.SCAL = dict(imports=["Work.XMut.ScalarZ"],
                   scal="GS.scalCheckZ raw ⟨%(d)d, %(sw)d, %(ux)d, %(us)d, %(uy)d⟩ %(lo)d %(n)d",
                   ychk="GS.yCheckZ raw ⟨%(d)d, %(sw)d, %(ux)d, %(us)d, %(uy)d⟩ 0 %(v)d")
    gb.spar = lambda D, path: dict(d=d, sw=sw, ux=ux, us=us, uy=uy,
                                   source="exact rational audit (scalar_audit_z) and exact packed mirror (zmirror)",
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
    # ---- the scalar half, sources mutable (m-gen)
    sc = base / ("%sScalOK.lean" % NAME)
    s = sc.read_text()
    for old, new in [("import Work.GCert.Scalar.Main", "import Work.XMut.ScalarZ"),
                     ("THE SCALAR HALF `scalOK : GS.ScalOK raw`",
                      "THE SCALAR HALF, SOURCES MUTABLE, `scalOKZ : GS.ScalOKZ raw` (x-mutable, Work/XMut/ScalarZ.lean)"),
                     ("theorem scalOK : GS.ScalOK raw where", "theorem scalOKZ : GS.ScalOKZ raw where"),
                     ("GS.XRun.of_check", "GS.XRunZ.of_check"), ("GS.YRun.of_check", "GS.YRunZ.of_check"),
                     ("#print axioms GXD.%s.scalOK\n" % NAME, "#print axioms GXD.%s.scalOKZ\n" % NAME)]:
        assert old in s, old
        s = s.replace(old, new)
    sc.write_text(s)
    # ---- gauge list, exit list, E7 checks (Shea's), with an exact mirror of each E7 batch
    gl = [r for r, f in c["ext"]]
    el = sorted(c["exits"])
    batches = [gl[i:i + 80] for i in range(0, len(gl), 80)]
    assert batches == eb0 and EPAR[1] == esw
    head = "namespace GXD.%s\nopen GXD\n" % NAME
    defs = [f"noncomputable def gl : List Nat := {gl}"]
    defs += [f"noncomputable def glb{i} : List Nat := {b}" for i, b in enumerate(batches)]
    defs += [f"/-- the exit registers (exit helpers of the hub word) -/\nnoncomputable def el : List Nat := {el}"]
    (base / ("Gen/%sGl.lean" % NAME)).write_text("import Work.GCert.Data.Gen.%s\n" % NAME + head + "\n".join(defs)
                                                 + "\nend GXD.%s\n" % NAME)
    checks = [f"theorem echk{i}_ok : GS.eCheckZ raw ⟨{d}, {esw}, {ux}, {us}, {uy}⟩ glb{i} = true := by decide +kernel"
              for i in range(len(batches))]
    (base / ("Gen/%sEChk.lean" % NAME)).write_text(
        "import Work.XMut.ECheckZ\nimport Work.GCert.Data.Gen.%sGl\nset_option maxRecDepth 100000\n" % NAME + head
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
                      f"GS.zCheckZ (GXD.nyRaw raw) {P} ({pre} ++ GXD.bitB lel {b} {'true' if val else 'false'}) (GXD.rowsXS (GXD.nyRaw raw)) = true"))
    for j, (e, k) in enumerate(zd["heavy_batches"]):
        pre = list(range(v, v + k))
        items.append(("zh%d" % j, f"GS.zCheckZ (GXD.nyRaw raw) {P} ({pre} ++ [{e}]) (GXD.rowsXS (GXD.nyRaw raw)) = true"))
    for j, ch in enumerate(zd["chunks"]):
        items.append(("zc%d" % j, f"GS.zCheckZ (GXD.nyRaw raw) {P} chunk{j} el = true"))
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
            "import Work.KLow.NoY\nimport Work.KMX.GlueXZ\nimport Work.GCert.Data.Gen.%sZL\nset_option maxRecDepth 100000\n" % NAME
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
/-! # (key: k-combine) the coverage facts of the zero checks of `GXD.{NAME}` and the three exit conditions
Generated by k-combine's generate_local_lean_xz.py (on k-lower's generate_local_lean_x.py; sources mutable).  Bit batches over the light exit list `lel` ({len(lel)} exits, {nb} bits),
one batch per heavy exit ({nh}), {nc} chunks of the other non-y registers; all on the y-free twin `GXD.nyRaw raw`. -/
namespace GXD.{NAME}
open GXD OAI.PowerSaving

theorem ny_v : (GXD.nyRaw raw).v = {v} := rfl

theorem hbit : ∀ b, b < {nb} → ∀ val : Bool, ∃ p pre, (∀ x ∈ pre, (GXD.nyRaw raw).v ≤ x ∧ x < 2 * (GXD.nyRaw raw).v) ∧
    GS.zCheckZ (GXD.nyRaw raw) p (pre ++ GXD.bitB lel b val) (GXD.rowsXS (GXD.nyRaw raw)) = true := by
  intro b hb val
  have h : {ors} := by omega
  rcases h with {rfls} <;> cases val
  exacts [{', '.join(bitcases)}]

theorem hhv : ∀ e ∈ hel, ∃ p pre, (∀ x ∈ pre, (GXD.nyRaw raw).v ≤ x ∧ x < 2 * (GXD.nyRaw raw).v) ∧
    GS.zCheckZ (GXD.nyRaw raw) p (pre ++ [e]) (GXD.rowsXS (GXD.nyRaw raw)) = true := by
  intro e he
  simp only [hel, List.mem_cons, List.not_mem_nil, or_false] at he
  rcases he with {' | '.join(['rfl'] * nh)}
  exacts [{hcases}]

theorem hch : ∀ l ∈ chunks, ∃ p, GS.zCheckZ (GXD.nyRaw raw) p l el = true := by
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

/-- **the three exit conditions** of the certificate (on its own total matrix; sources mutable) -/
theorem exitRC : OAI.PowerSaving.XEX.ExitRows (GLab.exitSet raw el) (GS.Mt raw scalOKZ.hv) ∧
    OAI.PowerSaving.XEX.ExitCols (GLab.exitSet raw el) (GS.Mt raw scalOKZ.hv) ∧
    OAI.PowerSaving.XEX.ExitColsX (GLab.exitSet raw el) (GS.Mt raw scalOKZ.hv) :=
  GXD.exit_of_twinZ raw scalOKZ.hv el (fun x hx => (scalOKZ.addB x hx).2.2.2)
    (GXD.exit_of_zchecksZ (GXD.nyRaw raw) scalOKZ.hv el
      (GXD.twin_hy raw scalOKZ.hv (fun x hx => (scalOKZ.addB x hx).2.2.2) scalOKZ.hy)
      lel hel {nb} chunks hlen (by decide) (fun r hr => hel2' r hr) hlel hhel hsplit hbit hhv hch
      (fun s hs h1 h2 => hcov' s hs h1 h2) hdisj)

end GXD.{NAME}
#print axioms GXD.{NAME}.exitRC
"""
    (base / ("Gen/%sZOK.lean" % NAME)).write_text(zok)
    result["zero_checks"] = dict(spar=zp, modules=mods, items=[nm for nm, _ in items])
    result["e7_sw"] = esw
    result["epar"] = list(EPAR)
    result["exact_packed_mirror"] = exact
    (a.out / ("logs/preflight-%s.json" % NAME)).write_text(json.dumps(result, indent=2, default=str) + "\n")
    print("ok", NAME, "E7 sw", esw, "zero checks", len(items), "modules", mods)


if __name__ == "__main__":
    main()
