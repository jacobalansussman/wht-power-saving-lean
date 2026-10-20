#!/usr/bin/env python3
"""m-gen (13 Oct 2026): generate the LOCAL Lean modules of the copy-free #588 certificate `C588M` (mutable sources).

Ryan Shea's research/complex-transfer-588-596/experiments/generate_local_lean.py (shea256/fourier-transform-below-nlogn
@ 5239484), which drives Jacob Sussman's gxgen.py / gxcore.py / gxbind.py UNCHANGED (copies in ../toolkit/tools/gx),
adapted to the mutable-source checks of x-mutable (XMut/ScalarZ.lean, XMut/ECheckZ.lean, XMut/Glue.lean):

* Python mirror (gxcore.mirror): the label side (E0 frames, E1 pairs and replay, E2 scatter, E6 histogram) runs
  unchanged; the scalar side's E3 kind test is the mutable-source kinds okAZK / okBZK (`chk_shape_z`), and E5
  asserts hid only, not the x restoration X3 (`chk_scalar_z`); E7 unchanged.
* SPar: (10, 21, 4, 8, 8) for the four source batches and the y check, (10, 8, 4, 8, 8) for the three E7 batches;
  both from the exact rational audit (scalar_audit_z.py, Shea's rule) and confirmed by the exact packed mirror
  of the kernel terms (zmirror.py); sw - 1 fails (zmirror --margin).
* Kernel checks named: GS.scalCheckZ / GS.yCheckZ (Scal<i>, YChk), GS.eCheckZ (EChk); the scalar half is
  `scalOKZ : GS.ScalOKZ raw` (GS.XRunZ.of_check / GS.YRunZ.of_check).
* Everything else as Shea: entrance-gauged ends (GLab.endKE with the gauge list gl), validE, the invoice patch
  (movement invoice N = R h + 2 v (h-1) + cst - sum of entrance ranks), the gauge list module.
* The XMut modules live at Work/XMut/*.lean in the overlay (module names Work.XMut.*), so the tree's
  `lean_lib Work` (globs Work.+) builds them without a lakefile change.

usage: generate_local_lean_z.py <cert.gcert.json.gz> --toolkit <toolkit/tools/gx> --out <dir> [--name C588M]
"""
import argparse
import gzip
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
SPAR = (10, 21, 4, 8, 8)
EPAR = (10, 8, 4, 8, 8)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificate", type=Path)
    ap.add_argument("--toolkit", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--name", default="C588M")
    a = ap.parse_args()
    NAME = a.name
    c = json.loads(gzip.decompress(a.certificate.read_bytes()))
    sys.path.insert(0, str(HERE))
    from scalar_audit_z import audit
    import zmirror as zm
    result = audit(c)
    assert tuple(result["SPar_by_Shea_rule"]) == SPAR, result["SPar_by_Shea_rule"]
    lc = zm.lean_data(c)
    lc["ext"] = c["ext"]
    exact = zm.suite(lc, SPAR, 4, EPAR, 80, True)
    assert exact["all_true"], exact
    os.environ["GX_OUT"] = str(a.out.resolve())
    sys.path.insert(0, str(a.toolkit.resolve()))
    import gxcore as gc
    import gxgen as gen
    import gxbind as gb

    gc.SC = 8                                     # Shea: entries in eighths

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
    gb.LAB["end_imports"] = ["Work.Gauge.Labels", "Work.GCert.Data.Gen.%sGl" % NAME]
    gb.LAB["end"] = "GLab.endKE par %(v)d ports (raw.ret.map Prod.fst) gl st0 st%(a)d st%(f)d"
    sys.argv = ["gxgen.py", str(a.certificate.resolve()), NAME, "gseg=2000", "batches=4", "frames=1"]
    gen.main()
    base = a.out / "Work/GCert/Data"
    # ---- the label half: validE over the gauge list (as Shea)
    cert = base / ("%sCert.lean" % NAME)
    s = cert.read_text()
    s = s.replace("import Work.GCert.Labels.Main", "import Work.Gauge.Labels\nimport Work.GCert.Data.Gen.%sGl" % NAME)
    assert "theorem valid : GLab.Valid raw aux where" in s
    s = s.replace("theorem valid : GLab.Valid raw aux where", "theorem validE : GLab.ValidE raw aux gl where")
    s = s.split("/-- price of the label half over ANY scalar half -/")[0]
    cert.write_text(s + "end GXD.%s\n#print axioms GXD.%s.validE\n" % (NAME, NAME))
    # ---- the invoice (as Shea): entrance gauges omit their initial ranks from the movement invoice
    price = base / ("%sPrice.lean" % NAME)
    s = price.read_text()
    old = f'{c["R"]} * {c["h"]} + 2 * {c["v"]} * ({c["h"]} - 1) + {c["cst"]}'
    assert s.count(old) == 1
    s = s.replace(old, str(c["N"]))
    s = s.replace("R h + 2 v (h - 1) + cst =", "R h + 2 v (h - 1) + cst - sum_q rank(ent q) =")
    price.write_text(s)
    # ---- the scalar half, sources mutable
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
    # ---- the gauge list and the E7 checks, sources mutable
    gl = [r for r, f in c["ext"]]
    batches = [gl[i:i + 80] for i in range(0, len(gl), 80)]
    head = "namespace GXD.%s\nopen GXD\n" % NAME
    defs = [f"noncomputable def gl : List Nat := {gl}"]
    defs += [f"noncomputable def glb{i} : List Nat := {b}" for i, b in enumerate(batches)]
    (base / ("Gen/%sGl.lean" % NAME)).write_text("import Work.GCert.Data.Gen.%s\n" % NAME + head + "\n".join(defs) +
                                                 "\nend GXD.%s\n" % NAME)
    ed, esw, eux, eus, euy = EPAR
    checks = [f"theorem echk{i}_ok : GS.eCheckZ raw ⟨{ed}, {esw}, {eux}, {eus}, {euy}⟩ glb{i} = true := by decide +kernel"
              for i in range(len(batches))]
    (base / ("Gen/%sEChk.lean" % NAME)).write_text(
        "import Work.XMut.ECheckZ\nimport Work.GCert.Data.Gen.%sGl\nset_option maxRecDepth 100000\n" % NAME + head +
        "\n".join(checks) + "\nend GXD.%s\n" % NAME +
        "\n".join(f"#print axioms GXD.{NAME}.echk{i}_ok" for i in range(len(batches))) + "\n")
    result["exact_packed_mirror"] = exact
    result["epar"] = list(EPAR)
    result["e_batches"] = [len(b) for b in batches]
    (a.out / ("logs/preflight-%s.json" % NAME)).write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
