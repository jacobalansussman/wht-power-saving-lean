#!/usr/bin/env python3
"""x-mutable: regenerate proj/XMut/Labels.lean from the vendored proj/Work/Gauge/Labels.lean (Nil Mamano's ValidE.glab,
via Ryan Shea's complex588 tree) by text substitution: GScal -> GScalZ, GLab -> GLabZ, glab -> glabZ.
usage: make_labelsZ.py <proj/Work/Gauge/Labels.lean> <proj/XMut/Labels.lean>"""
import sys
HDR = "import Work.Gauge.Labels\nimport XMut.GInv\n\n/-!\n# x-mutable (13 Oct 2026): the label half of a gauged certificate over a MUTABLE-SOURCE scalar half\n\n`ValidE.glabZ` / `ValidE.glab_ranksZ` = Nil Mamano's `ValidE.glab` / `ValidE.glab_ranks` (`Work/Gauge/Labels.lean`,\non Jacob Sussman's label check), word for word, with the scalar half `GScalZ` (no `hx`) in place of `GScal`.\nThe label check never used `hx`: it reads the scalar half only through `MA`, `MAi`, `MB`, `MBi`, `Cc` and\n`gauged`.  Generated from the vendored source by a text substitution (x-mutable/tools/make_labelsZ.sh).  No `sorry`.\n-/\n\nset_option linter.unusedVariables false\nset_option linter.unusedSectionVars false\nset_option linter.deprecated false\n\nnamespace GLab\nopen SSC GXD OAI.PowerSaving OAI.PowerSaving.Binary OAI.PowerSaving.GF OAI.PowerSaving.RAM\n  OAI.PowerSaving.CB Finset Matrix\n\n"
src = open(sys.argv[1]).read()
body = src[src.index("/-- **THE LABEL HALF of a gauged certificate**"):src.index("end GLab")]
for a, b in [("noncomputable def ValidE.glab ", "noncomputable def ValidE.glabZ "),
             ("theorem ValidE.glab_ranks ", "theorem ValidE.glab_ranksZ "),
             ("OAI.PowerSaving.GX.GScal (Fin c.v) (Fin c.R) C", "OAI.PowerSaving.GX.GScalZ (Fin c.v) (Fin c.R) C"),
             (": OAI.PowerSaving.GX.GLab (Fin a.p.h) S :=", ": OAI.PowerSaving.GX.GLabZ (Fin a.p.h) S :="),
             ("(V.glab S Ro", "(V.glabZ S Ro"),
             ("THE LABEL HALF of a gauged certificate**: `GX.GLab`",
              "THE LABEL HALF of a gauged certificate, sources mutable**: `GX.GLabZ`")]:
    body = body.replace(a, b)
open(sys.argv[2], "w").write(HDR + body + "end GLab\n\n#print axioms GLab.ValidE.glabZ\n#print axioms GLab.ValidE.glab_ranksZ\n")
