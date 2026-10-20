import Work.L3945F.EngineOfGeomW
import Work.L3945H.Prov

/-!
# fc-engine: THE engine theorem of the design-B FULL-cut lever instance L3945F (#3945 with the all-heads fill)

`MX.L3945F.engine : MX.L3945F.EngineStmt.{0}`: OpenAI's `hills_program` statement at the seam exponent
`1 - 1780693501/10^12` (DFT figure 1.7806935e-3), from the #3945 warp network with the all-heads fill of round
wht45/reopt (design B) (hub PR #3945, Dugongue's staged544 T360 terminal-gap fill, and its lineage; our port: p-lower,
p-wgen, p-warpchk, p-warpsem, p-a3, p-hostck, p-disp, p-hfirst, x-inst's warp heads-first, bB-data's design-B data)
with the FULL-cut lever (each of the 80,640 literal rank-48 heads removed entirely, target 25, and the other 57,312
idle rank-48 blocks split 24 + 24: every paid rank at most 24 < 25; e-r45m1's `cut_heads_full` and `splitWord`,
fc-engine's bridge, fc-rate's R5) and the mixed-target recursion (Theorem 1G: e-r12, e-r3prog, e-r3cost).
The providers, dispatch, window facts, residuals and unit price are design B's, unchanged (bB-data's
`Work.L3945H.Prov`: they do not depend on the cut).  This is the one theorem the generated chain
(`Work.FourierL3945F.Seam`) imports.
fc-engine (wht45/fullcut, 19-20 Oct 2026) on bB-data's Engine.lean; prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI
namespace PowerSaving
namespace MX
namespace L3945F

/-- **THE ENGINE THEOREM OF L3945F.** -/
theorem engine : EngineStmt.{0} :=
  engine_of_geomW L3945H.I0W L3945H.Ic3945W L3945H.hIc3945W (Fintype.card_fin 9) L3945H.hOk3945
    XHostT.Inst.sched3945H L3945H.hfit3945W L3945H.hpos3945W XHostT.Inst.sched3945H_heads L3945H.hU3945W

end L3945F
end MX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.MX.L3945F.engine
