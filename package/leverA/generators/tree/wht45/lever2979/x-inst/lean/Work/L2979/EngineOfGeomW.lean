import Work.L2979.Stmt
import Work.Mixed.Bridge.HFTest2979W

/-!
# x-inst: the adapter from the #2979 warp network to the seam's engine statement `MX.L2979.EngineStmt`

`MX.L2979.engine_of_geomW`: the geometric warp network of p-a3 (`xnetGeomW`, through x-inst's warp heads-first
certificate and p-bridge's bridge with the #2979 instance) gives OpenAI's `hills_program` statement at the seam
exponent `1 - 1755647301/10^12` (DFT figure 1.7556473e-3), i.e. `MX.L2979.EngineStmt.{0}` of the generated
`Work/L2979/Stmt.lean` (x-inst's tools/mk_l2979.py).  The start frames `base` are chosen as in p-a3's
`bridgeHostTW_gcert` (`base_ex`), the label `l` as in p-a3's `engine_of_ginvW_hostT` (any element of `H`).
What is left are the data of the #2979 network: the providers (x-words), the schedule `S` with `hfit`, positive
residuals, the head count and the unit price (x-tables).  x-final's `Work/L2979/Engine.lean` is ONE line:
`theorem engine : EngineStmt.{0} := engine_of_geomW I Ic hIc hH hOk S hfit hpos hheads hU` with those terms.
Written by x-inst (wht45/lever2979, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI
namespace PowerSaving
namespace MX
namespace L2979
open Binary Matrix Finset RAM RAM.Ty SS CB RF BR GF GX BC B401 XHost XHostR XHostD XHostT XHostTW

/-- **The engine statement of L2979 from the #2979 warp network and its data.** -/
theorem engine_of_geomW {H Sl C : Type} [Fintype H] [DecidableEq H] [Fintype Sl] [DecidableEq Sl]
    [Fintype C] [DecidableEq C] {Ok : ∀ {β : Type} [Fintype β] [DecidableEq β], GWin H β → Prop}
    (I : GInvW H (Fin 112) Sl C Ok) (Ic : Fin 74 → Fin 9 → GInvW H (Fin 112) Sl C Ok)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (hH : Fintype.card H = 9)
    (hOk : ∀ (l : H) (n : ℕ) (d : Orth (BC.LC (Fin 8) H l)), Ok (win401N l hH n d))
    (S : HostSchedT 74 Sl (Fin 112)) (hfit : S.fits (resIcW Ic))
    (hpos : ∀ ρ l' q, 0 < resIcW Ic ρ l' q)
    (hheads : (Finset.univ.filter (fun p : Fin 74 × Bool × Fin 112 => S.hrem p.1 4 p.2 = 48)).card = 16086)
    (hU : ∀ φ, unitPriceTW Ic S φ = listCost φ PB.unit2979) :
    EngineStmt.{0} := by
  intro cl m ι k v
  obtain ⟨l⟩ : Nonempty H := Fintype.card_pos_iff.mp (by omega)
  have hb : ∀ t, ∃ A : OBase H, I.tv t = A.v l := fun t => base_ex (I.tv t) (I.hunit t) (I.hgap t) l
  choose base hbase using hb
  exact PB.HFTest2979W.engine2979_geomW_unit I Ic hIc base l hbase hH (hOk l) S hfit hpos hheads hU m cl k v

end L2979
end MX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.MX.L2979.engine_of_geomW
