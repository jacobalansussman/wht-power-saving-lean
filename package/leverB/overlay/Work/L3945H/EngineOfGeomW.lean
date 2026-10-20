import Work.L3945H.Stmt
import Work.Mixed.Bridge.HFTest3945HW

/-!
# bB-data: the adapter from the design-B (#3945H) warp network to the seam's engine statement `MX.L3945H.EngineStmt`

`MX.L3945H.engine_of_geomW`: the geometric warp network of p-a3 (`xnetGeomW`, through the warp heads-first
certificate and p-bridge's bridge with the #3945 instance) gives OpenAI's `hills_program` statement at the seam
exponent `1 - 1758152001/10^12` (DFT figure 1.7581520e-3), i.e. `MX.L3945H.EngineStmt.{0}` of i-skel's generated
`Work/L3945H/Stmt.lean` (suffix L3945H).  The start frames `base` are chosen as in p-a3's `bridgeHostTW_gcert` (`base_ex`), the label
`l` as in p-a3's `engine_of_ginvW_hostT` (any element of `H`).  (i-final's L3945/EngineOfGeomW.lean with the design-B
numbers.)  bB-data (wht45/allheads, 19 Oct 2026, design B: #3945 with every child-8/10 bin a literal head); prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI
namespace PowerSaving
namespace MX
namespace L3945H
open Binary Matrix Finset RAM RAM.Ty SS CB RF BR GF GX BC B401 XHost XHostR XHostD XHostT XHostTW

/-- **The engine statement of L3945H from the design-B warp network and its data.** -/
theorem engine_of_geomW {H Sl C : Type} [Fintype H] [DecidableEq H] [Fintype Sl] [DecidableEq Sl]
    [Fintype C] [DecidableEq C] {Ok : ∀ {β : Type} [Fintype β] [DecidableEq β], GWin H β → Prop}
    (I : GInvW H (Fin 112) Sl C Ok) (Ic : Fin 360 → Fin 9 → GInvW H (Fin 112) Sl C Ok)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (hH : Fintype.card H = 9)
    (hOk : ∀ (l : H) (n : ℕ) (d : Orth (BC.LC (Fin 8) H l)), Ok (win401N l hH n d))
    (S : HostSchedT 360 Sl (Fin 112)) (hfit : S.fits (resIcW Ic))
    (hpos : ∀ ρ l' q, 0 < resIcW Ic ρ l' q)
    (hheads : (Finset.univ.filter (fun p : Fin 360 × Bool × Fin 112 => S.hrem p.1 4 p.2 = 48)).card = 80640)
    (hU : ∀ φ, unitPriceTW Ic S φ = listCost φ PB.unit3945H) :
    EngineStmt.{0} := by
  intro cl m ι k v
  obtain ⟨l⟩ : Nonempty H := Fintype.card_pos_iff.mp (by omega)
  have hb : ∀ t, ∃ A : OBase H, I.tv t = A.v l := fun t => base_ex (I.tv t) (I.hunit t) (I.hgap t) l
  choose base hbase using hb
  exact PB.HFTest3945HW.engine3945H_geomW_unit I Ic hIc base l hbase hH (hOk l) S hfit hpos hheads hU m cl k v

end L3945H
end MX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.MX.L3945H.engine_of_geomW
