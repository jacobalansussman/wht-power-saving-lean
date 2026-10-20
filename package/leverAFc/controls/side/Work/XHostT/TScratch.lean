import Work.XHostT.TDisp
import XExit.Pack

/-!
# (key: bridge-host-cross) PRIVATE SCRATCH of a call: the zero buffers of a borrow word

#3945's provider words (Dugongue's terminal-gap fill on the staged544 words; borrow and warp mechanism of PR #2226,
Cinna Davis) allocate 91 ZERO BUFFERS per call for the full-frame cleanup: allocated at zero, read, freed at zero
(the hub's `bw_emit.check_local`).  In our packages they are copies of the invocation: roles of the box
`Box T Sl C` outside the x, y and slot roles, which a route may start and end at ANY frame (erase gates) and whose
content every call erases (`shearF`, `shearB` are zero on the copies).  p-warpchk's `KLW.WValid.inv` puts them there
(`GInvX (Fin 9) (Fin 112) (Fin 408) (Fin nc)`, scratch registers `≥ 2v + R`).

This file gives every call such private scratch when a package lacks it, or gives two providers a common copies
type, WITHOUT changing anything the network or the host check reads:

* `fwd_reC`, `bwd_reC`: a route of an invocation on `Box T Sl C` is a route on `Box T Sl C'` for any embedding
  `e : C ↪ C'` of the copies; the copies outside the image of `e` (the private scratch) go from their start frame
  to their end frame by a free erase gate, untouched by the old route (`GRoute.lift`);
* `GInvE.reC`, `GInvX.reC` (`I.reC e`): the package with the copies re-indexed; SAME target vectors, entrance (and
  exit) labels and price (`rfl`); `GInvE.addScratch I Z := I.reC Function.Embedding.inl` adds `Z` private roles;
* `XNet.reC` (generated part, `TDispW.lean`): every call of the network with the private scratch; the host table,
  the banks, the host chains and the price are the SAME fields (`XNet.unit_reC : (N.reC e).unit = N.unit`, `rfl`);
* **`xnet_disp_scratch_certificate`**: two providers dispatched by a table, private scratch per call, price
  `|Γ| * (n0 * cost0 + n1 * cost1 + restU)`.

In the network the copies of a call are the roles `Sum.inr (g, c)` (`BRoleT`): not live, never a host (hosts are
data roles of a replica, `slotT_ne_copy`), never a bank, never paid (no term of the price mentions them); every call
finds them erased and leaves them erased, so the calls may share them in sequence.

Written by p-disp (wht45/lever, 19 Oct 2026) on Ryan Shea's complex588 tree (`GRoute.lift`, `xerase`) and p-a3's
network; prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace XHostT
open Binary Matrix Finset RAM SS CB RF BR GF GX XHost XHostD
noncomputable section

section Box
variable {T Sl C C' : Type} [Fintype T] [DecidableEq T] [Fintype Sl] [DecidableEq Sl]
  [Fintype C] [DecidableEq C] [Fintype C'] [DecidableEq C']

/-- the box of an invocation inside the box with the copies re-indexed by `e`. -/
def boxMap (e : C ↪ C') : Box T Sl C → Box T Sl C' := Sum.map id (Sum.map id e)

/-- the same as an embedding. -/
def boxEmb (e : C ↪ C') : Box T Sl C ↪ Box T Sl C' :=
  ⟨boxMap e, Function.injective_id.sumMap (Function.injective_id.sumMap e.injective)⟩

/-- the forward shear on the roles of the old box; the private scratch (copies outside the image of `e`) kept. -/
def shearFe (e : C ↪ C') (x : Box T Sl C' → ℂ) : Box T Sl C' → ℂ :=
  stamp (x ∘ bX) (x ∘ bY + x ∘ bX) (x ∘ bS) (fun c' => if c' ∈ Finset.univ.map e then 0 else x (bC c'))

/-- the backward shear on the roles of the old box; the private scratch kept. -/
def shearBe (e : C ↪ C') (x : Box T Sl C' → ℂ) : Box T Sl C' → ℂ :=
  stamp (x ∘ bX - x ∘ bY) (x ∘ bY) (x ∘ bS) (fun c' => if c' ∈ Finset.univ.map e then 0 else x (bC c'))

lemma img_mem (e : C ↪ C') (c : C) : e c ∈ Finset.univ.map e := Finset.mem_map_of_mem e (Finset.mem_univ c)

lemma not_img (e : C ↪ C') {c' : C'} (hi : (bC c' : Box T Sl C') ∉ covered (boxEmb (T := T) (Sl := Sl) e)) :
    c' ∉ Finset.univ.map e := fun h => by
  obtain ⟨c, -, rfl⟩ := Finset.mem_map.mp h
  exact hi ((cover_iff _ _).mpr ⟨bC c, rfl⟩)

lemma shearFe_in (e : C ↪ C') (x : Box T Sl C' → ℂ) :
    shearFe e x ∘ boxEmb e = shearF (x ∘ boxEmb e) := by
  funext i
  rcases i with ((t|t)|(q|c))
  · rfl
  · rfl
  · rfl
  · show (if e c ∈ Finset.univ.map e then (0:ℂ) else x (bC (e c))) = 0
    exact if_pos (img_mem e c)

lemma shearBe_in (e : C ↪ C') (x : Box T Sl C' → ℂ) :
    shearBe e x ∘ boxEmb e = shearB (x ∘ boxEmb e) := by
  funext i
  rcases i with ((t|t)|(q|c))
  · rfl
  · rfl
  · rfl
  · show (if e c ∈ Finset.univ.map e then (0:ℂ) else x (bC (e c))) = 0
    exact if_pos (img_mem e c)

lemma shearFe_out (e : C ↪ C') (x : Box T Sl C' → ℂ) (i : Box T Sl C')
    (hi : i ∉ covered (boxEmb (T := T) (Sl := Sl) e)) : shearFe e x i = x i := by
  rcases i with ((t|t)|(q|c'))
  · exact absurd ((cover_iff _ _).mpr ⟨bX t, rfl⟩) hi
  · exact absurd ((cover_iff _ _).mpr ⟨bY t, rfl⟩) hi
  · exact absurd ((cover_iff _ _).mpr ⟨bS q, rfl⟩) hi
  · show (if c' ∈ Finset.univ.map e then (0:ℂ) else x (bC c')) = x (bC c')
    exact if_neg (not_img e hi)

lemma shearBe_out (e : C ↪ C') (x : Box T Sl C' → ℂ) (i : Box T Sl C')
    (hi : i ∉ covered (boxEmb (T := T) (Sl := Sl) e)) : shearBe e x i = x i := by
  rcases i with ((t|t)|(q|c'))
  · exact absurd ((cover_iff _ _).mpr ⟨bX t, rfl⟩) hi
  · exact absurd ((cover_iff _ _).mpr ⟨bY t, rfl⟩) hi
  · exact absurd ((cover_iff _ _).mpr ⟨bS q, rfl⟩) hi
  · show (if c' ∈ Finset.univ.map e then (0:ℂ) else x (bC c')) = x (bC c')
    exact if_neg (not_img e hi)

/-- after the erase gate the private scratch is zero: the shear of the old box is the shear of the new one. -/
lemma shearFe_erase (e : C ↪ C') (x : Box T Sl C' → ℂ) : shearFe e (eraseC x) = shearF x := by
  funext i
  rcases i with ((t|t)|(q|c'))
  · rfl
  · rfl
  · rfl
  · show (if c' ∈ Finset.univ.map e then (0:ℂ) else (0:ℂ)) = 0
    exact ite_self _

lemma shearBe_erase (e : C ↪ C') (x : Box T Sl C' → ℂ) : shearBe e (eraseC x) = shearB x := by
  funext i
  rcases i with ((t|t)|(q|c'))
  · rfl
  · rfl
  · rfl
  · show (if c' ∈ Finset.univ.map e then (0:ℂ) else (0:ℂ)) = 0
    exact ite_self _

end Box

section Route
variable {H T Sl C C' : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype C'] [DecidableEq C']
variable {α : Type} [Fintype α] [DecidableEq α]

/-- the frames of the copies at the hinge: the old copies as at the start, the private scratch as at the end. -/
def cMid (e : C ↪ C') (c0 c1 : C' → CMat α) (c' : C') : CMat α :=
  if c' ∈ Finset.univ.map e then c0 c' else c1 c'

lemma iin_mid (e : C ↪ C') (Xm : XMap H α) (tv : T → Space H) (Z : Sl → CMat α) (c0 c1 : C' → CMat α) :
    iin Xm tv Z (cMid e c0 c1) ∘ boxEmb e = iin Xm tv Z (c0 ∘ e) := by
  funext i
  rcases i with ((t|t)|(q|c))
  · rfl
  · rfl
  · rfl
  · show cMid e c0 c1 (e c) = c0 (e c)
    exact if_pos (img_mem e c)

lemma iout_emb (e : C ↪ C') (Xm : XMap H α) (tv : T → Space H) (Z : Sl → CMat α) (c1 : C' → CMat α) :
    iout Xm tv Z c1 ∘ boxEmb e = iout Xm tv Z (c1 ∘ e) := by
  funext i
  rcases i with ((t|t)|(q|c)) <;> rfl

lemma mid_out (e : C ↪ C') (Xm : XMap H α) (tv : T → Space H) (Z0 Z1 : Sl → CMat α) (c0 c1 : C' → CMat α)
    (i : Box T Sl C') (hi : i ∉ covered (boxEmb (T := T) (Sl := Sl) e)) :
    iin Xm tv Z0 (cMid e c0 c1) i = iout Xm tv Z1 c1 i := by
  rcases i with ((t|t)|(q|c'))
  · exact absurd ((cover_iff _ _).mpr ⟨bX t, rfl⟩) hi
  · exact absurd ((cover_iff _ _).mpr ⟨bY t, rfl⟩) hi
  · exact absurd ((cover_iff _ _).mpr ⟨bS q, rfl⟩) hi
  · show cMid e c0 c1 c' = c1 c'
    exact if_neg (not_img e hi)

lemma erase_mid (e : C ↪ C') (Xm : XMap H α) (tv : T → Space H) (Z : Sl → CMat α) (c0 c1 : C' → CMat α) :
    GRoute (iin Xm tv Z c0) (iin Xm tv Z (cMid e c0 c1)) eraseC (fun _ => 0) :=
  GRoute.ofX (xerase (fun t => rfl) (fun t => rfl) (fun q => rfl))

/-- **a forward invocation route with more copies**: the old copies as before, the others private scratch. -/
theorem fwd_reC (e : C ↪ C') (Xm : XMap H α) (tv : T → Space H) (Z0 Z1 : Sl → CMat α) (c0 c1 : C' → CMat α)
    {cost : (ℕ → ℝ) → ℝ}
    (h : GRoute (iin Xm tv Z0 (c0 ∘ e)) (iout Xm tv Z1 (c1 ∘ e)) shearF cost) :
    GRoute (iin Xm tv Z0 c0) (iout Xm tv Z1 c1) shearF cost := by
  have e1 : GRoute (iin Xm tv Z0 (cMid e c0 c1)) (iout Xm tv Z1 c1) (shearFe e) cost := by
    refine GRoute.lift (boxEmb e) ?_ (fun x => shearFe_in e x) (fun x i hi => shearFe_out e x i hi)
      (fun i hi => mid_out e Xm tv Z0 Z1 c0 c1 i hi)
    rw [iin_mid, iout_emb]
    exact h
  exact (((erase_mid e Xm tv Z0 c0 c1).trans e1).castg (funext fun x => shearFe_erase e x)).cast
    (fun φ => by simp only [zero_add])

/-- **a backward invocation route with more copies**. -/
theorem bwd_reC (e : C ↪ C') (Xm : XMap H α) (tv : T → Space H) (Z0 Z1 : Sl → CMat α) (c0 c1 : C' → CMat α)
    {cost : (ℕ → ℝ) → ℝ}
    (h : GRoute (iin Xm tv Z0 (c0 ∘ e)) (iout Xm tv Z1 (c1 ∘ e)) shearB cost) :
    GRoute (iin Xm tv Z0 c0) (iout Xm tv Z1 c1) shearB cost := by
  have e1 : GRoute (iin Xm tv Z0 (cMid e c0 c1)) (iout Xm tv Z1 c1) (shearBe e) cost := by
    refine GRoute.lift (boxEmb e) ?_ (fun x => shearBe_in e x) (fun x i hi => shearBe_out e x i hi)
      (fun i hi => mid_out e Xm tv Z0 Z1 c0 c1 i hi)
    rw [iin_mid, iout_emb]
    exact h
  exact (((erase_mid e Xm tv Z0 c0 c1).trans e1).castg (funext fun x => shearBe_erase e x)).cast
    (fun φ => by simp only [zero_add])

end Route

section Pkg
variable {H T Sl C C' : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype C'] [DecidableEq C']

/-- **an invocation package with its copies re-indexed** by `e : C ↪ C'`: the copies outside the image of `e` are
private scratch of every call (any start and end frame, content erased); target vectors, entrance labels and
price unchanged. -/
def _root_.OAI.PowerSaving.GX.GInvE.reC (I : GInvE H T Sl C) (e : C ↪ C') : GInvE H T Sl C' where
  tv := I.tv
  hunit := I.hunit
  hgap := I.hgap
  cost := I.cost
  ent := I.ent
  fwd := fun W N Wn hN hW Z c0 c1 => fwd_reC e W.Xm I.tv Z _ c0 c1 (I.fwd W N Wn hN hW Z (c0 ∘ e) (c1 ∘ e))
  bwd := fun W N Wn hN hW Z c0 c1 => bwd_reC e W.Xm I.tv Z _ c0 c1 (I.bwd W N Wn hN hW Z (c0 ∘ e) (c1 ∘ e))

@[simp] lemma _root_.OAI.PowerSaving.GX.GInvE.reC_tv (I : GInvE H T Sl C) (e : C ↪ C') : (I.reC e).tv = I.tv := rfl
@[simp] lemma _root_.OAI.PowerSaving.GX.GInvE.reC_ent (I : GInvE H T Sl C) (e : C ↪ C') :
    (I.reC e).ent = I.ent := rfl
@[simp] lemma _root_.OAI.PowerSaving.GX.GInvE.reC_cost (I : GInvE H T Sl C) (e : C ↪ C') :
    (I.reC e).cost = I.cost := rfl

/-- **`Z` private scratch roles added to every call of `I`** (e.g. `Z = Fin 91`, the zero buffers). -/
def _root_.OAI.PowerSaving.GX.GInvE.addScratch (I : GInvE H T Sl C) (Z : Type) [Fintype Z] [DecidableEq Z] :
    GInvE H T Sl (C ⊕ Z) := I.reC Function.Embedding.inl

/-- the exit-frame package (x-exit-lean's `GInvX`, the form of p-warpchk's `WValid.inv`) with its copies
re-indexed. -/
def _root_.OAI.PowerSaving.XEX.GInvX.reC (I : XEX.GInvX H T Sl C) (e : C ↪ C') : XEX.GInvX H T Sl C' where
  tv := I.tv
  hunit := I.hunit
  hgap := I.hgap
  cost := I.cost
  ent := I.ent
  ex := I.ex
  hexG := I.hexG
  hexs := I.hexs
  fwd := fun W N Wn hN hW Z c0 c1 => fwd_reC e W.Xm I.tv Z _ c0 c1 (I.fwd W N Wn hN hW Z (c0 ∘ e) (c1 ∘ e))
  bwd := fun W N Wn hN hW Z c0 c1 => bwd_reC e W.Xm I.tv Z _ c0 c1 (I.bwd W N Wn hN hW Z (c0 ∘ e) (c1 ∘ e))

@[simp] lemma _root_.OAI.PowerSaving.XEX.GInvX.reC_cost (I : XEX.GInvX H T Sl C) (e : C ↪ C') :
    (I.reC e).cost = I.cost := rfl
@[simp] lemma _root_.OAI.PowerSaving.XEX.GInvX.reC_ent (I : XEX.GInvX H T Sl C) (e : C ↪ C') :
    (I.reC e).ent = I.ent := rfl
@[simp] lemma _root_.OAI.PowerSaving.XEX.GInvX.reC_ex (I : XEX.GInvX H T Sl C) (e : C ↪ C') :
    (I.reC e).ex = I.ex := rfl

/-- re-indexing the copies commutes with x-exit-lean's reduction to a gauged package (same labels, same price). -/
lemma _root_.OAI.PowerSaving.XEX.GInvX.reC_toE' (I : XEX.GInvX H T Sl C) (e : C ↪ C') :
    (I.reC e).toE'.ent = I.toE'.ent ∧ (I.reC e).toE'.cost = I.toE'.cost ∧ (I.reC e).toE'.tv = I.toE'.tv :=
  ⟨rfl, rfl, rfl⟩

end Pkg

section Private
variable {Γ T Sl C : Type} {nb nR : ℕ}

/-- the role of a slot of a call (a data role of a host replica, or a bank) is never a copy: the private scratch
is never hosted. -/
lemma slotT_ne_copy (hs : Sl → Option (Fin nR × Fin 9 × (Bool × T))) (u : (q : Sl) → hs q = none → Fin nb)
    (g : Γ) (q : Sl) (g' : Γ) (c : C) : slotT (C := C) hs u g q ≠ Sum.inr (g', c) := by
  rcases slotT_cases (C := C) hs u g q with ⟨_, e⟩ | ⟨_, _, _, _, _, e⟩ <;> rw [e] <;> exact Sum.inl_ne_inr

/-- the private scratch is not live (it is scratch of the certificate). -/
lemma copy_not_live (g : Γ) (c : C) (l : BLiveT Γ T nb nR) :
    (Sum.inl l : BRoleT Γ T nb C nR) ≠ Sum.inr (g, c) := Sum.inl_ne_inr

end Private
end
end XHostT
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHostT.fwd_reC
#print axioms OAI.PowerSaving.XHostT.bwd_reC
#print axioms OAI.PowerSaving.GX.GInvE.reC
#print axioms OAI.PowerSaving.GX.GInvE.addScratch
#print axioms OAI.PowerSaving.XEX.GInvX.reC
#print axioms OAI.PowerSaving.XEX.GInvX.reC_toE'
#print axioms OAI.PowerSaving.XHostT.slotT_ne_copy
