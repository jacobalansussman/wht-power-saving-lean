import Work.Mixed.Bridge.Inst3945H
import Work.Mixed.Bridge.HeadsFirst

/-!
# bB-data: the design-B (#3945H) bridge in the exact shape of p-hfirst's `XNet.heads_first` conclusion (round wht45/allheads)

x-inst's `HeadsFirst2979.lean` (p-bridge's `HeadsFirst.lean` with other numbers) with the design-B instance (`hfnetOf` is p-bridge's, imported).  `engine3945H_of_heads`
takes p-hfirst's existential (generic `e`, `L`, `G`, price `c0`), the unit profile (`c0 φ = listCost φ U`), the head
frames and line sets, and returns the engine theorem at the seam exponent `1 - 1758152001/10^12`
(DFT figure 1.7581520e-3).  Nothing here depends on the network type: any heads-first certificate of this shape
(plain providers or warp providers) plugs in.
Written by x-inst (wht45/lever2979) on p-bridge's file; instance 3945H by bB-data (wht45/allheads, 19 Oct 2026, design B: #3945 with every child-8/10 bin a literal head); prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.PB
open Binary Matrix Finset GF MX RAM

noncomputable section

universe U

section HF
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]

/-- **THE ENGINE THEOREM FOR DESIGN B (#3945H) FROM p-hfirst's CONCLUSION.**  `hnet` is `XNet.heads_first`'s conclusion
(`e = Sum.inl`, `L = headsL hb`, `G = card Γ`, `c0 = N.unit`); `hU` the unit profile of the network; `hcut` the
profile identity (`hcut3945H` for `unit3945H`); `hLlen`, `hσ` the counts. -/
theorem engine3945H_of_heads (e : σ → ρ) (he : Function.Injective e) (L : List σ) (hL : L.Nodup)
    (z : σ → Fin 48 → Space α) (ind : ∀ l, LinearIndependent F (z l)) (A : σ → OBase α) (d : σ → Finset α)
    (hd : ∀ l ∈ L, (d l).card = 48)
    (hzd : ∀ l (hl : l ∈ L), ∀ i, z l i = (A l).fam (d l) (Fin.cast (hd l hl).symm i))
    (G : ℕ) (c0 : (ℕ → ℝ) → ℝ) (U : List (ℕ × ℕ)) (hU : ∀ φ, c0 φ = listCost φ U)
    (hnet : ∃ body : GWord α ρ, body.Proper (Fintype.card α) ∧
      (∀ φ, body.costR φ = (G : ℝ) * c0 φ - (L.length : ℝ) * φ 48) ∧
      GLiveKernel e (headWord e L (fun _ => 48) z ind ++ body))
    (hα : Fintype.card α = 73)
    (hcut : ∀ r, (if r = 24 then 80640 else 0) + BlockAccounting.histFn U r
      = BlockAccounting.histFn cut3945H r + (if r = 24 + 24 then 80640 else 0))
    (hLlen : L.length = G * 80640) (hσ : Fintype.card σ = G * 322560) (hG : 1 ≤ G)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1758152001:ℝ)/(10:ℝ)^12) (k i)) := by
  obtain ⟨body, hP, hcost, hk⟩ := hnet
  have hcost' : ∀ φ, body.costR φ = (G : ℝ) * listCost φ U - (L.length : ℝ) * φ 48 := fun φ => by
    rw [hcost φ, hU φ]
  exact (hfnetOf e he L hL z ind A d hd hzd G U body hP hcost' hk).engine3945H hα rfl rfl hcut hLlen hσ hG m cl k v

/-- (f) for design B (#3945H) from p-hfirst's conclusion: the hypothesis bundle of `Eng.engine_seam_final`. -/
theorem bundle3945H_of_heads (e : σ → ρ) (he : Function.Injective e) (L : List σ) (hL : L.Nodup)
    (z : σ → Fin 48 → Space α) (ind : ∀ l, LinearIndependent F (z l)) (A : σ → OBase α) (d : σ → Finset α)
    (hd : ∀ l ∈ L, (d l).card = 48)
    (hzd : ∀ l (hl : l ∈ L), ∀ i, z l i = (A l).fam (d l) (Fin.cast (hd l hl).symm i))
    (G : ℕ) (c0 : (ℕ → ℝ) → ℝ) (U : List (ℕ × ℕ)) (hU : ∀ φ, c0 φ = listCost φ U)
    (hnet : ∃ body : GWord α ρ, body.Proper (Fintype.card α) ∧
      (∀ φ, body.costR φ = (G : ℝ) * c0 φ - (L.length : ℝ) * φ 48) ∧
      GLiveKernel e (headWord e L (fun _ => 48) z ind ++ body))
    (hα : Fintype.card α = 73)
    (hcut : ∀ r, (if r = 24 then 80640 else 0) + BlockAccounting.histFn U r
      = BlockAccounting.histFn cut3945H r + (if r = 24 + 24 then 80640 else 0))
    (hLlen : L.length = G * 80640) (hσ : Fintype.card σ = G * 322560) (hG : 1 ≤ G) :
    ∃ (E : Eng α ρ σ) (R : MX.Rule) (typ : ℕ → Fin R.K),
      E.e = e ∧ E.mx = 73 ∧ E.BridgeHyp R typ 49 ∧ R.RateQ zRate3945H 1 := by
  obtain ⟨body, hP, hcost, hk⟩ := hnet
  have hcost' : ∀ φ, body.costR φ = (G : ℝ) * listCost φ U - (L.length : ℝ) * φ 48 := fun φ => by
    rw [hcost φ, hU φ]
  obtain ⟨E, R, typ, -, h2, -, h4, hB, hR⟩ :=
    (hfnetOf e he L hL z ind A d hd hzd G U body hP hcost' hk).bundle3945H hα rfl rfl hcut hLlen hσ hG
  exact ⟨E, R, typ, h2, h4, hB, hR⟩

/-- **The engine theorem for design B (#3945H) with the cut data given directly** (p-hfirst's `cutA`, `cutS`, `bph_cut`: the
form `MX.cut_heads_half` takes, `a = b = 24`). -/
theorem engine3945H_of_cut (e : σ → ρ) (he : Function.Injective e) (L : List σ) (hL : L.Nodup)
    (z : σ → Fin 48 → Space α) (ind : ∀ l, LinearIndependent F (z l)) (A : σ → OBase α) (S : σ → Finset α)
    (hz : ∀ l ∈ L, ∀ x, bph (partA (a := 24) (b := 24) z l) x = bph ((A l).fam (S l)) x)
    (hS : ∀ l ∈ L, (S l).card = 24)
    (G : ℕ) (c0 : (ℕ → ℝ) → ℝ) (U : List (ℕ × ℕ)) (hU : ∀ φ, c0 φ = listCost φ U)
    (hnet : ∃ body : GWord α ρ, body.Proper (Fintype.card α) ∧
      (∀ φ, body.costR φ = (G : ℝ) * c0 φ - (L.length : ℝ) * φ 48) ∧
      GLiveKernel e (headWord e L (fun _ => 48) z ind ++ body))
    (hα : Fintype.card α = 73)
    (hcut : ∀ r, (if r = 24 then 80640 else 0) + BlockAccounting.histFn U r
      = BlockAccounting.histFn cut3945H r + (if r = 24 + 24 then 80640 else 0))
    (hLlen : L.length = G * 80640) (hσ : Fintype.card σ = G * 322560) (hG : 1 ≤ G)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1758152001:ℝ)/(10:ℝ)^12) (k i)) := by
  obtain ⟨body, hP, hcost, hk⟩ := hnet
  let N : HFNet α ρ σ :=
    { e := e, he := he, L := L, hL := hL, a := 24, b := 24, z := z, ind := ind, A := A, S := S, hz := hz,
      hS := hS, body := body, hP := hP, hk := hk, G := G, U := U,
      hcost := fun φ => by rw [hcost φ, hU φ] }
  exact N.engine3945H hα rfl rfl hcut hLlen hσ hG m cl k v

end HF

end
end PowerSaving.PB
end OAI

#print axioms OAI.PowerSaving.PB.engine3945H_of_heads
#print axioms OAI.PowerSaving.PB.bundle3945H_of_heads
#print axioms OAI.PowerSaving.PB.engine3945H_of_cut
