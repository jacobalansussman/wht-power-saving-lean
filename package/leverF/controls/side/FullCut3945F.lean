import Work.Mixed.HeadM1
import Work.Mixed.Hist3945F

/-!
# The mixed certificate of the FULL cut (M1 composed), and its separation on design B (suffix 3945F)

fc-rate (wht45/fullcut, 19 Oct 2026).  e-r45m1's M1 lemmas (`Work.Mixed.HeadM1`, unchanged) composed for the full cut:

* `fullCutWord e L A S body = splitWord 24 24 (body ++ shiftWord e L (cutShift A S))`: the heads removed, one free
  shift per head role, every rank-48 block that is left split 24 + 24;
* `cut_heads_full_split`: **the mixed certificate** of that word, from the unit certificate
  `GLiveKernel e (headWord e L rk z ind ++ body)` and the frame/lines hypothesis `hz` (M1's `cut_heads_full`, then
  `GLiveKernelM.splitWord`); target `blockMat (cutTgt L A S l).2` (the lines `(S l)ᶜ` on head roles, the kernel elsewhere);
* `cutTgt_rank_full`: with `card α = 73` and `|S l| = 48` the target ranks are 25 on head roles and 73 elsewhere;
* `hist_fullCutWord`: its histogram is `hist_splitWord 24 24` of the body's (the shifts are free);
* `proper_fullCutWord`, `proper25_of_hist3945F`: every paid rank of the word is in `[1, 24]`, below the cut target 25
  (separation `24 < 25`, here from the histogram, so for the actual word);
* `fullCut3945F_cert`: all of it for design B: if the unit word `headWord ++ body` has `G` copies of #3945's unit profile
  and `G * 80640` heads, then the full-cut word has the mixed certificate, target ranks {25, 73}, the histogram `G` copies
  of R5's `cut3945F`, and every paid rank in `[1, 24]`.
Prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false

namespace OAI
namespace PowerSaving.MX
open Binary GF Matrix Finset BlockAccounting

section Words
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]

/-- The full-cut word: the body (heads removed), one free shift per head role, every rank-48 block split 24 + 24. -/
noncomputable def fullCutWord (e : σ → ρ) (L : List σ) (A : σ → OBase α) (S : σ → Finset α) (body : GWord α ρ) : GWord α ρ :=
  splitWord 24 24 (body ++ shiftWord e L (cutShift A S))

/-- **The mixed certificate of the full cut.** -/
theorem cut_heads_full_split [DecidableEq σ] {e : σ → ρ} (he : Function.Injective e)
    (L : List σ) (hL : L.Nodup) (rk : σ → ℕ) (z : ∀ l, Fin (rk l) → Space α)
    (ind : ∀ l, LinearIndependent F (z l)) (A : σ → OBase α) (S : σ → Finset α)
    (hz : ∀ l ∈ L, ∀ x, bph (z l) x = bph ((A l).fam (S l)) x) (body : GWord α ρ)
    (hk : GLiveKernel e (headWord e L rk z ind ++ body)) :
    GLiveKernelM e (fun l => blockMat (cutTgt L A S l).2) (fullCutWord e L A S body) :=
  (cut_heads_full he L hL rk z ind A S hz body hk).splitWord 24 24

/-- Target ranks of the full cut: 25 on the head roles, 73 elsewhere. -/
theorem cutTgt_rank_full [DecidableEq σ] (L : List σ) (A : σ → OBase α) (S : σ → Finset α)
    (hα : Fintype.card α = 73) (hS : ∀ l ∈ L, (S l).card = 48) (l : σ) :
    (cutTgt L A S l).1 = if l ∈ L then 25 else 73 := by
  rw [cutTgt_rank]
  by_cases hl : l ∈ L
  · simp [hl, hα, hS l hl]
  · simp [hl, hα]

/-- Histogram of the full-cut word: `hist_splitWord 24 24` of the body's (the shifts are free). -/
theorem hist_fullCutWord (e : σ → ρ) (L : List σ) (A : σ → OBase α) (S : σ → Finset α) (body : GWord α ρ) (r : ℕ) :
    GWord.hist (fullCutWord e L A S body) r
      = (if r = 24 + 24 then 0 else GWord.hist body r) + (if r = 24 then GWord.hist body (24 + 24) else 0)
        + (if r = 24 then GWord.hist body (24 + 24) else 0) := by
  unfold fullCutWord
  have h1 := hist_splitWord 24 24 (by norm_num) (by norm_num) (body ++ shiftWord e L (cutShift A S)) r
  rw [h1]
  simp only [hist_append, hist_shiftWord, add_zero]

/-- Properness of the full-cut word from that of the body, for any bound `u > 24`. -/
theorem proper_fullCutWord (u : ℕ) (hu : 24 < u) (e : σ → ρ) (L : List σ) (A : σ → OBase α) (S : σ → Finset α)
    (body : GWord α ρ) (hbody : GWord.Proper u body) : GWord.Proper u (fullCutWord e L A S body) := by
  unfold fullCutWord
  apply proper_splitWord u 24 24 (by norm_num) (by norm_num) hu hu
  intro mv hmv
  rcases List.mem_append.mp hmv with h | h
  · exact hbody mv h
  · exact proper_shiftWord u e L (cutShift A S) mv h

/-- A paid letter of a word is counted by the word's histogram. -/
theorem ghist_pos_of_mem {w : GWord α ρ} {mv : GMove α ρ} (hmv : mv ∈ w) {r : ℕ} (hr : mv.shadow.rk? = some r) :
    0 < GWord.hist w r := by
  unfold GWord.hist BWord.hist GWord.shadow
  rw [List.countP_pos_iff]
  exact ⟨mv.shadow, List.mem_map.mpr ⟨mv, hmv, rfl⟩, by simp [hr]⟩

/-- **Separation on the actual word**: a word whose histogram is `G` copies of `cut3945F` has every paid rank in
`[1, 24]`, i.e. below the cut target 25. -/
theorem proper25_of_hist3945F {w : GWord α ρ} {G : ℕ} (h : ∀ r, GWord.hist w r = G * histFn cut3945F r) :
    GWord.Proper 25 w := by
  intro mv hmv
  cases mv with
  | old m =>
    cases m with
    | block l rk z ind =>
      have hpos := ghist_pos_of_mem hmv (r := rk) rfl
      rw [h rk] at hpos
      have hne : histFn cut3945F rk ≠ 0 := by
        intro h0
        rw [h0, mul_zero] at hpos
        exact lt_irrefl _ hpos
      have hk := mem_keys_of_histFn_ne_zero cut3945F rk hne
      unfold keys at hk
      rw [List.mem_toFinset, List.mem_map] at hk
      obtain ⟨b, hb, hbr⟩ := hk
      have hs := cut3945F_sep b hb
      show 1 ≤ rk ∧ rk < 25
      rw [← hbr]
      exact hs
    | localM g => trivial
    | shift l z => trivial
  | perm l Gp => trivial
  | phase l z c => trivial

/-- The body's histogram from "unit minus the heads". -/
theorem body_hist3945F {hb : ℕ → ℕ} {G : ℕ}
    (h : ∀ r, hb r + (if r = 24 + 24 then G * 80640 else 0) = G * histFn unit3945F r) :
    ∀ r, hb r = G * histFn body3945F r := by
  intro r
  have h1 := h r
  have h2 := hbody3945F r
  by_cases hr : r = 24 + 24
  · rw [if_pos hr] at h1 h2
    have h3 : G * histFn unit3945F r = G * histFn body3945F r + G * 80640 := by rw [← h2, mul_add]
    omega
  · rw [if_neg hr] at h1 h2
    rw [← h2] at h1
    simpa using h1

/-- The full-cut word's histogram on design B: `G` copies of R5's `cut3945F`. -/
theorem hist_fullCutWord3945F (e : σ → ρ) (L : List σ) (A : σ → OBase α) (S : σ → Finset α) (body : GWord α ρ)
    {G : ℕ} (hb : ∀ r, GWord.hist body r = G * histFn body3945F r) :
    ∀ r, GWord.hist (fullCutWord e L A S body) r = G * histFn cut3945F r := by
  intro r
  rw [hist_fullCutWord, hb r, hb (24 + 24), hsplit3945F r]
  split_ifs <;> ring

/-- **The full cut on design B, in one statement.**  From the unit certificate of `headWord ++ body` (heads of rank
48 = 24 + 24 on the lines `S l` of frames `A l`, `|S l| = 48`, `card α = 73`), the unit histogram (`G` copies of #3945's
profile) and the head count `G * 80640`: the full-cut word has the mixed certificate with targets of rank 25 (heads)
and 73 (elsewhere), its histogram is `G` copies of `cut3945F`, and every paid rank lies in `[1, 24]` (`24 < 25`). -/
theorem fullCut3945F_cert [DecidableEq σ] {e : σ → ρ} (he : Function.Injective e) (L : List σ) (hL : L.Nodup)
    (z : σ → Fin (24 + 24) → Space α) (ind : ∀ l, LinearIndependent F (z l)) (A : σ → OBase α) (S : σ → Finset α)
    (hα : Fintype.card α = 73) (hS : ∀ l ∈ L, (S l).card = 48)
    (hz : ∀ l ∈ L, ∀ x, bph (z l) x = bph ((A l).fam (S l)) x) (body : GWord α ρ)
    (hk : GLiveKernel e (headWord e L (fun _ => 24 + 24) z ind ++ body)) {G : ℕ}
    (hU : ∀ r, GWord.hist (headWord e L (fun _ => 24 + 24) z ind ++ body) r = G * histFn unit3945F r)
    (hLlen : L.length = G * 80640) :
    GLiveKernelM e (fun l => blockMat (cutTgt L A S l).2) (fullCutWord e L A S body) ∧
      (∀ l, (cutTgt L A S l).1 = if l ∈ L then 25 else 73) ∧
      (∀ r, GWord.hist (fullCutWord e L A S body) r = G * histFn cut3945F r) ∧
      GWord.Proper 25 (fullCutWord e L A S body) := by
  have hb : ∀ r, GWord.hist body r + (if r = 24 + 24 then G * 80640 else 0) = G * histFn unit3945F r := by
    intro r
    rw [← hU r, hist_append, hist_headWord, hLlen]
    exact Nat.add_comm _ _
  have hh := hist_fullCutWord3945F e L A S body (body_hist3945F hb)
  exact ⟨cut_heads_full_split he L hL (fun _ => 24 + 24) z ind A S hz body hk,
    cutTgt_rank_full L A S hα hS, hh, proper25_of_hist3945F hh⟩

end Words

end PowerSaving.MX
end OAI

#print axioms OAI.PowerSaving.MX.cut_heads_full_split
#print axioms OAI.PowerSaving.MX.cutTgt_rank_full
#print axioms OAI.PowerSaving.MX.hist_fullCutWord
#print axioms OAI.PowerSaving.MX.proper_fullCutWord
#print axioms OAI.PowerSaving.MX.ghist_pos_of_mem
#print axioms OAI.PowerSaving.MX.proper25_of_hist3945F
#print axioms OAI.PowerSaving.MX.body_hist3945F
#print axioms OAI.PowerSaving.MX.hist_fullCutWord3945F
#print axioms OAI.PowerSaving.MX.fullCut3945F_cert
