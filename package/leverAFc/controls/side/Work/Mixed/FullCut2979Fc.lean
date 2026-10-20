import Work.Mixed.HeadM1
import Work.Mixed.Hist2979Fc

/-!
# The mixed certificate of the FULL cut (M1 composed), and its separation on A-fc (suffix 2979Fc)

afc-cut (wht45/afc, 20 Oct 2026), fc-rate's `FullCut3945F.lean` with A-fc's numbers (the generic lemmas are repeated
here with the suffix `2979Fc`, so that this file does not import design B's data).  e-r45m1's M1 lemmas
(`Work.Mixed.HeadM1`, unchanged) composed for the full cut:

* `fullCutWord2979Fc e L A S body = splitWord 24 24 (body ++ shiftWord e L (cutShift A S))`: the heads removed, one
  free shift per head role, every rank-48 block that is left split 24 + 24 (on A-fc there is none);
* `cut_heads_full_split2979Fc`: **the mixed certificate** of that word, from the unit certificate
  `GLiveKernel e (headWord e L rk z ind ++ body)` and the frame/lines hypothesis `hz` (M1's `cut_heads_full`, then
  `GLiveKernelM.splitWord`); target `blockMat (cutTgt L A S l).2`;
* `cutTgt_rank_full2979Fc`: with `card α = 73` and `|S l| = 48` the target ranks are 25 on head roles, 73 elsewhere;
* `hist_fullCutWord2979Fc`: its histogram is `hist_splitWord 24 24` of the body's (the shifts are free), and
  `hist_fullCutWord_noSplit2979Fc`: **with no rank-48 block left in the body the split changes nothing** (the
  histogram of the word is the body's);
* `proper25_of_hist2979Fc`: every paid rank of a word with `G` copies of `cut2979Fc` is in `[1, 24]` (in fact
  `[1, 23]`), below the cut target 25;
* `fullCut2979Fc_cert`: all of it for A-fc: if the unit word `headWord ++ body` has `G` copies of A-fc's unit profile
  and `G * 16576` heads, then the full-cut word has the mixed certificate, target ranks {25, 73}, the histogram `G`
  copies of R5's `cut2979Fc` (= the body's), and every paid rank in `[1, 24]`.
Prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false

namespace OAI
namespace PowerSaving.MX
open Binary GF Matrix Finset BlockAccounting

section Words
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]

/-- The full-cut word: the body (heads removed), one free shift per head role, every rank-48 block split 24 + 24. -/
noncomputable def fullCutWord2979Fc (e : σ → ρ) (L : List σ) (A : σ → OBase α) (S : σ → Finset α)
    (body : GWord α ρ) : GWord α ρ :=
  splitWord 24 24 (body ++ shiftWord e L (cutShift A S))

/-- **The mixed certificate of the full cut.** -/
theorem cut_heads_full_split2979Fc [DecidableEq σ] {e : σ → ρ} (he : Function.Injective e)
    (L : List σ) (hL : L.Nodup) (rk : σ → ℕ) (z : ∀ l, Fin (rk l) → Space α)
    (ind : ∀ l, LinearIndependent F (z l)) (A : σ → OBase α) (S : σ → Finset α)
    (hz : ∀ l ∈ L, ∀ x, bph (z l) x = bph ((A l).fam (S l)) x) (body : GWord α ρ)
    (hk : GLiveKernel e (headWord e L rk z ind ++ body)) :
    GLiveKernelM e (fun l => blockMat (cutTgt L A S l).2) (fullCutWord2979Fc e L A S body) :=
  (cut_heads_full he L hL rk z ind A S hz body hk).splitWord 24 24

/-- Target ranks of the full cut: 25 on the head roles, 73 elsewhere. -/
theorem cutTgt_rank_full2979Fc [DecidableEq σ] (L : List σ) (A : σ → OBase α) (S : σ → Finset α)
    (hα : Fintype.card α = 73) (hS : ∀ l ∈ L, (S l).card = 48) (l : σ) :
    (cutTgt L A S l).1 = if l ∈ L then 25 else 73 := by
  rw [cutTgt_rank]
  by_cases hl : l ∈ L
  · simp [hl, hα, hS l hl]
  · simp [hl, hα]

/-- Histogram of the full-cut word: `hist_splitWord 24 24` of the body's (the shifts are free). -/
theorem hist_fullCutWord2979Fc (e : σ → ρ) (L : List σ) (A : σ → OBase α) (S : σ → Finset α) (body : GWord α ρ)
    (r : ℕ) :
    GWord.hist (fullCutWord2979Fc e L A S body) r
      = (if r = 24 + 24 then 0 else GWord.hist body r) + (if r = 24 then GWord.hist body (24 + 24) else 0)
        + (if r = 24 then GWord.hist body (24 + 24) else 0) := by
  unfold fullCutWord2979Fc
  have h1 := hist_splitWord 24 24 (by norm_num) (by norm_num) (body ++ shiftWord e L (cutShift A S)) r
  rw [h1]
  simp only [hist_append, hist_shiftWord, add_zero]

/-- **The split of nothing**: if the body has no rank-48 block, the full-cut word has the body's histogram. -/
theorem hist_fullCutWord_noSplit2979Fc (e : σ → ρ) (L : List σ) (A : σ → OBase α) (S : σ → Finset α)
    (body : GWord α ρ) (h0 : GWord.hist body (24 + 24) = 0) (r : ℕ) :
    GWord.hist (fullCutWord2979Fc e L A S body) r = GWord.hist body r := by
  have h := hist_fullCutWord2979Fc e L A S body r
  rw [h0] at h
  by_cases hr : r = 24 + 24
  · have h' : GWord.hist body r = 0 := by rw [hr]; exact h0
    have hr24 : r ≠ 24 := by omega
    rw [if_pos hr, if_neg hr24] at h
    omega
  · rw [if_neg hr] at h
    by_cases hr24 : r = 24
    · rw [if_pos hr24] at h
      omega
    · rw [if_neg hr24] at h
      omega

/-- Properness of the full-cut word from that of the body, for any bound `u > 24`. -/
theorem proper_fullCutWord2979Fc (u : ℕ) (hu : 24 < u) (e : σ → ρ) (L : List σ) (A : σ → OBase α)
    (S : σ → Finset α) (body : GWord α ρ) (hbody : GWord.Proper u body) :
    GWord.Proper u (fullCutWord2979Fc e L A S body) := by
  unfold fullCutWord2979Fc
  apply proper_splitWord u 24 24 (by norm_num) (by norm_num) hu hu
  intro mv hmv
  rcases List.mem_append.mp hmv with h | h
  · exact hbody mv h
  · exact proper_shiftWord u e L (cutShift A S) mv h

/-- A paid letter of a word is counted by the word's histogram. -/
theorem ghist_pos_of_mem2979Fc {w : GWord α ρ} {mv : GMove α ρ} (hmv : mv ∈ w) {r : ℕ}
    (hr : mv.shadow.rk? = some r) : 0 < GWord.hist w r := by
  unfold GWord.hist BWord.hist GWord.shadow
  rw [List.countP_pos_iff]
  exact ⟨mv.shadow, List.mem_map.mpr ⟨mv, hmv, rfl⟩, by simp [hr]⟩

/-- **Separation on the actual word**: a word whose histogram is `G` copies of `cut2979Fc` has every paid rank in
`[1, 24]` (in fact `[1, 23]`), i.e. below the cut target 25. -/
theorem proper25_of_hist2979Fc {w : GWord α ρ} {G : ℕ} (h : ∀ r, GWord.hist w r = G * histFn cut2979Fc r) :
    GWord.Proper 25 w := by
  intro mv hmv
  cases mv with
  | old m =>
    cases m with
    | block l rk z ind =>
      have hpos := ghist_pos_of_mem2979Fc hmv (r := rk) rfl
      rw [h rk] at hpos
      have hne : histFn cut2979Fc rk ≠ 0 := by
        intro h0
        rw [h0, mul_zero] at hpos
        exact lt_irrefl _ hpos
      have hk := mem_keys_of_histFn_ne_zero cut2979Fc rk hne
      unfold keys at hk
      rw [List.mem_toFinset, List.mem_map] at hk
      obtain ⟨b, hb, hbr⟩ := hk
      have hs := cut2979Fc_sep b hb
      show 1 ≤ rk ∧ rk < 25
      rw [← hbr]
      exact hs
    | localM g => trivial
    | shift l z => trivial
  | perm l Gp => trivial
  | phase l z c => trivial

/-- The body's histogram from "unit minus the heads". -/
theorem body_hist2979Fc {hb : ℕ → ℕ} {G : ℕ}
    (h : ∀ r, hb r + (if r = 24 + 24 then G * 16576 else 0) = G * histFn unit2979Fc r) :
    ∀ r, hb r = G * histFn body2979Fc r := by
  intro r
  have h1 := h r
  have h2 := hbody2979Fc r
  by_cases hr : r = 24 + 24
  · rw [if_pos hr] at h1 h2
    have h3 : G * histFn unit2979Fc r = G * histFn body2979Fc r + G * 16576 := by rw [← h2, mul_add]
    omega
  · rw [if_neg hr] at h1 h2
    rw [← h2] at h1
    simpa using h1

/-- The full-cut word's histogram on A-fc: `G` copies of R5's `cut2979Fc`. -/
theorem hist_fullCutWord2979Fc_cut (e : σ → ρ) (L : List σ) (A : σ → OBase α) (S : σ → Finset α)
    (body : GWord α ρ) {G : ℕ} (hb : ∀ r, GWord.hist body r = G * histFn body2979Fc r) :
    ∀ r, GWord.hist (fullCutWord2979Fc e L A S body) r = G * histFn cut2979Fc r := by
  intro r
  rw [hist_fullCutWord2979Fc, hb r, hb (24 + 24), hsplit2979Fc r]
  split_ifs <;> ring

/-- **The full cut on A-fc, in one statement.**  From the unit certificate of `headWord ++ body` (heads of rank
48 = 24 + 24 on the lines `S l` of frames `A l`, `|S l| = 48`, `card α = 73`), the unit histogram (`G` copies of A-fc's
profile) and the head count `G * 16576`: the full-cut word has the mixed certificate with targets of rank 25 (heads)
and 73 (elsewhere), the split changes nothing (its histogram is the body's), the histogram is `G` copies of
`cut2979Fc`, and every paid rank lies in `[1, 24]` (`23 < 25`). -/
theorem fullCut2979Fc_cert [DecidableEq σ] {e : σ → ρ} (he : Function.Injective e) (L : List σ) (hL : L.Nodup)
    (z : σ → Fin (24 + 24) → Space α) (ind : ∀ l, LinearIndependent F (z l)) (A : σ → OBase α) (S : σ → Finset α)
    (hα : Fintype.card α = 73) (hS : ∀ l ∈ L, (S l).card = 48)
    (hz : ∀ l ∈ L, ∀ x, bph (z l) x = bph ((A l).fam (S l)) x) (body : GWord α ρ)
    (hk : GLiveKernel e (headWord e L (fun _ => 24 + 24) z ind ++ body)) {G : ℕ}
    (hU : ∀ r, GWord.hist (headWord e L (fun _ => 24 + 24) z ind ++ body) r = G * histFn unit2979Fc r)
    (hLlen : L.length = G * 16576) :
    GLiveKernelM e (fun l => blockMat (cutTgt L A S l).2) (fullCutWord2979Fc e L A S body) ∧
      (∀ l, (cutTgt L A S l).1 = if l ∈ L then 25 else 73) ∧
      (∀ r, GWord.hist (fullCutWord2979Fc e L A S body) r = GWord.hist body r) ∧
      (∀ r, GWord.hist (fullCutWord2979Fc e L A S body) r = G * histFn cut2979Fc r) ∧
      GWord.Proper 25 (fullCutWord2979Fc e L A S body) := by
  have hb : ∀ r, GWord.hist body r + (if r = 24 + 24 then G * 16576 else 0) = G * histFn unit2979Fc r := by
    intro r
    rw [← hU r, hist_append, hist_headWord, hLlen]
    exact Nat.add_comm _ _
  have hbody := body_hist2979Fc hb
  have h0 : GWord.hist body (24 + 24) = 0 := by rw [hbody, body2979Fc_no48, mul_zero]
  have hh := hist_fullCutWord2979Fc_cut e L A S body hbody
  exact ⟨cut_heads_full_split2979Fc he L hL (fun _ => 24 + 24) z ind A S hz body hk,
    cutTgt_rank_full2979Fc L A S hα hS, hist_fullCutWord_noSplit2979Fc e L A S body h0, hh,
    proper25_of_hist2979Fc hh⟩

end Words

end PowerSaving.MX
end OAI

#print axioms OAI.PowerSaving.MX.cut_heads_full_split2979Fc
#print axioms OAI.PowerSaving.MX.cutTgt_rank_full2979Fc
#print axioms OAI.PowerSaving.MX.hist_fullCutWord2979Fc
#print axioms OAI.PowerSaving.MX.hist_fullCutWord_noSplit2979Fc
#print axioms OAI.PowerSaving.MX.proper_fullCutWord2979Fc
#print axioms OAI.PowerSaving.MX.ghist_pos_of_mem2979Fc
#print axioms OAI.PowerSaving.MX.proper25_of_hist2979Fc
#print axioms OAI.PowerSaving.MX.body_hist2979Fc
#print axioms OAI.PowerSaving.MX.hist_fullCutWord2979Fc_cut
#print axioms OAI.PowerSaving.MX.fullCut2979Fc_cert
