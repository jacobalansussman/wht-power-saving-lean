import Work.Mixed.Rate3945F

/-!
# The histogram bookkeeping of the FULL cut of design B (suffix 3945F): unit -> body -> cut

fc-rate (wht45/fullcut, 19 Oct 2026).  Design B (#3945 with the all-heads fill of wht45/reopt) has #3945's unit
profile `unit3945F` (= bB-data's `PB.unit3945H` = p-bridge's `PB.unit3945`, #3945's PROFILE.json); its 137,952
rank-48 blocks are 80,640 literal heads and 57,312 other idle blocks.
* M1's `cut_heads_full` removes the 80,640 heads from `headWord ++ body`: the body has the profile `body3945F`
  (`hbody3945F`);
* M1's `splitWord 24 24` splits each of the remaining 57,312 rank-48 blocks into two rank-24 blocks: by
  `hist_splitWord 24 24` the result has the profile `cut3945F` of R5 (`hsplit3945F` is stated in exactly the shape of
  the right-hand side of `hist_splitWord` with `a = b = 24` and `w` = the body);
* `hcut3945F`: unit -> cut in one identity; `mass3945F`: the rank mass drops by exactly 48 x 80,640;
* separation: `cut3945F_sep` (every paid rank < 25), `cut3945F_no48` (no rank-48 block is left);
* negatives (`decide`): one head more breaks `hbody3945F`, the unsplit body is not the cut profile, the half cut's
  remnants (80,640 rank-24 blocks) are not the full cut's 114,624.
Prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI.PowerSaving.MX
open BlockAccounting

/-- #3945's unit profile (design B changes only the fill, not the profile): `(rank, count)`. -/
def unit3945F : List (ℕ × ℕ) :=
  [(1, 8041680), (2, 2306880), (3, 729648), (4, 353160), (5, 69120), (8, 25920), (48, 137952)]

/-- The body after the full cut removes the 80,640 literal heads: the other 57,312 idle rank-48 blocks remain. -/
def body3945F : List (ℕ × ℕ) :=
  [(1, 8041680), (2, 2306880), (3, 729648), (4, 353160), (5, 69120), (8, 25920), (48, 57312)]

/-- the ranks that occur in the three profiles -/
def keys3945F : List ℕ := [1, 2, 3, 4, 5, 8, 24, 48]

theorem histFn_eq_zero_of_keys {l : List (ℕ × ℕ)} {K : List ℕ} (hl : ∀ q ∈ l, q.1 ∈ K) {r : ℕ} (hr : r ∉ K) :
    histFn l r = 0 := by
  by_contra hne
  have := mem_keys_of_histFn_ne_zero l r hne
  unfold keys at this
  rw [List.mem_toFinset, List.mem_map] at this
  obtain ⟨q, hq, rfl⟩ := this
  exact hr (hl q hq)

theorem unit3945F_keys : ∀ q ∈ unit3945F, q.1 ∈ keys3945F := by decide
theorem body3945F_keys : ∀ q ∈ body3945F, q.1 ∈ keys3945F := by decide
theorem cut3945F_keys : ∀ q ∈ cut3945F, q.1 ∈ keys3945F := by decide

/-- Removing the 80,640 literal heads (rank 48 = 24 + 24) from the unit profile leaves the body profile. -/
theorem hbody3945F : ∀ r, histFn body3945F r + (if r = 24 + 24 then 80640 else 0) = histFn unit3945F r := by
  intro r
  by_cases hr : r ∈ keys3945F
  · have h : ∀ r ∈ keys3945F, histFn body3945F r + (if r = 24 + 24 then 80640 else 0) = histFn unit3945F r := by
      decide
    exact h r hr
  · have h48 : r ≠ 24 + 24 := fun h => hr (by subst h; decide)
    rw [histFn_eq_zero_of_keys body3945F_keys hr, histFn_eq_zero_of_keys unit3945F_keys hr, if_neg h48]

/-- Splitting every rank-48 block of the body into 24 + 24 gives R5's cut profile, in the shape of the right-hand side
of M1's `hist_splitWord 24 24` (with `w` = the body). -/
theorem hsplit3945F : ∀ r, histFn cut3945F r
    = (if r = 24 + 24 then 0 else histFn body3945F r) + (if r = 24 then histFn body3945F (24 + 24) else 0)
      + (if r = 24 then histFn body3945F (24 + 24) else 0) := by
  intro r
  by_cases hr : r ∈ keys3945F
  · have h : ∀ r ∈ keys3945F, histFn cut3945F r
        = (if r = 24 + 24 then 0 else histFn body3945F r) + (if r = 24 then histFn body3945F (24 + 24) else 0)
          + (if r = 24 then histFn body3945F (24 + 24) else 0) := by
      decide
    exact h r hr
  · have h48 : r ≠ 24 + 24 := fun h => hr (by subst h; decide)
    have h24 : r ≠ 24 := fun h => hr (by subst h; decide)
    rw [histFn_eq_zero_of_keys cut3945F_keys hr, histFn_eq_zero_of_keys body3945F_keys hr, if_neg h48, if_neg h24]

/-- The same for any body whose histogram is the body profile (the form M1's lemmas produce). -/
theorem hsplit3945F_of {hb : ℕ → ℕ} (h : ∀ r, hb r + (if r = 24 + 24 then 80640 else 0) = histFn unit3945F r) :
    ∀ r, (if r = 24 + 24 then 0 else hb r) + (if r = 24 then hb (24 + 24) else 0)
      + (if r = 24 then hb (24 + 24) else 0) = histFn cut3945F r := by
  have hb' : ∀ r, hb r = histFn body3945F r := by
    intro r
    have h1 := h r
    have h2 := hbody3945F r
    omega
  intro r
  simp only [hb']
  exact (hsplit3945F r).symm

/-- Unit -> cut in one identity: the 137,952 rank-48 blocks go (80,640 removed, 57,312 split), 114,624 rank-24 come. -/
theorem hcut3945F : ∀ r, histFn cut3945F r + (if r = 24 + 24 then 80640 + 57312 else 0)
    = histFn unit3945F r + (if r = 24 then 2 * 57312 else 0) := by
  intro r
  by_cases hr : r ∈ keys3945F
  · have h : ∀ r ∈ keys3945F, histFn cut3945F r + (if r = 24 + 24 then 80640 + 57312 else 0)
        = histFn unit3945F r + (if r = 24 then 2 * 57312 else 0) := by
      decide
    exact h r hr
  · have h48 : r ≠ 24 + 24 := fun h => hr (by subst h; decide)
    have h24 : r ≠ 24 := fun h => hr (by subst h; decide)
    rw [histFn_eq_zero_of_keys cut3945F_keys hr, histFn_eq_zero_of_keys unit3945F_keys hr, if_neg h48, if_neg h24]

/-- The rank mass of one unit drops by exactly 48 x 80,640 (the heads; the split keeps the mass). -/
theorem mass3945F : (cut3945F.map fun b => b.1 * b.2).sum + 48 * 80640 = (unit3945F.map fun b => b.1 * b.2).sum := by
  decide

/-- The number of paid blocks: 11,641,032 = #3945's 11,664,360 - 80,640 heads + 57,312 extra halves. -/
theorem count3945F : (cut3945F.map Prod.snd).sum + 80640 = (unit3945F.map Prod.snd).sum + 57312 := by decide

/-- Separation of the full cut: every paid rank is below the cut target 25. -/
theorem cut3945F_sep : ∀ b ∈ cut3945F, 1 ≤ b.1 ∧ b.1 < 25 := by decide

/-- No paid rank-48 block is left after the full cut. -/
theorem cut3945F_no48 : histFn cut3945F 48 = 0 := by decide

theorem cut3945F_24 : histFn cut3945F 24 = 114624 := by decide

/-- Negative control: one literal head more (80,641) breaks the body identity at rank 48. -/
theorem hbody3945F_sharp : histFn body3945F 48 + 80641 ≠ histFn unit3945F 48 := by decide

/-- Negative control: the unsplit body is not the cut profile (rank 48). -/
theorem body3945F_ne_cut : histFn body3945F 48 ≠ histFn cut3945F 48 := by decide

/-- Negative control: the half cut's 80,640 rank-24 remnants are not the full cut's rank-24 count. -/
theorem half_ne_full3945F : histFn cut3945F 24 ≠ 80640 := by decide

end OAI.PowerSaving.MX

#print axioms OAI.PowerSaving.MX.hbody3945F
#print axioms OAI.PowerSaving.MX.hsplit3945F
#print axioms OAI.PowerSaving.MX.hsplit3945F_of
#print axioms OAI.PowerSaving.MX.hcut3945F
#print axioms OAI.PowerSaving.MX.mass3945F
#print axioms OAI.PowerSaving.MX.count3945F
#print axioms OAI.PowerSaving.MX.cut3945F_sep
#print axioms OAI.PowerSaving.MX.cut3945F_no48
#print axioms OAI.PowerSaving.MX.hbody3945F_sharp
#print axioms OAI.PowerSaving.MX.body3945F_ne_cut
#print axioms OAI.PowerSaving.MX.half_ne_full3945F
