import Work.XHost.HSched

/-!
# (key: bridge-host) A host schedule as data, checked by the kernel (`hostSchedK`)

Data (per certificate; `nT = 112` classes, host index `i = 224 c + 112 s + t` for window kind `c`, side `s`, class `t`):
* `hr`: 9 rows (stage `l`), entry `q` = the host index of slot `q` when it is hosted (`< 2016`), else `2016 + q`;
* `br`: 72 rows (`k = 9 ρ + l`), entry `q` = the bank of slot `q` (`< nb`) when it is NOT hosted, else `nb + q`;
* `pad` (per bank, `< 73`), `rem` (per host, `≤` the rank of its idle climb), `res` (per slot, the residual `9 - e_q`).
`hostSchedK` checks the lengths, entries distinct in every row (a `Nat` bitset, as Shea's / Mamano's `schedK`), host
row `k % 9` and bank row `k` agree on which slots are hosted, every host lies in its window (`inWin`), `pad < 73`,
`rem ≤ rank`, `72 * Σ res < 2^20`, and the two packed sums (digit-wise: every bank's load is `73 - pad`, every host's load
is `rank - rem`).  `HostSched.ofTab` and `HostSched.ofTab_fitsR`.  The list lemmas of the first section are copies of
those of Work/Bridge401/SchedK.lean (Shea's complex588 tree, "DFT New Sums").
Written by h-L3-build (wht43, 13 Oct 2026), prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace XHost
open Finset

/-! ## List lemmas (copies of Work/Bridge401/SchedK.lean) -/

/-- one row: every entry `< m`, no entry twice (`seen` = the bitset of the entries so far). -/
def hkInj (m : ℕ) : ℕ → List ℕ → Bool
  | _, [] => true
  | seen, x :: xs => (decide (x < m) && !(seen.testBit x)) && hkInj m (seen ||| (1 <<< x)) xs

/-- the entry of slot `q` in row `k`. -/
def tabG (rows : List (List ℕ)) (k q : ℕ) : ℕ := (rows.getD k []).getD q 0

lemma hkInj_sound (m : ℕ) : ∀ (seen : ℕ) (xs : List ℕ), hkInj m seen xs = true →
    (∀ x ∈ xs, x < m ∧ seen.testBit x = false) ∧ xs.Nodup := by
  intro seen xs
  induction xs generalizing seen with
  | nil => intro _; simp
  | cons x xs ih =>
    intro h
    simp only [hkInj, Bool.and_eq_true, decide_eq_true_eq, Bool.not_eq_true'] at h
    obtain ⟨⟨hx, hs⟩, hr⟩ := h
    obtain ⟨h1, h2⟩ := ih _ hr
    have hb : ∀ y, (seen ||| (1 <<< x)).testBit y = (seen.testBit y || decide (x = y)) := by
      intro y
      rw [Nat.testBit_or, Nat.one_shiftLeft, Nat.testBit_two_pow]
    refine ⟨fun y hy => ?_, List.nodup_cons.mpr ⟨fun hm => ?_, h2⟩⟩
    · rcases List.mem_cons.mp hy with rfl | hy
      · exact ⟨hx, hs⟩
      · obtain ⟨hy1, hy2⟩ := h1 y hy
        rw [hb, Bool.or_eq_false_iff] at hy2
        exact ⟨hy1, hy2.1⟩
    · have := (h1 x hm).2
      rw [hb] at this
      simp at this

lemma map_sum_eqH {α M : Type} [AddCommMonoid M] (f : α → M) (d : α) : ∀ l : List α,
    (l.map f).sum = ∑ k ∈ range l.length, f (l.getD k d) := by
  intro l
  induction l with
  | nil => simp
  | cons x xs ih =>
    rw [List.map_cons, List.sum_cons, ih, List.length_cons, Finset.sum_range_succ']
    simp only [List.getD_cons_succ, List.getD_cons_zero]
    exact add_comm _ _

lemma digits_boundH (B : ℕ) : ∀ (n : ℕ) (a : ℕ → ℕ), (∀ b < n, a b < B) →
    ∑ b ∈ range n, a b * B ^ b < B ^ n := by
  intro n
  induction n with
  | zero => intro a _; simp
  | succ n ih =>
    intro a ha
    rw [Finset.sum_range_succ, pow_succ]
    have h1 := ih a (fun b hb => ha b (by omega))
    have h2 : a n + 1 ≤ B := ha n (by omega)
    have h3 : (a n + 1) * B ^ n ≤ B * B ^ n := Nat.mul_le_mul_right _ h2
    nlinarith

lemma digits_uniqueH (B : ℕ) (hB : 0 < B) : ∀ (n : ℕ) (a c : ℕ → ℕ), (∀ b < n, a b < B) → (∀ b < n, c b < B) →
    ∑ b ∈ range n, a b * B ^ b = ∑ b ∈ range n, c b * B ^ b → ∀ b < n, a b = c b := by
  intro n
  induction n with
  | zero => intro a c _ _ _ b hb; omega
  | succ n ih =>
    intro a c ha hc h b hb
    rw [Finset.sum_range_succ, Finset.sum_range_succ] at h
    have hp : 0 < B ^ n := pow_pos hB n
    have ha' := digits_boundH B n a (fun b hb => ha b (by omega))
    have hc' := digits_boundH B n c (fun b hb => hc b (by omega))
    have hd : a n = c n := by
      have e1 : (∑ b ∈ range n, a b * B ^ b + a n * B ^ n) / B ^ n
          = (∑ b ∈ range n, c b * B ^ b + c n * B ^ n) / B ^ n := by rw [h]
      rw [Nat.add_mul_div_right _ _ hp, Nat.add_mul_div_right _ _ hp, Nat.div_eq_of_lt ha',
        Nat.div_eq_of_lt hc'] at e1
      omega
    have hs : ∑ b ∈ range n, a b * B ^ b = ∑ b ∈ range n, c b * B ^ b := by
      rw [hd] at h
      omega
    rcases Nat.lt_succ_iff_lt_or_eq.mp hb with hb | rfl
    · exact ih a c (fun b hb => ha b (by omega)) (fun b hb => hc b (by omega)) hs b hb
    · exact hd

lemma nodup_getD_injH {l : List ℕ} (hl : l.Nodup) {i j : ℕ} (hi : i < l.length) (hj : j < l.length)
    (h : l.getD i 0 = l.getD j 0) : i = j := by
  rw [List.getD_eq_getElem _ _ hi, List.getD_eq_getElem _ _ hj] at h
  exact hl.getElem_inj_iff.mp h

/-! ## Packed sums with a threshold -/

/-- the packed residuals of one row; only entries `< m` count (digit = the entry). -/
def hkRowF (m : ℕ) : List ℕ → List ℕ → ℕ
  | x :: xs, r :: rs => (if x < m then r <<< (20 * x) else 0) + hkRowF m xs rs
  | _, _ => 0

/-- the packed targets `top i - p_i`. -/
def hkTgtF (top : ℕ → ℕ) : ℕ → List ℕ → ℕ
  | _, [] => 0
  | i, p :: ps => ((top i - p) <<< (20 * i)) + hkTgtF top (i + 1) ps

lemma hkRowF_eq (m : ℕ) : ∀ (row res : List ℕ), row.length = res.length →
    hkRowF m row res = ∑ q ∈ range row.length,
      (if row.getD q 0 < m then res.getD q 0 * (2 ^ 20) ^ row.getD q 0 else 0) := by
  intro row
  induction row with
  | nil => intro res _; simp [hkRowF]
  | cons x xs ih =>
    intro res h
    cases res with
    | nil => simp at h
    | cons r rs =>
      simp only [List.length_cons, Nat.add_right_cancel_iff] at h
      rw [hkRowF, ih rs h, List.length_cons, Finset.sum_range_succ', Nat.shiftLeft_eq]
      simp only [List.getD_cons_succ, List.getD_cons_zero]
      rw [← pow_mul]
      ring

lemma hkTgtF_eq (top : ℕ → ℕ) : ∀ (i : ℕ) (ps : List ℕ),
    hkTgtF top i ps = ∑ b ∈ range ps.length, (top (i + b) - ps.getD b 0) * (2 ^ 20) ^ (i + b) := by
  intro i ps
  induction ps generalizing i with
  | nil => simp [hkTgtF]
  | cons p ps ih =>
    rw [hkTgtF, ih, List.length_cons, Finset.sum_range_succ', Nat.shiftLeft_eq]
    simp only [List.getD_cons_succ, List.getD_cons_zero, Nat.add_zero]
    have e : ∀ b, i + 1 + b = i + (b + 1) := fun b => by ring
    simp only [e]
    rw [← pow_mul]
    ring

lemma regroupF (n m : ℕ) (f g : ℕ → ℕ) (B : ℕ) :
    ∑ q ∈ range n, (if g q < m then f q * B ^ g q else 0)
      = ∑ b ∈ range m, (∑ q ∈ range n, if g q = b then f q else 0) * B ^ b := by
  simp_rw [Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun q _ => ?_)
  simp_rw [ite_mul, zero_mul]
  by_cases hg : g q < m
  · rw [if_pos hg, Finset.sum_ite_eq (range m) (g q) (fun b => f q * B ^ b),
      if_pos (Finset.mem_range.mpr hg)]
  · rw [if_neg hg]
    symm
    refine Finset.sum_eq_zero (fun b hb => ?_)
    rw [if_neg]
    intro e
    exact hg (e ▸ Finset.mem_range.mp hb)

/-- rows agree on which slots are hosted: host entry `< 2016` iff bank entry `≥ nb`. -/
def hkCons (nb : ℕ) : List ℕ → List ℕ → Bool
  | h :: hs, b :: bs => (decide (h < 2016) == !decide (b < nb)) && hkCons nb hs bs
  | [], [] => true
  | _, _ => false

lemma hkCons_sound (nb : ℕ) : ∀ (hs bs : List ℕ), hkCons nb hs bs = true →
    ∀ q, (hs.getD q 0 < 2016 ↔ ¬ bs.getD q 0 < nb) ∨ (hs.length ≤ q) := by
  intro hs
  induction hs with
  | nil => intro bs _ q; right; simp
  | cons h hs ih =>
    intro bs hc q
    cases bs with
    | nil => simp [hkCons] at hc
    | cons b bs =>
      simp only [hkCons, Bool.and_eq_true, beq_iff_eq] at hc
      obtain ⟨h1, h2⟩ := hc
      cases q with
      | zero =>
        left
        simp only [List.getD_cons_zero]
        by_cases e1 : h < 2016 <;> by_cases e2 : b < nb <;> simp_all
      | succ q =>
        rcases ih bs h2 q with h3 | h3
        · left
          simpa using h3
        · right
          simp only [List.length_cons]
          omega

/-- every host entry of row `l` lies in its window. -/
def hkWin (l : ℕ) (row : List ℕ) : Bool := row.all (fun e => !decide (e < 2016) || inWin (e / 224) l)

/-- `rem` fits under the ranks. -/
def hsRem : ℕ → List ℕ → Bool
  | _, [] => true
  | i, x :: xs => decide (x ≤ rankOf (i / 224)) && hsRem (i + 1) xs

lemma hsRem_sound : ∀ (i : ℕ) (xs : List ℕ), hsRem i xs = true →
    ∀ b < xs.length, xs.getD b 0 ≤ rankOf ((i + b) / 224) := by
  intro i xs
  induction xs generalizing i with
  | nil => intro _ b hb; simp at hb
  | cons x xs ih =>
    intro h b hb
    simp only [hsRem, Bool.and_eq_true, decide_eq_true_eq] at h
    cases b with
    | zero => simpa using h.1
    | succ b =>
      have := ih (i + 1) h.2 b (by simpa using hb)
      simp only [List.getD_cons_succ]
      have e : i + 1 + b = i + (b + 1) := by ring
      rw [e] at this
      exact this

/-! ## The check -/

/-- **The host-schedule check.** -/
def hostSchedK (R nb : ℕ) (hr br : List (List ℕ)) (pad rem res : List ℕ) : Bool :=
  decide (0 < nb) && (hr.length == 9) && hr.all (fun row => (row.length == R) && hkInj (2016 + R) 0 row) &&
  (br.length == 72) && br.all (fun row => (row.length == R) && hkInj (nb + R) 0 row) &&
  (List.range 72).all (fun k => hkCons nb (hr.getD (k % 9) []) (br.getD k [])) &&
  (List.range 9).all (fun l => hkWin l (hr.getD l [])) &&
  (pad.length == nb) && pad.all (fun p => decide (p < 73)) &&
  (rem.length == 2016) && hsRem 0 rem &&
  (res.length == R) && decide (72 * res.sum < 2 ^ 20) &&
  ((br.map (fun row => hkRowF nb row res)).sum == hkTgtF (fun _ => 73) 0 pad) &&
  ((hr.map (fun row => hkRowF 2016 row res)).sum == hkTgtF (fun i => rankOf (i / 224)) 0 rem)

section Sound
variable {R nb : ℕ} {hr br : List (List ℕ)} {pad rem res : List ℕ}

lemma hostSchedK_parts (h : hostSchedK R nb hr br pad rem res = true) :
    0 < nb ∧ hr.length = 9 ∧ (∀ row ∈ hr, row.length = R ∧ hkInj (2016 + R) 0 row = true) ∧
    br.length = 72 ∧ (∀ row ∈ br, row.length = R ∧ hkInj (nb + R) 0 row = true) ∧
    (∀ k < 72, hkCons nb (hr.getD (k % 9) []) (br.getD k []) = true) ∧
    (∀ l < 9, hkWin l (hr.getD l []) = true) ∧
    pad.length = nb ∧ (∀ p ∈ pad, p < 73) ∧ rem.length = 2016 ∧ hsRem 0 rem = true ∧
    res.length = R ∧ 72 * res.sum < 2 ^ 20 ∧
    (br.map (fun row => hkRowF nb row res)).sum = hkTgtF (fun _ => 73) 0 pad ∧
    (hr.map (fun row => hkRowF 2016 row res)).sum = hkTgtF (fun i => rankOf (i / 224)) 0 rem := by
  simp only [hostSchedK, Bool.and_eq_true, beq_iff_eq, List.all_eq_true, decide_eq_true_eq,
    List.mem_range] at h
  obtain ⟨⟨⟨⟨⟨⟨⟨⟨⟨⟨⟨⟨⟨⟨h1, h2⟩, h3⟩, h4⟩, h5⟩, h6⟩, h7⟩, h8⟩, h9⟩, h10⟩, h11⟩, h12⟩, h13⟩, h14⟩, h15⟩ := h
  exact ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15⟩

lemma hrow_mem (h : hostSchedK R nb hr br pad rem res = true) {l : ℕ} (hl : l < 9) : hr.getD l [] ∈ hr := by
  rw [List.getD_eq_getElem _ _ (by rw [(hostSchedK_parts h).2.1]; exact hl)]
  exact List.getElem_mem _

lemma brow_mem (h : hostSchedK R nb hr br pad rem res = true) {k : ℕ} (hk : k < 72) : br.getD k [] ∈ br := by
  rw [List.getD_eq_getElem _ _ (by rw [(hostSchedK_parts h).2.2.2.1]; exact hk)]
  exact List.getElem_mem _

lemma hrow_len (h : hostSchedK R nb hr br pad rem res = true) {l : ℕ} (hl : l < 9) : (hr.getD l []).length = R :=
  ((hostSchedK_parts h).2.2.1 _ (hrow_mem h hl)).1

lemma brow_len (h : hostSchedK R nb hr br pad rem res = true) {k : ℕ} (hk : k < 72) : (br.getD k []).length = R :=
  ((hostSchedK_parts h).2.2.2.2.1 _ (brow_mem h hk)).1

lemma hrow_nodup (h : hostSchedK R nb hr br pad rem res = true) {l : ℕ} (hl : l < 9) : (hr.getD l []).Nodup :=
  (hkInj_sound _ 0 _ ((hostSchedK_parts h).2.2.1 _ (hrow_mem h hl)).2).2

lemma brow_nodup (h : hostSchedK R nb hr br pad rem res = true) {k : ℕ} (hk : k < 72) : (br.getD k []).Nodup :=
  (hkInj_sound _ 0 _ ((hostSchedK_parts h).2.2.2.2.1 _ (brow_mem h hk)).2).2

lemma cons_at (h : hostSchedK R nb hr br pad rem res = true) {k q : ℕ} (hk : k < 72) (hq : q < R) :
    tabG hr (k % 9) q < 2016 ↔ ¬ tabG br k q < nb := by
  rcases hkCons_sound nb _ _ ((hostSchedK_parts h).2.2.2.2.2.1 k hk) q with h1 | h1
  · exact h1
  · have hl := hrow_len h (show k % 9 < 9 from Nat.mod_lt _ (by norm_num))
    omega

lemma win_at (h : hostSchedK R nb hr br pad rem res = true) {l q : ℕ} (hl : l < 9) (hq : q < R)
    (he : tabG hr l q < 2016) : inWin (tabG hr l q / 224) l = true := by
  have hw := (hostSchedK_parts h).2.2.2.2.2.2.1 l hl
  unfold hkWin at hw
  rw [List.all_eq_true] at hw
  have hm : tabG hr l q ∈ hr.getD l [] := by
    unfold tabG
    rw [List.getD_eq_getElem (hr.getD l []) 0 (by rw [hrow_len h hl]; exact hq)]
    exact List.getElem_mem _
  have := hw _ hm
  simp only [Bool.or_eq_true, Bool.not_eq_true', decide_eq_false_iff_not] at this
  rcases this with h1 | h1
  · exact absurd he h1
  · exact h1

end Sound

/-! ## Host indices -/

/-- decoding a host index `e < 2016`: window kind `e / 224`, side `e % 224 / 112`, class `e % 112`. -/
def decH (e : ℕ) : Fin 9 × (Bool × Fin 112) :=
  (⟨e / 224 % 9, Nat.mod_lt _ (by norm_num)⟩, (decide (e % 224 / 112 = 1), ⟨e % 112, Nat.mod_lt _ (by norm_num)⟩))

/-- the host index of `(c, s, t)`. -/
def encH (c : Fin 9) (j : Bool × Fin 112) : ℕ := c.val * 224 + (if j.1 then 112 else 0) + j.2.val

lemma encH_lt (c : Fin 9) (j : Bool × Fin 112) : encH c j < 2016 := by
  unfold encH
  have := c.isLt
  have := j.2.isLt
  split_ifs <;> omega

lemma encH_div (c : Fin 9) (j : Bool × Fin 112) : encH c j / 224 = c.val := by
  unfold encH
  have := j.2.isLt
  split_ifs <;> omega

lemma encH_decH {e : ℕ} (he : e < 2016) : encH (decH e).1 (decH e).2 = e := by
  unfold encH decH
  simp only
  by_cases h1 : e % 224 / 112 = 1
  · simp only [h1, decide_true, if_true]
    omega
  · simp only [h1, decide_false, Bool.false_eq_true, if_false]
    omega

lemma decH_encH (c : Fin 9) (j : Bool × Fin 112) : decH (encH c j) = (c, j) := by
  obtain ⟨s, t⟩ := j
  have hc := c.isLt
  have ht := t.isLt
  unfold decH encH
  cases s
  · simp only [Bool.false_eq_true, if_false]
    refine Prod.ext (Fin.ext ?_) (Prod.ext ?_ (Fin.ext ?_))
    · show (c.val * 224 + 0 + t.val) / 224 % 9 = c.val
      omega
    · show decide ((c.val * 224 + 0 + t.val) % 224 / 112 = 1) = false
      simp only [decide_eq_false_iff_not]
      omega
    · show (c.val * 224 + 0 + t.val) % 112 = t.val
      omega
  · simp only [if_true]
    refine Prod.ext (Fin.ext ?_) (Prod.ext ?_ (Fin.ext ?_))
    · show (c.val * 224 + 112 + t.val) / 224 % 9 = c.val
      omega
    · show decide ((c.val * 224 + 112 + t.val) % 224 / 112 = 1) = true
      simp only [decide_eq_true_eq]
      omega
    · show (c.val * 224 + 112 + t.val) % 112 = t.val
      omega

/-- the window of kind `c` is the stages `kz c ≤ l < kz c + kK c`. -/
lemma win_range : ∀ c l : Fin 9, inWin c.val l.val = true → kz c ≤ l.val ∧ l.val < kz c + kK c := by decide

lemma win_le : ∀ c : Fin 9, kz c + kK c ≤ 9 := by decide

/-! ## The host schedule of a checked table -/

/-- the host of slot `q` in stage `l`. -/
def hAtT (hr : List (List ℕ)) (l : Fin 9) (q : ℕ) : Option (Fin 9 × (Bool × Fin 112)) :=
  if tabG hr l.val q < 2016 then some (decH (tabG hr l.val q)) else none

lemma hAtT_some {hr : List (List ℕ)} {l : Fin 9} {q : ℕ} {c : Fin 9} {j : Bool × Fin 112} :
    hAtT hr l q = some (c, j) ↔ tabG hr l.val q = encH c j := by
  unfold hAtT
  constructor
  · intro h
    split_ifs at h with h1
    · have e := congrArg (fun p => encH p.1 p.2) (Option.some.inj h)
      simp only at e
      rw [encH_decH h1] at e
      exact e
  · intro h
    rw [h, if_pos (encH_lt c j), decH_encH]

lemma hAtT_none {hr : List (List ℕ)} {l : Fin 9} {q : ℕ} : hAtT hr l q = none ↔ ¬ tabG hr l.val q < 2016 := by
  unfold hAtT
  split_ifs with h1 <;> simp [h1]

/-- **A checked table as a host schedule** (`nT = 112`). -/
def HostSched.ofTab (R nb : ℕ) (hr br : List (List ℕ)) (pad rem res : List ℕ)
    (h : hostSchedK R nb hr br pad rem res = true) : HostSched (Fin R) (Fin 112) where
  nb := nb
  hostAt l q := hAtT hr l q.val
  hostInj l q q' j hq hq' := by
    obtain ⟨c, j⟩ := j
    have e1 := hAtT_some.mp hq
    have e2 := hAtT_some.mp hq'
    have hl := hrow_len h l.isLt
    exact Fin.ext (nodup_getD_injH (hrow_nodup h l.isLt) (by rw [hl]; exact q.isLt) (by rw [hl]; exact q'.isLt)
      (e1.trans e2.symm))
  hostWin l q c j hq := by
    have e1 := hAtT_some.mp hq
    have hw := win_at h l.isLt q.isLt (by rw [e1]; exact encH_lt c j)
    rw [e1, encH_div] at hw
    exact hw
  bk ρ l q := ⟨tabG br (9 * ρ.val + l.val) q.val % nb, Nat.mod_lt _ (hostSchedK_parts h).1⟩
  hinj ρ l q q' hq hq' he := by
    have hk : 9 * ρ.val + l.val < 72 := by have := ρ.isLt; have := l.isLt; omega
    have hm : (9 * ρ.val + l.val) % 9 = l.val := by have := l.isLt; omega
    have b1 : tabG br (9 * ρ.val + l.val) q.val < nb := by
      have := (cons_at h hk q.isLt).not
      rw [hm, not_not] at this
      exact this.mp (hAtT_none.mp hq)
    have b2 : tabG br (9 * ρ.val + l.val) q'.val < nb := by
      have := (cons_at h hk q'.isLt).not
      rw [hm, not_not] at this
      exact this.mp (hAtT_none.mp hq')
    have e := congrArg Fin.val he
    simp only [Nat.mod_eq_of_lt b1, Nat.mod_eq_of_lt b2] at e
    have hl := brow_len h hk
    exact Fin.ext (nodup_getD_injH (brow_nodup h hk) (by rw [hl]; exact q.isLt) (by rw [hl]; exact q'.isLt) e)
  pad b := pad.getD b.val 0
  hrem c j := rem.getD (encH c j) 0

section Fits
variable (R nb : ℕ) (hr br : List (List ℕ)) (pad rem res : List ℕ) (h : hostSchedK R nb hr br pad rem res = true)

/-- the load of a bank of a checked table, as a finite sum over rows and slots. -/
lemma ofTab_bload (b : Fin (HostSched.ofTab R nb hr br pad rem res h).nb) :
    LH (bUse (HostSched.ofTab R nb hr br pad rem res h).hostAt (HostSched.ofTab R nb hr br pad rem res h).bk)
        (fun q : Fin R => res.getD q.val 0) 72 b
      = ∑ k ∈ range 72, ∑ q ∈ range R, if tabG br k q = b.val then res.getD q 0 else 0 := by
  unfold LH rrH
  refine Finset.sum_congr rfl (fun k hk => ?_)
  have hk' : k < 72 := Finset.mem_range.mp hk
  rw [Finset.sum_filter]
  have ek : 9 * (k / 9 % 8) + k % 9 = k := by omega
  have e : ∀ q : Fin R, (bUse (HostSched.ofTab R nb hr br pad rem res h).hostAt
      (HostSched.ofTab R nb hr br pad rem res h).bk k q = some b) ↔ tabG br k q.val = b.val := by
    intro q
    have hc := cons_at h hk' q.isLt
    have hbl : b.val < nb := b.isLt
    unfold bUse
    split
    · rename_i hn
      have h1 : ¬ tabG hr (k % 9) q.val < 2016 := hAtT_none.mp hn
      have hb : tabG br k q.val < nb := by
        by_contra hb
        exact h1 (hc.mpr hb)
      constructor
      · intro e
        have e2 := congrArg Fin.val (Option.some.inj e)
        change tabG br (9 * (k / 9 % 8) + k % 9) q.val % nb = b.val at e2
        rw [ek, Nat.mod_eq_of_lt hb] at e2
        exact e2
      · intro e
        congr 1
        apply Fin.ext
        change tabG br (9 * (k / 9 % 8) + k % 9) q.val % nb = b.val
        rw [ek, Nat.mod_eq_of_lt hb]
        exact e
    · rename_i x hx
      have h1 : tabG hr (k % 9) q.val < 2016 := by
        by_contra h1
        have hn : hAtT hr ⟨k % 9, Nat.mod_lt _ (by norm_num)⟩ q.val = none := hAtT_none.mpr h1
        have hx' : hAtT hr ⟨k % 9, Nat.mod_lt _ (by norm_num)⟩ q.val = some x := hx
        rw [hn] at hx'
        exact absurd hx' (by simp)
      have hb : ¬ tabG br k q.val < nb := hc.mp h1
      constructor
      · intro e
        exact absurd e (by simp)
      · intro e
        exact absurd (by rw [e]; exact hbl) hb
  simp_rw [e]
  exact Fin.sum_univ_eq_sum_range (fun q => if tabG br k q = b.val then res.getD q 0 else 0) R

/-- the load of a host of a checked table, as a finite sum over the stages of its window. -/
lemma ofTab_hload (c : Fin 9) (j : Bool × Fin 112) :
    LH (hUse (HostSched.ofTab R nb hr br pad rem res h).hostAt c (kz c)) (fun q : Fin R => res.getD q.val 0) (kK c) j
      = ∑ k ∈ range (kK c), ∑ q ∈ range R, if tabG hr (kz c + k) q = encH c j then res.getD q 0 else 0 := by
  unfold LH rrH
  refine Finset.sum_congr rfl (fun k hk => ?_)
  have hk' : kz c + k < 9 := by have := win_le c; have := Finset.mem_range.mp hk; omega
  rw [Finset.sum_filter]
  have e : ∀ q : Fin R, (hUse (HostSched.ofTab R nb hr br pad rem res h).hostAt c (kz c) k q = some j) ↔
      tabG hr (kz c + k) q.val = encH c j := by
    intro q
    constructor
    · intro hq
      obtain ⟨hk2, e⟩ := hUse_to hq
      exact hAtT_some.mp e
    · intro hq
      exact hUse_of hk' (hAtT_some.mpr hq)
  simp_rw [e]
  exact Fin.sum_univ_eq_sum_range (fun q => if tabG hr (kz c + k) q = encH c j then res.getD q 0 else 0) R

include h in
/-- the window carries the whole load of a host. -/
lemma hload_window (c : Fin 9) (i : ℕ) (hi : i / 224 = c.val) (hi2 : i < 2016) :
    ∑ k ∈ range (kK c), ∑ q ∈ range R, (if tabG hr (kz c + k) q = i then res.getD q 0 else 0)
      = ∑ l ∈ range 9, ∑ q ∈ range R, (if tabG hr l q = i then res.getD q 0 else 0) := by
  rw [show (∑ k ∈ range (kK c), ∑ q ∈ range R, (if tabG hr (kz c + k) q = i then res.getD q 0 else 0))
      = ∑ l ∈ (range (kK c)).map (addLeftEmbedding (kz c)), ∑ q ∈ range R,
          (if tabG hr l q = i then res.getD q 0 else 0) from by rw [Finset.sum_map]; rfl]
  refine Finset.sum_subset (fun x hx => ?_) (fun x hx hn => ?_)
  · obtain ⟨a, ha, rfl⟩ := Finset.mem_map.mp hx
    have := win_le c
    have := Finset.mem_range.mp ha
    exact Finset.mem_range.mpr (by simp only [addLeftEmbedding_apply]; omega)
  · refine Finset.sum_eq_zero (fun q hq => ?_)
    rw [if_neg]
    intro e
    have hx' := Finset.mem_range.mp hx
    have hq' := Finset.mem_range.mp hq
    have hw := win_at h hx' hq' (by rw [e]; exact hi2)
    rw [e, hi] at hw
    obtain ⟨w1, w2⟩ := win_range c ⟨x, hx'⟩ hw
    exact hn (Finset.mem_map.mpr ⟨x - kz c, Finset.mem_range.mpr (by simp only at w1 w2; omega),
      by simp only [addLeftEmbedding_apply]; simp only at w1; omega⟩)

/-- **A checked table fits the residuals `res`.** -/
theorem HostSched.ofTab_fitsR : (HostSched.ofTab R nb hr br pad rem res h).fitsR (fun q => res.getD q.val 0) := by
  obtain ⟨hnb, hhl, hhrows, hbl, hbrows, hcons, hwin, hpl, hpd, hrl, hrem, hresl, hsum, hbeq, hheq⟩ :=
    hostSchedK_parts h
  have htot : ∀ (rows : List (List ℕ)) (n : ℕ), rows.length = n → n ≤ 72 → (∀ k < n, (rows.getD k []).length = R) →
      ∀ c, ∑ k ∈ range n, ∑ q ∈ range R, (if tabG rows k q = c then res.getD q 0 else 0) ≤ 72 * res.sum := by
    intro rows n hn hn72 hlen c
    have e1 : ∑ q ∈ range R, res.getD q 0 = res.sum := by
      rw [← hresl]
      have := map_sum_eqH (fun x : ℕ => x) 0 res
      simpa using this.symm
    calc ∑ k ∈ range n, ∑ q ∈ range R, (if tabG rows k q = c then res.getD q 0 else 0)
        ≤ ∑ k ∈ range n, ∑ q ∈ range R, res.getD q 0 := by
          refine Finset.sum_le_sum (fun k _ => Finset.sum_le_sum (fun q _ => ?_))
          split_ifs <;> omega
      _ = n * res.sum := by rw [e1, Finset.sum_const, Finset.card_range, smul_eq_mul]
      _ ≤ 72 * res.sum := Nat.mul_le_mul_right _ hn72
  -- the packed identities, digit by digit
  have hpacked : ∀ (rows : List (List ℕ)) (n m : ℕ), rows.length = n → (∀ k < n, (rows.getD k []).length = R) →
      (rows.map (fun row => hkRowF m row res)).sum
        = ∑ c ∈ range m, (∑ k ∈ range n, ∑ q ∈ range R, if tabG rows k q = c then res.getD q 0 else 0)
            * (2 ^ 20) ^ c := by
    intro rows n m hn hlen
    rw [map_sum_eqH (fun row => hkRowF m row res) [] rows, hn]
    have e2 : ∀ k ∈ range n, hkRowF m (rows.getD k []) res
        = ∑ c ∈ range m, (∑ q ∈ range R, if tabG rows k q = c then res.getD q 0 else 0) * (2 ^ 20) ^ c := by
      intro k hk
      have hl := hlen k (Finset.mem_range.mp hk)
      rw [hkRowF_eq _ _ _ (by rw [hl, hresl]), hl]
      exact regroupF R m (fun q => res.getD q 0) (fun q => tabG rows k q) (2 ^ 20)
    rw [Finset.sum_congr rfl e2, Finset.sum_comm]
    refine Finset.sum_congr rfl (fun c _ => ?_)
    rw [Finset.sum_mul]
  have hpad : ∀ b : Fin nb, pad.getD b.val 0 < 73 := fun b =>
    hpd _ (by rw [List.getD_eq_getElem _ _ (by omega)]; exact List.getElem_mem _)
  refine ⟨fun b => ?_, fun b => hpad b, fun c j => ?_⟩
  · -- banks
    show LH _ _ 72 b + pad.getD b.val 0 = 73
    rw [ofTab_bload R nb hr br pad rem res h b]
    have hT : hkTgtF (fun _ => 73) 0 pad = ∑ c ∈ range nb, (73 - pad.getD c 0) * (2 ^ 20) ^ c := by
      rw [hkTgtF_eq, hpl]
      simp only [Nat.zero_add]
    have hd := digits_uniqueH (2 ^ 20) (by norm_num) nb
      (fun c => ∑ k ∈ range 72, ∑ q ∈ range R, if tabG br k q = c then res.getD q 0 else 0)
      (fun c => 73 - pad.getD c 0)
      (fun c _ => lt_of_le_of_lt (htot br 72 hbl le_rfl (fun k hk => brow_len h hk) c) hsum)
      (fun c _ => by show 73 - pad.getD c 0 < 2 ^ 20; omega)
      (by rw [← hpacked br 72 nb hbl (fun k hk => brow_len h hk), ← hT]; exact hbeq) b.val b.isLt
    try simp only at hd
    rw [hd]
    have : pad.getD b.val 0 < 73 := hpad ⟨b.val, b.isLt⟩
    omega
  · -- hosts
    show LH _ _ (kK c) j + rem.getD (encH c j) 0 = rankOf c.val
    rw [ofTab_hload R nb hr br pad rem res h c j,
      hload_window R nb hr br pad rem res h c (encH c j) (encH_div c j) (encH_lt c j)]
    have hT : hkTgtF (fun i => rankOf (i / 224)) 0 rem
        = ∑ i ∈ range 2016, (rankOf (i / 224) - rem.getD i 0) * (2 ^ 20) ^ i := by
      rw [hkTgtF_eq, hrl]
      simp only [Nat.zero_add]
    have hrk : ∀ i, rankOf i ≤ 48 := by
      intro i
      unfold rankOf
      split <;> omega
    have hd := digits_uniqueH (2 ^ 20) (by norm_num) 2016
      (fun i => ∑ l ∈ range 9, ∑ q ∈ range R, if tabG hr l q = i then res.getD q 0 else 0)
      (fun i => rankOf (i / 224) - rem.getD i 0)
      (fun i _ => lt_of_le_of_lt (htot hr 9 hhl (by norm_num) (fun l hl => hrow_len h hl) i) hsum)
      (fun i _ => by have := hrk (i / 224); show rankOf (i / 224) - rem.getD i 0 < 2 ^ 20; omega)
      (by rw [← hpacked hr 9 2016 hhl (fun l hl => hrow_len h hl), ← hT]; exact hheq)
      (encH c j) (encH_lt c j)
    try simp only at hd
    rw [hd, encH_div]
    have := hsRem_sound 0 rem hrem (encH c j) (by rw [hrl]; exact encH_lt c j)
    rw [Nat.zero_add, encH_div] at this
    omega

end Fits

/-- the pad price of a checked table (as `rcost`-style list sum). -/
lemma HostSched.ofTab_pad (R nb : ℕ) (hr br : List (List ℕ)) (pad rem res : List ℕ)
    (h : hostSchedK R nb hr br pad rem res = true) (b : Fin nb) :
    (HostSched.ofTab R nb hr br pad rem res h).pad b = pad.getD b.val 0 := rfl

lemma HostSched.ofTab_hrem (R nb : ℕ) (hr br : List (List ℕ)) (pad rem res : List ℕ)
    (h : hostSchedK R nb hr br pad rem res = true) (c : Fin 9) (j : Bool × Fin 112) :
    (HostSched.ofTab R nb hr br pad rem res h).hrem c j = rem.getD (encH c j) 0 := rfl

end XHost
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHost.HostSched.ofTab_fitsR
