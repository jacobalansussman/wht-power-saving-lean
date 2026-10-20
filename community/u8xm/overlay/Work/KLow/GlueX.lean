import Work.Gauge.Glue
import Work.KLow.LabelsX
import Work.KLow.ZCheck
import Work.KLow.NoY

/-!
# (key: k-lower) Glue: a checked EXIT-FRAME certificate is an exit-frame invocation package (`XEX.GInvX`)

k-lower, 13 Oct 2026.  For ANY `c : GXD.Raw` (nothing is evaluated here):
* `exit_of_zchecks`  the two exit conditions on the total matrix `GS.Mt c hv` (Work/KLow/ExitEnd) from zero checks
                     (`GS.zCheck`, Work/KLow/ZCheck): BIT batches (sources = the exit registers whose index in the exit
                     list has bit `b` equal to `val`, rows = every slot) and CHUNKS (sources = the other non-y registers,
                     rows = the exit list), with the coverage facts `hcov`, `hdisj` (decidable);
* `glabXOf`          x-exit-lean's `GLabX` over DFT Bridge3's gauged scalar half `gscalOf` (Shea's Work/Gauge/Glue),
                     from `GLab.ValidX` (Work/KLow/LabelsX);
* `ginvXOf`          the exit-frame package `XEX.GInvX` (`GLabX.inv`); `ginvXOf_cost` its price (Shea's invoice:
                     the tail climbs are in `rB`, `rT = []`); `ginvXOf_res` the bank residuals
                     (`dim exit - dim start` for an exit slot, `h - dim start` otherwise).
Base: Shea's / Bridge3's `GXD.ginvOf` (Work/Gauge/Glue.lean) on Sussman's `Work/GCert/Data/Glue.lean`
(Apache-2.0); exit frames: x-exit-lean (XExit/Halves, XExit/Pack); the hub's tail and exit frames (#591/#593/#599,
Kothe's #799 lineage).  No `sorry`.
-/

set_option linter.unusedVariables false
set_option linter.unusedSectionVars false
set_option linter.deprecated false

namespace GXD
open SSC OAI.PowerSaving OAI.PowerSaving.CR OAI.PowerSaving.CB Matrix Finset

/-- the bit batch `(b, val)` of the exit list -/
def bitB (el : List Nat) (b : Nat) (val : Bool) : List Nat := el.filter (fun e => (el.idxOf e).testBit b == val)

/-- the slot registers -/
def slotsL (c : Raw) : List Nat := List.range' (2 * c.v) c.R

theorem memK_complete (r : Nat) : ∀ gl : List Nat, r ∈ gl → GLab.memK r gl = true := by
  intro gl
  induction gl with
  | nil => intro h; exact absurd h List.not_mem_nil
  | cons x gl ih =>
    intro h
    show Bool.rec (motive := fun _ => Bool) (GLab.memK r gl) true (Nat.beq x r) = true
    cases hb : Nat.beq x r with
    | true => rfl
    | false =>
      have hne := Nat.ne_of_beq_eq_false hb
      rcases List.mem_cons.mp h with e | e
      · exact absurd e.symm hne
      · exact ih e

theorem bit_sep (el : List Nat) (nb : Nat) (hlen : el.length ≤ 2 ^ nb) {r s : Nat} (hr : r ∈ el) (hs : s ∈ el)
    (hne : r ≠ s) :
    ∃ b, b < nb ∧ s ∈ bitB el b ((el.idxOf s).testBit b) ∧ r ∉ bitB el b ((el.idxOf s).testBit b) := by
  have hi : el.idxOf r ≠ el.idxOf s := fun h => hne ((List.idxOf_inj hr).mp h)
  have hr' : el.idxOf r < 2 ^ nb := lt_of_lt_of_le (List.idxOf_lt_length_of_mem hr) hlen
  have hs' : el.idxOf s < 2 ^ nb := lt_of_lt_of_le (List.idxOf_lt_length_of_mem hs) hlen
  obtain ⟨b, hb⟩ : ∃ b, (el.idxOf r).testBit b ≠ (el.idxOf s).testBit b := by
    by_contra h
    push_neg at h
    exact hi (Nat.eq_of_testBit_eq h)
  have hbn : b < nb := by
    by_contra h
    push_neg at h
    apply hb
    rw [Nat.testBit_eq_false_of_lt (lt_of_lt_of_le hr' (Nat.pow_le_pow_right (by norm_num) h)),
      Nat.testBit_eq_false_of_lt (lt_of_lt_of_le hs' (Nat.pow_le_pow_right (by norm_num) h))]
  refine ⟨b, hbn, ?_, ?_⟩
  · unfold bitB
    exact List.mem_filter.mpr ⟨hs, by simp⟩
  · unfold bitB
    intro h
    have h2 := (List.mem_filter.mp h).2
    simp only [beq_iff_eq] at h2
    exact hb h2

theorem mem_slotsL (c : Raw) (q : Fin c.R) : 2 * c.v + q.val ∈ slotsL c :=
  List.mem_range'_1.mpr ⟨by omega, by have := q.isLt; omega⟩

section
variable (c : Raw) (hv : 0 < c.v) (el : List Nat)

/-- **The two exit conditions from the zero checks.**  Exit batches (rows = every slot): bit batches over the LIGHT
exit list `lel`, one batch per HEAVY exit (`hel`), each after a prefix `pre` of y registers (fillers that move the
digits of the sources; y registers are never rows here).  Chunks (rows = the exit list): the other non-y sources. -/
theorem exit_of_zchecks
    (hy : ∀ (q : Fin c.R) (t : Fin c.v), GS.Mt c hv (lS q) (lY t) = 0)
    (lel hel : List Nat) (nb : Nat) (chunks : List (List Nat))
    (hlen : lel.length ≤ 2 ^ nb) (hnb : 1 ≤ nb)
    (hel2 : ∀ r ∈ el, 2 * c.v ≤ r)
    (hlel : ∀ r ∈ lel, GLab.memK r el = true) (hhel : ∀ r ∈ hel, r ∉ lel)
    (hsplit : ∀ r ∈ el, r ∈ lel ∨ r ∈ hel)
    (hbit : ∀ b, b < nb → ∀ val : Bool, ∃ p pre, (∀ x ∈ pre, x < 2 * c.v) ∧
      GS.zCheck c p (pre ++ bitB lel b val) (slotsL c) = true)
    (hhv : ∀ e ∈ hel, ∃ p pre, (∀ x ∈ pre, x < 2 * c.v) ∧ GS.zCheck c p (pre ++ [e]) (slotsL c) = true)
    (hch : ∀ l ∈ chunks, ∃ p, GS.zCheck c p l el = true)
    (hcov : ∀ s, s < 2 * c.v + c.R → (s < c.v ∨ 2 * c.v ≤ s) → GLab.memK s el = false → ∃ l ∈ chunks, s ∈ l)
    (hdisj : ∀ l ∈ chunks, ∀ s ∈ l, GLab.memK s el = false) :
    XEX.ExitRows (GLab.exitSet c el) (GS.Mt c hv) ∧ XEX.ExitCols (GLab.exitSet c el) (GS.Mt c hv) := by
  have hE : ∀ e ∈ GLab.exitSet c el, 2 * c.v + e.val ∈ el := fun e he =>
    GLab.memK_sound _ el (Finset.mem_filter.mp he).2
  have slot_ge : ∀ i : GS.Rl c, GS.emb c i ∈ slotsL c → 2 * c.v ≤ GS.emb c i := fun i h =>
    (List.mem_range'_1.mp h).1
  -- a non-exit source register: one chunk
  have viaChunk : ∀ (i j : GS.Rl c), GS.emb c i ∈ el → GLab.memK (GS.emb c j) el = false →
      GS.emb c j < 2 * c.v + c.R → (GS.emb c j < c.v ∨ 2 * c.v ≤ GS.emb c j) → GS.Mt c hv i j = 0 := by
    intro i j hi hj hlt hny
    obtain ⟨l, hl, hjl⟩ := hcov _ hlt hny hj
    obtain ⟨p, hp⟩ := hch l hl
    refine GS.zCheck_sound hv hp i j hjl hi (fun h => ?_)
    have h1 := hdisj l hl _ h
    rw [memK_complete _ el hi] at h1
    exact Bool.noConfusion h1
  -- an exit batch with source `j` and row `i` (a slot, not a source of the batch)
  have viaBatch : ∀ (i j : GS.Rl c) (pre l : List Nat) (p : GS.SPar), (∀ x ∈ pre, x < 2 * c.v) →
      GS.zCheck c p (pre ++ l) (slotsL c) = true → GS.emb c j ∈ l → GS.emb c i ∈ slotsL c → GS.emb c i ∉ l →
      GS.Mt c hv i j = 0 := by
    intro i j pre l p hpre hp hj hi hni
    refine GS.zCheck_sound hv hp i j (List.mem_append_right _ hj) hi (fun h => ?_)
    rcases List.mem_append.mp h with h | h
    · have := hpre _ h; have := slot_ge i hi; omega
    · exact hni h
  -- a source exit `s` and a row slot `r ≠ s`
  have sep : ∀ (i j : GS.Rl c), GS.emb c j ∈ el → GS.emb c i ∈ slotsL c → GS.emb c i ≠ GS.emb c j →
      GS.Mt c hv i j = 0 := by
    intro i j hj hi hne
    rcases hsplit _ hj with hjl | hjh
    · by_cases hil : GS.emb c i ∈ lel
      · obtain ⟨b, hb, hjb, hib⟩ := bit_sep lel nb hlen hil hjl hne
        obtain ⟨p, pre, hpre, hp⟩ := hbit b hb _
        exact viaBatch i j pre _ p hpre hp hjb hi hib
      · obtain ⟨p, pre, hpre, hp⟩ := hbit 0 hnb ((lel.idxOf (GS.emb c j)).testBit 0)
        refine viaBatch i j pre _ p hpre hp ?_ hi (fun h => hil (List.mem_filter.mp h).1)
        unfold bitB; exact List.mem_filter.mpr ⟨hjl, by simp⟩
    · obtain ⟨p, pre, hpre, hp⟩ := hhv _ hjh
      exact viaBatch i j pre [GS.emb c j] p hpre hp (List.mem_singleton_self _) hi
        (fun h => hne (List.mem_singleton.mp h))
  refine ⟨fun e he j hj => ?_, fun e he q hq => ?_⟩
  · have hr := hE e he
    rcases j with ((t|t)|q)
    · have ht : GLab.memK t.val el = false := by
        cases h : GLab.memK t.val el with
        | false => rfl
        | true => have := hel2 _ (GLab.memK_sound _ el h); have := t.isLt; omega
      exact viaChunk (lS e) (lX t) hr ht (show t.val < 2 * c.v + c.R by have := t.isLt; omega)
        (Or.inl (show t.val < c.v from t.isLt))
    · exact hy e t
    · have hqe : q ≠ e := fun h => hj (by rw [h])
      cases h : GLab.memK (2 * c.v + q.val) el with
      | false =>
        exact viaChunk (lS e) (lS q) hr h (show 2 * c.v + q.val < 2 * c.v + c.R by have := q.isLt; omega)
          (Or.inr (show 2 * c.v ≤ 2 * c.v + q.val by omega))
      | true =>
        have hne : 2 * c.v + e.val ≠ 2 * c.v + q.val := fun h' => hqe (Fin.ext (by omega))
        exact sep (lS e) (lS q) (GLab.memK_sound _ el h) (mem_slotsL c e) hne
  · have hr := hE e he
    have hqe : 2 * c.v + q.val ≠ 2 * c.v + e.val := fun h => hq (Fin.ext (by omega))
    exact sep (lS q) (lS e) hr (mem_slotsL c q) hqe

end

section
variable (c : Raw) (hv : 0 < c.v)
    (hAi : GS.MAi c hv * GS.MA c hv = 1) (hBi : GS.MBi c hv * GS.MB c hv = 1)
    (hx : ∀ (t : Fin c.v) (j : GS.Rl c), (GS.MB c hv * scatM (GS.Jrm c * GS.Ccm c) * GS.MA c hv)
      (lX t) j = if lX t = j then 1 else 0)
    (hy : ∀ (i : GS.Rl c) (t : Fin c.v), (GS.MB c hv * scatM (GS.Jrm c * GS.Ccm c) * GS.MA c hv)
      i (lY t) = if i = lY t then 1 else 0)
    (hid : ∀ S t : Fin c.v, (GS.MB c hv * scatM (GS.Jrm c * GS.Ccm c) * GS.MA c hv)
      (lY S) (lX t) = if S = t then 1 else 0)
    (es : List (GS.SPar × List Nat)) (he : ∀ e ∈ es, GS.eCheck c e.1 e.2 = true)
    (a : GLab.Aux) (el : List Nat) (V : GLab.ValidX c a (es.flatMap Prod.snd) el)
    (hR : XEX.ExitRows (GLab.exitSet c el) (GS.Mt c hv)) (hC : XEX.ExitCols (GLab.exitSet c el) (GS.Mt c hv))

/-- the exit label half of a checked exit-frame certificate over its gauged scalar half -/
noncomputable def glabXOf : XEX.GLabX (Fin a.p.h) (gscalOf c hv hAi hBi hx hy hid es he) :=
  V.glabX (gscalOf c hv hAi hBi hx hy hid es he) (roles c hv) GS.toCoef
    (scalOf_MA c hv hAi hBi hx hy hid) (scalOf_MAi c hv hAi hBi hx hy hid)
    (scalOf_MB c hv hAi hBi hx hy hid) (scalOf_MBi c hv hAi hBi hx hy hid)
    (reg c) (reg_slot c hv hAi hBi hx hy hid) (reg_mem c)
    (fun q hq => (gscalOf_gauged c hv hAi hBi hx hy hid es he q).mpr hq) hR hC

/-- **a checked exit-frame certificate is an exit-frame invocation package** (`GLabX.inv`) -/
noncomputable def ginvXOf : XEX.GInvX (Fin a.p.h) (Fin c.v) (Fin c.R) (Fin c.ret.length) :=
  (glabXOf c hv hAi hBi hx hy hid es he a el V hR hC).inv

/-- its price: Shea's invoice (the tail climbs are part of `rB`) -/
theorem ginvXOf_cost (φ : ℕ → ℝ) :
    (ginvXOf c hv hAi hBi hx hy hid es he a el V hR hC).cost φ
      = rcost φ a.rA + (∑ k : Fin c.ret.length, bcost φ (GLab.cdim (a.p.tab.get (a.SA.get (reg c k)))))
        + rcost φ a.rB + rcost φ a.rF := by
  obtain ⟨h1, h2, h3, h4, h5⟩ := V.glabX_ranks (gscalOf c hv hAi hBi hx hy hid es he) (roles c hv) GS.toCoef
    (scalOf_MA c hv hAi hBi hx hy hid) (scalOf_MAi c hv hAi hBi hx hy hid)
    (scalOf_MB c hv hAi hBi hx hy hid) (scalOf_MBi c hv hAi hBi hx hy hid)
    (reg c) (reg_slot c hv hAi hBi hx hy hid) (reg_mem c)
    (fun q hq => (gscalOf_gauged c hv hAi hBi hx hy hid es he q).mpr hq) hR hC
  show XEX.GLabX.cost (glabXOf c hv hAi hBi hx hy hid es he a el V hR hC) φ = _
  unfold XEX.GLabX.cost
  unfold glabXOf
  rw [h1, h2, h3, h4, rcost_nil]
  simp only [h5, add_zero]

/-- its bank residuals -/
theorem ginvXOf_res (q : Fin c.R) :
    ((ginvXOf c hv hAi hBi hx hy hid es he a el V hR hC).ex q).dim
      - ((ginvXOf c hv hAi hBi hx hy hid es he a el V hR hC).ent q).dim
      = Bool.rec (motive := fun _ => Nat) (a.p.h - GLab.cdim (a.p.tab.get (a.S0.get (2 * c.v + q.val))))
          (GLab.cdim (a.p.tab.get (a.SF.get (2 * c.v + q.val)))
            - GLab.cdim (a.p.tab.get (a.S0.get (2 * c.v + q.val))))
          (GLab.memK (2 * c.v + q.val) el) :=
  V.glabX_res (gscalOf c hv hAi hBi hx hy hid es he) (roles c hv) GS.toCoef
    (scalOf_MA c hv hAi hBi hx hy hid) (scalOf_MAi c hv hAi hBi hx hy hid)
    (scalOf_MB c hv hAi hBi hx hy hid) (scalOf_MBi c hv hAi hBi hx hy hid)
    (reg c) (reg_slot c hv hAi hBi hx hy hid) (reg_mem c)
    (fun q hq => (gscalOf_gauged c hv hAi hBi hx hy hid es he q).mpr hq) hR hC q

end

end GXD

#print axioms GXD.exit_of_zchecks
#print axioms GXD.ginvXOf
#print axioms GXD.ginvXOf_cost
#print axioms GXD.ginvXOf_res
