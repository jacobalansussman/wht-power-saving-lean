import Work.XMut.Glue
import Work.KLow.GlueX
import Work.KMX.ZCheckZ
import Work.KMX.LabelsXZ

/-!
# (key: k-combine) Glue: a checked EXIT-FRAME certificate with MUTABLE sources is an exit-frame invocation package

k-combine, 13 Oct 2026.  For ANY `c : GXD.Raw` (nothing is evaluated here):
* `exit_of_zchecksZ`  k-lower's `exit_of_zchecks` (Work/KLow/GlueX.lean) with the mutable-source zero check
                      `GS.zCheckZ` (Work/KMX/ZCheckZ.lean) and the rows of the exit batches extended by the x registers
                      (`rowsXS`), so that it also gives the third exit condition `ExitColsX` (no x row of `Mt` reads an
                      exit helper); the y fillers of a batch are y registers (`v ≤ x < 2 v`);
* `exit_of_twinZ`     the three conditions of `c` from those of k-lower's y-free twin `nyRaw c` (`nyRaw_Mt`);
* `glabXOfZ`          `XEX.GLabXZ` (Work/KMX/HalvesZ.lean) over x-mutable's gauged mutable-source scalar half
                      `gscalOfZ` (XMut/Glue.lean), from k-lower's label check `GLab.ValidX`;
* `ginvXOfZ`          the exit-frame package `XEX.GInvX`; `ginvXOfZ_cost` (Shea's invoice, the tail climbs in `rB`,
                      `rT = []`), `ginvXOfZ_res` (bank residuals as k-lower's `ginvXOf_res`).
Credits: k-lower (exit conditions, zero checks, ValidX, ginvXOf), x-mutable (ScalOKZ, eCheckZ, gscalOfZ, ginvOfZ),
x-exit-lean (GLabX, GInvX), Ryan Shea / DFT Bridge3 / Nil Mamano (ginvOf, the gauged route), Jacob Sussman (engine).
No `sorry`.
-/

set_option linter.unusedVariables false
set_option linter.unusedSectionVars false
set_option linter.deprecated false

namespace GXD
open SSC OAI.PowerSaving OAI.PowerSaving.CR OAI.PowerSaving.CB Matrix Finset

/-- the x registers -/
def xsL (c : Raw) : List Nat := List.range c.v

/-- the rows of an exit batch: the x registers and the slot registers -/
def rowsXS (c : Raw) : List Nat := xsL c ++ slotsL c

theorem mem_rowsXS_slot (c : Raw) {r : Nat} (h : r ∈ slotsL c) : r ∈ rowsXS c := List.mem_append_right _ h
theorem mem_rowsXS_x (c : Raw) (t : Fin c.v) : t.val ∈ rowsXS c :=
  List.mem_append_left _ (List.mem_range.mpr t.isLt)

section
variable (c : Raw) (hv : 0 < c.v) (el : List Nat)

/-- **The three exit conditions from the zero checks, sources mutable.** -/
theorem exit_of_zchecksZ
    (hy : ∀ (q : Fin c.R) (t : Fin c.v), GS.Mt c hv (lS q) (lY t) = 0)
    (lel hel : List Nat) (nb : Nat) (chunks : List (List Nat))
    (hlen : lel.length ≤ 2 ^ nb) (hnb : 1 ≤ nb)
    (hel2 : ∀ r ∈ el, 2 * c.v ≤ r)
    (hlel : ∀ r ∈ lel, GLab.memK r el = true) (hhel : ∀ r ∈ hel, r ∉ lel)
    (hsplit : ∀ r ∈ el, r ∈ lel ∨ r ∈ hel)
    (hbit : ∀ b, b < nb → ∀ val : Bool, ∃ p pre, (∀ x ∈ pre, c.v ≤ x ∧ x < 2 * c.v) ∧
      GS.zCheckZ c p (pre ++ bitB lel b val) (rowsXS c) = true)
    (hhv : ∀ e ∈ hel, ∃ p pre, (∀ x ∈ pre, c.v ≤ x ∧ x < 2 * c.v) ∧
      GS.zCheckZ c p (pre ++ [e]) (rowsXS c) = true)
    (hch : ∀ l ∈ chunks, ∃ p, GS.zCheckZ c p l el = true)
    (hcov : ∀ s, s < 2 * c.v + c.R → (s < c.v ∨ 2 * c.v ≤ s) → GLab.memK s el = false → ∃ l ∈ chunks, s ∈ l)
    (hdisj : ∀ l ∈ chunks, ∀ s ∈ l, GLab.memK s el = false) :
    XEX.ExitRows (GLab.exitSet c el) (GS.Mt c hv) ∧ XEX.ExitCols (GLab.exitSet c el) (GS.Mt c hv) ∧
      XEX.ExitColsX (GLab.exitSet c el) (GS.Mt c hv) := by
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
    refine GS.zCheckZ_sound hv hp i j hjl hi (fun h => ?_)
    have h1 := hdisj l hl _ h
    rw [memK_complete _ el hi] at h1
    exact Bool.noConfusion h1
  -- an exit batch with source `j` and a row `i` (an x register or a slot, not a source of the batch)
  have viaBatch : ∀ (i j : GS.Rl c) (pre l : List Nat) (p : GS.SPar), (∀ x ∈ pre, c.v ≤ x ∧ x < 2 * c.v) →
      GS.zCheckZ c p (pre ++ l) (rowsXS c) = true → GS.emb c j ∈ l → GS.emb c i ∈ rowsXS c →
      (GS.emb c i < c.v ∨ 2 * c.v ≤ GS.emb c i) → GS.emb c i ∉ l → GS.Mt c hv i j = 0 := by
    intro i j pre l p hpre hp hj hi hxs hni
    refine GS.zCheckZ_sound hv hp i j (List.mem_append_right _ hj) hi (fun h => ?_)
    rcases List.mem_append.mp h with h | h
    · have := hpre _ h; omega
    · exact hni h
  -- a source exit `s` and a row `i ≠ s` that is an x register or a slot, and not an exit unless `rows` separate
  have sepAny : ∀ (i j : GS.Rl c), GS.emb c j ∈ el → GS.emb c i ∈ rowsXS c →
      (GS.emb c i < c.v ∨ 2 * c.v ≤ GS.emb c i) → GS.emb c i ≠ GS.emb c j → GS.Mt c hv i j = 0 := by
    intro i j hj hi hxs hne
    rcases hsplit _ hj with hjl | hjh
    · by_cases hil : GS.emb c i ∈ lel
      · obtain ⟨b, hb, hjb, hib⟩ := bit_sep lel nb hlen hil hjl hne
        obtain ⟨p, pre, hpre, hp⟩ := hbit b hb _
        exact viaBatch i j pre _ p hpre hp hjb hi hxs hib
      · obtain ⟨p, pre, hpre, hp⟩ := hbit 0 hnb ((lel.idxOf (GS.emb c j)).testBit 0)
        refine viaBatch i j pre _ p hpre hp ?_ hi hxs (fun h => hil (List.mem_filter.mp h).1)
        unfold bitB; exact List.mem_filter.mpr ⟨hjl, by simp⟩
    · obtain ⟨p, pre, hpre, hp⟩ := hhv _ hjh
      exact viaBatch i j pre [GS.emb c j] p hpre hp (List.mem_singleton_self _) hi hxs
        (fun h => hne (List.mem_singleton.mp h))
  refine ⟨fun e he j hj => ?_, fun e he q hq => ?_, fun e he t => ?_⟩
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
        exact sepAny (lS e) (lS q) (GLab.memK_sound _ el h) (mem_rowsXS_slot c (mem_slotsL c e))
          (Or.inr (show 2 * c.v ≤ 2 * c.v + e.val by omega)) hne
  · have hr := hE e he
    have hqe : 2 * c.v + q.val ≠ 2 * c.v + e.val := fun h => hq (Fin.ext (by omega))
    exact sepAny (lS q) (lS e) hr (mem_rowsXS_slot c (mem_slotsL c q))
      (Or.inr (show 2 * c.v ≤ 2 * c.v + q.val by omega)) hqe
  · have hr := hE e he
    have hne : t.val ≠ 2 * c.v + e.val := by have := t.isLt; omega
    exact sepAny (lX t) (lS e) hr (mem_rowsXS_x c t) (Or.inl t.isLt) hne

/-- the three exit conditions of `c` from those of its y-free twin (exit rows, slot rows and x rows are not y rows) -/
theorem exit_of_twinZ
    (hB : ∀ x ∈ GS.addsB c, c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v)
    (h : XEX.ExitRows (GLab.exitSet (nyRaw c) el) (GS.Mt (nyRaw c) hv) ∧
      XEX.ExitCols (GLab.exitSet (nyRaw c) el) (GS.Mt (nyRaw c) hv) ∧
      XEX.ExitColsX (GLab.exitSet (nyRaw c) el) (GS.Mt (nyRaw c) hv)) :
    XEX.ExitRows (GLab.exitSet c el) (GS.Mt c hv) ∧ XEX.ExitCols (GLab.exitSet c el) (GS.Mt c hv) ∧
      XEX.ExitColsX (GLab.exitSet c el) (GS.Mt c hv) := by
  refine ⟨fun e he j hj => ?_, fun e he q hq => ?_, fun e he t => ?_⟩
  · rw [nyRaw_Mt c hv hB (lS e) j (by simp [iY, lS])]
    exact h.1 e he j hj
  · rw [nyRaw_Mt c hv hB (lS q) (lS e) (by simp [iY, lS])]
    exact h.2.1 e he q hq
  · rw [nyRaw_Mt c hv hB (lX t) (lS e) (by simp [iY, lX])]
    exact h.2.2 e he t

end

section
variable (c : Raw) (K : GS.ScalOKZ c)
    (es : List (GS.SPar × List Nat)) (he : ∀ e ∈ es, GS.eCheckZ c e.1 e.2 = true)
    (a : GLab.Aux) (el : List Nat) (V : GLab.ValidX c a (es.flatMap Prod.snd) el)
    (hR : XEX.ExitRows (GLab.exitSet c el) (GS.Mt c K.hv)) (hC : XEX.ExitCols (GLab.exitSet c el) (GS.Mt c K.hv))
    (hX : XEX.ExitColsX (GLab.exitSet c el) (GS.Mt c K.hv))

/-- the exit label half of a checked exit-frame certificate over its gauged mutable-source scalar half -/
noncomputable def glabXOfZ : XEX.GLabXZ (Fin a.p.h) (gscalOfZ c K es he) :=
  V.glabXZ (gscalOfZ c K es he) (roles c K.hv) GS.toCoef rfl rfl rfl rfl
    (reg c) (reg_slotZ c K es he) (reg_mem c)
    (fun q hq => (gscalOfZ_gauged c K es he q).mpr hq) hR hC hX

/-- **a checked exit-frame certificate with mutable sources is an exit-frame invocation package** -/
noncomputable def ginvXOfZ : XEX.GInvX (Fin a.p.h) (Fin c.v) (Fin c.R) (Fin c.ret.length) :=
  (glabXOfZ c K es he a el V hR hC hX).inv

/-- its price: Shea's invoice (the tail climbs are part of `rB`) -/
theorem ginvXOfZ_cost (φ : ℕ → ℝ) :
    (ginvXOfZ c K es he a el V hR hC hX).cost φ
      = rcost φ a.rA + (∑ k : Fin c.ret.length, bcost φ (GLab.cdim (a.p.tab.get (a.SA.get (reg c k)))))
        + rcost φ a.rB + rcost φ a.rF := by
  obtain ⟨h1, h2, h3, h4, h5⟩ := V.glabXZ_ranks (gscalOfZ c K es he) (roles c K.hv) GS.toCoef rfl rfl rfl rfl
    (reg c) (reg_slotZ c K es he) (reg_mem c)
    (fun q hq => (gscalOfZ_gauged c K es he q).mpr hq) hR hC hX
  show XEX.GLabXZ.cost (glabXOfZ c K es he a el V hR hC hX) φ = _
  unfold XEX.GLabXZ.cost
  unfold glabXOfZ
  rw [h1, h2, h3, h4, rcost_nil]
  simp only [h5, add_zero]

/-- its bank residuals -/
theorem ginvXOfZ_res (q : Fin c.R) :
    ((ginvXOfZ c K es he a el V hR hC hX).ex q).dim - ((ginvXOfZ c K es he a el V hR hC hX).ent q).dim
      = Bool.rec (motive := fun _ => Nat) (a.p.h - GLab.cdim (a.p.tab.get (a.S0.get (2 * c.v + q.val))))
          (GLab.cdim (a.p.tab.get (a.SF.get (2 * c.v + q.val)))
            - GLab.cdim (a.p.tab.get (a.S0.get (2 * c.v + q.val))))
          (GLab.memK (2 * c.v + q.val) el) :=
  V.glabXZ_res (gscalOfZ c K es he) (roles c K.hv) GS.toCoef rfl rfl rfl rfl
    (reg c) (reg_slotZ c K es he) (reg_mem c)
    (fun q hq => (gscalOfZ_gauged c K es he q).mpr hq) hR hC hX q

end

end GXD

#print axioms GXD.exit_of_zchecksZ
#print axioms GXD.exit_of_twinZ
#print axioms GXD.ginvXOfZ
#print axioms GXD.ginvXOfZ_cost
#print axioms GXD.ginvXOfZ_res
