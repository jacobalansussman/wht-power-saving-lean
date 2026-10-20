import Work.XHost.HStage

/-!
# (key: bridge-host) The certificate of the HOSTED nine-stage word

Copy of Work/Bridge401/GNetCert.lean (Ryan Shea's complex588 tree; written by Nil Mamano's gen401.py) for `GHost401`:
the scalar map of a run is Shea's `chainQ` of the bank map `bk` (a hosted stage has the scalar map of a banked one), so
the sweep, the exchange and the scalar identity are Shea's, word for word.  Written by h-L3-build/gen/gen_hnet.py,
prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace XHost
open Binary Matrix Finset RAM SS CB RF BR GF GX
noncomputable section

section Net
variable {H T Sl C α Γ : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α]
  [Fintype Γ] [DecidableEq Γ]
variable {nb : ℕ} (N : GHost401 H T Sl C α Γ nb)

/-- the whole run `ρ` before the exchange, on the roles of a run. -/
lemma GHost401.run (ρ : Fin 8) :
    GRoute (N.S0 ρ) (N.S25 ρ) (chainQ (N.bk ρ)) (fun φ => (Fintype.card Γ : ℝ) * N.runUnit φ) := by
  have R := ((((((((((((((((((((((((((((((((((N.adapt0 ρ)).trans (N.stage0 ρ)).trans (N.adapt1 ρ)).trans (N.stage1 ρ)).trans ((N.climbA2 ρ))).trans (GRoute.ofX (N.gateA2 ρ))).trans (N.adapt2 ρ)).trans (N.stage2 ρ)).trans ((N.climbB2 ρ))).trans (GRoute.ofX (N.gateB2 ρ))).trans (N.adapt3 ρ)).trans (N.stage3 ρ)).trans ((N.climbA3 ρ))).trans (GRoute.ofX (N.gateA3 ρ))).trans (N.adapt4 ρ)).trans (N.stage4 ρ)).trans ((N.climbB3 ρ))).trans (GRoute.ofX (N.gateB3 ρ))).trans (N.adapt5 ρ)).trans (N.stage5 ρ)).trans ((N.climbA4 ρ))).trans (GRoute.ofX (N.gateA4 ρ))).trans (N.adapt6 ρ)).trans (N.stage6 ρ)).trans ((N.climbB4 ρ))).trans (GRoute.ofX (N.gateB4 ρ))).trans (N.adapt7 ρ)).trans (N.stage7 ρ)).trans (N.adapt8 ρ)).trans (N.stage8 ρ)).trans ((N.fin ρ))).trans (GRoute.ofX (N.gateC4 ρ))).trans (GRoute.ofX (N.gateC3 ρ))).trans (GRoute.ofX (N.gateC2 ρ))
  refine R.cast (fun φ => ?_)
  simp only [GHost401.runUnit]
  ring

/-- run `ρ` on all roles: the other replicas wait (before their run at `InitD`, after it at `FinD`). -/
lemma GHost401.runFull (ρ : Fin 8) :
    GRoute (N.Sfull ρ.castSucc) (N.Sfull ρ.succ) (liftRun4 ρ (chainQ (N.bk ρ)))
      (fun φ => (Fintype.card Γ : ℝ) * N.runUnit φ) := by
  refine GRoute.lift (emRun4 ρ) (h := chainQ (N.bk ρ)) ?_ (fun x => ?_) (fun x i hi => ?_) (fun i hi => ?_)
  · have e1 : N.Sfull ρ.castSucc ∘ (emRun4 ρ) = N.S0 ρ := by
      funext r
      rcases r with ((d|⟨g,b⟩)|⟨g,c⟩)
      · exact if_neg (lt_irrefl _)
      · rfl
      · rfl
    have e2 : N.Sfull ρ.succ ∘ (emRun4 ρ) = N.S25 ρ := by
      funext r
      rcases r with ((d|⟨g,b⟩)|⟨g,c⟩)
      · exact if_pos (Nat.lt_succ_self _)
      · have e : 9 * ρ.succ.val = 9 * ρ.val + 9 := by
          rw [Fin.val_succ]
          ring
        show N.Zb (9 * ρ.succ.val) g b = N.Zb (9 * ρ.val + 9) g b
        rw [e]
      · rfl
    rw [e1, e2]
    exact N.run ρ
  · funext r
    rcases r with ((d|⟨g,b⟩)|⟨g,c⟩)
    · exact if_pos rfl
    · rfl
    · rfl
  · rcases i with ((⟨ρ',d⟩|⟨g,b⟩)|⟨g,c⟩)
    · by_cases e : ρ' = ρ
      · subst e
        exact absurd ((cover_iff _ _).mpr ⟨Sum.inl (Sum.inl d), rfl⟩) hi
      · exact if_neg e
    · exact absurd ((cover_iff _ _).mpr ⟨Sum.inl (Sum.inr (g, b)), rfl⟩) hi
    · exact absurd ((cover_iff _ _).mpr ⟨Sum.inr (g, c), rfl⟩) hi
  · rcases i with ((⟨ρ',d⟩|⟨g,b⟩)|⟨g,c⟩)
    · by_cases e : ρ' = ρ
      · subst e
        exact absurd ((cover_iff _ _).mpr ⟨Sum.inl (Sum.inl d), rfl⟩) hi
      · have e' : ρ'.val ≠ ρ.val := fun h => e (Fin.ext h)
        show (if ρ'.val < ρ.castSucc.val then N.FinD d else N.InitD d)
          = (if ρ'.val < ρ.succ.val then N.FinD d else N.InitD d)
        simp only [Fin.val_castSucc, Fin.val_succ]
        split_ifs <;> first | rfl | (exfalso; omega)
    · exact absurd ((cover_iff _ _).mpr ⟨Sum.inl (Sum.inr (g, b)), rfl⟩) hi
    · exact absurd ((cover_iff _ _).mpr ⟨Sum.inr (g, c), rfl⟩) hi

/-- **the 8 runs one after the other.** -/
lemma GHost401.sweep :
    GRoute (N.Sfull 0) (N.Sfull (Fin.last 8)) (sweepMap (fun ρ : Fin 8 => liftRun4 ρ (chainQ (N.bk ρ))))
      (fun φ => ((8 : ℕ) : ℝ) * ((Fintype.card Γ : ℝ) * N.runUnit φ)) :=
  GRoute.sweep_const 8 N.Sfull _ (fun _ φ => (Fintype.card Γ : ℝ) * N.runUnit φ) _ (fun _ _ => rfl)
    (fun ρ => N.runFull ρ)

/-- **the padding blocks**: ONE block of rank `pad b` per bank. -/
lemma GHost401.padR :
    GRoute (N.Sfull (Fin.last 8)) N.Sfin id
      (fun φ => (Fintype.card Γ : ℝ) * ∑ b, bcost φ (N.pad b)) := by
  have H1 := GRoute.reach_all (N.Sfull (Fin.last 8)) N.Sfin
    (mk401 (fun _ _ _ => 0) (fun _ b φ => bcost φ (N.pad b)) (fun _ _ _ => 0)) (by
      rintro ((⟨ρ',d⟩|⟨g,b⟩)|⟨g,c⟩)
      · have e : N.Sfull (Fin.last 8) (Sum.inl (Sum.inl (ρ', d))) = N.FinD d := if_pos ρ'.isLt
        show GReach (N.Sfull (Fin.last 8) (Sum.inl (Sum.inl (ρ', d)))) (N.FinD d) (fun _ => 0)
        rw [e]
        exact GReach.refl _
      · exact N.bpad g b
      · exact GReach.refl _)
  refine H1.cast (fun φ => ?_)
  simp_rw [mk401_map (fun c : (ℕ → ℝ) → ℝ => c φ)]
  rw [sum_mk401]
  simp only [sum_const, card_univ, nsmul_eq_mul, mul_zero, add_zero, zero_add, Finset.sum_const_zero]

/-- the free shifts at the end. -/
def GHost401.corr : BRole401 Γ T nb C → Space α :=
  mk401 (fun _ => mkD4 (fun _ _ => 0) (fun g t => N.u g t) (fun _ _ => 0) (fun g t => N.u g t) (fun _ _ => 0) (fun g t => N.u g t) (fun _ _ => 0) (fun g t => N.u g t)) (fun _ _ => 0) (fun _ _ => 0)

def GHost401.Sz : BRole401 Γ T nb C → CMat α := fun r => shift (N.corr r) * N.Sfin r

/-- **The scratch certificate of the nine-stage word with entrance gauges and residual banks (PR #401)**, in whole
blocks, with its exact price for every price list. -/
theorem ghost_certificate :
    GCert α (Sum.inl : BLive401 Γ T nb → BRole401 Γ T nb C)
      (fun φ => (Fintype.card Γ : ℝ) * N.unit φ) := by
  have R1 := N.sweep.trans N.padR
  have R2 : GRoute N.Sfin N.Sz id (fun _ => 0) := GRoute.shifts N.hm N.Sfin N.corr
  have R3 : GRoute N.Sz (fun r => N.Sz (bpull401 r)) (actPoint brec401) (fun _ => 0) :=
    GRoute.gate brec401 N.Sz _ (fun i j h => by
      have hj : j = bpull401 i := by
        by_contra hne
        exact h (by simp [brec401, hne])
      rw [hj])
  have Rt := ((R1.trans R2).trans R3).cast (c' := fun φ => (Fintype.card Γ : ℝ) * N.unit φ)
    (fun φ => by simp only [GHost401.unit]; push_cast; ring)
  have h0 : ∀ (ρ' : Fin 8) (d : BDat4 Γ T), N.Sfull 0 (Sum.inl (Sum.inl (ρ', d))) = N.InitD d :=
    fun ρ' d => by simp [GHost401.Sfull, mk401]
  have hS : ∀ r, ∃ M, N.Sfull 0 r * M = 1 := by
    rintro ((⟨ρ',d⟩|⟨g,b⟩)|⟨g,c⟩)
    · rw [h0]
      rcases d with ((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))
      · exact N.inv0 g t
      · refine ⟨1, ?_⟩
        show (N.G0 g).Xm.Φ 0 * 1 = 1
        rw [N.z0, Matrix.mul_one]
      · exact N.inv0 g t
      · refine ⟨1, ?_⟩
        show (N.G0 g).Xm.Φ 0 * 1 = 1
        rw [N.z0, Matrix.mul_one]
      · exact N.inv0 g t
      · refine ⟨1, ?_⟩
        show (N.G0 g).Xm.Φ 0 * 1 = 1
        rw [N.z0, Matrix.mul_one]
      · exact N.inv0 g t
      · refine ⟨1, ?_⟩
        show (N.G0 g).Xm.Φ 0 * 1 = 1
        rw [N.z0, Matrix.mul_one]
    · refine ⟨1, ?_⟩
      show N.Zb 0 g b * 1 = 1
      rw [N.z0b, Matrix.mul_one]
    · exact ⟨1, Matrix.mul_one _⟩
  have hT : ∀ l : BLive401 Γ T nb, N.Sz (bpull401 (Sum.inl l)) = kernel α * N.Sfull 0 (Sum.inl l) := by
    rintro (⟨ρ', d⟩|⟨g,b⟩)
    · rw [h0]
      rcases d with ((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))
      · show shift (N.u g t) * N.YF g t = kernel α * (N.G0 g).Xm.Φ (tt (N.I.tv t))
        exact (N.term g t).symm
      · show shift (0 : Space α) * kernel α = kernel α * (N.G0 g).Xm.Φ 0
        rw [shift_zero, Matrix.one_mul, N.z0, Matrix.mul_one]
      · show shift (N.u g t) * N.YF g t = kernel α * (N.G0 g).Xm.Φ (tt (N.I.tv t))
        exact (N.term g t).symm
      · show shift (0 : Space α) * kernel α = kernel α * (N.G0 g).Xm.Φ 0
        rw [shift_zero, Matrix.one_mul, N.z0, Matrix.mul_one]
      · show shift (N.u g t) * N.YF g t = kernel α * (N.G0 g).Xm.Φ (tt (N.I.tv t))
        exact (N.term g t).symm
      · show shift (0 : Space α) * kernel α = kernel α * (N.G0 g).Xm.Φ 0
        rw [shift_zero, Matrix.one_mul, N.z0, Matrix.mul_one]
      · show shift (N.u g t) * N.YF g t = kernel α * (N.G0 g).Xm.Φ (tt (N.I.tv t))
        exact (N.term g t).symm
      · show shift (0 : Space α) * kernel α = kernel α * (N.G0 g).Xm.Φ 0
        rw [shift_zero, Matrix.one_mul, N.z0, Matrix.mul_one]
    · show shift (0 : Space α) * kernel α = kernel α * N.Zb 0 g b
      rw [shift_zero, Matrix.one_mul, N.z0b, Matrix.mul_one]
  have hg : ∀ x : BRole401 Γ T nb C → ℂ,
      (∀ s, (∀ l : BLive401 Γ T nb, (Sum.inl l : BRole401 Γ T nb C) ≠ s) → x s = 0) →
      ∀ l : BLive401 Γ T nb,
        (actPoint brec401 ∘ (id ∘ (id ∘ sweepMap (fun ρ : Fin 8 => liftRun4 ρ (chainQ (N.bk ρ)))))) x (Sum.inl l)
          = x (Sum.inl l) := by
    intro x _ l
    rw [Function.comp_apply, Function.comp_apply, Function.comp_apply, id_eq, id_eq, actPoint_brec401]
    rcases l with (⟨ρ, d⟩|⟨g,b⟩)
    · rcases d with ((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))
      · show ((1:ℚ):ℂ) * sweepMap (fun i : Fin 8 => liftRun4 i (chainQ (N.bk i))) x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inl (Sum.inr ((g,t)))))))) = _
        rw [sweep_data4 N.bk, chainQ_eq]
        show ((1:ℚ):ℂ) * x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inl (Sum.inl ((g,t)))))))) = x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inl (Sum.inl ((g,t))))))))
        first | (push_cast; ring) | simp
      · show ((-1:ℚ):ℂ) * sweepMap (fun i : Fin 8 => liftRun4 i (chainQ (N.bk i))) x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inl (Sum.inl ((g,t)))))))) = _
        rw [sweep_data4 N.bk, chainQ_eq]
        show ((-1:ℚ):ℂ) * (-x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inl (Sum.inr ((g,t))))))))) = x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inl (Sum.inr ((g,t))))))))
        first | (push_cast; ring) | simp
      · show ((1:ℚ):ℂ) * sweepMap (fun i : Fin 8 => liftRun4 i (chainQ (N.bk i))) x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inr (Sum.inr ((g,t)))))))) = _
        rw [sweep_data4 N.bk, chainQ_eq]
        show ((1:ℚ):ℂ) * x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inr (Sum.inl ((g,t)))))))) = x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inr (Sum.inl ((g,t))))))))
        first | (push_cast; ring) | simp
      · show ((-1:ℚ):ℂ) * sweepMap (fun i : Fin 8 => liftRun4 i (chainQ (N.bk i))) x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inr (Sum.inl ((g,t)))))))) = _
        rw [sweep_data4 N.bk, chainQ_eq]
        show ((-1:ℚ):ℂ) * (-x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inr (Sum.inr ((g,t))))))))) = x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inl (Sum.inr (Sum.inr ((g,t))))))))
        first | (push_cast; ring) | simp
      · show ((1:ℚ):ℂ) * sweepMap (fun i : Fin 8 => liftRun4 i (chainQ (N.bk i))) x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inr (Sum.inr ((g,t))))))) = _
        rw [sweep_data4 N.bk, chainQ_eq]
        show ((1:ℚ):ℂ) * x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inr (Sum.inl ((g,t))))))) = x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inr (Sum.inl ((g,t)))))))
        first | (push_cast; ring) | simp
      · show ((-1:ℚ):ℂ) * sweepMap (fun i : Fin 8 => liftRun4 i (chainQ (N.bk i))) x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inr (Sum.inl ((g,t))))))) = _
        rw [sweep_data4 N.bk, chainQ_eq]
        show ((-1:ℚ):ℂ) * (-x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inr (Sum.inr ((g,t)))))))) = x (Sum.inl (Sum.inl (ρ, Sum.inl (Sum.inr (Sum.inr ((g,t)))))))
        first | (push_cast; ring) | simp
      · show ((1:ℚ):ℂ) * sweepMap (fun i : Fin 8 => liftRun4 i (chainQ (N.bk i))) x (Sum.inl (Sum.inl (ρ, Sum.inr (Sum.inr ((g,t)))))) = _
        rw [sweep_data4 N.bk, chainQ_eq]
        show ((1:ℚ):ℂ) * x (Sum.inl (Sum.inl (ρ, Sum.inr (Sum.inl ((g,t)))))) = x (Sum.inl (Sum.inl (ρ, Sum.inr (Sum.inl ((g,t))))))
        first | (push_cast; ring) | simp
      · show ((-1:ℚ):ℂ) * sweepMap (fun i : Fin 8 => liftRun4 i (chainQ (N.bk i))) x (Sum.inl (Sum.inl (ρ, Sum.inr (Sum.inl ((g,t)))))) = _
        rw [sweep_data4 N.bk, chainQ_eq]
        show ((-1:ℚ):ℂ) * (-x (Sum.inl (Sum.inl (ρ, Sum.inr (Sum.inr ((g,t))))))) = x (Sum.inl (Sum.inl (ρ, Sum.inr (Sum.inr ((g,t))))))
        first | (push_cast; ring) | simp
    · show ((1:ℚ):ℂ) * sweepMap (fun i : Fin 8 => liftRun4 i (chainQ (N.bk i))) x (Sum.inl (Sum.inr (g, b))) = _
      rw [sweep_bank4 N.bk]
      first | (push_cast; ring) | simp
  exact gcert_of_groute (Sum.inl : BLive401 Γ T nb → BRole401 Γ T nb C) (N.Sfull 0)
    (fun r => N.Sz (bpull401 r)) _ _ Rt hS hg hT

end Net
end
end XHost
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHost.GHost401.run
#print axioms OAI.PowerSaving.XHost.GHost401.runFull
#print axioms OAI.PowerSaving.XHost.ghost_certificate
