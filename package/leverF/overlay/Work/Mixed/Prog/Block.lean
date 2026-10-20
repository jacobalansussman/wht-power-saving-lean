import Work.Mixed.Prog.Route

/-!
# Mixed engine, R3 program part 4: one block letter of a class = ONE call (agent e-r3prog)

The handler of the lockstep program has ONE port for top jobs and classes:

* `portIn cl`  = `(n, I, kind, μ, o, own jobs, supply)` (`kind = 0`: top job; `kind = h+1`: class of
                 height `h`; `μ` its multiplier; `o` the supply offset of its subtree);
* `portOut cl` = `(own jobs answered, consumed supply answered, sizes of the consumed supply)`.

`CCtx.block_letter`: gather the child's input, call, read the answers back.  Work linear in the
class and in the child's input; charge = the child's cost.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF
noncomputable section

universe U

/-- A tape of jobs (each job a tape of `2^n` entries). -/
abbrev jobsT (cl : Paint) : Ty := a (a (c cl))
/-- The input of every call: `(n, I, kind, μ, o, own jobs, supply)`. -/
abbrev portIn (cl : Paint) : Ty := p w (p sc (p w (p w (p w (p (jobsT cl) (jobsT cl))))))
/-- The output of every call: `(own jobs answered, consumed supply answered, consumed sizes)`. -/
abbrev portOut (cl : Paint) : Ty := p (jobsT cl) (p (jobsT cl) (a w))
abbrev portM (cl : Paint) : Port := some (portIn cl, portOut cl)

/-- An input value of the port. -/
def mkIn (cl : Paint) (n kind μ o : ℕ) (own sup : Tape (Tape ℂ)) : T (portIn cl) :=
  (n, Complex.I, kind, μ, o, own, sup)

/-- A tape of jobs from a family of job functions. -/
def jobsTape (n N : ℕ) (J : ℕ → Bits n → ℂ) : Tape (Tape ℂ) :=
  (Layout.std N).tape (fun t => (bits n).tape (J t))

section Block
variable {α β ρ : Type*} [Fintype α] [DecidableEq α] [Fintype β] [DecidableEq β]
  [Fintype ρ] [DecidableEq ρ]
variable {k : ℕ}

lemma blockJob_congr (X : Layout α) (S : Split α β (Fin k)) (X' : Layout β) (hX' : X'.n + k = X.n)
    (l : ρ) {NU n : ℕ} {col ib : ℕ → ℕ} (hfit : ∀ p : Fin NU, X.n * col p + ib p = n)
    (V : ∀ p : Fin NU, Bits (ib p) → Grain α ρ 0 (col p)) (p p' : Fin NU) (hp : p = p') :
    blockJob S X' (col p) (ib p) n (fitk X' hX' (hfit p)) l (V p) =
      blockJob S X' (col p') (ib p') n (fitk X' hX' (hfit p')) l (V p') := by
  subst hp; rfl

/-- A fed slot whose parent unit is `p` carries the block job of `p`. -/
lemma cjob_fed (X : Layout α) (S : Split α β (Fin k)) (X' : Layout β) (hX' : X'.n + k = X.n)
    (l : ρ) {NU n : ℕ} {col ib : ℕ → ℕ} (hfit : ∀ p : Fin NU, X.n * col p + ib p = n)
    (V : ∀ p : Fin NU, Bits (ib p) → Grain α ρ 0 (col p))
    (fl pp si : ℕ → ℕ) (hpp : ∀ t, fl t ≠ 0 → pp t < NU) (sup : ℕ → Bits n → ℂ) (t : ℕ)
    (p : Fin NU) (hf : fl t ≠ 0) (hpt : pp t = p) :
    cjob X S X' hX' l hfit V fl pp si hpp sup t =
      blockJob S X' (col p) (ib p) n (fitk X' hX' (hfit p)) l (V p) := by
  unfold cjob
  rw [dif_pos hf]
  exact blockJob_congr X S X' hX' l hfit V ⟨pp t, hpp t hf⟩ p (Fin.ext hpt)

end Block

section BlockRAM
variable {ι : Type U} {m : Bool} {cl : Paint} {b : Bank (portM cl) ι} {s : Ty} {x : ι → T s}
  {n NU : ι → ℕ} {col ib : ι → ℕ → ℕ}
variable {α β ρ : Type*} [Fintype α] [DecidableEq α] [Fintype β] [DecidableEq β]
  [Fintype ρ] [DecidableEq ρ]
variable {X : Layout α} {R : Layout ρ} {k : ℕ}

/-- **One block letter of a class: gather, ONE call, read back.** -/
theorem CCtx.block_letter (h : CCtx m b.B s x X.n n NU col ib)
    {V : ∀ i (p : Fin (NU i)), Bits (ib i p) → Grain α ρ 0 (col i p)}
    (ct : Has m b.B s (a (a (a (c cl)))) x (fun i => ctab X R (NU i) (col i) (ib i) (V i)))
    (S : Split α β (Fin k)) (X' : Layout β) (hX' : X'.n + k = X.n) (l : ρ)
    {NJ BS kd μ o : ι → ℕ} {fl pp si slot sz : ι → ℕ → ℕ} {sup : ∀ i, ℕ → Bits (n i) → ℂ}
    (hpp : ∀ i t, fl i t ≠ 0 → pp i t < NU i)
    (hsi : ∀ i t, t < NJ i → fl i t = 0 → si i t < BS i)
    (hslot : ∀ i (p : Fin (NU i)), slot i p < NJ i ∧ fl i (slot i p) ≠ 0 ∧ pp i (slot i p) = p ∧
      sz i (slot i p) = k * col i p)
    (hNJ : Knows m b.B s x NJ) (hNJs : Small b.B NJ)
    (tfl : Has m b.B s (a w) x (fun i => (Layout.std (NJ i)).tape (fun t => fl i t)))
    (tpp : Has m b.B s (a w) x (fun i => (Layout.std (NJ i)).tape (fun t => pp i t)))
    (tsi : Has m b.B s (a w) x (fun i => (Layout.std (NJ i)).tape (fun t => si i t)))
    (tslot : Has m b.B s (a w) x (fun i => (Layout.std (NU i)).tape (fun p => slot i p)))
    (tsup : Has m b.B s (jobsT cl) x (fun i => jobsTape (n i) (BS i) (sup i)))
    (hkd : Knows m b.B s x kd) (hμ : Knows m b.B s x μ) (ho : Knows m b.B s x o)
    {ACC : ι → Tape (Tape ℂ)} {SZ : ι → Tape ℕ} {τ : ι → ℕ}
    (hk : ∀ i, (b.use i (mkIn cl (n i) (kd i) (μ i) (o i)
        (jobsTape (n i) (NJ i) (cjob X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i)))
        (jobsTape (n i) (BS i) (sup i)))).OK
      (jobsTape (n i) (NJ i) (fun t => hlow (n i) (sz i t)
          (cjob X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i) t)),
        ACC i, SZ i) (τ i) (b.ext i)) :
    Able m b s (p (a (a (a (c cl)))) (portOut cl)) x
      (fun i => (ctab X R (NU i) (col i) (ib i)
          (fun p q => ferryM 0 (col i p) l (S.lift (kernel (Fin k))) (V i p q)),
        (jobsTape (n i) (NJ i) (fun t => hlow (n i) (sz i t)
          (cjob X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i) t)),
        ACC i, SZ i)))
      (fun i => NJ i * 2^(n i) + NU i * 2^(n i)) τ := by
  -- the child's input, then the call
  have HG := h.gather ct S X' hX' l hpp hsi hNJ hNJs tfl tpp tsi tsup
  have Hin : Can m b.B s (portIn cl) x
      (fun i => mkIn cl (n i) (kd i) (μ i) (o i)
        (jobsTape (n i) (NJ i) (cjob X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i)))
        (jobsTape (n i) (BS i) (sup i)))
      (fun i => NJ i * 2^(n i)) :=
    (h.nK.mono (fun _ => Nat.zero_le _)).pair ((h.imag.mono (fun _ => Nat.zero_le _)).pair
      ((hkd.mono (fun _ => Nat.zero_le _)).pair ((hμ.mono (fun _ => Nat.zero_le _)).pair
        ((ho.mono (fun _ => Nat.zero_le _)).pair (HG.pair (tsup.mono (fun _ => Nat.zero_le _)))))))
  have Hcall := (Able.of_can (b:=b) Hin).comp (Able.use_call (m:=m) (b:=b) hk)
  -- read back
  let outp (i : ι) : T (portOut cl) :=
    (jobsTape (n i) (NJ i) (fun t => hlow (n i) (sz i t)
      (cjob X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i) t)), ACC i, SZ i)
  have h' := h.ext outp
  have F : ∀ {u : Ty} {z : ι → T u}, Has m b.B s u x z → Has m b.B (p s (portOut cl)) u (fun i => (x i, outp i)) z :=
    fun g => Can.first.then_do g
  have hY : ∀ i (p : Fin (NU i)) (q : Bits (ib i p)) (y : Sky α (col i p)),
      hlow (n i) (sz i (slot i p))
        (cjob X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i) (slot i p))
        (placeS S X' (col i p) (ib i p) (n i) (fitk X' hX' (h'.fit i p)) (q, y)) =
      ferryM 0 (col i p) l (S.lift (kernel (Fin k))) (V i p q) (l, (0, y)) := by
    intro i p q y
    obtain ⟨_, hf, hpt, hsz⟩ := hslot i p
    rw [cjob_fed X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i) (slot i p) p hf hpt,
      hsz]
    exact blockJob_hlow S X' _ _ _ _ l (V i p) q 0 y
  have tout : Has m b.B (p s (portOut cl)) (jobsT cl) (fun i => (x i, outp i))
      (fun i => jobsTape (n i) (NJ i) (fun t => hlow (n i) (sz i t)
        (cjob X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i) t))) :=
    Can.fst Can.second
  have HS := h'.scatter (Y:=fun i t => hlow (n i) (sz i t)
      (cjob X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i) t))
    (F ct) S X' hX' l (fun i p => (hslot i p).1) (F tslot) tout hY
  have Hout := Able.of_can (b:=b) (HS.pair (Can.second.mono (fun _ => Nat.zero_le _)))
  refine (Hcall.save Hout).weaken (Dom.of_le (fun i => ?_)) (fun i => by simp)
  simp

end BlockRAM

end
end PowerSaving.RAM
end OAI
