import Work.PfLean.Bounds
import Work.PfLean.Dem3Defs

/-!
# The demand of the block-3 hook, and the block-3 bound `14 W ≤ 5 · cost` (px-lean, 14 October 2026)

Block size 3 (`α = Fin 3`), `f = 5` digit columns, certificate `hookN 3` (dimension 280).
We prove `dem3 : 224 ≤ raise (hookN 3) (matL 5 (kernel (Fin 3)))` and hence, with pf-lean's
`block3_bound_of_dem`, the unconditional `block3_bound`.

Method (no rank is trusted):
* `entryF`: the value of `kernel` applied to `Alt (δ σ)` at a tuple `τ` is
  `ca^15 · Pr(τ0, σ0) · det [Pr(τ_(b+1), σ_(a+1))]`, `Pr(x, y) = ∏ᵢ (1 if xᵢ = yᵢ else -i)`.
* Translations commute with the kernel and keep `N`; character sums of translates
  (`nv`, functionals `Ψ`) make the `224 × 224` system block diagonal: 8 blocks of `28 × 28`
  (`dbl_sum`, from the orthogonality of the characters, itself from `wal * walI = 1`).
* Each block is a matrix over `ℤ[i]` (`GenR`, written for any commutative ring and `ι`).
  Its determinant is nonzero because its image under `ℤ[i] → ZMod 65537`, `i ↦ 256`, has an
  explicit inverse: `Dem3Comp.blkOK` checks `B · G ≡ 1 (mod 65537)` by kernel evaluation
  (`decide +kernel`, no `native_decide`), where `G` is recomputed from the formula.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
set_option maxHeartbeats 8000000

namespace OAI
namespace PowerSaving.PfLean
open Binary GF Matrix Module Submodule Complex
open Dem3Data Dem3Comp

namespace Dem3

abbrev S3 := Space (Fin 3)

/-! ### 1. Rings: the block entry for any commutative ring -/

/-- `∏ᵢ (1 if xᵢ = yᵢ else -ι)`. -/
def PrR {R : Type*} [CommRing R] (ι : R) (x y : S3) : R := ∏ i, (if x i = y i then 1 else -ι)

/-- The character `(-1)^(s·u)`. -/
def chiR {R : Type*} [CommRing R] (s u : S3) : R := if dot s u = 0 then 1 else -1

/-- The block entry for the character `s`, row tuple `τ`, column tuple `σ`. -/
noncomputable def GenR {R : Type*} [CommRing R] (ι : R) (s : S3) (τ σ : Sky (Fin 3) 5) : R :=
  ∑ u : S3, chiR s u * (PrR ι (τ 0 + u) (σ 0) *
    Matrix.det (Matrix.of fun b a : Fin 4 => PrR ι (τ b.succ + u) (σ a.succ)))

lemma PrR_map {R R' : Type*} [CommRing R] [CommRing R'] (φ : R →+* R') (ι : R) (x y : S3) :
    φ (PrR ι x y) = PrR (φ ι) x y := by
  unfold PrR; rw [map_prod]; refine Finset.prod_congr rfl (fun i _ => ?_); split_ifs <;> simp

lemma chiR_map {R R' : Type*} [CommRing R] [CommRing R'] (φ : R →+* R') (s u : S3) :
    φ (chiR s u) = chiR s u := by
  unfold chiR; split_ifs <;> simp

lemma GenR_map {R R' : Type*} [CommRing R] [CommRing R'] (φ : R →+* R') (ι : R) (s : S3)
    (τ σ : Sky (Fin 3) 5) : φ (GenR ι s τ σ) = GenR (φ ι) s τ σ := by
  unfold GenR
  rw [map_sum]
  refine Finset.sum_congr rfl (fun u _ => ?_)
  rw [map_mul, map_mul, RingHom.map_det, chiR_map, PrR_map]
  congr 3
  ext b a
  simp [RingHom.mapMatrix_apply, PrR_map]

/-! ### 2. The kernel entry formula -/

lemma matAct_del5 (M : CMat (Fin 3)) (σ τ : Sky (Fin 3) 5) :
    matAct 5 M (del σ) τ = ∏ l, M (τ l) (σ l) := by
  simp [matAct, del, Matrix.mulVec_single_one, Matrix.col_apply, digitProd_apply]

lemma littleC_eq (a b : F) : littleC a b = ca * (if a = b then 1 else -I) := by
  unfold littleC
  split_ifs
  · ring
  · simp only [ca, cb]; linear_combination (1/2 : ℂ) * I_sq

lemma kernel3 (x y : S3) : kernel (Fin 3) x y = ca ^ 3 * PrR I x y := by
  simp only [kernel, digitProd_apply, littleC_eq, Finset.prod_mul_distrib, Finset.prod_const,
    Finset.card_univ, Fintype.card_fin, PrR]

/-- **The entry formula.** -/
lemma entryF (τ σ : Sky (Fin 3) 5) :
    matAct 5 (kernel (Fin 3)) (Alt 3 (del σ)) τ =
      ca ^ 15 * (PrR I (τ 0) (σ 0) *
        Matrix.det (Matrix.of fun b a : Fin 4 => PrR I (τ b.succ) (σ a.succ))) := by
  rw [matAct_Alt, Alt_apply, Finset.sum_apply, Matrix.det_apply, Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun π _ => ?_)
  rw [Pi.smul_apply, permV_apply, matAct_del5, Fin.prod_univ_succ]
  simp only [Function.comp_apply, lift1_zero, lift1_succ, kernel3, Matrix.of_apply,
    Finset.prod_mul_distrib, Finset.prod_const, Finset.card_univ, Fintype.card_fin,
    Units.smul_def, zsmul_eq_mul, smul_eq_mul]
  ring

/-! ### 3. Translations and characters -/

/-- Add `g` to every digit column. -/
def tr (g : S3) (τ : Sky (Fin 3) 5) : Sky (Fin 3) 5 := fun l => τ l + g

lemma tr_tr (g h : S3) (τ : Sky (Fin 3) 5) : tr h (tr g τ) = tr (g + h) τ := by
  funext l; simp [tr, add_assoc]

lemma shift_act (h : S3) (v : Sky (Fin 3) 5 → ℂ) (τ : Sky (Fin 3) 5) :
    matAct 5 (shift h) v τ = v (tr h τ) := by
  rw [matAct_mono (shift h) (· + h) (fun _ => 1) (fun x y => by simp [Binary.shift])]
  simp only [Finset.prod_const_one, one_mul]; rfl

lemma K3_shift (h : S3) (v : Sky (Fin 3) 5 → ℂ) :
    matAct 5 (kernel (Fin 3)) (matAct 5 (shift h) v) =
      matAct 5 (shift h) (matAct 5 (kernel (Fin 3)) v) := by
  rw [matAct_mul, matAct_mul]
  congr 1
  rw [← wrap_lum, shift_phase]; exact wrap_comm _ _

lemma walI3 (x y : S3) : walI (Fin 3) x y = sign (dot x y) / 8 := by
  simp only [walI, digitProd_apply, littleH', Finset.prod_div_distrib, Finset.prod_const,
    Finset.card_univ, Fintype.card_fin, dot, sign_sum]
  norm_num

lemma chi_orth (s s' : S3) :
    ∑ h : S3, sign (dot s h) * sign (dot s' h) = if s = s' then 8 else 0 := by
  have h1 := congrFun (congrFun (wal_right (Fin 3)) s) s'
  rw [Matrix.mul_apply, Matrix.one_apply] at h1
  simp only [wal_apply, walI3] at h1
  have h2 : ∑ h : S3, sign (dot s h) * sign (dot s' h) =
      8 * ∑ h : S3, sign (dot s h) * (sign (dot h s') / 8) := by
    rw [Finset.mul_sum]; refine Finset.sum_congr rfl (fun h _ => ?_); rw [dot_symm h s']; ring
  rw [h2, h1]; split_ifs <;> norm_num

lemma inner_shift (s h : S3) (F : S3 → ℂ) :
    ∑ g : S3, sign (dot s g) * F (g + h) = sign (dot s h) * ∑ u : S3, sign (dot s u) * F u := by
  rw [Finset.mul_sum]
  refine Fintype.sum_equiv (Equiv.addRight h) _ _ (fun g => ?_)
  simp only [Equiv.coe_addRight, dot_add_right, sign_add]
  have := sign_sq (dot s h)
  linear_combination (-(sign (dot s g) * F (g + h))) * this

/-- **Block diagonality.** -/
lemma dbl_sum (s s' : S3) (F : S3 → ℂ) :
    ∑ g : S3, sign (dot s g) * ∑ h : S3, sign (dot s' h) * F (g + h) =
      if s = s' then 8 * ∑ u : S3, sign (dot s u) * F u else 0 := by
  have e1 : ∑ g : S3, sign (dot s g) * ∑ h : S3, sign (dot s' h) * F (g + h) =
      ∑ h : S3, sign (dot s' h) * ∑ g : S3, sign (dot s g) * F (g + h) := by
    simp_rw [Finset.mul_sum]
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl (fun h _ => Finset.sum_congr rfl (fun g _ => by ring))
  rw [e1]
  simp_rw [inner_shift s _ F]
  have e2 : ∑ h : S3, sign (dot s' h) * (sign (dot s h) * ∑ u : S3, sign (dot s u) * F u) =
      (∑ h : S3, sign (dot s h) * sign (dot s' h)) * ∑ u : S3, sign (dot s u) * F u := by
    rw [Finset.sum_mul]; exact Finset.sum_congr rfl (fun h _ => by ring)
  rw [e2, chi_orth]; split_ifs <;> simp

/-! ### 4. Addresses, rows and columns from the data -/

/-- The address with code `k` (bit `i` of `k`). -/
def ad (k : ℕ) : S3 := fun i => (((k >>> i.val) % 2 : ℕ) : ZMod 2)

lemma ad_zero : ad 0 = 0 := by funext i; simp [ad]

def usD (s r : ℕ) : List ℕ := (rowsD.getD s []).getD r []
def tsD (s c : ℕ) : List ℕ := (colsD.getD s []).getD c []

/-- Row tuple `(0, u1, u2, u3, u4)`. -/
def rowT (s r : ℕ) : Sky (Fin 3) 5 := Fin.cons 0 (fun b : Fin 4 => ad ((usD s r).getD b.val 0))
/-- Column tuple `(0, 0, t1, t2, t3)`. -/
def colT (s c : ℕ) : Sky (Fin 3) 5 := tup 0 (fun a : Fin 3 => ad ((tsD s c).getD a.val 0))

lemma rowT_zero (s r : ℕ) : rowT s r 0 = 0 := rfl
lemma rowT_succ (s r : ℕ) (b : Fin 4) : rowT s r b.succ = ad ((usD s r).getD b.val 0) := rfl
lemma colT_zero (s c : ℕ) : colT s c 0 = 0 := rfl
lemma colT_succ (s c : ℕ) (a : Fin 4) :
    colT s c a.succ = ad ((0 :: tsD s c).getD a.val 0) := by
  cases a using Fin.cases with
  | zero => simp [colT, ad_zero]
  | succ a => simp [colT]

lemma ad_inj3 : ∀ a b : Fin 8, ad a.val 0 = ad b.val 0 → ad a.val 1 = ad b.val 1 →
    ad a.val 2 = ad b.val 2 → a = b := by decide

lemma ad_bij : Function.Bijective (fun k : Fin 8 => ad k.val) := by
  rw [Fintype.bijective_iff_injective_and_card]
  refine ⟨fun a b h => ad_inj3 a b (congrFun h 0) (congrFun h 1) (congrFun h 2), ?_⟩
  simp [S3, Space, Fintype.card_fun, ZMod.card]

lemma sum_S3 {M : Type*} [AddCommMonoid M] (f : S3 → M) :
    ∑ u : S3, f u = sum8 (fun k => f (ad k)) := by
  rw [← ad_bij.sum_comp f, Fin.sum_univ_eight]; rfl

lemma sum8_congr {M : Type*} [AddCommMonoid M] (f g : ℕ → M) (h : ∀ k < 8, f k = g k) :
    sum8 f = sum8 g := by
  unfold sum8
  rw [h 0 (by norm_num), h 1 (by norm_num), h 2 (by norm_num), h 3 (by norm_num),
    h 4 (by norm_num), h 5 (by norm_num), h 6 (by norm_num), h 7 (by norm_num)]

lemma det_eq_det4 {R : Type*} [CommRing R] (M : Matrix (Fin 4) (Fin 4) R) :
    M.det = det4 (fun b a => M b a) := by
  rw [Matrix.det_succ_row_zero]
  simp [Fin.sum_univ_succ, Matrix.det_fin_three, Fin.succAbove, det4, Matrix.submatrix_apply]
  ring

/-! ### 5. The link to the kernel computation -/

lemma pw_link0 : ∀ a u b : Fin 8,
    ((if (ad a.val + ad u.val) 0 = ad b.val 0 then 1 else -256) *
      (if (ad a.val + ad u.val) 1 = ad b.val 1 then 1 else -256) *
      (if (ad a.val + ad u.val) 2 = ad b.val 2 then 1 else (-256 : ℤ))) =
      pwZ (a.val ^^^ u.val ^^^ b.val) := by decide

lemma pw_link (a u b : ℕ) (ha : a < 8) (hu : u < 8) (hb : b < 8) :
    PrR (256 : ℤ) (ad a + ad u) (ad b) = pwZ (a ^^^ u ^^^ b) := by
  have := pw_link0 ⟨a, ha⟩ ⟨u, hu⟩ ⟨b, hb⟩
  simp only [PrR, Fin.prod_univ_three]
  exact this

lemma chi_link0 : ∀ s u : Fin 8,
    (if (ad s.val 0 * ad u.val 0 + ad s.val 1 * ad u.val 1 + ad s.val 2 * ad u.val 2) = 0
      then (1 : ℤ) else -1) = chiZ s.val u.val := by decide

lemma chi_link (s u : ℕ) (hs : s < 8) (hu : u < 8) : chiR (ad s) (ad u) = (chiZ s u : ℤ) := by
  have := chi_link0 ⟨s, hs⟩ ⟨u, hu⟩
  simp only [chiR, dot, Fin.sum_univ_three]
  exact this

/-- Shape facts about the data (checked by the kernel). -/
lemma data_ok : ∀ s < 8, ∀ r < 28, ∀ b < 4,
    (usD s r).getD b 0 < 8 ∧ (usD s r).getD b 0 ≠ 0 := by decide +kernel

lemma data_ok' : ∀ s < 8, ∀ c < 28, ∀ a < 4, (0 :: tsD s c).getD a 0 < 8 := by decide +kernel

lemma GenR_link (s r c : ℕ) (hs : s < 8) (hr : r < 28) (hc : c < 28) :
    GenR (256 : ℤ) (ad s) (rowT s r) (colT s c) = GenZ s (usD s r) (tsD s c) := by
  unfold GenR GenZ
  rw [sum_S3]
  apply sum8_congr
  intro k hk
  rw [chi_link s k hs hk, rowT_zero, colT_zero, ← ad_zero, pw_link 0 k 0 (by norm_num) hk
    (by norm_num), det_eq_det4]
  simp only [Nat.zero_xor, Nat.xor_zero, Matrix.of_apply, rowT_succ, colT_succ]
  rw [mul_assoc]
  congr 3
  funext b a
  exact pw_link _ _ _ (data_ok s hs r hr b.val b.isLt).1 hk (data_ok' s hs c hc a.val a.isLt)

/-! ### 6. Nonsingularity of the blocks -/

abbrev P : ℕ := 65537

lemma sq256 : (256 : ZMod P) * 256 = ((-1 : ℤ) : ZMod P) := by decide

/-- `ℤ[i] → ZMod 65537`, `i ↦ 256`. -/
def phiP : GaussianInt →+* ZMod P := Zsqrtd.lift ⟨256, sq256⟩

lemma phiP_i : phiP (⟨0, 1⟩ : GaussianInt) = 256 := by simp [phiP]

lemma toC_i : GaussianInt.toComplex (⟨0, 1⟩ : GaussianInt) = I := by
  simp [GaussianInt.toComplex_def']

lemma cast_dotN (f g : ℕ → ℤ) (n : ℕ) :
    ((dotN f g n : ℤ) : ZMod P) = ∑ k : Fin n, ((f k : ℤ) : ZMod P) * ((g k : ℤ) : ZMod P) := by
  induction n with
  | zero => simp [dotN]
  | succ n ih => rw [dotN, Fin.sum_univ_castSucc]; push_cast; rw [ih]; rfl

/-- The certificate gives an inverse mod `P`. -/
lemma block_inv (s : ℕ) (hchk : BlockCheck s) :
    (Matrix.of fun i k : Fin 28 => ((Bf s i k : ℤ) : ZMod P)) *
      (Matrix.of fun k j : Fin 28 => ((Gf s k j : ℤ) : ZMod P)) = 1 := by
  ext i j
  have h := hchk i i.isLt j j.isLt
  rw [Matrix.mul_apply, Matrix.one_apply]
  have h2 : ((dotN (Bf s i) (fun k => Gf s k j) 28 - if (i : ℕ) = j then 1 else 0 : ℤ) : ZMod P)
      = 0 := by
    rw [ZMod.intCast_zmod_eq_zero_iff_dvd]; exact Int.dvd_of_emod_eq_zero h
  rw [Int.cast_sub, cast_dotN, sub_eq_zero] at h2
  simp only [Matrix.of_apply]
  rw [h2]
  by_cases hij : i = j <;> simp [hij, Fin.val_inj]

instance factP : Fact (1 < P) := ⟨by norm_num⟩

lemma block_det (s : ℕ) (hs : s < 8) (hchk : BlockCheck s) :
    (Matrix.of fun r c : Fin 28 => GenR I (ad s) (rowT s r) (colT s c)).det ≠ 0 := by
  set ZG : Matrix (Fin 28) (Fin 28) GaussianInt :=
    Matrix.of fun r c => GenR (⟨0, 1⟩ : GaussianInt) (ad s) (rowT s r) (colT s c) with hZG
  have hC : (Matrix.of fun r c : Fin 28 => GenR I (ad s) (rowT s r) (colT s c)) =
      GaussianInt.toComplex.mapMatrix ZG := by
    ext r c; simp [hZG, RingHom.mapMatrix_apply, GenR_map, toC_i]
  have hP : phiP.mapMatrix ZG = Matrix.of fun k j : Fin 28 => ((Gf s k j : ℤ) : ZMod P) := by
    ext r c
    simp only [hZG, RingHom.mapMatrix_apply, Matrix.map_apply, Matrix.of_apply, GenR_map, phiP_i]
    have := GenR_map (Int.castRingHom (ZMod P)) (256 : ℤ) (ad s) (rowT s r) (colT s c)
    rw [GenR_link s r c hs r.isLt c.isLt] at this
    simp only [eq_intCast, Int.cast_ofNat] at this
    rw [← this]; rfl
  have hdet : (phiP.mapMatrix ZG).det ≠ 0 := by
    rw [hP]
    have := congrArg Matrix.det (block_inv s hchk)
    rw [Matrix.det_mul, Matrix.det_one] at this
    exact right_ne_zero_of_mul_eq_one this
  rw [hC, ← RingHom.map_det]
  intro h0
  apply hdet
  rw [← RingHom.map_det]
  have : ZG.det = 0 := by
    have h1 : GaussianInt.toComplex ZG.det = GaussianInt.toComplex 0 := by rw [h0, map_zero]
    exact GaussianInt.toComplex_inj.mp h1
  rw [this, map_zero]

/-! ### 7. The demand -/

/-- Vector of the hook: a character sum of translates of `Alt δ(column tuple)`. -/
noncomputable def nv (p : Fin 8 × Fin 28) : Sky (Fin 3) 5 → ℂ :=
  ∑ h : S3, sign (dot (ad p.1) h) • matAct 5 (shift h) (Alt 3 (del (colT p.1 p.2)))

lemma nv_mem (p : Fin 8 × Fin 28) : nv p ∈ hookN (α := Fin 3) 3 := by
  apply Submodule.sum_mem; intro h _
  apply Submodule.smul_mem
  apply (hook_kept (α := Fin 3) 3).shift
  exact ⟨_, del_mem_Dsub 3 _ (tup_01 _ _), rfl⟩

/-- The value of the kernel on a vector `nv`, at a translated row tuple. -/
lemma M_nv (p : Fin 8 × Fin 28) (τ : Sky (Fin 3) 5) :
    matAct 5 (kernel (Fin 3)) (nv p) τ =
      ∑ h : S3, sign (dot (ad p.1) h) *
        matAct 5 (kernel (Fin 3)) (Alt 3 (del (colT p.1 p.2))) (tr h τ) := by
  rw [nv, ← matL_apply, map_sum, Finset.sum_apply]
  refine Finset.sum_congr rfl (fun h _ => ?_)
  rw [map_smul, Pi.smul_apply, matL_apply, K3_shift, shift_act, smul_eq_mul]

lemma char_val (s : S3) (τ σ : Sky (Fin 3) 5) :
    ∑ u : S3, sign (dot s u) * matAct 5 (kernel (Fin 3)) (Alt 3 (del σ)) (tr u τ) =
      ca ^ 15 * GenR I s τ σ := by
  rw [GenR, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun u _ => ?_)
  rw [entryF]
  simp only [tr, chiR, sign]
  ring

lemma ca_ne : (ca : ℂ) ≠ 0 := by
  simp [ca, Complex.ext_iff]

lemma rowT_ne (s r : ℕ) (hs : s < 8) (hr : r < 28) (b : Fin 4) :
    rowT s r b.succ ≠ rowT s r 0 := by
  rw [rowT_succ, rowT_zero]
  have h := data_ok s hs r hr b.val b.isLt
  intro he
  have h0 : ad ((usD s r).getD b.val 0) = ad 0 := by rw [he, ad_zero]
  have key : ∀ x : Fin 8, ad x.val 0 = ad 0 0 → ad x.val 1 = ad 0 1 → ad x.val 2 = ad 0 2 →
      x.val = 0 := by decide
  exact h.2 (key ⟨_, h.1⟩ (congrFun h0 0) (congrFun h0 1) (congrFun h0 2))

/-- **The demand of the block-3 hook is at least 224**, given the eight kernel checks. -/
theorem dem3_of_check (hchk : ∀ s < 8, BlockCheck s) :
    224 ≤ raise (hookN (α := Fin 3) 3) (matL 5 (kernel (Fin 3))) := by
  let y : Fin 8 × Fin 28 → (Xf (Fin 3) 5 ⧸ hookN (α := Fin 3) 3) := fun p =>
    (hookN (α := Fin 3) 3).mkQ (matL 5 (kernel (Fin 3)) (nv p))
  have hli : LinearIndependent ℂ y := by
    rw [Fintype.linearIndependent_iff]
    intro a ha
    have hmem : matL 5 (kernel (Fin 3)) (∑ p, a p • nv p) ∈ hookN (α := Fin 3) 3 := by
      rw [← Submodule.Quotient.mk_eq_zero, ← Submodule.mkQ_apply, map_sum, map_sum]
      simpa [y, map_smul] using ha
    -- the functional `Ψ s r` at block `s`, row `r`
    have hblock : ∀ (s : Fin 8) (r : Fin 28),
        ∑ c : Fin 28, (Matrix.of fun r c : Fin 28 => GenR I (ad s) (rowT s r) (colT s c)) r c *
          a (s, c) = 0 := by
      intro s r
      have hz : ∑ g : S3, sign (dot (ad s) g) *
          matL 5 (kernel (Fin 3)) (∑ p, a p • nv p) (tr g (rowT s r)) = 0 := by
        apply Finset.sum_eq_zero; intro g _
        rw [hookN_vanish 3 hmem (tr g (rowT s r)) (fun i h => rowT_ne s r s.isLt r.isLt i
          (add_right_cancel h)), mul_zero]
      simp only [map_sum, map_smul, Finset.sum_apply, Pi.smul_apply, smul_eq_mul, matL_apply,
        M_nv, tr_tr] at hz
      have hz2 : ∑ p : Fin 8 × Fin 28, a p * (if ad s = ad p.1 then
          8 * (ca ^ 15 * GenR I (ad s) (rowT s r) (colT p.1 p.2)) else 0) = 0 := by
        have e : ∑ g : S3, sign (dot (ad s) g) * ∑ p : Fin 8 × Fin 28, a p *
            ∑ h : S3, sign (dot (ad p.1) h) * matAct 5 (kernel (Fin 3))
              (Alt 3 (del (colT p.1 p.2))) (tr (g + h) (rowT s r)) =
            ∑ p : Fin 8 × Fin 28, a p * (if ad s = ad p.1 then
              8 * (ca ^ 15 * GenR I (ad s) (rowT s r) (colT p.1 p.2)) else 0) := by
          calc _ = ∑ p : Fin 8 × Fin 28, a p * ∑ g : S3, sign (dot (ad s) g) *
                ∑ h : S3, sign (dot (ad p.1) h) * matAct 5 (kernel (Fin 3))
                  (Alt 3 (del (colT p.1 p.2))) (tr (g + h) (rowT s r)) := by
                simp_rw [Finset.mul_sum]
                rw [Finset.sum_comm]
                exact Finset.sum_congr rfl (fun p _ => Finset.sum_congr rfl
                  (fun g _ => Finset.sum_congr rfl (fun h _ => by ring)))
            _ = _ := by
                refine Finset.sum_congr rfl (fun p _ => ?_)
                rw [← char_val]
                exact congrArg (fun t => a p * t) (dbl_sum (ad s) (ad p.1)
                  (fun w => matAct 5 (kernel (Fin 3)) (Alt 3 (del (colT p.1 p.2)))
                    (tr w (rowT s r))))
        rw [← e]; exact hz
      rw [Fintype.sum_prod_type] at hz2
      have hs' : ∀ s' : Fin 8, (ad s = ad s'.val ↔ s = s') := fun s' =>
        ⟨fun h => ad_bij.1 h, fun h => by rw [h]⟩
      simp only [hs'] at hz2
      rw [Finset.sum_eq_single s (fun s' _ hne => by
          simp [Ne.symm hne]) (fun h => absurd (Finset.mem_univ s) h)] at hz2
      simp only [ite_true] at hz2
      have hk : (8 * ca ^ 15 : ℂ) ≠ 0 := mul_ne_zero (by norm_num) (pow_ne_zero _ ca_ne)
      have : (8 * ca ^ 15) * ∑ c : Fin 28,
          GenR I (ad s) (rowT s r) (colT s c) * a (s, c) = 0 := by
        rw [Finset.mul_sum, ← hz2]
        refine Finset.sum_congr rfl (fun c _ => by ring)
      simpa [Matrix.of_apply, hk] using this
    intro p
    obtain ⟨s, c⟩ := p
    have hv : (Matrix.of fun r c : Fin 28 => GenR I (ad s) (rowT s r) (colT s c)) *ᵥ
        (fun c => a (s, c)) = 0 := by
      funext r; simp only [Matrix.mulVec, dotProduct, Pi.zero_apply]; exact hblock s r
    have := Matrix.eq_zero_of_mulVec_eq_zero (block_det s s.isLt (hchk s s.isLt)) hv
    exact congrFun this c
  have hspan : span ℂ (Set.range y) ≤
      (hookN (α := Fin 3) 3).map ((hookN (α := Fin 3) 3).mkQ ∘ₗ matL 5 (kernel (Fin 3))) := by
    rw [Submodule.span_le]; rintro _ ⟨p, rfl⟩; exact ⟨_, nv_mem p, rfl⟩
  have h := Submodule.finrank_mono hspan
  rw [finrank_span_eq_card hli] at h
  have hc : Fintype.card (Fin 8 × Fin 28) = 224 := by simp
  rw [hc] at h
  exact h

end Dem3

end PowerSaving.PfLean
end OAI

#print axioms OAI.PowerSaving.PfLean.Dem3.dem3_of_check
