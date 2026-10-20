import Work.GCert.Chain.Phased
import XMut.Dirt

/-!
# x-mutable (13 Oct 2026): the scalar half WITHOUT read-only sources

`ScalZ` = Jacob Sussman's `GX.Scal` (`Work/GCert/Chain/Phased.lean`) with the field `hx` (the x rows of the
product are unit rows: no gate writes a source) REMOVED.  The forward word ends with the derived gate `g8z`
(x AND slot rows of the inverse) instead of `g8` (slot rows only); the backward word with `h8z`.  The two scalar
identities are word for word those of `Scal.fwd_scalar` / `Scal.bwd_scalar`: `y += x` (forward), `x -= y`
(backward), slots and sources restored, for ARBITRARY slot contents.

The proofs are Sussman's, with `XMut.hFz` / `XMut.hBz` for `Scal.hF` / `Scal.hB`.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace GX
open Binary Matrix Finset RAM SS CB BR CR GF
noncomputable section

/-- **The scalar half, sources mutable** (`Scal` without `hx`). -/
structure ScalZ (T Sl C : Type) [Fintype T] [DecidableEq T]
    [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] where
  Cc : Matrix C Sl ℚ
  Jr : Matrix T C ℚ
  MA : Matrix (L3 T Sl) (L3 T Sl) ℚ
  MAi : Matrix (L3 T Sl) (L3 T Sl) ℚ
  MB : Matrix (L3 T Sl) (L3 T Sl) ℚ
  MBi : Matrix (L3 T Sl) (L3 T Sl) ℚ
  hAi : MAi * MA = 1
  hBi : MBi * MB = 1
  hy : ∀ i t, (MB * scatM (Jr * Cc) * MA) i (lY t) = if i = lY t then 1 else 0
  hid : ∀ S t, (MB * scatM (Jr * Cc) * MA) (lY S) (lX t) = if S = t then 1 else 0

/-- a classic scalar half is a mutable-source one (forget `hx`). -/
def Scal.toZ {T Sl C : Type} [Fintype T] [DecidableEq T]
    [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] (K : Scal T Sl C) : ScalZ T Sl C :=
  ⟨K.Cc, K.Jr, K.MA, K.MAi, K.MB, K.MBi, K.hAi, K.hBi, K.hy, K.hid⟩

section Alg
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]

namespace ScalZ
variable (K : ScalZ T Sl C)

def Mt : Matrix (L3 T Sl) (L3 T Sl) ℚ := K.MB * scatM (K.Jr * K.Cc) * K.MA
def Mti : Matrix (L3 T Sl) (L3 T Sl) ℚ := K.MAi * scatMi (K.Jr * K.Cc) * K.MBi

lemma im : K.Mti * K.Mt = 1 := by
  unfold Mti Mt
  have hB : ∀ X : Matrix (L3 T Sl) (L3 T Sl) ℚ, K.MBi * (K.MB * X) = X := fun X => by
    rw [← Matrix.mul_assoc, K.hBi, Matrix.one_mul]
  have hS : ∀ X : Matrix (L3 T Sl) (L3 T Sl) ℚ,
      scatMi (K.Jr * K.Cc) * (scatM (K.Jr * K.Cc) * X) = X := fun X => by
    rw [← Matrix.mul_assoc, scat_inv, Matrix.one_mul]
  rw [Matrix.mul_assoc, Matrix.mul_assoc, Matrix.mul_assoc, hB, hS, K.hAi]

lemma mi : K.Mt * K.Mti = 1 := mul_eq_one_comm.mp K.im

lemma hmy : K.Mt * PY = PY := mul_PY_of (M := K.Mt) K.hy
lemma hiy : K.Mti * PY = PY := dirt_iy K.hmy K.im

lemma hF : XMut.g8z K.Mti * K.Mt * g1 K.Mt = 1 + PY * K.Mt * PX := XMut.hFz K.hmy K.im
lemma hB : g1i K.Mt * K.Mti * XMut.h8z K.Mt = 1 - PY * K.Mt * PX := XMut.hBz K.hmy K.im

/-- live part of the forward word before the scatter: `y -= K s`, then phase A. -/
def f1 (u : L3 T Sl → ℂ) : L3 T Sl → ℂ := actPoint K.MA (actPoint (g1 K.Mt) u)
/-- after the scatter: phase B, then the x AND slot rows of the inverse. -/
def f2 (u : L3 T Sl → ℂ) : L3 T Sl → ℂ := actPoint (XMut.g8z K.Mti) (actPoint K.MB u)
def b1 (u : L3 T Sl → ℂ) : L3 T Sl → ℂ := actPoint K.MAiᵀ (actPoint (g1i K.Mt)ᵀ u)
def b2 (u : L3 T Sl → ℂ) : L3 T Sl → ℂ := actPoint (XMut.h8z K.Mt)ᵀ (actPoint K.MBiᵀ u)

lemma fwd_live (u : L3 T Sl → ℂ) :
    K.f2 (actPoint (scatM (K.Jr * K.Cc)) (K.f1 u)) = actPoint (1 + PY * K.Mt * PX) u := by
  rw [← K.hF]
  show actPoint (XMut.g8z K.Mti) (actPoint K.MB (actPoint (scatM (K.Jr * K.Cc))
      (actPoint K.MA (actPoint (g1 K.Mt) u))))
    = actPoint (XMut.g8z K.Mti * (K.MB * scatM (K.Jr * K.Cc) * K.MA) * g1 K.Mt) u
  rw [actPoint_mul, actPoint_mul, actPoint_mul, actPoint_mul]

lemma bwd_live (u : L3 T Sl → ℂ) :
    K.b2 (actPoint ((scatMi (K.Jr * K.Cc))ᵀ) (K.b1 u))
      = actPoint ((1 - PY * K.Mt * PX)ᵀ) u := by
  rw [← K.hB]
  show actPoint (XMut.h8z K.Mt)ᵀ (actPoint K.MBiᵀ (actPoint (scatMi (K.Jr * K.Cc))ᵀ
      (actPoint K.MAiᵀ (actPoint (g1i K.Mt)ᵀ u))))
    = actPoint ((g1i K.Mt * (K.MAi * scatMi (K.Jr * K.Cc) * K.MBi) * XMut.h8z K.Mt)ᵀ) u
  rw [Matrix.transpose_mul, Matrix.transpose_mul, Matrix.transpose_mul, Matrix.transpose_mul,
    actPoint_mul, actPoint_mul, actPoint_mul, actPoint_mul]

lemma liveAct_glueZ (φ : (L3 T Sl → ℂ) → (L3 T Sl → ℂ)) (u : L3 T Sl → ℂ) (c : C → ℂ) :
    liveAct φ (glue u c) = glue (φ u) c :=
  congrArg (fun v => glue (φ v) c) (glue_up u c)

/-- **Scalar identity of the forward word, sources mutable**: `y += x`, sources and slots restored. -/
theorem fwd_scalar (z : Box T Sl C → ℂ) :
    eraseC (liveAct K.f2 (gateS bY bC K.Jr (gateS bC bS K.Cc (liveAct K.f1 (eraseC z)))))
      = shearF z := by
  have h1 : liveAct K.f1 (eraseC z) = glue (K.f1 (z ∘ upE)) (fun _ => 0) := by
    rw [erase_glue z, liveAct_glueZ]
  have hs : gateS bY bC K.Jr (gateS bC bS K.Cc (glue (K.f1 (z ∘ upE)) (fun _ => 0))) ∘ upE
      = actPoint (scatM (K.Jr * K.Cc)) (K.f1 (z ∘ upE)) := by
    have h := scat_glue K.Jr K.Cc (glue (K.f1 (z ∘ upE)) (fun _ => 0)) rfl
    rw [glue_up] at h
    exact h
  rw [h1]
  unfold liveAct
  rw [erase_of_glue, hs, K.fwd_live, act_shear K.Mt K.hid]
  exact glue_shearF z

/-- **Scalar identity of the backward word, sources mutable**: `x -= y`, sources and slots restored. -/
theorem bwd_scalar (z : Box T Sl C → ℂ) :
    eraseC (liveAct K.b2 (gateS bS bC (-K.Ccᵀ) (gateS bC bY K.Jrᵀ (liveAct K.b1 (eraseC z)))))
      = shearB z := by
  have h1 : liveAct K.b1 (eraseC z) = glue (K.b1 (z ∘ upE)) (fun _ => 0) := by
    rw [erase_glue z, liveAct_glueZ]
  have hs : gateS bS bC (-K.Ccᵀ) (gateS bC bY K.Jrᵀ (glue (K.b1 (z ∘ upE)) (fun _ => 0))) ∘ upE
      = actPoint ((scatMi (K.Jr * K.Cc))ᵀ) (K.b1 (z ∘ upE)) := by
    have h := scat_glueT K.Jr K.Cc (glue (K.b1 (z ∘ upE)) (fun _ => 0)) rfl
    rw [glue_up] at h
    exact h
  rw [h1]
  unfold liveAct
  rw [erase_of_glue, hs, K.bwd_live, act_shearT K.Mt K.hid]
  exact glue_shearB z

end ScalZ
end Alg
end
end GX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.GX.ScalZ.fwd_scalar
#print axioms OAI.PowerSaving.GX.ScalZ.bwd_scalar
