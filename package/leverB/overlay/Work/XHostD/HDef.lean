import Work.XHost.HRoles
import Work.XHostD.Prov

/-!
# (key: bridge-host) The nine-stage word of PR #401 with SAME-RUN HOSTING: data and frames

Copy of Work/Bridge401/GNetDef.lean (Ryan Shea's complex588 tree; `G401` by Nil Mamano, written by his gen401.py, on
Dugongue #401 / jacklightChen #377 / Lydakis #275 / DaysSky #363; engine and E8 unit by Jacob Sussman) with idle-axis
hosting (DaysSky, CrocSwap/integer-mult-bounds PR #442): in stage `l` the slot `q` goes to the idle data role
`hostAt l q = some (c, s, t)` (window kind `c`, side `s`, class `(g, t)` of the run's replica; the same table in every
run) or, when `hostAt l q = none`, to the bank `bk ρ l q`.  The idle climbs A2..A4, B2..B4 and the finals F1..F3 become
host chains `Hb`/`Hp`: free adapters, the hosted residuals, ONE remaining block of rank `hrem c j` at the old place.
Written by h-L3-build/gen/gen_hnet.py (wht43, 13 Oct 2026), prepared with assistance from Claude.  PER-RUN version (one host table per run) written by fsp-perrun/tools/mk_perrun.py (wht43, 14 Oct 2026) from h-L3-build's file by checked text substitution; prepared with assistance from Claude.  TWO-PROVIDER version (a per-call invocation `Ic ρ l` with the target vectors of `I`, per-call entrance labels and residuals) written by lg-a2a4/tools/mk_disp.py (wht43, 14 Oct 2026) from fsp-perrun's file by checked text substitution; prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace XHostD
open Binary Matrix Finset RAM SS CB RF BR GF GX XHost
noncomputable section

structure GHost401 (H T Sl C α Γ : Type) (nb : ℕ) [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
    [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α] where
  /-- the gauged invocation (`GInvE`). -/
  I : GInvE H T Sl C
  /-- the invocation of call `(ρ, l)` (its provider word), with the target vectors of `I`. -/
  Ic : Fin 8 → Fin 9 → GInvE H T Sl C
  hIc : ∀ ρ l, (Ic ρ l).tv = I.tv
  G0 : Γ → GWin H α
  G1 : Γ → GWin H α
  G2 : Γ → GWin H α
  G3 : Γ → GWin H α
  G4 : Γ → GWin H α
  G5 : Γ → GWin H α
  G6 : Γ → GWin H α
  G7 : Γ → GWin H α
  G8 : Γ → GWin H α
  r1 : T → Γ ≃ Γ
  r2 : T → Γ ≃ Γ
  r3 : T → Γ ≃ Γ
  r4 : T → Γ ≃ Γ
  r5 : T → Γ ≃ Γ
  r6 : T → Γ ≃ Γ
  r7 : T → Γ ≃ Γ
  r8 : T → Γ ≃ Γ
  N0 : Fin 8 → Γ → Sl → CMat α
  Wn0 : Fin 8 → Γ → Sl → CMat α
  N1 : Fin 8 → Γ → Sl → CMat α
  Wn1 : Fin 8 → Γ → Sl → CMat α
  N2 : Fin 8 → Γ → Sl → CMat α
  Wn2 : Fin 8 → Γ → Sl → CMat α
  N3 : Fin 8 → Γ → Sl → CMat α
  Wn3 : Fin 8 → Γ → Sl → CMat α
  N4 : Fin 8 → Γ → Sl → CMat α
  Wn4 : Fin 8 → Γ → Sl → CMat α
  N5 : Fin 8 → Γ → Sl → CMat α
  Wn5 : Fin 8 → Γ → Sl → CMat α
  N6 : Fin 8 → Γ → Sl → CMat α
  Wn6 : Fin 8 → Γ → Sl → CMat α
  N7 : Fin 8 → Γ → Sl → CMat α
  Wn7 : Fin 8 → Γ → Sl → CMat α
  N8 : Fin 8 → Γ → Sl → CMat α
  Wn8 : Fin 8 → Γ → Sl → CMat α
  /-- the host of slot `q` in stage `l`, the SAME in every run: (window kind, side, class); `none`: a bank. -/
  hostAt : Fin 8 → Fin 9 → Sl → Option (Fin 9 × (Bool × T))
  hostInj : ∀ ρ l q q' j, hostAt ρ l q = some j → hostAt ρ l q' = some j → q = q'
  hostWin : ∀ ρ l q c j, hostAt ρ l q = some (c, j) → inWin c.val l.val = true
  /-- the bank of a NON-hosted slot `q` in stage `l` of run `ρ` (junk on hosted slots). -/
  bk : Fin 8 → Fin 9 → Sl → Fin nb
  hinj : ∀ ρ l q q', hostAt ρ l q = none → hostAt ρ l q' = none → bk ρ l q = bk ρ l q' → q = q'
  Zb : ℕ → Γ → Fin nb → CMat α
  Zp : ℕ → Γ → Fin nb → CMat α
  pad : Fin nb → ℕ
  /-- the host chains: frame of host `j` of kind `c` before (`Hb`) and after (`Hp`) the adapter of stage `k`. -/
  Hb : Fin 8 → Fin 9 → ℕ → Γ → Bool × T → CMat α
  Hp : Fin 8 → Fin 9 → ℕ → Γ → Bool × T → CMat α
  hrem : Fin 8 → Fin 9 → Bool × T → ℕ
  u : Γ → T → Space α
  YF : Γ → T → CMat α
  hm : 2 ≤ Fintype.card α
  cX4 : (ℕ → ℝ) → ℝ
  cY4 : (ℕ → ℝ) → ℝ
  inv0 : ∀ g t, ∃ N : CMat α, (G0 g).Xm.Φ (tt (I.tv t)) * N = 1
  z0 : ∀ g, (G0 g).Xm.Φ 0 = 1
  hN0 : ∀ (ρ : Fin 8) d q, (G0 d).X.Φ ((Ic ρ 0).ent q) * N0 ρ d q = 1
  hW0 : ∀ (ρ : Fin 8) d q, (G0 d).X.Φ ((Ic ρ 0).ent q).top = Wn0 ρ d q * (G0 d).X.Φ ((Ic ρ 0).ent q)
  hN1 : ∀ (ρ : Fin 8) d q, (G1 d).X.Φ ((Ic ρ 1).ent q) * N1 ρ d q = 1
  hW1 : ∀ (ρ : Fin 8) d q, (G1 d).X.Φ ((Ic ρ 1).ent q).top = Wn1 ρ d q * (G1 d).X.Φ ((Ic ρ 1).ent q)
  hN2 : ∀ (ρ : Fin 8) d q, (G2 d).X.Φ ((Ic ρ 2).ent q) * N2 ρ d q = 1
  hW2 : ∀ (ρ : Fin 8) d q, (G2 d).X.Φ ((Ic ρ 2).ent q).top = Wn2 ρ d q * (G2 d).X.Φ ((Ic ρ 2).ent q)
  hN3 : ∀ (ρ : Fin 8) d q, (G3 d).X.Φ ((Ic ρ 3).ent q) * N3 ρ d q = 1
  hW3 : ∀ (ρ : Fin 8) d q, (G3 d).X.Φ ((Ic ρ 3).ent q).top = Wn3 ρ d q * (G3 d).X.Φ ((Ic ρ 3).ent q)
  hN4 : ∀ (ρ : Fin 8) d q, (G4 d).X.Φ ((Ic ρ 4).ent q) * N4 ρ d q = 1
  hW4 : ∀ (ρ : Fin 8) d q, (G4 d).X.Φ ((Ic ρ 4).ent q).top = Wn4 ρ d q * (G4 d).X.Φ ((Ic ρ 4).ent q)
  hN5 : ∀ (ρ : Fin 8) d q, (G5 d).X.Φ ((Ic ρ 5).ent q) * N5 ρ d q = 1
  hW5 : ∀ (ρ : Fin 8) d q, (G5 d).X.Φ ((Ic ρ 5).ent q).top = Wn5 ρ d q * (G5 d).X.Φ ((Ic ρ 5).ent q)
  hN6 : ∀ (ρ : Fin 8) d q, (G6 d).X.Φ ((Ic ρ 6).ent q) * N6 ρ d q = 1
  hW6 : ∀ (ρ : Fin 8) d q, (G6 d).X.Φ ((Ic ρ 6).ent q).top = Wn6 ρ d q * (G6 d).X.Φ ((Ic ρ 6).ent q)
  hN7 : ∀ (ρ : Fin 8) d q, (G7 d).X.Φ ((Ic ρ 7).ent q) * N7 ρ d q = 1
  hW7 : ∀ (ρ : Fin 8) d q, (G7 d).X.Φ ((Ic ρ 7).ent q).top = Wn7 ρ d q * (G7 d).X.Φ ((Ic ρ 7).ent q)
  hN8 : ∀ (ρ : Fin 8) d q, (G8 d).X.Φ ((Ic ρ 8).ent q) * N8 ρ d q = 1
  hW8 : ∀ (ρ : Fin 8) d q, (G8 d).X.Φ ((Ic ρ 8).ent q).top = Wn8 ρ d q * (G8 d).X.Φ ((Ic ρ 8).ent q)
  x01 : ∀ g t, (G0 g).Xm.Φ 1 = (G1 (r1 t g)).Xm.Φ (tt (I.tv t))
  y01 : ∀ g t, (G0 g).Xm.Φ (1 + tt (I.tv t)) = (G1 (r1 t g)).Xm.Φ 0
  x12 : ∀ g t, (G1 (r1 t g)).Xm.Φ 1 = (G2 (r2 t g)).Xm.Φ (tt (I.tv t))
  y12 : ∀ g t, (G1 (r1 t g)).Xm.Φ (1 + tt (I.tv t)) = (G2 (r2 t g)).Xm.Φ 0
  x23 : ∀ g t, (G2 (r2 t g)).Xm.Φ 1 = (G3 (r3 t g)).Xm.Φ (tt (I.tv t))
  y23 : ∀ g t, (G2 (r2 t g)).Xm.Φ (1 + tt (I.tv t)) = (G3 (r3 t g)).Xm.Φ 0
  x34 : ∀ g t, (G3 (r3 t g)).Xm.Φ 1 = (G4 (r4 t g)).Xm.Φ (tt (I.tv t))
  y34 : ∀ g t, (G3 (r3 t g)).Xm.Φ (1 + tt (I.tv t)) = (G4 (r4 t g)).Xm.Φ 0
  x45 : ∀ g t, (G4 (r4 t g)).Xm.Φ 1 = (G5 (r5 t g)).Xm.Φ (tt (I.tv t))
  y45 : ∀ g t, (G4 (r4 t g)).Xm.Φ (1 + tt (I.tv t)) = (G5 (r5 t g)).Xm.Φ 0
  x56 : ∀ g t, (G5 (r5 t g)).Xm.Φ 1 = (G6 (r6 t g)).Xm.Φ (tt (I.tv t))
  y56 : ∀ g t, (G5 (r5 t g)).Xm.Φ (1 + tt (I.tv t)) = (G6 (r6 t g)).Xm.Φ 0
  x67 : ∀ g t, (G6 (r6 t g)).Xm.Φ 1 = (G7 (r7 t g)).Xm.Φ (tt (I.tv t))
  y67 : ∀ g t, (G6 (r6 t g)).Xm.Φ (1 + tt (I.tv t)) = (G7 (r7 t g)).Xm.Φ 0
  x78 : ∀ g t, (G7 (r7 t g)).Xm.Φ 1 = (G8 (r8 t g)).Xm.Φ (tt (I.tv t))
  y78 : ∀ g t, (G7 (r7 t g)).Xm.Φ (1 + tt (I.tv t)) = (G8 (r8 t g)).Xm.Φ 0
  /-- the banks (non-hosted slots only). -/
  z0b : ∀ g b, Zb 0 g b = 1
  hP0 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 0) g b) (Zp (9 * ρ.val + 0) g b) (fun _ => 0)
  hZ0 : ∀ (ρ : Fin 8) g q, hostAt ρ 0 q = none → Wn0 ρ g q * Zp (9 * ρ.val + 0) g (bk ρ 0 q) = Zb (9 * ρ.val + 1) g (bk ρ 0 q)
  hZn0 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt ρ 0 q = none → bk ρ 0 q ≠ b) → Zb (9 * ρ.val + 1) g b = Zp (9 * ρ.val + 0) g b
  hP1 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 1) g b) (Zp (9 * ρ.val + 1) g b) (fun _ => 0)
  hZ1 : ∀ (ρ : Fin 8) g q, hostAt ρ 1 q = none → Wn1 ρ g q * Zp (9 * ρ.val + 1) g (bk ρ 1 q) = Zb (9 * ρ.val + 2) g (bk ρ 1 q)
  hZn1 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt ρ 1 q = none → bk ρ 1 q ≠ b) → Zb (9 * ρ.val + 2) g b = Zp (9 * ρ.val + 1) g b
  hP2 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 2) g b) (Zp (9 * ρ.val + 2) g b) (fun _ => 0)
  hZ2 : ∀ (ρ : Fin 8) g q, hostAt ρ 2 q = none → Wn2 ρ g q * Zp (9 * ρ.val + 2) g (bk ρ 2 q) = Zb (9 * ρ.val + 3) g (bk ρ 2 q)
  hZn2 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt ρ 2 q = none → bk ρ 2 q ≠ b) → Zb (9 * ρ.val + 3) g b = Zp (9 * ρ.val + 2) g b
  hP3 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 3) g b) (Zp (9 * ρ.val + 3) g b) (fun _ => 0)
  hZ3 : ∀ (ρ : Fin 8) g q, hostAt ρ 3 q = none → Wn3 ρ g q * Zp (9 * ρ.val + 3) g (bk ρ 3 q) = Zb (9 * ρ.val + 4) g (bk ρ 3 q)
  hZn3 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt ρ 3 q = none → bk ρ 3 q ≠ b) → Zb (9 * ρ.val + 4) g b = Zp (9 * ρ.val + 3) g b
  hP4 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 4) g b) (Zp (9 * ρ.val + 4) g b) (fun _ => 0)
  hZ4 : ∀ (ρ : Fin 8) g q, hostAt ρ 4 q = none → Wn4 ρ g q * Zp (9 * ρ.val + 4) g (bk ρ 4 q) = Zb (9 * ρ.val + 5) g (bk ρ 4 q)
  hZn4 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt ρ 4 q = none → bk ρ 4 q ≠ b) → Zb (9 * ρ.val + 5) g b = Zp (9 * ρ.val + 4) g b
  hP5 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 5) g b) (Zp (9 * ρ.val + 5) g b) (fun _ => 0)
  hZ5 : ∀ (ρ : Fin 8) g q, hostAt ρ 5 q = none → Wn5 ρ g q * Zp (9 * ρ.val + 5) g (bk ρ 5 q) = Zb (9 * ρ.val + 6) g (bk ρ 5 q)
  hZn5 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt ρ 5 q = none → bk ρ 5 q ≠ b) → Zb (9 * ρ.val + 6) g b = Zp (9 * ρ.val + 5) g b
  hP6 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 6) g b) (Zp (9 * ρ.val + 6) g b) (fun _ => 0)
  hZ6 : ∀ (ρ : Fin 8) g q, hostAt ρ 6 q = none → Wn6 ρ g q * Zp (9 * ρ.val + 6) g (bk ρ 6 q) = Zb (9 * ρ.val + 7) g (bk ρ 6 q)
  hZn6 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt ρ 6 q = none → bk ρ 6 q ≠ b) → Zb (9 * ρ.val + 7) g b = Zp (9 * ρ.val + 6) g b
  hP7 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 7) g b) (Zp (9 * ρ.val + 7) g b) (fun _ => 0)
  hZ7 : ∀ (ρ : Fin 8) g q, hostAt ρ 7 q = none → Wn7 ρ g q * Zp (9 * ρ.val + 7) g (bk ρ 7 q) = Zb (9 * ρ.val + 8) g (bk ρ 7 q)
  hZn7 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt ρ 7 q = none → bk ρ 7 q ≠ b) → Zb (9 * ρ.val + 8) g b = Zp (9 * ρ.val + 7) g b
  hP8 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 8) g b) (Zp (9 * ρ.val + 8) g b) (fun _ => 0)
  hZ8 : ∀ (ρ : Fin 8) g q, hostAt ρ 8 q = none → Wn8 ρ g q * Zp (9 * ρ.val + 8) g (bk ρ 8 q) = Zb (9 * ρ.val + 9) g (bk ρ 8 q)
  hZn8 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt ρ 8 q = none → bk ρ 8 q ≠ b) → Zb (9 * ρ.val + 9) g b = Zp (9 * ρ.val + 8) g b
  bpad : ∀ g b, GReach (Zb 72 g b) (kernel α) (fun φ => bcost φ (pad b))
  /-- the host chains (XHost.host_frames): start, free adapters, hosted residuals, ONE remaining block. -/
  hH0_A2 : ∀ (ρ : Fin 8) g t, Hb ρ 0 0 g (false, t) = (G0 g).Xm.Φ (tt (I.tv t)) ∧ Hb ρ 0 0 g (true, t) = (G0 g).Xm.Φ 0
  hHP_A2_0 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 0 0 g j) (Hp ρ 0 0 g j) (fun _ => 0)
  hHZ_A2_0 : ∀ (ρ : Fin 8) g q j, hostAt ρ 0 q = some (0, j) → Wn0 ρ g q * Hp ρ 0 0 g j = Hb ρ 0 1 g j
  hHZn_A2_0 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 0 q ≠ some (0, j)) → Hb ρ 0 1 g j = Hp ρ 0 0 g j
  hHP_A2_1 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 0 1 g j) (Hp ρ 0 1 g j) (fun _ => 0)
  hHZ_A2_1 : ∀ (ρ : Fin 8) g q j, hostAt ρ 1 q = some (0, j) → Wn1 ρ g q * Hp ρ 0 1 g j = Hb ρ 0 2 g j
  hHZn_A2_1 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 1 q ≠ some (0, j)) → Hb ρ 0 2 g j = Hp ρ 0 1 g j
  hHE_A2 : ∀ (ρ : Fin 8) g t, GReach (Hb ρ 0 2 g (false, t)) ((G1 (r1 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem ρ 0 (false, t))) ∧
    GReach (Hb ρ 0 2 g (true, t)) ((G1 (r1 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem ρ 0 (true, t)))
  hH0_B2 : ∀ (ρ : Fin 8) g t, Hb ρ 1 2 g (false, t) = (G2 (r2 t g)).Xm.Φ (tt (I.tv t)) ∧ Hb ρ 1 2 g (true, t) = (G2 (r2 t g)).Xm.Φ 0
  hHP_B2_2 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 1 2 g j) (Hp ρ 1 2 g j) (fun _ => 0)
  hHZ_B2_2 : ∀ (ρ : Fin 8) g q j, hostAt ρ 2 q = some (1, j) → Wn2 ρ g q * Hp ρ 1 2 g j = Hb ρ 1 3 g j
  hHZn_B2_2 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 2 q ≠ some (1, j)) → Hb ρ 1 3 g j = Hp ρ 1 2 g j
  hHE_B2 : ∀ (ρ : Fin 8) g t, GReach (Hb ρ 1 3 g (false, t)) ((G2 (r2 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem ρ 1 (false, t))) ∧
    GReach (Hb ρ 1 3 g (true, t)) ((G2 (r2 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem ρ 1 (true, t)))
  hH0_A3 : ∀ (ρ : Fin 8) g t, Hb ρ 2 0 g (false, t) = (G0 g).Xm.Φ (tt (I.tv t)) ∧ Hb ρ 2 0 g (true, t) = (G0 g).Xm.Φ 0
  hHP_A3_0 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 2 0 g j) (Hp ρ 2 0 g j) (fun _ => 0)
  hHZ_A3_0 : ∀ (ρ : Fin 8) g q j, hostAt ρ 0 q = some (2, j) → Wn0 ρ g q * Hp ρ 2 0 g j = Hb ρ 2 1 g j
  hHZn_A3_0 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 0 q ≠ some (2, j)) → Hb ρ 2 1 g j = Hp ρ 2 0 g j
  hHP_A3_1 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 2 1 g j) (Hp ρ 2 1 g j) (fun _ => 0)
  hHZ_A3_1 : ∀ (ρ : Fin 8) g q j, hostAt ρ 1 q = some (2, j) → Wn1 ρ g q * Hp ρ 2 1 g j = Hb ρ 2 2 g j
  hHZn_A3_1 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 1 q ≠ some (2, j)) → Hb ρ 2 2 g j = Hp ρ 2 1 g j
  hHP_A3_2 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 2 2 g j) (Hp ρ 2 2 g j) (fun _ => 0)
  hHZ_A3_2 : ∀ (ρ : Fin 8) g q j, hostAt ρ 2 q = some (2, j) → Wn2 ρ g q * Hp ρ 2 2 g j = Hb ρ 2 3 g j
  hHZn_A3_2 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 2 q ≠ some (2, j)) → Hb ρ 2 3 g j = Hp ρ 2 2 g j
  hHP_A3_3 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 2 3 g j) (Hp ρ 2 3 g j) (fun _ => 0)
  hHZ_A3_3 : ∀ (ρ : Fin 8) g q j, hostAt ρ 3 q = some (2, j) → Wn3 ρ g q * Hp ρ 2 3 g j = Hb ρ 2 4 g j
  hHZn_A3_3 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 3 q ≠ some (2, j)) → Hb ρ 2 4 g j = Hp ρ 2 3 g j
  hHE_A3 : ∀ (ρ : Fin 8) g t, GReach (Hb ρ 2 4 g (false, t)) ((G3 (r3 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem ρ 2 (false, t))) ∧
    GReach (Hb ρ 2 4 g (true, t)) ((G3 (r3 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem ρ 2 (true, t)))
  hH0_B3 : ∀ (ρ : Fin 8) g t, Hb ρ 3 4 g (false, t) = (G4 (r4 t g)).Xm.Φ (tt (I.tv t)) ∧ Hb ρ 3 4 g (true, t) = (G4 (r4 t g)).Xm.Φ 0
  hHP_B3_4 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 3 4 g j) (Hp ρ 3 4 g j) (fun _ => 0)
  hHZ_B3_4 : ∀ (ρ : Fin 8) g q j, hostAt ρ 4 q = some (3, j) → Wn4 ρ g q * Hp ρ 3 4 g j = Hb ρ 3 5 g j
  hHZn_B3_4 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 4 q ≠ some (3, j)) → Hb ρ 3 5 g j = Hp ρ 3 4 g j
  hHE_B3 : ∀ (ρ : Fin 8) g t, GReach (Hb ρ 3 5 g (false, t)) ((G4 (r4 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem ρ 3 (false, t))) ∧
    GReach (Hb ρ 3 5 g (true, t)) ((G4 (r4 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem ρ 3 (true, t)))
  hH0_A4 : ∀ (ρ : Fin 8) g t, Hb ρ 4 0 g (false, t) = (G0 g).Xm.Φ (tt (I.tv t)) ∧ Hb ρ 4 0 g (true, t) = (G0 g).Xm.Φ 0
  hHP_A4_0 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 4 0 g j) (Hp ρ 4 0 g j) (fun _ => 0)
  hHZ_A4_0 : ∀ (ρ : Fin 8) g q j, hostAt ρ 0 q = some (4, j) → Wn0 ρ g q * Hp ρ 4 0 g j = Hb ρ 4 1 g j
  hHZn_A4_0 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 0 q ≠ some (4, j)) → Hb ρ 4 1 g j = Hp ρ 4 0 g j
  hHP_A4_1 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 4 1 g j) (Hp ρ 4 1 g j) (fun _ => 0)
  hHZ_A4_1 : ∀ (ρ : Fin 8) g q j, hostAt ρ 1 q = some (4, j) → Wn1 ρ g q * Hp ρ 4 1 g j = Hb ρ 4 2 g j
  hHZn_A4_1 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 1 q ≠ some (4, j)) → Hb ρ 4 2 g j = Hp ρ 4 1 g j
  hHP_A4_2 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 4 2 g j) (Hp ρ 4 2 g j) (fun _ => 0)
  hHZ_A4_2 : ∀ (ρ : Fin 8) g q j, hostAt ρ 2 q = some (4, j) → Wn2 ρ g q * Hp ρ 4 2 g j = Hb ρ 4 3 g j
  hHZn_A4_2 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 2 q ≠ some (4, j)) → Hb ρ 4 3 g j = Hp ρ 4 2 g j
  hHP_A4_3 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 4 3 g j) (Hp ρ 4 3 g j) (fun _ => 0)
  hHZ_A4_3 : ∀ (ρ : Fin 8) g q j, hostAt ρ 3 q = some (4, j) → Wn3 ρ g q * Hp ρ 4 3 g j = Hb ρ 4 4 g j
  hHZn_A4_3 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 3 q ≠ some (4, j)) → Hb ρ 4 4 g j = Hp ρ 4 3 g j
  hHP_A4_4 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 4 4 g j) (Hp ρ 4 4 g j) (fun _ => 0)
  hHZ_A4_4 : ∀ (ρ : Fin 8) g q j, hostAt ρ 4 q = some (4, j) → Wn4 ρ g q * Hp ρ 4 4 g j = Hb ρ 4 5 g j
  hHZn_A4_4 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 4 q ≠ some (4, j)) → Hb ρ 4 5 g j = Hp ρ 4 4 g j
  hHP_A4_5 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 4 5 g j) (Hp ρ 4 5 g j) (fun _ => 0)
  hHZ_A4_5 : ∀ (ρ : Fin 8) g q j, hostAt ρ 5 q = some (4, j) → Wn5 ρ g q * Hp ρ 4 5 g j = Hb ρ 4 6 g j
  hHZn_A4_5 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 5 q ≠ some (4, j)) → Hb ρ 4 6 g j = Hp ρ 4 5 g j
  hHE_A4 : ∀ (ρ : Fin 8) g t, GReach (Hb ρ 4 6 g (false, t)) ((G5 (r5 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem ρ 4 (false, t))) ∧
    GReach (Hb ρ 4 6 g (true, t)) ((G5 (r5 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem ρ 4 (true, t)))
  hH0_B4 : ∀ (ρ : Fin 8) g t, Hb ρ 5 6 g (false, t) = (G6 (r6 t g)).Xm.Φ (tt (I.tv t)) ∧ Hb ρ 5 6 g (true, t) = (G6 (r6 t g)).Xm.Φ 0
  hHP_B4_6 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 5 6 g j) (Hp ρ 5 6 g j) (fun _ => 0)
  hHZ_B4_6 : ∀ (ρ : Fin 8) g q j, hostAt ρ 6 q = some (5, j) → Wn6 ρ g q * Hp ρ 5 6 g j = Hb ρ 5 7 g j
  hHZn_B4_6 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 6 q ≠ some (5, j)) → Hb ρ 5 7 g j = Hp ρ 5 6 g j
  hHE_B4 : ∀ (ρ : Fin 8) g t, GReach (Hb ρ 5 7 g (false, t)) ((G6 (r6 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem ρ 5 (false, t))) ∧
    GReach (Hb ρ 5 7 g (true, t)) ((G6 (r6 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem ρ 5 (true, t)))
  hH0_F1 : ∀ (ρ : Fin 8) g t, Hb ρ 6 3 g (false, t) = (G2 (r2 t g)).Xm.Φ 1 ∧ Hb ρ 6 3 g (true, t) = (G2 (r2 t g)).Xm.Φ (1 + tt (I.tv t))
  hHP_F1_3 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 6 3 g j) (Hp ρ 6 3 g j) (fun _ => 0)
  hHZ_F1_3 : ∀ (ρ : Fin 8) g q j, hostAt ρ 3 q = some (6, j) → Wn3 ρ g q * Hp ρ 6 3 g j = Hb ρ 6 4 g j
  hHZn_F1_3 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 3 q ≠ some (6, j)) → Hb ρ 6 4 g j = Hp ρ 6 3 g j
  hHP_F1_4 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 6 4 g j) (Hp ρ 6 4 g j) (fun _ => 0)
  hHZ_F1_4 : ∀ (ρ : Fin 8) g q j, hostAt ρ 4 q = some (6, j) → Wn4 ρ g q * Hp ρ 6 4 g j = Hb ρ 6 5 g j
  hHZn_F1_4 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 4 q ≠ some (6, j)) → Hb ρ 6 5 g j = Hp ρ 6 4 g j
  hHP_F1_5 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 6 5 g j) (Hp ρ 6 5 g j) (fun _ => 0)
  hHZ_F1_5 : ∀ (ρ : Fin 8) g q j, hostAt ρ 5 q = some (6, j) → Wn5 ρ g q * Hp ρ 6 5 g j = Hb ρ 6 6 g j
  hHZn_F1_5 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 5 q ≠ some (6, j)) → Hb ρ 6 6 g j = Hp ρ 6 5 g j
  hHP_F1_6 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 6 6 g j) (Hp ρ 6 6 g j) (fun _ => 0)
  hHZ_F1_6 : ∀ (ρ : Fin 8) g q j, hostAt ρ 6 q = some (6, j) → Wn6 ρ g q * Hp ρ 6 6 g j = Hb ρ 6 7 g j
  hHZn_F1_6 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 6 q ≠ some (6, j)) → Hb ρ 6 7 g j = Hp ρ 6 6 g j
  hHP_F1_7 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 6 7 g j) (Hp ρ 6 7 g j) (fun _ => 0)
  hHZ_F1_7 : ∀ (ρ : Fin 8) g q j, hostAt ρ 7 q = some (6, j) → Wn7 ρ g q * Hp ρ 6 7 g j = Hb ρ 6 8 g j
  hHZn_F1_7 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 7 q ≠ some (6, j)) → Hb ρ 6 8 g j = Hp ρ 6 7 g j
  hHP_F1_8 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 6 8 g j) (Hp ρ 6 8 g j) (fun _ => 0)
  hHZ_F1_8 : ∀ (ρ : Fin 8) g q j, hostAt ρ 8 q = some (6, j) → Wn8 ρ g q * Hp ρ 6 8 g j = Hb ρ 6 9 g j
  hHZn_F1_8 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 8 q ≠ some (6, j)) → Hb ρ 6 9 g j = Hp ρ 6 8 g j
  hHE_F1 : ∀ (ρ : Fin 8) g t, GReach (Hb ρ 6 9 g (false, t)) (kernel α) (fun φ => bcost φ (hrem ρ 6 (false, t))) ∧
    GReach (Hb ρ 6 9 g (true, t)) (YF g t) (fun φ => bcost φ (hrem ρ 6 (true, t)))
  hH0_F2 : ∀ (ρ : Fin 8) g t, Hb ρ 7 5 g (false, t) = (G4 (r4 t g)).Xm.Φ 1 ∧ Hb ρ 7 5 g (true, t) = (G4 (r4 t g)).Xm.Φ (1 + tt (I.tv t))
  hHP_F2_5 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 7 5 g j) (Hp ρ 7 5 g j) (fun _ => 0)
  hHZ_F2_5 : ∀ (ρ : Fin 8) g q j, hostAt ρ 5 q = some (7, j) → Wn5 ρ g q * Hp ρ 7 5 g j = Hb ρ 7 6 g j
  hHZn_F2_5 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 5 q ≠ some (7, j)) → Hb ρ 7 6 g j = Hp ρ 7 5 g j
  hHP_F2_6 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 7 6 g j) (Hp ρ 7 6 g j) (fun _ => 0)
  hHZ_F2_6 : ∀ (ρ : Fin 8) g q j, hostAt ρ 6 q = some (7, j) → Wn6 ρ g q * Hp ρ 7 6 g j = Hb ρ 7 7 g j
  hHZn_F2_6 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 6 q ≠ some (7, j)) → Hb ρ 7 7 g j = Hp ρ 7 6 g j
  hHP_F2_7 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 7 7 g j) (Hp ρ 7 7 g j) (fun _ => 0)
  hHZ_F2_7 : ∀ (ρ : Fin 8) g q j, hostAt ρ 7 q = some (7, j) → Wn7 ρ g q * Hp ρ 7 7 g j = Hb ρ 7 8 g j
  hHZn_F2_7 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 7 q ≠ some (7, j)) → Hb ρ 7 8 g j = Hp ρ 7 7 g j
  hHP_F2_8 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 7 8 g j) (Hp ρ 7 8 g j) (fun _ => 0)
  hHZ_F2_8 : ∀ (ρ : Fin 8) g q j, hostAt ρ 8 q = some (7, j) → Wn8 ρ g q * Hp ρ 7 8 g j = Hb ρ 7 9 g j
  hHZn_F2_8 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 8 q ≠ some (7, j)) → Hb ρ 7 9 g j = Hp ρ 7 8 g j
  hHE_F2 : ∀ (ρ : Fin 8) g t, GReach (Hb ρ 7 9 g (false, t)) (kernel α) (fun φ => bcost φ (hrem ρ 7 (false, t))) ∧
    GReach (Hb ρ 7 9 g (true, t)) (YF g t) (fun φ => bcost φ (hrem ρ 7 (true, t)))
  hH0_F3 : ∀ (ρ : Fin 8) g t, Hb ρ 8 7 g (false, t) = (G6 (r6 t g)).Xm.Φ 1 ∧ Hb ρ 8 7 g (true, t) = (G6 (r6 t g)).Xm.Φ (1 + tt (I.tv t))
  hHP_F3_7 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 8 7 g j) (Hp ρ 8 7 g j) (fun _ => 0)
  hHZ_F3_7 : ∀ (ρ : Fin 8) g q j, hostAt ρ 7 q = some (8, j) → Wn7 ρ g q * Hp ρ 8 7 g j = Hb ρ 8 8 g j
  hHZn_F3_7 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 7 q ≠ some (8, j)) → Hb ρ 8 8 g j = Hp ρ 8 7 g j
  hHP_F3_8 : ∀ (ρ : Fin 8) g j, GReach (Hb ρ 8 8 g j) (Hp ρ 8 8 g j) (fun _ => 0)
  hHZ_F3_8 : ∀ (ρ : Fin 8) g q j, hostAt ρ 8 q = some (8, j) → Wn8 ρ g q * Hp ρ 8 8 g j = Hb ρ 8 9 g j
  hHZn_F3_8 : ∀ (ρ : Fin 8) g j, (∀ q, hostAt ρ 8 q ≠ some (8, j)) → Hb ρ 8 9 g j = Hp ρ 8 8 g j
  hHE_F3 : ∀ (ρ : Fin 8) g t, GReach (Hb ρ 8 9 g (false, t)) (kernel α) (fun φ => bcost φ (hrem ρ 8 (false, t))) ∧
    GReach (Hb ρ 8 9 g (true, t)) (YF g t) (fun φ => bcost φ (hrem ρ 8 (true, t)))
  /-- pair 4 after stage 8: ONE final block (not hosted). -/
  fX4 : ∀ g t, XReach ((G8 (r8 t g)).Xm.Φ 1) (kernel α) cX4
  fY4 : ∀ g t, XReach ((G8 (r8 t g)).Xm.Φ (1 + tt (I.tv t))) (YF g t) cY4
  term : ∀ g t, kernel α * (G0 g).Xm.Φ (tt (I.tv t)) = shift (u g t) * YF g t

section Net
variable {H T Sl C α Γ : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α]
  [Fintype Γ] [DecidableEq Γ]
variable {nb : ℕ} (N : GHost401 H T Sl C α Γ nb)

/-- price of the remaining climbs of the host chains of kind `c` (both sides, all classes), per element of `Γ`. -/
def GHost401.hcost (ρ : Fin 8) (c : Fin 9) (φ : ℕ → ℝ) : ℝ :=
  ∑ t : T, (bcost φ (N.hrem ρ c (false, t)) + bcost φ (N.hrem ρ c (true, t)))

/-- price of run `ρ` per element of `Γ`: its 9 calls (each its own provider's price), the remaining climbs of run `ρ`, pair 4's final. -/
def GHost401.runUnit (ρ : Fin 8) (φ : ℕ → ℝ) : ℝ :=
  (N.Ic ρ 0).cost φ + (N.Ic ρ 1).cost φ + (N.Ic ρ 2).cost φ + (N.Ic ρ 3).cost φ + (N.Ic ρ 4).cost φ + (N.Ic ρ 5).cost φ + (N.Ic ρ 6).cost φ + (N.Ic ρ 7).cost φ + (N.Ic ρ 8).cost φ
    + (N.hcost ρ 0 φ + N.hcost ρ 1 φ + N.hcost ρ 2 φ + N.hcost ρ 3 φ + N.hcost ρ 4 φ + N.hcost ρ 5 φ
    + N.hcost ρ 6 φ + N.hcost ρ 7 φ + N.hcost ρ 8 φ) + (Fintype.card T : ℝ) * (N.cX4 φ + N.cY4 φ)

/-- price of one unit (one element of `Γ`): the 8 runs (each with its own remaining climbs) and the padding blocks. -/
def GHost401.unit (φ : ℕ → ℝ) : ℝ :=
  (∑ ρ : Fin 8, N.runUnit ρ φ) + ∑ b, bcost φ (N.pad b)

/-- the invocation of call `(ρ, l)` with the target vectors of `I` DEFINITIONALLY (`retagTv`). -/
def GHost401.Jc (ρ : Fin 8) (l : Fin 9) : GInvE H T Sl C := XHostD.retagTv (N.Ic ρ l) N.I.tv (N.hIc ρ l)

/-- the conditions of the role map of stage `l`. -/
lemma GHost401.hok (ρ : Fin 8) (l : Fin 9) : HostOK (pairL l.val) (N.hostAt ρ l) (N.bk ρ l) :=
  ⟨N.hinj ρ l, N.hostInj ρ l, fun q c j h => win_pair c l (N.hostWin ρ l q c j h),
    fun q q' c c' j j' h h' e => win_uniq c c' l (N.hostWin ρ l q c j h) (N.hostWin ρ l q' c' j' h') e⟩

def GHost401.xi0 (g : Γ) (t : T) : CMat α := (N.G0 g).Xm.Φ (tt (N.I.tv t))
def GHost401.yi0 (g : Γ) (t : T) : CMat α := (N.G0 g).Xm.Φ 0
def GHost401.xi1 (g : Γ) (t : T) : CMat α := (N.G1 (N.r1 t g)).Xm.Φ (tt (N.I.tv t))
def GHost401.yi1 (g : Γ) (t : T) : CMat α := (N.G1 (N.r1 t g)).Xm.Φ 0
def GHost401.xi2 (g : Γ) (t : T) : CMat α := (N.G2 (N.r2 t g)).Xm.Φ (tt (N.I.tv t))
def GHost401.yi2 (g : Γ) (t : T) : CMat α := (N.G2 (N.r2 t g)).Xm.Φ 0
def GHost401.xi3 (g : Γ) (t : T) : CMat α := (N.G3 (N.r3 t g)).Xm.Φ (tt (N.I.tv t))
def GHost401.yi3 (g : Γ) (t : T) : CMat α := (N.G3 (N.r3 t g)).Xm.Φ 0
def GHost401.xi4 (g : Γ) (t : T) : CMat α := (N.G4 (N.r4 t g)).Xm.Φ (tt (N.I.tv t))
def GHost401.yi4 (g : Γ) (t : T) : CMat α := (N.G4 (N.r4 t g)).Xm.Φ 0
def GHost401.xi5 (g : Γ) (t : T) : CMat α := (N.G5 (N.r5 t g)).Xm.Φ (tt (N.I.tv t))
def GHost401.yi5 (g : Γ) (t : T) : CMat α := (N.G5 (N.r5 t g)).Xm.Φ 0
def GHost401.xi6 (g : Γ) (t : T) : CMat α := (N.G6 (N.r6 t g)).Xm.Φ (tt (N.I.tv t))
def GHost401.yi6 (g : Γ) (t : T) : CMat α := (N.G6 (N.r6 t g)).Xm.Φ 0
def GHost401.xi7 (g : Γ) (t : T) : CMat α := (N.G7 (N.r7 t g)).Xm.Φ (tt (N.I.tv t))
def GHost401.yi7 (g : Γ) (t : T) : CMat α := (N.G7 (N.r7 t g)).Xm.Φ 0
def GHost401.xi8 (g : Γ) (t : T) : CMat α := (N.G8 (N.r8 t g)).Xm.Φ (tt (N.I.tv t))
def GHost401.yi8 (g : Γ) (t : T) : CMat α := (N.G8 (N.r8 t g)).Xm.Φ 0
def GHost401.xo2 (g : Γ) (t : T) : CMat α := (N.G2 (N.r2 t g)).Xm.Φ 1
def GHost401.yo2 (g : Γ) (t : T) : CMat α := (N.G2 (N.r2 t g)).Xm.Φ (1 + tt (N.I.tv t))
def GHost401.xo4 (g : Γ) (t : T) : CMat α := (N.G4 (N.r4 t g)).Xm.Φ 1
def GHost401.yo4 (g : Γ) (t : T) : CMat α := (N.G4 (N.r4 t g)).Xm.Φ (1 + tt (N.I.tv t))
def GHost401.xo6 (g : Γ) (t : T) : CMat α := (N.G6 (N.r6 t g)).Xm.Φ 1
def GHost401.yo6 (g : Γ) (t : T) : CMat α := (N.G6 (N.r6 t g)).Xm.Φ (1 + tt (N.I.tv t))
def GHost401.xo8 (g : Γ) (t : T) : CMat α := (N.G8 (N.r8 t g)).Xm.Φ 1
def GHost401.yo8 (g : Γ) (t : T) : CMat α := (N.G8 (N.r8 t g)).Xm.Φ (1 + tt (N.I.tv t))

/-- the data of a replica before its run and after it. -/
def GHost401.InitD : BDat4 Γ T → CMat α := mkD4 N.xi0 N.yi0 N.xi0 N.yi0 N.xi0 N.yi0 N.xi0 N.yi0
def GHost401.FinD : BDat4 Γ T → CMat α := mkD4 (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF

def GHost401.S0 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi0 N.yi0 N.xi0 N.yi0 N.xi0 N.yi0 N.xi0 N.yi0 (fun g b => N.Zb (9 * ρ.val + 0) g b) (cOne Γ C α)
def GHost401.S1 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi0 N.yi0 (fun g t => N.Hp ρ 0 0 g (false, t)) (fun g t => N.Hp ρ 0 0 g (true, t)) (fun g t => N.Hp ρ 2 0 g (false, t)) (fun g t => N.Hp ρ 2 0 g (true, t)) (fun g t => N.Hp ρ 4 0 g (false, t)) (fun g t => N.Hp ρ 4 0 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 0) g b) (cOne Γ C α)
def GHost401.S2 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi1 N.yi1 (fun g t => N.Hb ρ 0 1 g (false, t)) (fun g t => N.Hb ρ 0 1 g (true, t)) (fun g t => N.Hb ρ 2 1 g (false, t)) (fun g t => N.Hb ρ 2 1 g (true, t)) (fun g t => N.Hb ρ 4 1 g (false, t)) (fun g t => N.Hb ρ 4 1 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 1) g b) (cOne Γ C α)
def GHost401.S3 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi1 N.yi1 (fun g t => N.Hp ρ 0 1 g (false, t)) (fun g t => N.Hp ρ 0 1 g (true, t)) (fun g t => N.Hp ρ 2 1 g (false, t)) (fun g t => N.Hp ρ 2 1 g (true, t)) (fun g t => N.Hp ρ 4 1 g (false, t)) (fun g t => N.Hp ρ 4 1 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 1) g b) (cOne Γ C α)
def GHost401.S4 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi2 N.yi2 (fun g t => N.Hb ρ 0 2 g (false, t)) (fun g t => N.Hb ρ 0 2 g (true, t)) (fun g t => N.Hb ρ 2 2 g (false, t)) (fun g t => N.Hb ρ 2 2 g (true, t)) (fun g t => N.Hb ρ 4 2 g (false, t)) (fun g t => N.Hb ρ 4 2 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 2) g b) (cOne Γ C α)
def GHost401.S5 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi2 N.yi2 N.xi2 N.yi2 (fun g t => N.Hb ρ 2 2 g (false, t)) (fun g t => N.Hb ρ 2 2 g (true, t)) (fun g t => N.Hb ρ 4 2 g (false, t)) (fun g t => N.Hb ρ 4 2 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 2) g b) (cOne Γ C α)
def GHost401.S6 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi2 N.yi2 (fun g t => N.Hp ρ 1 2 g (false, t)) (fun g t => N.Hp ρ 1 2 g (true, t)) (fun g t => N.Hp ρ 2 2 g (false, t)) (fun g t => N.Hp ρ 2 2 g (true, t)) (fun g t => N.Hp ρ 4 2 g (false, t)) (fun g t => N.Hp ρ 4 2 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 2) g b) (cOne Γ C α)
def GHost401.S7 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xo2 N.yo2 (fun g t => N.Hb ρ 1 3 g (false, t)) (fun g t => N.Hb ρ 1 3 g (true, t)) (fun g t => N.Hb ρ 2 3 g (false, t)) (fun g t => N.Hb ρ 2 3 g (true, t)) (fun g t => N.Hb ρ 4 3 g (false, t)) (fun g t => N.Hb ρ 4 3 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 3) g b) (cOne Γ C α)
def GHost401.S8 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xo2 N.yo2 N.xi3 N.yi3 (fun g t => N.Hb ρ 2 3 g (false, t)) (fun g t => N.Hb ρ 2 3 g (true, t)) (fun g t => N.Hb ρ 4 3 g (false, t)) (fun g t => N.Hb ρ 4 3 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 3) g b) (cOne Γ C α)
def GHost401.S9 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp ρ 6 3 g (false, t)) (fun g t => N.Hp ρ 6 3 g (true, t)) N.xi3 N.yi3 (fun g t => N.Hp ρ 2 3 g (false, t)) (fun g t => N.Hp ρ 2 3 g (true, t)) (fun g t => N.Hp ρ 4 3 g (false, t)) (fun g t => N.Hp ρ 4 3 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 3) g b) (cOne Γ C α)
def GHost401.S10 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb ρ 6 4 g (false, t)) (fun g t => N.Hb ρ 6 4 g (true, t)) N.xi4 N.yi4 (fun g t => N.Hb ρ 2 4 g (false, t)) (fun g t => N.Hb ρ 2 4 g (true, t)) (fun g t => N.Hb ρ 4 4 g (false, t)) (fun g t => N.Hb ρ 4 4 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 4) g b) (cOne Γ C α)
def GHost401.S11 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb ρ 6 4 g (false, t)) (fun g t => N.Hb ρ 6 4 g (true, t)) N.xi4 N.yi4 N.xi4 N.yi4 (fun g t => N.Hb ρ 4 4 g (false, t)) (fun g t => N.Hb ρ 4 4 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 4) g b) (cOne Γ C α)
def GHost401.S12 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp ρ 6 4 g (false, t)) (fun g t => N.Hp ρ 6 4 g (true, t)) N.xi4 N.yi4 (fun g t => N.Hp ρ 3 4 g (false, t)) (fun g t => N.Hp ρ 3 4 g (true, t)) (fun g t => N.Hp ρ 4 4 g (false, t)) (fun g t => N.Hp ρ 4 4 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 4) g b) (cOne Γ C α)
def GHost401.S13 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb ρ 6 5 g (false, t)) (fun g t => N.Hb ρ 6 5 g (true, t)) N.xo4 N.yo4 (fun g t => N.Hb ρ 3 5 g (false, t)) (fun g t => N.Hb ρ 3 5 g (true, t)) (fun g t => N.Hb ρ 4 5 g (false, t)) (fun g t => N.Hb ρ 4 5 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 5) g b) (cOne Γ C α)
def GHost401.S14 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb ρ 6 5 g (false, t)) (fun g t => N.Hb ρ 6 5 g (true, t)) N.xo4 N.yo4 N.xi5 N.yi5 (fun g t => N.Hb ρ 4 5 g (false, t)) (fun g t => N.Hb ρ 4 5 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 5) g b) (cOne Γ C α)
def GHost401.S15 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp ρ 6 5 g (false, t)) (fun g t => N.Hp ρ 6 5 g (true, t)) (fun g t => N.Hp ρ 7 5 g (false, t)) (fun g t => N.Hp ρ 7 5 g (true, t)) N.xi5 N.yi5 (fun g t => N.Hp ρ 4 5 g (false, t)) (fun g t => N.Hp ρ 4 5 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 5) g b) (cOne Γ C α)
def GHost401.S16 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb ρ 6 6 g (false, t)) (fun g t => N.Hb ρ 6 6 g (true, t)) (fun g t => N.Hb ρ 7 6 g (false, t)) (fun g t => N.Hb ρ 7 6 g (true, t)) N.xi6 N.yi6 (fun g t => N.Hb ρ 4 6 g (false, t)) (fun g t => N.Hb ρ 4 6 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 6) g b) (cOne Γ C α)
def GHost401.S17 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb ρ 6 6 g (false, t)) (fun g t => N.Hb ρ 6 6 g (true, t)) (fun g t => N.Hb ρ 7 6 g (false, t)) (fun g t => N.Hb ρ 7 6 g (true, t)) N.xi6 N.yi6 N.xi6 N.yi6 (fun g b => N.Zb (9 * ρ.val + 6) g b) (cOne Γ C α)
def GHost401.S18 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp ρ 6 6 g (false, t)) (fun g t => N.Hp ρ 6 6 g (true, t)) (fun g t => N.Hp ρ 7 6 g (false, t)) (fun g t => N.Hp ρ 7 6 g (true, t)) N.xi6 N.yi6 (fun g t => N.Hp ρ 5 6 g (false, t)) (fun g t => N.Hp ρ 5 6 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 6) g b) (cOne Γ C α)
def GHost401.S19 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb ρ 6 7 g (false, t)) (fun g t => N.Hb ρ 6 7 g (true, t)) (fun g t => N.Hb ρ 7 7 g (false, t)) (fun g t => N.Hb ρ 7 7 g (true, t)) N.xo6 N.yo6 (fun g t => N.Hb ρ 5 7 g (false, t)) (fun g t => N.Hb ρ 5 7 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 7) g b) (cOne Γ C α)
def GHost401.S20 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb ρ 6 7 g (false, t)) (fun g t => N.Hb ρ 6 7 g (true, t)) (fun g t => N.Hb ρ 7 7 g (false, t)) (fun g t => N.Hb ρ 7 7 g (true, t)) N.xo6 N.yo6 N.xi7 N.yi7 (fun g b => N.Zb (9 * ρ.val + 7) g b) (cOne Γ C α)
def GHost401.S21 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp ρ 6 7 g (false, t)) (fun g t => N.Hp ρ 6 7 g (true, t)) (fun g t => N.Hp ρ 7 7 g (false, t)) (fun g t => N.Hp ρ 7 7 g (true, t)) (fun g t => N.Hp ρ 8 7 g (false, t)) (fun g t => N.Hp ρ 8 7 g (true, t)) N.xi7 N.yi7 (fun g b => N.Zp (9 * ρ.val + 7) g b) (cOne Γ C α)
def GHost401.S22 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb ρ 6 8 g (false, t)) (fun g t => N.Hb ρ 6 8 g (true, t)) (fun g t => N.Hb ρ 7 8 g (false, t)) (fun g t => N.Hb ρ 7 8 g (true, t)) (fun g t => N.Hb ρ 8 8 g (false, t)) (fun g t => N.Hb ρ 8 8 g (true, t)) N.xi8 N.yi8 (fun g b => N.Zb (9 * ρ.val + 8) g b) (cOne Γ C α)
def GHost401.S23 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp ρ 6 8 g (false, t)) (fun g t => N.Hp ρ 6 8 g (true, t)) (fun g t => N.Hp ρ 7 8 g (false, t)) (fun g t => N.Hp ρ 7 8 g (true, t)) (fun g t => N.Hp ρ 8 8 g (false, t)) (fun g t => N.Hp ρ 8 8 g (true, t)) N.xi8 N.yi8 (fun g b => N.Zp (9 * ρ.val + 8) g b) (cOne Γ C α)
def GHost401.S24 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb ρ 6 9 g (false, t)) (fun g t => N.Hb ρ 6 9 g (true, t)) (fun g t => N.Hb ρ 7 9 g (false, t)) (fun g t => N.Hb ρ 7 9 g (true, t)) (fun g t => N.Hb ρ 8 9 g (false, t)) (fun g t => N.Hb ρ 8 9 g (true, t)) N.xo8 N.yo8 (fun g b => N.Zb (9 * ρ.val + 9) g b) (cOne Γ C α)
def GHost401.S25 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF (fun g b => N.Zb (9 * ρ.val + 9) g b) (cOne Γ C α)

/-- frames of all roles at the start of run `k` (`k = 8`: after the last run). -/
def GHost401.Sfull (k : Fin 9) : BRole401 Γ T nb C → CMat α :=
  mk401 (fun ρ' d => if ρ'.val < k.val then N.FinD d else N.InitD d) (fun g b => N.Zb (9 * k.val) g b) (cOne Γ C α)

/-- frames of all roles after the padding climbs. -/
def GHost401.Sfin : BRole401 Γ T nb C → CMat α :=
  mk401 (fun _ d => N.FinD d) (fun _ _ => kernel α) (cOne Γ C α)

end Net
end
end XHostD
end PowerSaving
end OAI

