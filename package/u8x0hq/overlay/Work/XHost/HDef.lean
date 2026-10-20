import Work.XHost.HRoles

/-!
# (key: bridge-host) The nine-stage word of PR #401 with SAME-RUN HOSTING: data and frames

Copy of Work/Bridge401/GNetDef.lean (Ryan Shea's complex588 tree; `G401` by Nil Mamano, written by his gen401.py, on
Dugongue #401 / jacklightChen #377 / Lydakis #275 / DaysSky #363; engine and E8 unit by Jacob Sussman) with idle-axis
hosting (DaysSky, CrocSwap/integer-mult-bounds PR #442): in stage `l` the slot `q` goes to the idle data role
`hostAt l q = some (c, s, t)` (window kind `c`, side `s`, class `(g, t)` of the run's replica; the same table in every
run) or, when `hostAt l q = none`, to the bank `bk ρ l q`.  The idle climbs A2..A4, B2..B4 and the finals F1..F3 become
host chains `Hb`/`Hp`: free adapters, the hosted residuals, ONE remaining block of rank `hrem c j` at the old place.
Written by h-L3-build/gen/gen_hnet.py (wht43, 13 Oct 2026), prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace XHost
open Binary Matrix Finset RAM SS CB RF BR GF GX
noncomputable section

structure GHost401 (H T Sl C α Γ : Type) (nb : ℕ) [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
    [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α] where
  /-- the gauged invocation (`GInvE`). -/
  I : GInvE H T Sl C
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
  N0 : Γ → Sl → CMat α
  Wn0 : Γ → Sl → CMat α
  N1 : Γ → Sl → CMat α
  Wn1 : Γ → Sl → CMat α
  N2 : Γ → Sl → CMat α
  Wn2 : Γ → Sl → CMat α
  N3 : Γ → Sl → CMat α
  Wn3 : Γ → Sl → CMat α
  N4 : Γ → Sl → CMat α
  Wn4 : Γ → Sl → CMat α
  N5 : Γ → Sl → CMat α
  Wn5 : Γ → Sl → CMat α
  N6 : Γ → Sl → CMat α
  Wn6 : Γ → Sl → CMat α
  N7 : Γ → Sl → CMat α
  Wn7 : Γ → Sl → CMat α
  N8 : Γ → Sl → CMat α
  Wn8 : Γ → Sl → CMat α
  /-- the host of slot `q` in stage `l`, the SAME in every run: (window kind, side, class); `none`: a bank. -/
  hostAt : Fin 9 → Sl → Option (Fin 9 × (Bool × T))
  hostInj : ∀ l q q' j, hostAt l q = some j → hostAt l q' = some j → q = q'
  hostWin : ∀ l q c j, hostAt l q = some (c, j) → inWin c.val l.val = true
  /-- the bank of a NON-hosted slot `q` in stage `l` of run `ρ` (junk on hosted slots). -/
  bk : Fin 8 → Fin 9 → Sl → Fin nb
  hinj : ∀ ρ l q q', hostAt l q = none → hostAt l q' = none → bk ρ l q = bk ρ l q' → q = q'
  Zb : ℕ → Γ → Fin nb → CMat α
  Zp : ℕ → Γ → Fin nb → CMat α
  pad : Fin nb → ℕ
  /-- the host chains: frame of host `j` of kind `c` before (`Hb`) and after (`Hp`) the adapter of stage `k`. -/
  Hb : Fin 9 → ℕ → Γ → Bool × T → CMat α
  Hp : Fin 9 → ℕ → Γ → Bool × T → CMat α
  hrem : Fin 9 → Bool × T → ℕ
  u : Γ → T → Space α
  YF : Γ → T → CMat α
  hm : 2 ≤ Fintype.card α
  cX4 : (ℕ → ℝ) → ℝ
  cY4 : (ℕ → ℝ) → ℝ
  inv0 : ∀ g t, ∃ N : CMat α, (G0 g).Xm.Φ (tt (I.tv t)) * N = 1
  z0 : ∀ g, (G0 g).Xm.Φ 0 = 1
  hN0 : ∀ d q, (G0 d).X.Φ (I.ent q) * N0 d q = 1
  hW0 : ∀ d q, (G0 d).X.Φ (I.ent q).top = Wn0 d q * (G0 d).X.Φ (I.ent q)
  hN1 : ∀ d q, (G1 d).X.Φ (I.ent q) * N1 d q = 1
  hW1 : ∀ d q, (G1 d).X.Φ (I.ent q).top = Wn1 d q * (G1 d).X.Φ (I.ent q)
  hN2 : ∀ d q, (G2 d).X.Φ (I.ent q) * N2 d q = 1
  hW2 : ∀ d q, (G2 d).X.Φ (I.ent q).top = Wn2 d q * (G2 d).X.Φ (I.ent q)
  hN3 : ∀ d q, (G3 d).X.Φ (I.ent q) * N3 d q = 1
  hW3 : ∀ d q, (G3 d).X.Φ (I.ent q).top = Wn3 d q * (G3 d).X.Φ (I.ent q)
  hN4 : ∀ d q, (G4 d).X.Φ (I.ent q) * N4 d q = 1
  hW4 : ∀ d q, (G4 d).X.Φ (I.ent q).top = Wn4 d q * (G4 d).X.Φ (I.ent q)
  hN5 : ∀ d q, (G5 d).X.Φ (I.ent q) * N5 d q = 1
  hW5 : ∀ d q, (G5 d).X.Φ (I.ent q).top = Wn5 d q * (G5 d).X.Φ (I.ent q)
  hN6 : ∀ d q, (G6 d).X.Φ (I.ent q) * N6 d q = 1
  hW6 : ∀ d q, (G6 d).X.Φ (I.ent q).top = Wn6 d q * (G6 d).X.Φ (I.ent q)
  hN7 : ∀ d q, (G7 d).X.Φ (I.ent q) * N7 d q = 1
  hW7 : ∀ d q, (G7 d).X.Φ (I.ent q).top = Wn7 d q * (G7 d).X.Φ (I.ent q)
  hN8 : ∀ d q, (G8 d).X.Φ (I.ent q) * N8 d q = 1
  hW8 : ∀ d q, (G8 d).X.Φ (I.ent q).top = Wn8 d q * (G8 d).X.Φ (I.ent q)
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
  hZ0 : ∀ (ρ : Fin 8) g q, hostAt 0 q = none → Wn0 g q * Zp (9 * ρ.val + 0) g (bk ρ 0 q) = Zb (9 * ρ.val + 1) g (bk ρ 0 q)
  hZn0 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt 0 q = none → bk ρ 0 q ≠ b) → Zb (9 * ρ.val + 1) g b = Zp (9 * ρ.val + 0) g b
  hP1 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 1) g b) (Zp (9 * ρ.val + 1) g b) (fun _ => 0)
  hZ1 : ∀ (ρ : Fin 8) g q, hostAt 1 q = none → Wn1 g q * Zp (9 * ρ.val + 1) g (bk ρ 1 q) = Zb (9 * ρ.val + 2) g (bk ρ 1 q)
  hZn1 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt 1 q = none → bk ρ 1 q ≠ b) → Zb (9 * ρ.val + 2) g b = Zp (9 * ρ.val + 1) g b
  hP2 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 2) g b) (Zp (9 * ρ.val + 2) g b) (fun _ => 0)
  hZ2 : ∀ (ρ : Fin 8) g q, hostAt 2 q = none → Wn2 g q * Zp (9 * ρ.val + 2) g (bk ρ 2 q) = Zb (9 * ρ.val + 3) g (bk ρ 2 q)
  hZn2 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt 2 q = none → bk ρ 2 q ≠ b) → Zb (9 * ρ.val + 3) g b = Zp (9 * ρ.val + 2) g b
  hP3 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 3) g b) (Zp (9 * ρ.val + 3) g b) (fun _ => 0)
  hZ3 : ∀ (ρ : Fin 8) g q, hostAt 3 q = none → Wn3 g q * Zp (9 * ρ.val + 3) g (bk ρ 3 q) = Zb (9 * ρ.val + 4) g (bk ρ 3 q)
  hZn3 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt 3 q = none → bk ρ 3 q ≠ b) → Zb (9 * ρ.val + 4) g b = Zp (9 * ρ.val + 3) g b
  hP4 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 4) g b) (Zp (9 * ρ.val + 4) g b) (fun _ => 0)
  hZ4 : ∀ (ρ : Fin 8) g q, hostAt 4 q = none → Wn4 g q * Zp (9 * ρ.val + 4) g (bk ρ 4 q) = Zb (9 * ρ.val + 5) g (bk ρ 4 q)
  hZn4 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt 4 q = none → bk ρ 4 q ≠ b) → Zb (9 * ρ.val + 5) g b = Zp (9 * ρ.val + 4) g b
  hP5 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 5) g b) (Zp (9 * ρ.val + 5) g b) (fun _ => 0)
  hZ5 : ∀ (ρ : Fin 8) g q, hostAt 5 q = none → Wn5 g q * Zp (9 * ρ.val + 5) g (bk ρ 5 q) = Zb (9 * ρ.val + 6) g (bk ρ 5 q)
  hZn5 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt 5 q = none → bk ρ 5 q ≠ b) → Zb (9 * ρ.val + 6) g b = Zp (9 * ρ.val + 5) g b
  hP6 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 6) g b) (Zp (9 * ρ.val + 6) g b) (fun _ => 0)
  hZ6 : ∀ (ρ : Fin 8) g q, hostAt 6 q = none → Wn6 g q * Zp (9 * ρ.val + 6) g (bk ρ 6 q) = Zb (9 * ρ.val + 7) g (bk ρ 6 q)
  hZn6 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt 6 q = none → bk ρ 6 q ≠ b) → Zb (9 * ρ.val + 7) g b = Zp (9 * ρ.val + 6) g b
  hP7 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 7) g b) (Zp (9 * ρ.val + 7) g b) (fun _ => 0)
  hZ7 : ∀ (ρ : Fin 8) g q, hostAt 7 q = none → Wn7 g q * Zp (9 * ρ.val + 7) g (bk ρ 7 q) = Zb (9 * ρ.val + 8) g (bk ρ 7 q)
  hZn7 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt 7 q = none → bk ρ 7 q ≠ b) → Zb (9 * ρ.val + 8) g b = Zp (9 * ρ.val + 7) g b
  hP8 : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + 8) g b) (Zp (9 * ρ.val + 8) g b) (fun _ => 0)
  hZ8 : ∀ (ρ : Fin 8) g q, hostAt 8 q = none → Wn8 g q * Zp (9 * ρ.val + 8) g (bk ρ 8 q) = Zb (9 * ρ.val + 9) g (bk ρ 8 q)
  hZn8 : ∀ (ρ : Fin 8) g b, (∀ q, hostAt 8 q = none → bk ρ 8 q ≠ b) → Zb (9 * ρ.val + 9) g b = Zp (9 * ρ.val + 8) g b
  bpad : ∀ g b, GReach (Zb 72 g b) (kernel α) (fun φ => bcost φ (pad b))
  /-- the host chains (XHost.host_frames): start, free adapters, hosted residuals, ONE remaining block. -/
  hH0_A2 : ∀ g t, Hb 0 0 g (false, t) = (G0 g).Xm.Φ (tt (I.tv t)) ∧ Hb 0 0 g (true, t) = (G0 g).Xm.Φ 0
  hHP_A2_0 : ∀ g j, GReach (Hb 0 0 g j) (Hp 0 0 g j) (fun _ => 0)
  hHZ_A2_0 : ∀ g q j, hostAt 0 q = some (0, j) → Wn0 g q * Hp 0 0 g j = Hb 0 1 g j
  hHZn_A2_0 : ∀ g j, (∀ q, hostAt 0 q ≠ some (0, j)) → Hb 0 1 g j = Hp 0 0 g j
  hHP_A2_1 : ∀ g j, GReach (Hb 0 1 g j) (Hp 0 1 g j) (fun _ => 0)
  hHZ_A2_1 : ∀ g q j, hostAt 1 q = some (0, j) → Wn1 g q * Hp 0 1 g j = Hb 0 2 g j
  hHZn_A2_1 : ∀ g j, (∀ q, hostAt 1 q ≠ some (0, j)) → Hb 0 2 g j = Hp 0 1 g j
  hHE_A2 : ∀ g t, GReach (Hb 0 2 g (false, t)) ((G1 (r1 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem 0 (false, t))) ∧
    GReach (Hb 0 2 g (true, t)) ((G1 (r1 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem 0 (true, t)))
  hH0_B2 : ∀ g t, Hb 1 2 g (false, t) = (G2 (r2 t g)).Xm.Φ (tt (I.tv t)) ∧ Hb 1 2 g (true, t) = (G2 (r2 t g)).Xm.Φ 0
  hHP_B2_2 : ∀ g j, GReach (Hb 1 2 g j) (Hp 1 2 g j) (fun _ => 0)
  hHZ_B2_2 : ∀ g q j, hostAt 2 q = some (1, j) → Wn2 g q * Hp 1 2 g j = Hb 1 3 g j
  hHZn_B2_2 : ∀ g j, (∀ q, hostAt 2 q ≠ some (1, j)) → Hb 1 3 g j = Hp 1 2 g j
  hHE_B2 : ∀ g t, GReach (Hb 1 3 g (false, t)) ((G2 (r2 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem 1 (false, t))) ∧
    GReach (Hb 1 3 g (true, t)) ((G2 (r2 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem 1 (true, t)))
  hH0_A3 : ∀ g t, Hb 2 0 g (false, t) = (G0 g).Xm.Φ (tt (I.tv t)) ∧ Hb 2 0 g (true, t) = (G0 g).Xm.Φ 0
  hHP_A3_0 : ∀ g j, GReach (Hb 2 0 g j) (Hp 2 0 g j) (fun _ => 0)
  hHZ_A3_0 : ∀ g q j, hostAt 0 q = some (2, j) → Wn0 g q * Hp 2 0 g j = Hb 2 1 g j
  hHZn_A3_0 : ∀ g j, (∀ q, hostAt 0 q ≠ some (2, j)) → Hb 2 1 g j = Hp 2 0 g j
  hHP_A3_1 : ∀ g j, GReach (Hb 2 1 g j) (Hp 2 1 g j) (fun _ => 0)
  hHZ_A3_1 : ∀ g q j, hostAt 1 q = some (2, j) → Wn1 g q * Hp 2 1 g j = Hb 2 2 g j
  hHZn_A3_1 : ∀ g j, (∀ q, hostAt 1 q ≠ some (2, j)) → Hb 2 2 g j = Hp 2 1 g j
  hHP_A3_2 : ∀ g j, GReach (Hb 2 2 g j) (Hp 2 2 g j) (fun _ => 0)
  hHZ_A3_2 : ∀ g q j, hostAt 2 q = some (2, j) → Wn2 g q * Hp 2 2 g j = Hb 2 3 g j
  hHZn_A3_2 : ∀ g j, (∀ q, hostAt 2 q ≠ some (2, j)) → Hb 2 3 g j = Hp 2 2 g j
  hHP_A3_3 : ∀ g j, GReach (Hb 2 3 g j) (Hp 2 3 g j) (fun _ => 0)
  hHZ_A3_3 : ∀ g q j, hostAt 3 q = some (2, j) → Wn3 g q * Hp 2 3 g j = Hb 2 4 g j
  hHZn_A3_3 : ∀ g j, (∀ q, hostAt 3 q ≠ some (2, j)) → Hb 2 4 g j = Hp 2 3 g j
  hHE_A3 : ∀ g t, GReach (Hb 2 4 g (false, t)) ((G3 (r3 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem 2 (false, t))) ∧
    GReach (Hb 2 4 g (true, t)) ((G3 (r3 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem 2 (true, t)))
  hH0_B3 : ∀ g t, Hb 3 4 g (false, t) = (G4 (r4 t g)).Xm.Φ (tt (I.tv t)) ∧ Hb 3 4 g (true, t) = (G4 (r4 t g)).Xm.Φ 0
  hHP_B3_4 : ∀ g j, GReach (Hb 3 4 g j) (Hp 3 4 g j) (fun _ => 0)
  hHZ_B3_4 : ∀ g q j, hostAt 4 q = some (3, j) → Wn4 g q * Hp 3 4 g j = Hb 3 5 g j
  hHZn_B3_4 : ∀ g j, (∀ q, hostAt 4 q ≠ some (3, j)) → Hb 3 5 g j = Hp 3 4 g j
  hHE_B3 : ∀ g t, GReach (Hb 3 5 g (false, t)) ((G4 (r4 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem 3 (false, t))) ∧
    GReach (Hb 3 5 g (true, t)) ((G4 (r4 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem 3 (true, t)))
  hH0_A4 : ∀ g t, Hb 4 0 g (false, t) = (G0 g).Xm.Φ (tt (I.tv t)) ∧ Hb 4 0 g (true, t) = (G0 g).Xm.Φ 0
  hHP_A4_0 : ∀ g j, GReach (Hb 4 0 g j) (Hp 4 0 g j) (fun _ => 0)
  hHZ_A4_0 : ∀ g q j, hostAt 0 q = some (4, j) → Wn0 g q * Hp 4 0 g j = Hb 4 1 g j
  hHZn_A4_0 : ∀ g j, (∀ q, hostAt 0 q ≠ some (4, j)) → Hb 4 1 g j = Hp 4 0 g j
  hHP_A4_1 : ∀ g j, GReach (Hb 4 1 g j) (Hp 4 1 g j) (fun _ => 0)
  hHZ_A4_1 : ∀ g q j, hostAt 1 q = some (4, j) → Wn1 g q * Hp 4 1 g j = Hb 4 2 g j
  hHZn_A4_1 : ∀ g j, (∀ q, hostAt 1 q ≠ some (4, j)) → Hb 4 2 g j = Hp 4 1 g j
  hHP_A4_2 : ∀ g j, GReach (Hb 4 2 g j) (Hp 4 2 g j) (fun _ => 0)
  hHZ_A4_2 : ∀ g q j, hostAt 2 q = some (4, j) → Wn2 g q * Hp 4 2 g j = Hb 4 3 g j
  hHZn_A4_2 : ∀ g j, (∀ q, hostAt 2 q ≠ some (4, j)) → Hb 4 3 g j = Hp 4 2 g j
  hHP_A4_3 : ∀ g j, GReach (Hb 4 3 g j) (Hp 4 3 g j) (fun _ => 0)
  hHZ_A4_3 : ∀ g q j, hostAt 3 q = some (4, j) → Wn3 g q * Hp 4 3 g j = Hb 4 4 g j
  hHZn_A4_3 : ∀ g j, (∀ q, hostAt 3 q ≠ some (4, j)) → Hb 4 4 g j = Hp 4 3 g j
  hHP_A4_4 : ∀ g j, GReach (Hb 4 4 g j) (Hp 4 4 g j) (fun _ => 0)
  hHZ_A4_4 : ∀ g q j, hostAt 4 q = some (4, j) → Wn4 g q * Hp 4 4 g j = Hb 4 5 g j
  hHZn_A4_4 : ∀ g j, (∀ q, hostAt 4 q ≠ some (4, j)) → Hb 4 5 g j = Hp 4 4 g j
  hHP_A4_5 : ∀ g j, GReach (Hb 4 5 g j) (Hp 4 5 g j) (fun _ => 0)
  hHZ_A4_5 : ∀ g q j, hostAt 5 q = some (4, j) → Wn5 g q * Hp 4 5 g j = Hb 4 6 g j
  hHZn_A4_5 : ∀ g j, (∀ q, hostAt 5 q ≠ some (4, j)) → Hb 4 6 g j = Hp 4 5 g j
  hHE_A4 : ∀ g t, GReach (Hb 4 6 g (false, t)) ((G5 (r5 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem 4 (false, t))) ∧
    GReach (Hb 4 6 g (true, t)) ((G5 (r5 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem 4 (true, t)))
  hH0_B4 : ∀ g t, Hb 5 6 g (false, t) = (G6 (r6 t g)).Xm.Φ (tt (I.tv t)) ∧ Hb 5 6 g (true, t) = (G6 (r6 t g)).Xm.Φ 0
  hHP_B4_6 : ∀ g j, GReach (Hb 5 6 g j) (Hp 5 6 g j) (fun _ => 0)
  hHZ_B4_6 : ∀ g q j, hostAt 6 q = some (5, j) → Wn6 g q * Hp 5 6 g j = Hb 5 7 g j
  hHZn_B4_6 : ∀ g j, (∀ q, hostAt 6 q ≠ some (5, j)) → Hb 5 7 g j = Hp 5 6 g j
  hHE_B4 : ∀ g t, GReach (Hb 5 7 g (false, t)) ((G6 (r6 t g)).Xm.Φ 1) (fun φ => bcost φ (hrem 5 (false, t))) ∧
    GReach (Hb 5 7 g (true, t)) ((G6 (r6 t g)).Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (hrem 5 (true, t)))
  hH0_F1 : ∀ g t, Hb 6 3 g (false, t) = (G2 (r2 t g)).Xm.Φ 1 ∧ Hb 6 3 g (true, t) = (G2 (r2 t g)).Xm.Φ (1 + tt (I.tv t))
  hHP_F1_3 : ∀ g j, GReach (Hb 6 3 g j) (Hp 6 3 g j) (fun _ => 0)
  hHZ_F1_3 : ∀ g q j, hostAt 3 q = some (6, j) → Wn3 g q * Hp 6 3 g j = Hb 6 4 g j
  hHZn_F1_3 : ∀ g j, (∀ q, hostAt 3 q ≠ some (6, j)) → Hb 6 4 g j = Hp 6 3 g j
  hHP_F1_4 : ∀ g j, GReach (Hb 6 4 g j) (Hp 6 4 g j) (fun _ => 0)
  hHZ_F1_4 : ∀ g q j, hostAt 4 q = some (6, j) → Wn4 g q * Hp 6 4 g j = Hb 6 5 g j
  hHZn_F1_4 : ∀ g j, (∀ q, hostAt 4 q ≠ some (6, j)) → Hb 6 5 g j = Hp 6 4 g j
  hHP_F1_5 : ∀ g j, GReach (Hb 6 5 g j) (Hp 6 5 g j) (fun _ => 0)
  hHZ_F1_5 : ∀ g q j, hostAt 5 q = some (6, j) → Wn5 g q * Hp 6 5 g j = Hb 6 6 g j
  hHZn_F1_5 : ∀ g j, (∀ q, hostAt 5 q ≠ some (6, j)) → Hb 6 6 g j = Hp 6 5 g j
  hHP_F1_6 : ∀ g j, GReach (Hb 6 6 g j) (Hp 6 6 g j) (fun _ => 0)
  hHZ_F1_6 : ∀ g q j, hostAt 6 q = some (6, j) → Wn6 g q * Hp 6 6 g j = Hb 6 7 g j
  hHZn_F1_6 : ∀ g j, (∀ q, hostAt 6 q ≠ some (6, j)) → Hb 6 7 g j = Hp 6 6 g j
  hHP_F1_7 : ∀ g j, GReach (Hb 6 7 g j) (Hp 6 7 g j) (fun _ => 0)
  hHZ_F1_7 : ∀ g q j, hostAt 7 q = some (6, j) → Wn7 g q * Hp 6 7 g j = Hb 6 8 g j
  hHZn_F1_7 : ∀ g j, (∀ q, hostAt 7 q ≠ some (6, j)) → Hb 6 8 g j = Hp 6 7 g j
  hHP_F1_8 : ∀ g j, GReach (Hb 6 8 g j) (Hp 6 8 g j) (fun _ => 0)
  hHZ_F1_8 : ∀ g q j, hostAt 8 q = some (6, j) → Wn8 g q * Hp 6 8 g j = Hb 6 9 g j
  hHZn_F1_8 : ∀ g j, (∀ q, hostAt 8 q ≠ some (6, j)) → Hb 6 9 g j = Hp 6 8 g j
  hHE_F1 : ∀ g t, GReach (Hb 6 9 g (false, t)) (kernel α) (fun φ => bcost φ (hrem 6 (false, t))) ∧
    GReach (Hb 6 9 g (true, t)) (YF g t) (fun φ => bcost φ (hrem 6 (true, t)))
  hH0_F2 : ∀ g t, Hb 7 5 g (false, t) = (G4 (r4 t g)).Xm.Φ 1 ∧ Hb 7 5 g (true, t) = (G4 (r4 t g)).Xm.Φ (1 + tt (I.tv t))
  hHP_F2_5 : ∀ g j, GReach (Hb 7 5 g j) (Hp 7 5 g j) (fun _ => 0)
  hHZ_F2_5 : ∀ g q j, hostAt 5 q = some (7, j) → Wn5 g q * Hp 7 5 g j = Hb 7 6 g j
  hHZn_F2_5 : ∀ g j, (∀ q, hostAt 5 q ≠ some (7, j)) → Hb 7 6 g j = Hp 7 5 g j
  hHP_F2_6 : ∀ g j, GReach (Hb 7 6 g j) (Hp 7 6 g j) (fun _ => 0)
  hHZ_F2_6 : ∀ g q j, hostAt 6 q = some (7, j) → Wn6 g q * Hp 7 6 g j = Hb 7 7 g j
  hHZn_F2_6 : ∀ g j, (∀ q, hostAt 6 q ≠ some (7, j)) → Hb 7 7 g j = Hp 7 6 g j
  hHP_F2_7 : ∀ g j, GReach (Hb 7 7 g j) (Hp 7 7 g j) (fun _ => 0)
  hHZ_F2_7 : ∀ g q j, hostAt 7 q = some (7, j) → Wn7 g q * Hp 7 7 g j = Hb 7 8 g j
  hHZn_F2_7 : ∀ g j, (∀ q, hostAt 7 q ≠ some (7, j)) → Hb 7 8 g j = Hp 7 7 g j
  hHP_F2_8 : ∀ g j, GReach (Hb 7 8 g j) (Hp 7 8 g j) (fun _ => 0)
  hHZ_F2_8 : ∀ g q j, hostAt 8 q = some (7, j) → Wn8 g q * Hp 7 8 g j = Hb 7 9 g j
  hHZn_F2_8 : ∀ g j, (∀ q, hostAt 8 q ≠ some (7, j)) → Hb 7 9 g j = Hp 7 8 g j
  hHE_F2 : ∀ g t, GReach (Hb 7 9 g (false, t)) (kernel α) (fun φ => bcost φ (hrem 7 (false, t))) ∧
    GReach (Hb 7 9 g (true, t)) (YF g t) (fun φ => bcost φ (hrem 7 (true, t)))
  hH0_F3 : ∀ g t, Hb 8 7 g (false, t) = (G6 (r6 t g)).Xm.Φ 1 ∧ Hb 8 7 g (true, t) = (G6 (r6 t g)).Xm.Φ (1 + tt (I.tv t))
  hHP_F3_7 : ∀ g j, GReach (Hb 8 7 g j) (Hp 8 7 g j) (fun _ => 0)
  hHZ_F3_7 : ∀ g q j, hostAt 7 q = some (8, j) → Wn7 g q * Hp 8 7 g j = Hb 8 8 g j
  hHZn_F3_7 : ∀ g j, (∀ q, hostAt 7 q ≠ some (8, j)) → Hb 8 8 g j = Hp 8 7 g j
  hHP_F3_8 : ∀ g j, GReach (Hb 8 8 g j) (Hp 8 8 g j) (fun _ => 0)
  hHZ_F3_8 : ∀ g q j, hostAt 8 q = some (8, j) → Wn8 g q * Hp 8 8 g j = Hb 8 9 g j
  hHZn_F3_8 : ∀ g j, (∀ q, hostAt 8 q ≠ some (8, j)) → Hb 8 9 g j = Hp 8 8 g j
  hHE_F3 : ∀ g t, GReach (Hb 8 9 g (false, t)) (kernel α) (fun φ => bcost φ (hrem 8 (false, t))) ∧
    GReach (Hb 8 9 g (true, t)) (YF g t) (fun φ => bcost φ (hrem 8 (true, t)))
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
def GHost401.hcost (c : Fin 9) (φ : ℕ → ℝ) : ℝ :=
  ∑ t : T, (bcost φ (N.hrem c (false, t)) + bcost φ (N.hrem c (true, t)))

/-- price of ONE run per element of `Γ`: 9 invocations, the remaining climbs, pair 4's final. -/
def GHost401.runUnit (φ : ℕ → ℝ) : ℝ :=
  9 * N.I.cost φ + (N.hcost 0 φ + N.hcost 1 φ + N.hcost 2 φ + N.hcost 3 φ + N.hcost 4 φ + N.hcost 5 φ
    + N.hcost 6 φ + N.hcost 7 φ + N.hcost 8 φ) + (Fintype.card T : ℝ) * (N.cX4 φ + N.cY4 φ)

/-- price of one unit (one element of `Γ`): 8 runs and the padding blocks of the banks. -/
def GHost401.unit (φ : ℕ → ℝ) : ℝ :=
  8 * N.runUnit φ + ∑ b, bcost φ (N.pad b)

/-- the conditions of the role map of stage `l`. -/
lemma GHost401.hok (ρ : Fin 8) (l : Fin 9) : HostOK (pairL l.val) (N.hostAt l) (N.bk ρ l) :=
  ⟨N.hinj ρ l, N.hostInj l, fun q c j h => win_pair c l (N.hostWin l q c j h),
    fun q q' c c' j j' h h' e => win_uniq c c' l (N.hostWin l q c j h) (N.hostWin l q' c' j' h') e⟩

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
def GHost401.S1 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi0 N.yi0 (fun g t => N.Hp 0 0 g (false, t)) (fun g t => N.Hp 0 0 g (true, t)) (fun g t => N.Hp 2 0 g (false, t)) (fun g t => N.Hp 2 0 g (true, t)) (fun g t => N.Hp 4 0 g (false, t)) (fun g t => N.Hp 4 0 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 0) g b) (cOne Γ C α)
def GHost401.S2 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi1 N.yi1 (fun g t => N.Hb 0 1 g (false, t)) (fun g t => N.Hb 0 1 g (true, t)) (fun g t => N.Hb 2 1 g (false, t)) (fun g t => N.Hb 2 1 g (true, t)) (fun g t => N.Hb 4 1 g (false, t)) (fun g t => N.Hb 4 1 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 1) g b) (cOne Γ C α)
def GHost401.S3 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi1 N.yi1 (fun g t => N.Hp 0 1 g (false, t)) (fun g t => N.Hp 0 1 g (true, t)) (fun g t => N.Hp 2 1 g (false, t)) (fun g t => N.Hp 2 1 g (true, t)) (fun g t => N.Hp 4 1 g (false, t)) (fun g t => N.Hp 4 1 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 1) g b) (cOne Γ C α)
def GHost401.S4 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi2 N.yi2 (fun g t => N.Hb 0 2 g (false, t)) (fun g t => N.Hb 0 2 g (true, t)) (fun g t => N.Hb 2 2 g (false, t)) (fun g t => N.Hb 2 2 g (true, t)) (fun g t => N.Hb 4 2 g (false, t)) (fun g t => N.Hb 4 2 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 2) g b) (cOne Γ C α)
def GHost401.S5 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi2 N.yi2 N.xi2 N.yi2 (fun g t => N.Hb 2 2 g (false, t)) (fun g t => N.Hb 2 2 g (true, t)) (fun g t => N.Hb 4 2 g (false, t)) (fun g t => N.Hb 4 2 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 2) g b) (cOne Γ C α)
def GHost401.S6 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xi2 N.yi2 (fun g t => N.Hp 1 2 g (false, t)) (fun g t => N.Hp 1 2 g (true, t)) (fun g t => N.Hp 2 2 g (false, t)) (fun g t => N.Hp 2 2 g (true, t)) (fun g t => N.Hp 4 2 g (false, t)) (fun g t => N.Hp 4 2 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 2) g b) (cOne Γ C α)
def GHost401.S7 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xo2 N.yo2 (fun g t => N.Hb 1 3 g (false, t)) (fun g t => N.Hb 1 3 g (true, t)) (fun g t => N.Hb 2 3 g (false, t)) (fun g t => N.Hb 2 3 g (true, t)) (fun g t => N.Hb 4 3 g (false, t)) (fun g t => N.Hb 4 3 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 3) g b) (cOne Γ C α)
def GHost401.S8 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 N.xo2 N.yo2 N.xi3 N.yi3 (fun g t => N.Hb 2 3 g (false, t)) (fun g t => N.Hb 2 3 g (true, t)) (fun g t => N.Hb 4 3 g (false, t)) (fun g t => N.Hb 4 3 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 3) g b) (cOne Γ C α)
def GHost401.S9 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp 6 3 g (false, t)) (fun g t => N.Hp 6 3 g (true, t)) N.xi3 N.yi3 (fun g t => N.Hp 2 3 g (false, t)) (fun g t => N.Hp 2 3 g (true, t)) (fun g t => N.Hp 4 3 g (false, t)) (fun g t => N.Hp 4 3 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 3) g b) (cOne Γ C α)
def GHost401.S10 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb 6 4 g (false, t)) (fun g t => N.Hb 6 4 g (true, t)) N.xi4 N.yi4 (fun g t => N.Hb 2 4 g (false, t)) (fun g t => N.Hb 2 4 g (true, t)) (fun g t => N.Hb 4 4 g (false, t)) (fun g t => N.Hb 4 4 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 4) g b) (cOne Γ C α)
def GHost401.S11 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb 6 4 g (false, t)) (fun g t => N.Hb 6 4 g (true, t)) N.xi4 N.yi4 N.xi4 N.yi4 (fun g t => N.Hb 4 4 g (false, t)) (fun g t => N.Hb 4 4 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 4) g b) (cOne Γ C α)
def GHost401.S12 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp 6 4 g (false, t)) (fun g t => N.Hp 6 4 g (true, t)) N.xi4 N.yi4 (fun g t => N.Hp 3 4 g (false, t)) (fun g t => N.Hp 3 4 g (true, t)) (fun g t => N.Hp 4 4 g (false, t)) (fun g t => N.Hp 4 4 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 4) g b) (cOne Γ C α)
def GHost401.S13 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb 6 5 g (false, t)) (fun g t => N.Hb 6 5 g (true, t)) N.xo4 N.yo4 (fun g t => N.Hb 3 5 g (false, t)) (fun g t => N.Hb 3 5 g (true, t)) (fun g t => N.Hb 4 5 g (false, t)) (fun g t => N.Hb 4 5 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 5) g b) (cOne Γ C α)
def GHost401.S14 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb 6 5 g (false, t)) (fun g t => N.Hb 6 5 g (true, t)) N.xo4 N.yo4 N.xi5 N.yi5 (fun g t => N.Hb 4 5 g (false, t)) (fun g t => N.Hb 4 5 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 5) g b) (cOne Γ C α)
def GHost401.S15 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp 6 5 g (false, t)) (fun g t => N.Hp 6 5 g (true, t)) (fun g t => N.Hp 7 5 g (false, t)) (fun g t => N.Hp 7 5 g (true, t)) N.xi5 N.yi5 (fun g t => N.Hp 4 5 g (false, t)) (fun g t => N.Hp 4 5 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 5) g b) (cOne Γ C α)
def GHost401.S16 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb 6 6 g (false, t)) (fun g t => N.Hb 6 6 g (true, t)) (fun g t => N.Hb 7 6 g (false, t)) (fun g t => N.Hb 7 6 g (true, t)) N.xi6 N.yi6 (fun g t => N.Hb 4 6 g (false, t)) (fun g t => N.Hb 4 6 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 6) g b) (cOne Γ C α)
def GHost401.S17 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb 6 6 g (false, t)) (fun g t => N.Hb 6 6 g (true, t)) (fun g t => N.Hb 7 6 g (false, t)) (fun g t => N.Hb 7 6 g (true, t)) N.xi6 N.yi6 N.xi6 N.yi6 (fun g b => N.Zb (9 * ρ.val + 6) g b) (cOne Γ C α)
def GHost401.S18 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp 6 6 g (false, t)) (fun g t => N.Hp 6 6 g (true, t)) (fun g t => N.Hp 7 6 g (false, t)) (fun g t => N.Hp 7 6 g (true, t)) N.xi6 N.yi6 (fun g t => N.Hp 5 6 g (false, t)) (fun g t => N.Hp 5 6 g (true, t)) (fun g b => N.Zp (9 * ρ.val + 6) g b) (cOne Γ C α)
def GHost401.S19 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb 6 7 g (false, t)) (fun g t => N.Hb 6 7 g (true, t)) (fun g t => N.Hb 7 7 g (false, t)) (fun g t => N.Hb 7 7 g (true, t)) N.xo6 N.yo6 (fun g t => N.Hb 5 7 g (false, t)) (fun g t => N.Hb 5 7 g (true, t)) (fun g b => N.Zb (9 * ρ.val + 7) g b) (cOne Γ C α)
def GHost401.S20 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb 6 7 g (false, t)) (fun g t => N.Hb 6 7 g (true, t)) (fun g t => N.Hb 7 7 g (false, t)) (fun g t => N.Hb 7 7 g (true, t)) N.xo6 N.yo6 N.xi7 N.yi7 (fun g b => N.Zb (9 * ρ.val + 7) g b) (cOne Γ C α)
def GHost401.S21 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp 6 7 g (false, t)) (fun g t => N.Hp 6 7 g (true, t)) (fun g t => N.Hp 7 7 g (false, t)) (fun g t => N.Hp 7 7 g (true, t)) (fun g t => N.Hp 8 7 g (false, t)) (fun g t => N.Hp 8 7 g (true, t)) N.xi7 N.yi7 (fun g b => N.Zp (9 * ρ.val + 7) g b) (cOne Γ C α)
def GHost401.S22 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb 6 8 g (false, t)) (fun g t => N.Hb 6 8 g (true, t)) (fun g t => N.Hb 7 8 g (false, t)) (fun g t => N.Hb 7 8 g (true, t)) (fun g t => N.Hb 8 8 g (false, t)) (fun g t => N.Hb 8 8 g (true, t)) N.xi8 N.yi8 (fun g b => N.Zb (9 * ρ.val + 8) g b) (cOne Γ C α)
def GHost401.S23 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hp 6 8 g (false, t)) (fun g t => N.Hp 6 8 g (true, t)) (fun g t => N.Hp 7 8 g (false, t)) (fun g t => N.Hp 7 8 g (true, t)) (fun g t => N.Hp 8 8 g (false, t)) (fun g t => N.Hp 8 8 g (true, t)) N.xi8 N.yi8 (fun g b => N.Zp (9 * ρ.val + 8) g b) (cOne Γ C α)
def GHost401.S24 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun g t => N.Hb 6 9 g (false, t)) (fun g t => N.Hb 6 9 g (true, t)) (fun g t => N.Hb 7 9 g (false, t)) (fun g t => N.Hb 7 9 g (true, t)) (fun g t => N.Hb 8 9 g (false, t)) (fun g t => N.Hb 8 9 g (true, t)) N.xo8 N.yo8 (fun g b => N.Zb (9 * ρ.val + 9) g b) (cOne Γ C α)
def GHost401.S25 (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF (fun g b => N.Zb (9 * ρ.val + 9) g b) (cOne Γ C α)

/-- frames of all roles at the start of run `k` (`k = 8`: after the last run). -/
def GHost401.Sfull (k : Fin 9) : BRole401 Γ T nb C → CMat α :=
  mk401 (fun ρ' d => if ρ'.val < k.val then N.FinD d else N.InitD d) (fun g b => N.Zb (9 * k.val) g b) (cOne Γ C α)

/-- frames of all roles after the padding climbs. -/
def GHost401.Sfin : BRole401 Γ T nb C → CMat α :=
  mk401 (fun _ d => N.FinD d) (fun _ _ => kernel α) (cOne Γ C α)

end Net
end
end XHost
end PowerSaving
end OAI

