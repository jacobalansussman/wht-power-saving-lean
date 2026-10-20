import Work.PfLean.Dem3Data

/-!
# Definitions for the kernel-checked computation of the block-3 hook demand (px-lean, 14 October 2026)

For each of the 8 characters `s` of `F₂³`, the `28 × 28` integer matrix `Gf s` (entries computed here
from the closed formula, with `i ↦ 256`, a square root of `-1` mod `65537`) times the certificate
`Bf s` (data from `Dem3Data`, untrusted) is the identity mod `65537`.  The checks `blkOK*` are
evaluated by the kernel (`decide +kernel`); no `native_decide`, no extra axiom.
-/

namespace OAI.PowerSaving.PfLean.Dem3Comp
open Dem3Data

/-- `(-256)^(popcount z)` for `z < 8`. -/
def pwZ (z : ℕ) : ℤ := [1, -256, -256, 65536, -256, 65536, 65536, -16777216].getD z 0
/-- `(-1)^(popcount (s &&& u))` for `s, u < 8`. -/
def chiZ (s u : ℕ) : ℤ := [1, -1, -1, 1, -1, 1, 1, -1].getD (s &&& u) 0

/-- The `4 × 4` determinant, by Laplace expansion on rows `0, 1`. -/
def det4 {R : Type*} [CommRing R] (m : Fin 4 → Fin 4 → R) : R :=
  (m 0 0 * m 1 1 - m 0 1 * m 1 0) * (m 2 2 * m 3 3 - m 2 3 * m 3 2)
  - (m 0 0 * m 1 2 - m 0 2 * m 1 0) * (m 2 1 * m 3 3 - m 2 3 * m 3 1)
  + (m 0 0 * m 1 3 - m 0 3 * m 1 0) * (m 2 1 * m 3 2 - m 2 2 * m 3 1)
  + (m 0 1 * m 1 2 - m 0 2 * m 1 1) * (m 2 0 * m 3 3 - m 2 3 * m 3 0)
  - (m 0 1 * m 1 3 - m 0 3 * m 1 1) * (m 2 0 * m 3 2 - m 2 2 * m 3 0)
  + (m 0 2 * m 1 3 - m 0 3 * m 1 2) * (m 2 0 * m 3 1 - m 2 1 * m 3 0)

def sum8 {R : Type*} [AddCommMonoid R] (f : ℕ → R) : R :=
  f 0 + f 1 + f 2 + f 3 + f 4 + f 5 + f 6 + f 7

/-- The block entry for character `s`, row `(0, us)`, column `(0, 0, ts)`. -/
def GenZ (s : ℕ) (us ts : List ℕ) : ℤ :=
  sum8 fun u => chiZ s u * pwZ u *
    det4 (fun b a => pwZ (us.getD b.val 0 ^^^ u ^^^ (0 :: ts).getD a.val 0))

/-- Certificate entry `(i, k)` of block `s`. -/
def Bf (s i k : ℕ) : ℤ := (((invD.getD s []).getD i []).getD k 0 : ℕ)
/-- Matrix entry `(k, j)` of block `s`, computed from the formula. -/
def Gf (s k j : ℕ) : ℤ := GenZ s ((rowsD.getD s []).getD k []) ((colsD.getD s []).getD j [])

def dotN (f g : ℕ → ℤ) : ℕ → ℤ
  | 0 => 0
  | n + 1 => dotN f g n + f n * g n

/-- The check for one block. -/
abbrev BlockCheck (s : ℕ) : Prop :=
  ∀ i < 28, ∀ j < 28, (dotN (Bf s i) (fun k => Gf s k j) 28 - if i = j then 1 else 0) % 65537 = 0

end OAI.PowerSaving.PfLean.Dem3Comp
