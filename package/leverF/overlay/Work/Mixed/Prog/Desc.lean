import Work.Mixed.Prog.Exit

/-!
# Mixed engine, R3 program part 6: the engine data, the descriptors and their calls (agent e-r3prog)

* `Eng`        the fixed data: the word `wd` (a `GWord` with targets `blockMat (tz l)` on the live roles,
               rl-prove Theorem 1G's rule), layouts, the separation hypothesis;
* `MDesc`       ONE descriptor type for the handler: top jobs `Do(n, B)` and classes of the lockstep tree;
* `Eng.enc`, `Eng.out`   the input and the answer of every call (port of `Block.lean`);
* `Eng.V0`     the class table loaded from a class's own jobs (live roles through the target splits,
               scratch roles zero);
* `Eng.kidsC`  the calls of a class: one child class per block letter in word order (one exit for a leaf
               class), each with its input as the program builds it after the letters before it;
* `Eng.srcSizes`  the sizes of the supply vectors a class's subtree consumes, in consumption order.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF MXR
noncomputable section

/-- One descriptor type for all calls of the lockstep program. -/
inductive MDesc
  | top (n B : ℕ) (v : ℕ → Bits n → ℂ)
  | cls (n h μ o BS : ℕ) (J sup : ℕ → Bits n → ℂ)

/-- The fixed data of the mixed engine. -/
structure Eng (α ρ σ : Type*) [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
    [Fintype σ] [DecidableEq σ] where
  X : Layout α
  R : Layout ρ
  Q : Layout σ
  e : σ → ρ
  he : Function.Injective e
  tr : σ → ℕ
  tz : ∀ l, Fin (tr l) → Space α
  tind : ∀ l, LinearIndependent F (tz l)
  tr_le : ∀ l, tr l ≤ X.n
  tr_pos : ∀ l, 1 ≤ tr l
  wd : GWord α ρ
  hP : wd.Proper X.n
  /-- separation (rl-prove (S)): every paid rank is below every target rank -/
  sep : ∀ mv ∈ wd, ∀ r, mv.rk? = some r → ∀ l, r < tr l
  W_pos : 0 < Q.n
  three : 3 ≤ X.n
  /-- the mixed certificate (`GLiveKernelM e (fun l => blockMat (tz l)) wd` of rb-r45m1, same body) -/
  hw : ∀ (f : ℕ) (v : Data α ρ f), (∀ s, (∀ l, e l ≠ s) → v s = 0) →
    ∀ l, gwalk f wd v (e l) = matAct f (blockMat (tz l)) (v (e l))

namespace Eng
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]
variable (E : Eng α ρ σ)

/-- Number of live roles. -/
abbrev W : ℕ := E.Q.n
/-- Target rank of the live role at position `a`. -/
def tq (a : ℕ) : ℕ := if h : a < E.Q.n then E.tr (E.Q.atFin ⟨a, h⟩) else 0
/-- Axes of an array. -/
abbrev mx : ℕ := E.X.n

/-! ### the units of a class -/
abbrev NUc (n h μ : ℕ) : ℕ := MXR.NU E.W E.tq E.mx n h μ
abbrev Uc (n h μ p : ℕ) : ℕ := MXR.U E.W E.tq E.mx n h μ p
abbrev PCc (n h μ t : ℕ) : ℕ := MXR.PC E.W E.tq E.mx n h μ t
abbrev cUc (n h μ p : ℕ) : ℕ := MXR.cU E.W E.tq E.mx n h μ p
abbrev iUc (n h μ p : ℕ) : ℕ := MXR.iU E.W E.tq E.mx n h μ p

lemma fitc (n h μ : ℕ) (p : Fin (E.NUc n h μ)) : E.X.n * E.cUc n h μ p + E.iUc n h μ p = n :=
  MXR.iU_add p.2

/-- Size of own slot `t = p*W + a` of a class. -/
def slotSize (n h μ t : ℕ) : ℕ := E.tq (t % E.W) * E.cUc n h μ (t / E.W)

/-- The class state: element data of every unit and instance. -/
abbrev CState (n h μ : ℕ) := ∀ p : Fin (E.NUc n h μ), Bits (E.iUc n h μ p) → Grain α ρ 0 (E.cUc n h μ p)

/-! ### target splits of the live roles -/
/-- The target splitting of live role `l`. -/
def St (l : σ) : Split α (Fin (Fintype.card α - E.tr l)) (Fin (E.tr l)) :=
  Classical.choose (exists_split (E.tz l) (E.tind l))
lemma St_spec (l : σ) (i : Fin (E.tr l)) : (E.St l).e.symm (0, eu i) = E.tz l i :=
  Classical.choose_spec (exists_split (E.tz l) (E.tind l)) i
/-- Its rest layout. -/
def Xt (l : σ) : Layout (Fin (Fintype.card α - E.tr l)) := Layout.someLayout _
lemma Xt_n (l : σ) : (E.Xt l).n + E.tr l = E.X.n := by
  have e1 : (E.Xt l).n = Fintype.card α - E.tr l := by
    rw [← Layout.card (E.Xt l)]; simp
  have e2 : Fintype.card α = E.X.n := Layout.card E.X
  have := E.tr_le l
  omega

/-- The live role behind a role, if any. -/
def liveOf (l' : ρ) : Option σ := if h : ∃ s0, E.e s0 = l' then some (Classical.choose h) else none

lemma liveOf_e (s0 : σ) : E.liveOf (E.e s0) = some s0 := by
  have h : ∃ s1, E.e s1 = E.e s0 := ⟨s0, rfl⟩
  unfold liveOf
  rw [dif_pos h]
  exact congrArg some (E.he (Classical.choose_spec h))

lemma liveOf_none (l' : ρ) (h : ∀ s0, E.e s0 ≠ l') : E.liveOf l' = none := by
  unfold liveOf; rw [dif_neg (fun ⟨s0, hs⟩ => h s0 hs)]

/-- The placement of a live role's job in the instances of unit `p`. -/
def placeT (n h μ : ℕ) (p : Fin (E.NUc n h μ)) (l : σ) :
    Bits (E.iUc n h μ p) × Sky α (E.cUc n h μ p) ≃ Bits n :=
  placeS (E.St l) (E.Xt l) (E.cUc n h μ p) (E.iUc n h μ p) n (fitk (E.Xt l) (E.Xt_n l) (E.fitc n h μ p))

/-- **The class table loaded from a class's own jobs**: live role `e l` of instance `q` of unit
`p` reads job `p*W + loc l` through the target placement; scratch roles are zero. -/
def V0 (n h μ : ℕ) (J : ℕ → Bits n → ℂ) : E.CState n h μ :=
  fun p q x => match E.liveOf x.1 with
    | some l => J (p.val * E.W + E.Q.loc l) (E.placeT n h μ p l (q, x.2.2))
    | none => 0

/-- One letter on every element of a class. -/
def stepV {n h μ : ℕ} (V : E.CState n h μ) (mv : GMove α ρ) : E.CState n h μ :=
  fun p q => gferry 0 (E.cUc n h μ p) mv (V p q)

/-- The whole word on every element. -/
def runV {n h μ : ℕ} (V : E.CState n h μ) (w : GWord α ρ) : E.CState n h μ :=
  fun p q => gcoast 0 (E.cUc n h μ p) w (V p q)

/-! ### routing of a block letter of rank `r` from class `(n, h+1, μ)` to its child `(n, h, μ*r)` -/

/-- Code of the parent unit feeding child slot `t`. -/
def parCode (n h μ r t : ℕ) : ℕ := t % E.W + E.W * E.Uc n h (μ * r) (t / E.W)
/-- Child slot `t` is fed (by a parent unit). -/
def fedP (n h μ r t : ℕ) : Prop :=
  t < E.NUc n h (μ * r) * E.W ∧ inRS E.W E.tq E.mx n (h+1) μ (E.parCode n h μ r t)
instance (n h μ r : ℕ) : DecidablePred (E.fedP n h μ r) := fun t => by unfold fedP; infer_instance
/-- Child slot `t` is a source. -/
def srcP (n h μ r t : ℕ) : Prop := t < E.NUc n h (μ * r) * E.W ∧ ¬ inRS E.W E.tq E.mx n (h+1) μ (E.parCode n h μ r t)
instance (n h μ r : ℕ) : DecidablePred (E.srcP n h μ r) := fun t => by unfold srcP; infer_instance

def flR (n h μ r t : ℕ) : ℕ := if E.fedP n h μ r t then 1 else 0
def ppR (n h μ r t : ℕ) : ℕ := E.PCc n (h+1) μ (E.parCode n h μ r t)
/-- Supply index of a source slot: running offset plus the number of source slots before it. -/
def siR (n h μ r o t : ℕ) : ℕ := o + Nat.count (E.srcP n h μ r) t
/-- Child slot of a parent unit. -/
def slotR (n h μ r p : ℕ) : ℕ :=
  E.PCc n h (μ * r) (E.Uc n (h+1) μ p / E.W) * E.W + E.Uc n (h+1) μ p % E.W

lemma ppR_lt (n h μ r t : ℕ) (hf : E.flR n h μ r t ≠ 0) : E.ppR n h μ r t < E.NUc n (h+1) μ := by
  unfold flR at hf
  split_ifs at hf with h1
  · exact PC_lt h1.2
  · exact absurd rfl hf

/-! ### sizes of consumed sources -/

/-- Ranks of the block letters of a word, in order. -/
def blockRanks (w : GWord α ρ) : List ℕ := w.filterMap GMove.rk?

/-- Sizes of the source slots of the child `(n, h, μ*r)` of a block letter of rank `r`, in slot order. -/
def ownSrc (n h μ r : ℕ) : List ℕ :=
  (List.range (E.NUc n h (μ * r) * E.W)).filterMap
    (fun t => if E.srcP n h μ r t then some (E.slotSize n h (μ * r) t) else none)

/-- Sizes of the supply vectors consumed by the subtree of a class (not its own slots), in order. -/
def srcSizes (n : ℕ) : ℕ → ℕ → List ℕ
  | 0, _ => []
  | h+1, μ => (blockRanks E.wd).flatMap (fun r => E.ownSrc n h μ r ++ srcSizes n h (μ * r))

def cons (n h μ : ℕ) : ℕ := (E.srcSizes n h μ).length
def csz (n h μ k : ℕ) : ℕ := (E.srcSizes n h μ).getD k 0

/-! ### the calls of a class -/

/-- The block letter's split (as in `open_moveB`). -/
def bS {r : ℕ} (z : Fin r → Space α) (ind : LinearIndependent F z) :
    Split α (Fin (Fintype.card α - r)) (Fin r) := Classical.choose (exists_split z ind)
def bX (r : ℕ) : Layout (Fin (Fintype.card α - r)) := Layout.someLayout _
lemma bX_n {r : ℕ} (hr : r ≤ E.X.n) : (bX (α:=α) r).n + r = E.X.n := by
  have e1 : (bX (α:=α) r).n = Fintype.card α - r := by
    rw [← Layout.card (bX (α:=α) r)]; simp
  have e2 : Fintype.card α = E.X.n := Layout.card E.X
  omega

/-- The child class of a block letter (class `(n, h+1, μ)`, state `V`, running offset `o`). -/
def childC (n h μ BS : ℕ) (sup : ℕ → Bits n → ℂ) (V : E.CState n (h+1) μ) (o : ℕ)
    (l : ρ) (r : ℕ) (z : Fin r → Space α) (ind : LinearIndependent F z) (hr : r ≤ E.X.n) : MDesc :=
  MDesc.cls n h (μ * r) (o + (E.ownSrc n h μ r).length) BS
    (cjob E.X (bS z ind) (bX (α:=α) r) (E.bX_n hr) l (E.fitc n (h+1) μ) V (E.flR n h μ r) (E.ppR n h μ r)
      (E.siR n h μ r o) (E.ppR_lt n h μ r) sup) sup

/-- The exit of a block letter of a leaf class (class `(n, 0, μ)`, state `V`). -/
def exitC (n μ : ℕ) (V : E.CState n 0 μ) (l : ρ) (r : ℕ) (z : Fin r → Space α)
    (ind : LinearIndependent F z) (hr : r ≤ E.X.n) : MDesc :=
  if hk : r * μ ≤ n then
    MDesc.top (r * μ) (E.NUc n 0 μ * 2^(n - r * μ))
      (toVec n (r * μ) hk (exJobs E.X (bS z ind) (bX (α:=α) r) (E.bX_n hr) l (E.fitc n 0 μ) V))
  else MDesc.top 0 0 (fun _ _ => 0)

/-- The calls of a class of height `h+1` along a word: state and running offset are threaded. -/
def kidsUp (n h μ BS : ℕ) (sup : ℕ → Bits n → ℂ) : E.CState n (h+1) μ → ℕ → GWord α ρ → List MDesc
  | _, _, [] => []
  | V, o, mv :: ws => match mv with
    | .old (.block l r z ind) =>
      (if hr : r ≤ E.X.n then [E.childC n h μ BS sup V o l r z ind hr] else []) ++
        kidsUp n h μ BS sup (E.stepV V mv) (o + (E.ownSrc n h μ r).length + E.cons n h (μ * r)) ws
    | _ => kidsUp n h μ BS sup (E.stepV V mv) o ws

/-- The calls of a leaf class along a word. -/
def kidsLeaf (n μ : ℕ) : E.CState n 0 μ → GWord α ρ → List MDesc
  | _, [] => []
  | V, mv :: ws => match mv with
    | .old (.block l r z ind) =>
      (if hr : r ≤ E.X.n then [E.exitC n μ V l r z ind hr] else []) ++ kidsLeaf n μ (E.stepV V mv) ws
    | _ => kidsLeaf n μ (E.stepV V mv) ws

/-- **The calls of a class.** -/
def kidsC : MDesc → List MDesc
  | MDesc.top _ _ _ => []
  | MDesc.cls n 0 μ _ _ J _ => E.kidsLeaf n μ (E.V0 n 0 μ J) E.wd
  | MDesc.cls n (h+1) μ o BS J sup => E.kidsUp n h μ BS sup (E.V0 n (h+1) μ J) o E.wd

/-! ### input and answer of every call -/

/-- Input of a call. -/
def enc (cl : Paint) : MDesc → T (portIn cl)
  | MDesc.top n B v => mkIn cl n 0 0 0 (jobsTape n B v) (jobsTape n 0 (fun _ _ => 0))
  | MDesc.cls n h μ o BS J sup => mkIn cl n (h+1) μ o (jobsTape n (E.NUc n h μ * E.W) J) (jobsTape n BS sup)

/-- Answer of a call. -/
def out (cl : Paint) : MDesc → T (portOut cl)
  | MDesc.top n B v => (jobsTape n B (fun t => ripple n (v t)), jobsTape n 0 (fun _ _ => 0),
      (Layout.std 0).tape (fun _ => 0))
  | MDesc.cls n h μ o BS J sup =>
      (jobsTape n (E.NUc n h μ * E.W) (fun t => hlow n (E.slotSize n h μ t) (J t)),
       jobsTape n (E.cons n h μ) (fun k => hlow n (E.csz n h μ k) (sup (o + k))),
       (Layout.std (E.cons n h μ)).tape (fun k => E.csz n h μ k))

end Eng
end
end PowerSaving.RAM
end OAI
