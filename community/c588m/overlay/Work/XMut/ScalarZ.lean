import Work.GCert.Scalar.Main
import Work.XMut.Scal

/-!
# x-mutable (13 Oct 2026): THE CONCRETE SCALAR CHECK WITH MUTABLE SOURCES (`ScalOKZ`)

Jacob Sussman's kernel scalar check (`Work/GCert/Scalar/{Seg,Runs,Check,Cols,Main}.lean`) with two changes only:
* the gate kinds: `okAZK` (phase A: no target named; sources MAY be written, from sources or slots) and `okBZK`
  (phase B: a target is read only by a target; sources MAY be written) instead of `okAK` / `okBK`;
* the final test of the x-column replay checks the y rows only (`hid`); the x rows are not required to come
  back (no `hx`): the mutable-source invocation theorem (`XMut/GInv.lean`) restores them by its derived gate.
The kernel replay machinery (`gGates`, `gScat`, `gChunks`, `gFin`, their soundness lemmas) is reused unchanged:
it is generic in the kind predicate.  Output: `ScalOKZ.toScalZ : GX.ScalZ` (`hy`, `hid`, `hAi`, `hBi`).
No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
set_option linter.unusedSimpArgs false
set_option linter.deprecated false

namespace GS
open SSC GXD OAI.PowerSaving OAI.PowerSaving.Binary OAI.PowerSaving.RAM OAI.PowerSaving.SS
  OAI.PowerSaving.CB OAI.PowerSaving.CR Finset Matrix

/-- kinds of phase A, sources mutable: both registers in range, neither a target -/
def okAZK (v v2 nr t s : Nat) : Bool :=
  Nat.ble (Nat.succ t) nr && (Nat.ble (Nat.succ s) nr &&
    ((Nat.ble (Nat.succ t) v || Nat.ble v2 t) && (Nat.ble (Nat.succ s) v || Nat.ble v2 s)))

/-- kinds of phase B, sources mutable: both registers in range, a target is read only by a target -/
def okBZK (v v2 nr t s : Nat) : Bool :=
  Nat.ble (Nat.succ t) nr && (Nat.ble (Nat.succ s) nr &&
    (Nat.ble (Nat.succ s) v || (Nat.ble v2 s || (Nat.ble v t && Nat.ble (Nat.succ t) v2))))

theorem okAZ_spec {v v2 nr t s : Nat} (h : okAZK v v2 nr t s = true) :
    t < nr ∧ s < nr ∧ (t < v ∨ v2 ≤ t) ∧ (s < v ∨ v2 ≤ s) := by
  simp only [okAZK, Bool.and_eq_true, Bool.or_eq_true, Nat.ble_eq] at h
  omega

theorem okBZ_spec {v v2 nr t s : Nat} (h : okBZK v v2 nr t s = true) :
    t < nr ∧ s < nr ∧ (s < v ∨ v2 ≤ s ∨ (v ≤ t ∧ t < v2)) := by
  simp only [okBZK, Bool.and_eq_true, Bool.or_eq_true, Nat.ble_eq] at h
  omega

abbrev okAZ (c : Raw) : Nat → Nat → Bool := okAZK c.v (2 * c.v) (2 * c.v + c.R)
abbrev okBZ (c : Raw) : Nat → Nat → Bool := okBZK c.v (2 * c.v) (2 * c.v + c.R)

/-! ### kernel terms -/
section
variable (c : Raw) (p : SPar) (v2 nr tm K m2 rl : Nat)

noncomputable def gItemZ (it : Item) (tr : Trie) (k : Trie → Bool) : Bool :=
  Item.rec (motive := fun _ => Bool)
    (fun gs => gGates tm K m2 p.ux p.us p.uy c.v v2 (okAZK c.v v2 nr) (fun _ => true) gs tr k)
    (gScat tm K m2 p.ux p.us p.uy c.v v2 (okSK c.v v2 nr) c.ret rl (scatRows c) c.v tr k)
    (fun gs => gGates tm K m2 p.ux p.us p.uy c.v v2 (okBZK c.v v2 nr) (fun _ => true) gs tr k) it

noncomputable def gItemsZ : List Item → Trie → (Trie → Bool) → Bool :=
  List.rec (motive := fun _ => Trie → (Trie → Bool) → Bool)
    (fun tr k => k tr)
    (fun it _ ih tr k => gItemZ c p v2 nr tm K m2 rl it tr fun tr' => ih tr' k)
end

/-- **a segment of the scalar check, sources mutable**: as `segCheck`, the final test on the y rows only. -/
noncomputable def segCheckZ (c : Raw) (p : SPar) (lo n : Nat) (pre : Option Trie) (items : List Item)
    (post : Option Trie) : Bool :=
  gPre c p lo n fun v2 nr hi tm K m2 => forceN (List.length c.ret) fun rl =>
    Option.rec (motive := fun _ => (Trie → Bool) → Bool)
      (fun k => gInit p.sw p.ux lo n 0 (Trie.mk p.d) k) (fun t k => k t) pre fun tr0 =>
        gItemsZ c p v2 nr tm K m2 rl items tr0 fun tr3 =>
          Option.rec (motive := fun _ => Bool)
            (gFin p.sw K m2 p.uy c.v lo hi tr3 c.v 0)
            (fun t => beqK tr3 t) post

noncomputable def scalCheckZ (c : Raw) (p : SPar) (lo n : Nat) : Bool :=
  segCheckZ c p lo n none (prog c) none

/-- the y check with the phase-B kinds of mutable sources -/
noncomputable def yCheckZ (c : Raw) (p : SPar) (lo n : Nat) : Bool :=
  gPre c p lo n fun v2 nr hi tm K m2 =>
    gInit p.sw p.uy (Nat.add c.v lo) n 0 (Trie.mk p.d) fun tr0 =>
      gChunks tm K m2 p.ux p.us p.uy c.v v2 (okBZK c.v v2 nr) (selY c.v v2) c.B tr0 fun tr3 =>
        gFin p.sw K m2 p.uy c.v lo hi tr3 c.v 0

/-! ### propositions -/
section
variable (c : Raw) (p : SPar) (lo n : Nat)

def ItemRunZ : Item → SSt → SSt → Prop
  | .cA gs, a, b => GRun p.sw n (U c p) (okAZ c) (fun _ => true) a (gs.flatMap Gate.adds) b
  | .scat, a, b => GRun p.sw n (U c p) (okS c) (fun _ => true) a (scatAdds c) b ∧ ScatShape c
  | .cB gs, a, b => GRun p.sw n (U c p) (okBZ c) (fun _ => true) a (gs.flatMap Gate.adds) b

def ItemsRunZ : List Item → SSt → SSt → Prop
  | [], a, b => b = a
  | it :: l, a, b => ∃ m, ItemRunZ c p n it a m ∧ ItemsRunZ l m b

def PostStZ : Option Trie → SSt → Prop
  | none, b => SmallSt p.sw n b ∧ FinSt c p lo n c.v p.uy b
  | some t, b => GOk p.d p.sw n t ∧ b = absG p.d (64 + p.sw * n) t

def SegOKZ (pre : Option Trie) (items : List Item) (post : Option Trie) : Prop :=
  PGood c p lo n ∧
    match pre with
    | none => ∃ a b, InitSt p n lo p.ux a ∧ ItemsRunZ c p n items a b ∧ PostStZ c p lo n post b
    | some t => GOk p.d p.sw n t →
        ∃ b, ItemsRunZ c p n items (absG p.d (64 + p.sw * n) t) b ∧ PostStZ c p lo n post b

theorem ItemsRunZ_append (l1 l2 : List Item) : ∀ (a m b : SSt),
    ItemsRunZ c p n l1 a m → ItemsRunZ c p n l2 m b → ItemsRunZ c p n (l1 ++ l2) a b := by
  induction l1 with
  | nil =>
    intro a m b h1 h2
    have e : m = a := h1
    rw [e] at h2
    exact h2
  | cons x l1 ih =>
    intro a m b h1 h2
    obtain ⟨m1, hs, hr⟩ := h1
    exact ⟨m1, hs, ih m1 m b hr h2⟩

theorem ItemsRunZ_split (l1 l2 : List Item) : ∀ (a b : SSt),
    ItemsRunZ c p n (l1 ++ l2) a b → ∃ m, ItemsRunZ c p n l1 a m ∧ ItemsRunZ c p n l2 m b := by
  induction l1 with
  | nil => intro a b h; exact ⟨a, rfl, h⟩
  | cons x l1 ih =>
    intro a b h
    obtain ⟨m1, hs, hr⟩ := h
    obtain ⟨m, r1, r2⟩ := ih m1 b hr
    exact ⟨m, ⟨m1, hs, r1⟩, r2⟩

theorem SegOKZ.append' {pre post : Option Trie} {l1 l2 : List Item} {T : Trie}
    (h1 : SegOKZ c p lo n pre l1 (some T)) (h2 : SegOKZ c p lo n (some T) l2 post) :
    SegOKZ c p lo n pre (l1 ++ l2) post := by
  obtain ⟨hp, h1⟩ := h1
  obtain ⟨_, h2⟩ := h2
  refine ⟨hp, ?_⟩
  cases pre with
  | none =>
    obtain ⟨a, b, hi, hr, hT, eb⟩ := h1
    obtain ⟨b', hr', hpost⟩ := h2 hT
    rw [← eb] at hr'
    exact ⟨a, b', hi, ItemsRunZ_append c p n l1 l2 a b b' hr hr', hpost⟩
  | some t =>
    intro ht
    obtain ⟨b, hr, hT, eb⟩ := h1 ht
    obtain ⟨b', hr', hpost⟩ := h2 hT
    rw [← eb] at hr'
    exact ⟨b', ItemsRunZ_append c p n l1 l2 _ b b' hr hr', hpost⟩

theorem ItemsRunZ_cA (cs : List (List Gate)) : ∀ (a b : SSt),
    ItemsRunZ c p n (cs.map Item.cA) a b →
    GRun p.sw n (U c p) (okAZ c) (fun _ => true) a (cs.flatten.flatMap Gate.adds) b := by
  induction cs with
  | nil => intro a b h; exact h
  | cons gs cs ih =>
    intro a b h
    obtain ⟨m, h1, h2⟩ := h
    rw [List.flatten_cons, List.flatMap_append]
    exact GRun_append _ _ _ _ _ _ _ _ _ _ h1 (ih m b h2)

theorem ItemsRunZ_cB (cs : List (List Gate)) : ∀ (a b : SSt),
    ItemsRunZ c p n (cs.map Item.cB) a b →
    GRun p.sw n (U c p) (okBZ c) (fun _ => true) a (cs.flatten.flatMap Gate.adds) b := by
  induction cs with
  | nil => intro a b h; exact h
  | cons gs cs ih =>
    intro a b h
    obtain ⟨m, h1, h2⟩ := h
    rw [List.flatten_cons, List.flatMap_append]
    exact GRun_append _ _ _ _ _ _ _ _ _ _ h1 (ih m b h2)
end

theorem SegOKZ.append {c : Raw} {p : SPar} {lo n : Nat} {pre post : Option Trie}
    {l1 l2 : List Item} {T : Trie} (h1 : SegOKZ c p lo n pre l1 (some T))
    (h2 : SegOKZ c p lo n (some T) l2 post) : SegOKZ c p lo n pre (l1 ++ l2) post :=
  SegOKZ.append' c p lo n h1 h2

/-- **an accepted replay of the main phase, sources mutable**, for the sources `lo .. lo+n-1` -/
def XRunZ (c : Raw) (lo n : Nat) : Prop :=
  ∃ (p : SPar) (a0 a1 a2 a3 : SSt), PGood c p lo n ∧ InitSt p n lo p.ux a0 ∧
    GRun p.sw n (U c p) (okAZ c) (fun _ => true) a0 (addsA c) a1 ∧
    GRun p.sw n (U c p) (okS c) (fun _ => true) a1 (scatAdds c) a2 ∧ ScatShape c ∧
    GRun p.sw n (U c p) (okBZ c) (fun _ => true) a2 (addsB c) a3 ∧
    SmallSt p.sw n a3 ∧ FinSt c p lo n c.v p.uy a3

def YRunZ (c : Raw) (lo n : Nat) : Prop :=
  ∃ (p : SPar) (a0 a3 : SSt), PGood c p lo n ∧ InitSt p n (c.v + lo) p.uy a0 ∧
    GRun p.sw n (U c p) (okBZ c) (selY c.v (2 * c.v)) a0 (addsB c) a3 ∧
    SmallSt p.sw n a3 ∧ FinSt c p lo n c.v p.uy a3

theorem XRunZ.of_seg {c : Raw} {p : SPar} {lo n : Nat} (h : SegOKZ c p lo n none (prog c) none) :
    XRunZ c lo n := by
  obtain ⟨hp, a, b, hi, hr, hsm, hfy⟩ := h
  obtain ⟨m1, r1, r2⟩ := ItemsRunZ_split c p n _ _ a b hr
  obtain ⟨m2, ⟨rs, hsh⟩, r3⟩ := r2
  exact ⟨p, a, m1, m2, b, hp, hi, ItemsRunZ_cA c p n _ a m1 r1, rs, hsh,
    ItemsRunZ_cB c p n _ m2 b r3, hsm, hfy⟩

/-! ### the kernel checks give the replays -/
section
variable (c : Raw) (p : SPar) (lo n : Nat) (hp : PGood c p lo n)
include hp

theorem gItemZ_sound (it : Item) (tr : Trie) (k : Trie → Bool) (hI : GOk p.d p.sw n tr)
    (h : gItemZ c p (2 * c.v) (2 * c.v + c.R) (maskW p.sw n) (64 + p.sw * n) (2^(64 + p.sw * n))
      c.ret.length it tr k = true) :
    ∃ tr', GOk p.d p.sw n tr' ∧ k tr' = true ∧
      ItemRunZ c p n it (absG p.d (64 + p.sw * n) tr) (absG p.d (64 + p.sw * n) tr') := by
  have hu := U_pos c p lo n hp
  cases it with
  | cA gs =>
    exact gGates_sound p.d p.sw n _ _ p.ux p.us p.uy c.v (2 * c.v) (okAZ c) (fun _ => true) hp.hsw rfl
      rfl hu gs tr k hI h
  | scat =>
    obtain ⟨tr', hI', hk, hlen, hb, hr⟩ := gScat_sound p.d p.sw n _ _ p.ux p.us p.uy c.v (2 * c.v)
      (okS c) hp.hsw rfl rfl hu c.ret c.ret.length (scatRows c) c.v tr k hI h
    exact ⟨tr', hI', hk, hr, hlen, hb⟩
  | cB gs =>
    exact gGates_sound p.d p.sw n _ _ p.ux p.us p.uy c.v (2 * c.v) (okBZ c) (fun _ => true) hp.hsw rfl
      rfl hu gs tr k hI h

theorem gItemsZ_sound (items : List Item) : ∀ (tr : Trie) (k : Trie → Bool), GOk p.d p.sw n tr →
    gItemsZ c p (2 * c.v) (2 * c.v + c.R) (maskW p.sw n) (64 + p.sw * n) (2^(64 + p.sw * n))
      c.ret.length items tr k = true →
    ∃ tr', GOk p.d p.sw n tr' ∧ k tr' = true ∧
      ItemsRunZ c p n items (absG p.d (64 + p.sw * n) tr) (absG p.d (64 + p.sw * n) tr') := by
  induction items with
  | nil => intro tr k hI h; exact ⟨tr, hI, h, rfl⟩
  | cons it items ih =>
    intro tr k hI h
    obtain ⟨tr1, hI1, hk1, r1⟩ := gItemZ_sound c p lo n hp it tr _ hI h
    obtain ⟨tr2, hI2, hk2, r2⟩ := ih tr1 k hI1 hk1
    exact ⟨tr2, hI2, hk2, _, r1, r2⟩
end

theorem SegOKZ.of_check {c : Raw} {p : SPar} {lo n : Nat} {pre post : Option Trie}
    {items : List Item} (h : segCheckZ c p lo n pre items post = true) :
    SegOKZ c p lo n pre items post := by
  unfold segCheckZ at h
  obtain ⟨hp, h1⟩ := pgood_of (gPre_sound c p lo n _ h)
  rw [forceN_eq] at h1
  refine ⟨hp, ?_⟩
  have hpost : ∀ tr3, GOk p.d p.sw n tr3 →
      Option.rec (motive := fun _ => Bool)
        (gFin p.sw (64 + p.sw * n) (2^(64 + p.sw * n)) p.uy c.v lo (lo + n) tr3 c.v 0)
        (fun t => beqK tr3 t) post = true →
      PostStZ c p lo n post (absG p.d (64 + p.sw * n) tr3) := by
    intro tr3 hI3 hk
    cases post with
    | none =>
      exact ⟨hI3.small, fin_of c p lo n hp c.v p.uy tr3 hI3 hk⟩
    | some t =>
      have e : tr3 = t := beqK_eq tr3 t hk
      rw [← e]
      exact ⟨hI3, rfl⟩
  cases pre with
  | none =>
    have h2 : gInit p.sw p.ux lo n 0 (Trie.mk p.d) (fun tr0 =>
        gItemsZ c p (2 * c.v) (2 * c.v + c.R) (maskW p.sw n) (64 + p.sw * n) (2^(64 + p.sw * n))
          c.ret.length items tr0 fun tr3 => Option.rec (motive := fun _ => Bool)
            (gFin p.sw (64 + p.sw * n) (2^(64 + p.sw * n)) p.uy c.v lo (lo + n) tr3 c.v 0)
            (fun t => beqK tr3 t) post) = true := h1
    have hle : lo + 0 + n ≤ 2^p.d := by have := hp.hd; have := hp.hn; omega
    obtain ⟨t0, hwf, hk, hval⟩ := gInit_sound p.d p.sw p.ux lo n 0 (Trie.mk p.d) _ (Trie.wf_mk p.d)
      hle h2
    obtain ⟨hI0, hinit⟩ := init_ok p.d p.sw n p.ux lo hp.hsw hp.hxs hp.hx2 t0 hwf hval
    obtain ⟨tr3, hI3, hk3, hr⟩ := gItemsZ_sound c p lo n hp items t0 _ hI0 hk
    exact ⟨_, _, hinit, hr, hpost tr3 hI3 hk3⟩
  | some t =>
    intro hI0
    have h2 : gItemsZ c p (2 * c.v) (2 * c.v + c.R) (maskW p.sw n) (64 + p.sw * n)
        (2^(64 + p.sw * n)) c.ret.length items t (fun tr3 => Option.rec (motive := fun _ => Bool)
          (gFin p.sw (64 + p.sw * n) (2^(64 + p.sw * n)) p.uy c.v lo (lo + n) tr3 c.v 0)
          (fun t => beqK tr3 t) post) = true := h1
    obtain ⟨tr3, hI3, hk3, hr⟩ := gItemsZ_sound c p lo n hp items t _ hI0 h2
    exact ⟨_, hr, hpost tr3 hI3 hk3⟩

theorem XRunZ.of_check {c : Raw} {p : SPar} {lo n : Nat} (h : scalCheckZ c p lo n = true) :
    XRunZ c lo n := XRunZ.of_seg (SegOKZ.of_check h)

theorem YRunZ.of_check {c : Raw} {p : SPar} {lo n : Nat} (h : yCheckZ c p lo n = true) :
    YRunZ c lo n := by
  unfold yCheckZ at h
  obtain ⟨hp, h1⟩ := pgood_of (gPre_sound c p lo n _ h)
  have hle : c.v + lo + 0 + n ≤ 2^p.d := by have := hp.hd; have := hp.hn; omega
  obtain ⟨t0, hwf, hk, hval⟩ := gInit_sound p.d p.sw p.uy (c.v + lo) n 0 (Trie.mk p.d) _
    (Trie.wf_mk p.d) hle h1
  obtain ⟨hI0, hinit⟩ := init_ok p.d p.sw n p.uy (c.v + lo) hp.hsw hp.hys hp.hy2 t0 hwf hval
  obtain ⟨tr3, hI3, hk3, hr⟩ := gChunks_sound p.d p.sw n _ _ p.ux p.us p.uy c.v (2 * c.v) (okBZ c)
    (selY c.v (2 * c.v)) hp.hsw rfl rfl (U_pos c p lo n hp) c.B t0 _ hI0 hk
  exact ⟨p, _, _, hp, hinit, hr, hI3.small, fin_of c p lo n hp c.v p.uy tr3 hI3 hk3⟩

/-! ### the y rows of the source columns (`hid`) -/
variable {c : Raw}

theorem XRunZ.addA {lo n : Nat} (h : XRunZ c lo n) : ∀ x ∈ addsA c,
    x.1 < 2 * c.v + c.R ∧ x.2.1 < 2 * c.v + c.R ∧ x.1 ≠ x.2.1 ∧ (x.1 < c.v ∨ 2 * c.v ≤ x.1) ∧
      (x.2.1 < c.v ∨ 2 * c.v ≤ x.2.1) := by
  obtain ⟨p, a0, a1, a2, a3, hp, hi, rA, rS, hsh, rB, hsm, fy⟩ := h
  intro x hx
  obtain ⟨ho, hne⟩ := GRun_all _ _ _ _ _ _ _ _ rA x hx rfl
  obtain ⟨h1, h2, h3, h4⟩ := okAZ_spec ho
  exact ⟨h1, h2, hne, h3, h4⟩

theorem XRunZ.addB {lo n : Nat} (h : XRunZ c lo n) : ∀ x ∈ addsB c,
    x.1 < 2 * c.v + c.R ∧ x.2.1 < 2 * c.v + c.R ∧ x.1 ≠ x.2.1 ∧
      (c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v) := by
  obtain ⟨p, a0, a1, a2, a3, hp, hi, rA, rS, hsh, rB, hsm, fy⟩ := h
  intro x hx
  obtain ⟨ho, hne⟩ := GRun_all _ _ _ _ _ _ _ _ rB x hx rfl
  obtain ⟨h1, h2, h3⟩ := okBZ_spec ho
  exact ⟨h1, h2, hne, by omega⟩

theorem XRunZ.ycol (hv : 0 < c.v) {lo n : Nat} (h : XRunZ c lo n) (t : Fin c.v) (hlo : lo ≤ t.val)
    (hhi : t.val < lo + n) :
    ∀ S : Fin c.v, Mt c hv (lY S) (lX t) = if S = t then 1 else 0 := by
  obtain ⟨p, a0, a1, a2, a3, hp, hi, rA, rS, hsh, rB, hsm, fy⟩ := h
  have hokA : ∀ t s, okAZ c t s = true → emb c (unemb c hv t) = t ∧ emb c (unemb c hv s) = s :=
    fun t s ho => ⟨emb_unemb hv t (okAZ_spec ho).1, emb_unemb hv s (okAZ_spec ho).2.1⟩
  have hokB : ∀ t s, okBZ c t s = true → emb c (unemb c hv t) = t ∧ emb c (unemb c hv s) = s :=
    fun t s ho => ⟨emb_unemb hv t (okBZ_spec ho).1, emb_unemb hv s (okBZ_spec ho).2.1⟩
  have hCA : Cm p.sw n (U c p) (emb c) a1 = MA c hv * Cm p.sw n (U c p) (emb c) a0 := by
    have := GRun_Cm p.sw n (U c p) (okAZ c) (fun _ => true) (emb c) (unemb c hv) (unemb_emb hv) hokA
      (addsA c) a0 a1 rA
    rw [micF_true] at this
    exact this
  have hCB : Cm p.sw n (U c p) (emb c) a3 = MB c hv * Cm p.sw n (U c p) (emb c) a2 := by
    have := GRun_Cm p.sw n (U c p) (okBZ c) (fun _ => true) (emb c) (unemb c hv) (unemb_emb hv) hokB
      (addsB c) a2 a3 rB
    rw [micF_true] at this
    exact this
  have hCS := cm_scat p lo n hp a1 a2 rS hsh
  have hus : p.ux < 2^p.sw :=
    lt_of_lt_of_le hp.hxs (Nat.pow_le_pow_right (by norm_num) (by omega))
  have hT : t.val - lo < n := by omega
  have hC0 : ∀ j : Rl c, Cm p.sw n (U c p) (emb c) a0 j ⟨t.val - lo, hT⟩
      = if j = lX t then 1 else 0 := by
    intro j
    have h1 := cont_start p.sw n (U c p) a0 lo p.ux hp.hx hus hi (emb c j) (t.val - lo) hT
      (fun hr => U_x c p _ (by have := hp.hn; omega))
    show gcont p.sw (U c p) a0 (emb c j) (t.val - lo) = _
    rw [h1]
    have e : lo + (t.val - lo) = t.val := by omega
    by_cases hj : j = lX t
    · rw [if_pos hj, hj, if_pos (by rw [e]; rfl)]
    · have hne : ¬ emb c j = lo + (t.val - lo) := by
        intro he
        apply hj
        rw [← unemb_emb hv j, he, e, unemb_x hv t.val t.isLt]
      rw [if_neg hj, if_neg hne]
  have key : ∀ i : Rl c, Mt c hv i (lX t) = Cm p.sw n (U c p) (emb c) a3 i ⟨t.val - lo, hT⟩ := by
    intro i
    have e1 := mul_unitCol (Mt c hv) (Cm p.sw n (U c p) (emb c) a0) (lX t) ⟨t.val - lo, hT⟩ hC0 i
    rw [← e1]
    unfold Mt
    rw [Matrix.mul_assoc, Matrix.mul_assoc, ← hCA, ← hCS, ← hCB]
  intro S
  rw [key]
  have hf := fy S.val S.isLt
  have h1 := cont_end p.sw n (U c p) a3 (c.v + S.val) p.uy (S.val - lo)
    (lo ≤ S.val ∧ S.val < lo + n) hp.hsw hp.hy hp.hys
    (U_y c p _ (by omega) (by have := S.isLt; omega)) (hsm (c.v + S.val)).2 hf (t.val - lo) hT
  show gcont p.sw (U c p) a3 (c.v + S.val) (t.val - lo) = _
  rw [h1]
  by_cases e : S = t
  · rw [if_pos e, if_pos ⟨by rw [e]; exact ⟨hlo, hhi⟩, by rw [e]⟩]
  · have hne : ¬ ((lo ≤ S.val ∧ S.val < lo + n) ∧ t.val - lo = S.val - lo) := by
      intro hh
      exact e (Fin.ext (by omega))
    rw [if_neg e, if_neg hne]

/-! ### THE SCALAR HALF, SOURCES MUTABLE -/

structure ScalOKZ (c : Raw) : Prop where
  hv : 0 < c.v
  hs : ∀ t, t < c.v → ∃ lo n, lo ≤ t ∧ t < lo + n ∧ XRunZ c lo n
  hyc : ∀ t, t < c.v → ∃ lo n, lo ≤ t ∧ t < lo + n ∧ YRunZ c lo n

namespace ScalOKZ
variable (K : ScalOKZ c)
include K

theorem addA : ∀ x ∈ addsA c, x.1 < 2 * c.v + c.R ∧ x.2.1 < 2 * c.v + c.R ∧ x.1 ≠ x.2.1 ∧
    (x.1 < c.v ∨ 2 * c.v ≤ x.1) ∧ (x.2.1 < c.v ∨ 2 * c.v ≤ x.2.1) := by
  obtain ⟨lo, n, _, _, h⟩ := K.hs 0 K.hv
  exact h.addA

theorem addB : ∀ x ∈ addsB c, x.1 < 2 * c.v + c.R ∧ x.2.1 < 2 * c.v + c.R ∧ x.1 ≠ x.2.1 ∧
    (c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v) := by
  obtain ⟨lo, n, _, _, h⟩ := K.hs 0 K.hv
  exact h.addB

theorem un_ne (t s : Nat) (ht : t < 2 * c.v + c.R) (hs : s < 2 * c.v + c.R) (hne : t ≠ s) :
    unemb c K.hv t ≠ unemb c K.hv s := fun e => hne (by
  rw [← emb_unemb K.hv t ht, e, emb_unemb K.hv s hs])

theorem hAi : MAi c K.hv * MA c K.hv = 1 := by
  apply mul_eq_one_comm.mp
  exact matP_inv _ _ (gatesDistinct_of _ _ (fun t s cf hm => by
    obtain ⟨x, hx, rfl, rfl⟩ := mem_mic _ _ _ _ hm
    obtain ⟨h1, h2, h3, _⟩ := K.addA x hx
    exact K.un_ne _ _ h1 h2 h3))

theorem hBi : MBi c K.hv * MB c K.hv = 1 := by
  apply mul_eq_one_comm.mp
  exact matP_inv _ _ (gatesDistinct_of _ _ (fun t s cf hm => by
    obtain ⟨x, hx, rfl, rfl⟩ := mem_mic _ _ _ _ hm
    obtain ⟨h1, h2, h3, _⟩ := K.addB x hx
    exact K.un_ne _ _ h1 h2 h3))

/-- **the block `y ← x` is the identity** -/
theorem hid (S t : Fin c.v) :
    (MB c K.hv * scatM (Jrm c * Ccm c) * MA c K.hv) (lY S) (lX t) = if S = t then 1 else 0 := by
  obtain ⟨lo, n, hlo, hhi, h⟩ := K.hs t.val t.isLt
  exact h.ycol K.hv t hlo hhi S

theorem mb_ycol (i : Rl c) (t : Fin c.v) : MB c K.hv i (lY t) = if i = lY t then 1 else 0 := by
  have hf := ycol_filter (unemb c K.hv) (selY c.v (2 * c.v)) (addsB c)
    (fun x hx hsel τ hτ => by
      obtain ⟨_, h2, _⟩ := K.addB x hx
      have := selY_spec.mpr (unemb_of_y K.hv x.2.1 h2 τ hτ)
      rw [hsel] at this
      exact Bool.false_ne_true this)
    (fun x hx hsel => by
      obtain ⟨_, _, _, h5⟩ := K.addB x hx
      obtain ⟨g1, g2⟩ := selY_spec.mp hsel
      obtain ⟨g3, g4⟩ := h5 g1 g2
      exact ⟨_, unemb_y K.hv x.1 g3 g4⟩) i t
  show matP (unemb c K.hv) (mic (addsB c)) i (lY t) = _
  rw [hf]
  have htg : ∀ t' s cf, Micro.add t' s cf ∈ micF (selY c.v (2 * c.v)) (addsB c) →
      ∃ τ : Fin c.v, unemb c K.hv t' = lY τ := by
    intro t' s cf hm
    obtain ⟨x, hx, rfl, rfl, hsel⟩ := mem_micF _ _ _ _ _ hm
    obtain ⟨_, _, _, h5⟩ := K.addB x hx
    obtain ⟨g1, g2⟩ := selY_spec.mp hsel
    obtain ⟨g3, g4⟩ := h5 g1 g2
    exact ⟨_, unemb_y K.hv x.1 g3 g4⟩
  rcases i with (τ | S) | q
  · exact matP_row _ _ _ (fun t' s cf hm hh => by
      obtain ⟨τ', hτ'⟩ := htg t' s cf hm
      rw [hτ'] at hh; cases hh) _
  · obtain ⟨lo, n, hlo, hhi, p, a0, a3, hp, hi, hr, hsm, fy⟩ := K.hyc t.val t.isLt
    have hokB : ∀ t s, okBZ c t s = true →
        emb c (unemb c K.hv t) = t ∧ emb c (unemb c K.hv s) = s :=
      fun t s ho => ⟨emb_unemb K.hv t (okBZ_spec ho).1, emb_unemb K.hv s (okBZ_spec ho).2.1⟩
    have hC := GRun_Cm p.sw n (U c p) (okBZ c) (selY c.v (2 * c.v)) (emb c) (unemb c K.hv)
      (unemb_emb K.hv) hokB (addsB c) a0 a3 hr
    have hus : p.uy < 2^p.sw :=
      lt_of_lt_of_le hp.hys (Nat.pow_le_pow_right (by norm_num) (by omega))
    have hT : t.val - lo < n := by omega
    have hC0 : ∀ j : Rl c, Cm p.sw n (U c p) (emb c) a0 j ⟨t.val - lo, hT⟩
        = if j = lY t then 1 else 0 := by
      intro j
      have h1 := cont_start p.sw n (U c p) a0 (c.v + lo) p.uy hp.hy hus hi (emb c j) (t.val - lo) hT
        (fun hr => U_y c p _ (by omega) (by have := hp.hn; omega))
      show gcont p.sw (U c p) a0 (emb c j) (t.val - lo) = _
      rw [h1]
      have e : c.v + lo + (t.val - lo) = c.v + t.val := by omega
      by_cases hj : j = lY t
      · rw [if_pos hj, hj, if_pos (by rw [e]; rfl)]
      · have hne : ¬ emb c j = c.v + lo + (t.val - lo) := by
          intro he
          apply hj
          rw [← unemb_emb K.hv j, he, e, unemb_y K.hv (c.v + t.val) (by omega)
            (by have := t.isLt; omega)]
          exact congrArg (fun x => (Sum.inl (Sum.inr x) : Rl c)) (Fin.ext (by simp))
        rw [if_neg hj, if_neg hne]
    have e1 := mul_unitCol (matP (unemb c K.hv) (micF (selY c.v (2 * c.v)) (addsB c)))
      (Cm p.sw n (U c p) (emb c) a0) (lY t) ⟨t.val - lo, hT⟩ hC0 (lY S)
    rw [← e1, ← hC]
    have h1 := cont_end p.sw n (U c p) a3 (c.v + S.val) p.uy (S.val - lo)
      (lo ≤ S.val ∧ S.val < lo + n) hp.hsw hp.hy hp.hys
      (U_y c p _ (by omega) (by have := S.isLt; omega)) (hsm (c.v + S.val)).2 (fy S.val S.isLt)
      (t.val - lo) hT
    show gcont p.sw (U c p) a3 (c.v + S.val) (t.val - lo) = _
    rw [h1]
    by_cases e : S = t
    · rw [e, if_pos ⟨⟨hlo, hhi⟩, rfl⟩, if_pos rfl]
    · have hne : ¬ ((lo ≤ S.val ∧ S.val < lo + n) ∧ t.val - lo = S.val - lo) := by
        intro hh
        exact e (Fin.ext (by omega))
      rw [if_neg hne, if_neg (fun hh => e (Sum.inr.inj (Sum.inl.inj hh)))]
  · exact matP_row _ _ _ (fun t' s cf hm hh => by
      obtain ⟨τ', hτ'⟩ := htg t' s cf hm
      rw [hτ'] at hh; cases hh) _

/-- **the y columns are unit columns** -/
theorem hy (i : Rl c) (t : Fin c.v) :
    (MB c K.hv * scatM (Jrm c * Ccm c) * MA c K.hv) i (lY t) = if i = lY t then 1 else 0 := by
  have hA : ∀ i : Rl c, MA c K.hv i (lY t) = if i = lY t then 1 else 0 :=
    matP_col (unemb c K.hv) (lY t) (mic (addsA c)) (fun t' s cf hm hh => by
    obtain ⟨x, hx, rfl, rfl⟩ := mem_mic _ _ _ _ hm
    obtain ⟨_, h2, _, _, h5⟩ := K.addA x hx
    have := unemb_of_y K.hv x.2.1 h2 t hh
    omega)
  have hS : ∀ i : Rl c, scatM (Jrm c * Ccm c) i (lY t) = if i = lY t then 1 else 0 := by
    intro i
    unfold scatM
    rw [Matrix.add_apply, Matrix.one_apply, embYS_colY, add_zero]
  have hcol := unitCol_mul _ _ (lY t) hS hA
  rw [Matrix.mul_assoc, mul_unitCol (MB c K.hv) _ (lY t) (lY t) hcol i]
  exact K.mb_ycol i t

end ScalOKZ

/-- **a checked mutable-source certificate gives the mutable-source scalar half** -/
noncomputable def ScalOKZ.toScalZ (K : ScalOKZ c) :
    OAI.PowerSaving.GX.ScalZ (Fin c.v) (Fin c.R) (Fin c.ret.length) :=
  ⟨Ccm c, Jrm c, MA c K.hv, MAi c K.hv, MB c K.hv, MBi c K.hv, K.hAi, K.hBi, K.hy, K.hid⟩

end GS

#print axioms GS.XRunZ.of_check
#print axioms GS.YRunZ.of_check
#print axioms GS.ScalOKZ.hy
#print axioms GS.ScalOKZ.hid
#print axioms GS.ScalOKZ.toScalZ
