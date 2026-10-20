import Work.FourierL2979H.Main
import Lean
import Std.Data.HashMap
import Std.Data.HashSet
/-!
BosDump (h-bos, round wht45/release, 19 Oct 2026): fingerprint every constant of a submission's closure, so that two
packagings of the same proof (the packager's stage and the merged stage of `bos_merge.py`) can be compared constant
by constant.  Run inside a BUILT stage:  `BOSDUMP_OUT=<file> lake env lean BosDump.lean`.
Tree mode (the packager's source tree, already built): replace `import Submission` by `import <root module>`, set
BOSDUMP_MODULES to a file listing the closure's modules (one per line), and run lean with LEAN_PATH (hbos-check-fast.sh).

One line per constant whose module is a `Submission` module (the closure and the adapter):
  kind <TAB> normalised name <TAB> universe parameters <TAB> type hash <TAB> value hash <TAB> flag <TAB> module <TAB> name
Hashes are 64-bit structural hashes (binder names and metadata ignored, binder info kept) of the expressions after
normalising the names that embed the module they were made in, or that Lean may share between declarations of one
module (its caches): a private name `_private.<Module>.<k>.x` becomes `x.‹private›`; an auxiliary constant (a name
component starting with `_auxLemma`, `_aux_`, `_proof_`, `proof_`, `_simp_` or `match_`, or a name with macro scopes)
is named by the hash of its kind, normalised type and, unless it is a theorem or a metaprogramming object,
normalised value.  The value hash is that of the definition body (definitions, opaques), of the
proof (theorems: reported, not required to agree), of the constructor list (inductives), of the rules (recursors).
Flag column: 1 = auxiliary, m = metaprogramming object (type head in `Lean`/IO: notation parsers, macros), else 0.
Universe parameter names lose their macro scopes.  Last lines: `AXIOMS <name> <axioms...>` for each name in BOSDUMP_AXIOMS (default `Race.bound`).
-/
open Lean Elab Command

namespace BosDump

@[inline] def mix (h x : UInt64) : UInt64 :=
  let z : UInt64 := h ^^^ (x + (0x9e3779b97f4a7c15 : UInt64) + (h <<< (6 : UInt64)) + (h >>> (2 : UInt64)))
  let z : UInt64 := (z ^^^ (z >>> (30 : UInt64))) * (0xbf58476d1ce4e5b9 : UInt64)
  let z : UInt64 := (z ^^^ (z >>> (27 : UInt64))) * (0x94d049bb133111eb : UInt64)
  z ^^^ (z >>> (31 : UInt64))

def hashStr (s : String) : UInt64 := s.foldl (fun h c => mix h c.toNat.toUInt64) 0x51ed27

partial def hashNat (n : Nat) (h : UInt64 := 0x7f4a) : UInt64 :=
  if n < 18446744073709551616 then mix h n.toUInt64
  else hashNat (n >>> 64) (mix h (n % 18446744073709551616).toUInt64)

def auxPrefixes : List String := ["_auxLemma", "_aux_", "_proof_", "proof_", "_simp_", "match_"]

def isAuxName (n : Name) : Bool :=
  n.hasMacroScopes || n.components.any fun c =>
    match c with
    | .str _ s => auxPrefixes.any (s.startsWith ·)
    | _ => false

/-- private names lose their module part; everything else unchanged (aux names are handled by the caller) -/
def normPrivate (n : Name) : Name :=
  if isPrivateName n then
    match privateToUserName? n with
    | some u => Name.str u "‹private›"
    | none => n
  else n

abbrev Memo := Std.HashMap USize UInt64

/-- universe levels; parameter names lose their macro scopes (hygienic names embed the module and a counter) -/
def hLevel : Level → UInt64
  | .zero => 21
  | .succ l => mix 22 (hLevel l)
  | .max a b => mix 23 (mix (hLevel a) (hLevel b))
  | .imax a b => mix 24 (mix (hLevel a) (hLevel b))
  | .param n => mix 25 n.eraseMacroScopes.hash
  | .mvar m => mix 26 m.name.hash

def hLevels (us : List Level) : UInt64 := us.foldl (fun a u => mix a (hLevel u)) 27

def biCode : BinderInfo → UInt64
  | .default => 0 | .implicit => 1 | .strictImplicit => 2 | .instImplicit => 3

unsafe def hExprU (nm : Name → UInt64) (e : Expr) : StateM Memo UInt64 := do
  let key := ptrAddrUnsafe e
  if let some h := (← get)[key]? then return h
  let h ← match e with
    | .bvar i => pure (mix 1 i.toUInt64)
    | .fvar f => pure (mix 2 f.name.hash)
    | .mvar m => pure (mix 3 m.name.hash)
    | .sort l => pure (mix 4 (hLevel l))
    | .const n us => pure (mix 5 (mix (nm n) (hLevels us)))
    | .app f a => do
        let x ← hExprU nm f; let y ← hExprU nm a; pure (mix 6 (mix x y))
    | .lam _ t b bi => do
        let x ← hExprU nm t; let y ← hExprU nm b; pure (mix 7 (mix (mix x y) (biCode bi)))
    | .forallE _ t b bi => do
        let x ← hExprU nm t; let y ← hExprU nm b; pure (mix 8 (mix (mix x y) (biCode bi)))
    | .letE _ t v b nd => do
        let x ← hExprU nm t; let y ← hExprU nm v; let z ← hExprU nm b
        pure (mix 9 (mix (mix x y) (mix z (if nd then 1 else 0))))
    | .lit (.natVal n) => pure (mix 10 (hashNat n))
    | .lit (.strVal s) => pure (mix 11 (hashStr s))
    | .mdata _ b => hExprU nm b
    | .proj s i b => do
        let x ← hExprU nm b; pure (mix 12 (mix (nm s) (mix i.toUInt64 x)))
  modify (·.insert key h)
  return h

@[implemented_by hExprU] opaque hExpr (nm : Name → UInt64) (e : Expr) : StateM Memo UInt64

def runH (nm : Name → UInt64) (e : Expr) (memo : Memo) : UInt64 × Memo := Id.run ((hExpr nm e).run memo)

/-- metaprogramming objects (parser descriptors of notations, macros, elaborators, initialisers): their values embed
the module's name in name literals; they are compared by type only -/
def isMeta (ci : ConstantInfo) : Bool :=
  match ci with
  | .defnInfo _ | .opaqueInfo _ =>
    match ci.type.getForallBody.getAppFn.constName? with
    | some n => (`Lean).isPrefixOf n || n == ``IO || n == ``EIO || n == ``BaseIO
    | none => false
  | _ => false

def kindOf : ConstantInfo → String
  | .axiomInfo _ => "axiom" | .defnInfo _ => "def" | .thmInfo _ => "thm" | .opaqueInfo _ => "opaque"
  | .quotInfo _ => "quot" | .inductInfo _ => "induct" | .ctorInfo _ => "ctor" | .recInfo _ => "rec"

end BosDump

open BosDump in
#eval show CommandElabM Unit from do
  let env ← getEnv
  let out := (← IO.getEnv "BOSDUMP_OUT").getD "bosdump.tsv"
  let axNames := ((← IO.getEnv "BOSDUMP_AXIOMS").getD "Race.bound").splitOn " " |>.filter (· ≠ "")
  let modNames := env.header.moduleNames
  let modOf (n : Name) : Option Name := (env.getModuleIdxFor? n).map fun i => modNames[i.toNat]!
  -- the closure: the `Submission` modules, or (tree mode) the modules listed in the file BOSDUMP_MODULES
  let listed : Std.HashSet Name ← match (← IO.getEnv "BOSDUMP_MODULES") with
    | some f => do
        let ls ← IO.FS.lines f
        pure (ls.foldl (fun s l => if l.trim.isEmpty then s else s.insert l.trim.toName) {})
    | none => pure {}
  let inClosure (n : Name) : Bool :=
    match modOf n with
    | some m => if listed.isEmpty then (`Submission).isPrefixOf m else listed.contains m
    | none => false
  let mine := env.constants.fold (fun acc n ci => if inClosure n then acc.push (n, ci) else acc) #[]
  -- pass 1: each aux constant is named by its content (kind, type, value unless a theorem), aux-blind
  let nm0 : Name → UInt64 := fun n => if isAuxName n then hashStr "‹aux›" else (normPrivate n).hash
  let mut memo0 : Memo := {}
  let mut auxH : Std.HashMap Name UInt64 := {}
  for (n, ci) in mine do
    if isAuxName n then
      let (th, m1) := runH nm0 ci.type memo0
      let (vh, m2) := if isMeta ci then (0, m1) else match ci with
        | .defnInfo d => runH nm0 d.value m1
        | .opaqueInfo d => runH nm0 d.value m1
        | _ => (0, m1)
      memo0 := m2; auxH := auxH.insert n (mix (hashStr (kindOf ci)) (mix th vh))
  -- names outside the closure keep their own name (Mathlib, Lean, ...)
  let nm : Name → UInt64 := fun n =>
    if isAuxName n then (auxH.get? n).getD (mix (hashStr "‹aux-out›") (normPrivate n).hash)
    else (normPrivate n).hash
  let h ← IO.FS.Handle.mk out .write
  let mut memo : Memo := {}
  let mut count := 0
  for (n, ci) in mine do
    let (th, m1) := runH nm ci.type memo
    memo := m1
    let (vh, m2) : UInt64 × Memo := match ci with
      | .defnInfo d => runH nm d.value memo
      | .thmInfo d => runH nm d.value memo
      | .opaqueInfo d => runH nm d.value memo
      | .inductInfo d =>
          (d.ctors.foldl (fun a c => mix a (nm c))
            (mix (mix d.numParams.toUInt64 d.numIndices.toUInt64)
              (mix (if d.isRec then 1 else 0) (d.all.foldl (fun a c => mix a (nm c)) 0))), memo)
      | .ctorInfo d =>
          (mix (nm d.induct) (mix d.cidx.toUInt64 (mix d.numParams.toUInt64 d.numFields.toUInt64)), memo)
      | .recInfo d => Id.run do
          let mut mm := memo
          let mut a : UInt64 := mix d.numMotives.toUInt64 (mix d.numMinors.toUInt64 (if d.k then 1 else 0))
          for r in d.rules do
            let (rh, m') := runH nm r.rhs mm
            mm := m'; a := mix a (mix (nm r.ctor) (mix r.nfields.toUInt64 rh))
          return (a, mm)
      | .quotInfo _ => (5, memo)
      | .axiomInfo _ => (0, memo)
    memo := m2
    let aux := isAuxName n
    let shown := if aux then s!"‹aux›.{(auxH.get? n).getD 0}" else toString (normPrivate n)
    let safety := match ci with
      | .defnInfo d => match d.safety with | .safe => "" | .unsafe => "-unsafe" | .partial => "-partial"
      | .opaqueInfo d => if d.isUnsafe then "-unsafe" else ""
      | .axiomInfo d => if d.isUnsafe then "-unsafe" else ""
      | _ => ""
    let lps := ci.levelParams.map Name.eraseMacroScopes
    let flag := if aux then "1" else if isMeta ci then "m" else "0"
    h.putStrLn s!"{kindOf ci}{safety}\t{shown}\t{lps}\t{th}\t{vh}\t{flag}\t{(modOf n).getD .anonymous}\t{n}"
    count := count + 1
  for a in axNames do
    let n := a.toName
    if env.contains n then
      let axs ← Lean.collectAxioms n
      let axs := axs.qsort (fun x y => x.toString < y.toString)
      h.putStrLn s!"AXIOMS\t{a}\t{axs.toList}"
    else
      h.putStrLn s!"AXIOMS\t{a}\tMISSING"
  h.flush
  logInfo m!"BosDump: {count} constants of the closure written to {out}"
