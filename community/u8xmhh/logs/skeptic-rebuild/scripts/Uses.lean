import Lean
open Lean

/-! v-fs-rebuild (skeptic): which modules does the PROOF of the final theorems really reach?
Transitive closure through types, definition bodies, theorem PROOFS, inductives, constructors, recursors
(the same walk as tools/Compare.lean with withProofs := true), grouped by defining module.
usage: lake env lean --run Uses.lean <module> <theorem>... -/

namespace VU
partial def implicitGo (e : Expr) : StateM (Std.HashSet Expr × Array Name) Unit := do
  if (← get).1.contains e then return
  modify fun (seen, out) => (seen.insert e, out)
  match e with
  | .proj s _ b => modify (fun (seen, out) => (seen, out.push s)); implicitGo b
  | .lit (.natVal _) => modify fun (seen, out) => (seen, out.push ``Nat)
  | .lit (.strVal _) => modify fun (seen, out) => (seen, out.push ``String)
  | .app f a => implicitGo f; implicitGo a
  | .lam _ t b _ => implicitGo t; implicitGo b
  | .forallE _ t b _ => implicitGo t; implicitGo b
  | .letE _ t v b _ => implicitGo t; implicitGo v; implicitGo b
  | .mdata _ b => implicitGo b
  | _ => pure ()

def exprNames (e : Expr) : Array Name := e.getUsedConstants ++ ((implicitGo e).run ({}, #[])).2.2

def deps (ci : ConstantInfo) : Array Name := Id.run do
  let mut out := exprNames ci.type
  match ci with
  | .defnInfo v => out := out ++ exprNames v.value
  | .opaqueInfo v => out := out ++ exprNames v.value
  | .thmInfo v => out := out ++ exprNames v.value
  | .inductInfo v => out := out ++ v.ctors.toArray ++ v.all.toArray
  | .ctorInfo v => out := out.push v.induct
  | .recInfo v =>
    out := out ++ v.all.toArray
    for r in v.rules do
      out := (out ++ exprNames r.rhs).push r.ctor
  | _ => pure ()
  return out

def closure (env : Environment) (roots : Array Name) : NameSet × Array Name := Id.run do
  let mut seen : NameSet := {}
  let mut missing : Array Name := #[]
  let mut stack := roots
  while !stack.isEmpty do
    let n := stack.back!
    stack := stack.pop
    if seen.contains n then continue
    seen := seen.insert n
    match env.find? n with
    | none => missing := missing.push n
    | some ci => stack := stack ++ deps ci
  return (seen, missing)

def moduleOf (env : Environment) (n : Name) : String :=
  match env.getModuleIdxFor? n with
  | some i => toString (env.header.moduleNames[i.toNat]!)
  | none => "<current>"

/-- `Work.XHost.Data.HU8XMSched` ↦ `Work.XHost` ; `Mathlib.X.Y` ↦ `Mathlib` ; `OAI.Computability.X` ↦ `OAI.Computability` -/
def group (m : String) : String :=
  let parts := m.splitOn "."
  match parts with
  | "Work" :: d :: _ => "Work." ++ d
  | "OAI" :: d :: _ => "OAI." ++ d
  | p :: _ => p
  | [] => m
end VU

open VU in
unsafe def main (args : List String) : IO UInt32 := do
  let modl :: thms := args | IO.eprintln "usage"; return 2
  initSearchPath (← findSysroot)
  let env ← importModules #[{ module := modl.toName }] {} (trustLevel := 0)
  IO.println s!"module {modl}: {env.header.moduleNames.size} modules loaded"
  let mut rc : UInt32 := 0
  for t in thms do
    let T := t.toName
    let (cl, miss) := closure env #[T]
    IO.println s!"\n=== root {T}: proof closure {cl.size} constants, unresolved {miss.size}"
    -- by group
    let mut cnt : Std.HashMap String Nat := {}
    let mut byMod : Std.HashMap String (Array Name) := {}
    for n in cl.toList do
      let m := moduleOf env n
      cnt := cnt.insert (group m) (cnt.getD (group m) 0 + 1)
      byMod := byMod.insert m ((byMod.getD m #[]).push n)
    let groups := (cnt.toList.map (fun (g, c) => (g, c))).toArray.qsort (fun a b => a.1 < b.1)
    for (g, c) in groups do
      if g.startsWith "Work" || g.startsWith "OAI" || g.startsWith "WHT" || g.startsWith "XExit" then
        IO.println s!"  group {g}: {c}"
    let others := groups.filter (fun (g, _) => !(g.startsWith "Work" || g.startsWith "OAI" || g.startsWith "WHT" || g.startsWith "XExit"))
    IO.println s!"  other groups: {others.toList}"
    -- modules of interest: every module of XHost, KMX, XMut, and the U8XM/U8XMHH data
    let mods := (byMod.toList.map (·.1)).toArray.qsort (· < ·)
    for m in mods do
      if m.startsWith "Work.XHost" || m.startsWith "Work.KMX" || m.startsWith "Work.XMut" || m.startsWith "Work.FourierU8XMHH"
         || m.startsWith "Work.FourierU8XM." || m.startsWith "Work.FourierC588" || m.startsWith "Work.GCert.Data.U8XM"
         || m.startsWith "Work.GCert.Data.InstU8XM" || m.startsWith "Work.GCert.Data.Gen.U8XM" || m.startsWith "Work.KLow" then
        let ns := (byMod.getD m #[]).qsort (fun a b => a.toString < b.toString)
        IO.println s!"  module {m}: {ns.size} constants; e.g. {(ns.toList.take 6)}"
    -- named checkpoints
    let want : List Name := [`OAI.PowerSaving.XHost.engine_of_ginv_host, `OAI.PowerSaving.XHost.HostSched.ofTab_fitsR,
      `OAI.PowerSaving.XHost.hostSchedK, `GXD.U8XMHH.hs_ok, `GXD.U8XMHH.hsched, `GXD.U8XMHH.hfits, `GXD.U8XMHH.hprice,
      `GXD.U8XMHH.hwidth, `GXD.U8XM.ginv, `GXD.ginvXOfZ, `GXD.U8XM.scalOKZ, `GXD.U8XM.validX, `GXD.exit_of_zchecksZ,
      `OAI.PowerSaving.BlockAccounting.fold_U8XMHHx, `OAI.PowerSaving.BlockAccounting.blocksU8XMHHx,
      `OAI.PowerSaving.BridgeRate.U8XMHH.unit_cost_x, `OAI.PowerSaving.SeamU8XMHH.engine_U8XMHHx,
      `OAI.PowerSaving.RAM.hillsU8XMHH_program, `OAI.PowerSaving.BlockAccounting.fold_U8XMx]
    for w in want do
      IO.println s!"  contains {w}: {cl.contains w}  [{if (env.find? w).isSome then moduleOf env w else "not declared"}]"
    let axs := cl.toList.filter (fun n => match env.find? n with | some (.axiomInfo _) => true | _ => false)
    IO.println s!"  axioms reached: {axs}"
    if !(cl.contains `OAI.PowerSaving.XHost.engine_of_ginv_host) then rc := 1
  return rc
