import Work.FourierC588M.Main
open Lean

/-! v-m-rebuild (13 Oct 2026): is the XMut engine actually on the proof path of the final theorems?
Walk the constant closure (types and values, inductives with constructors/recursors) of each root and report
which of the named engine/certificate constants it reaches, and how many constants of each overlay module family. -/

partial def closureOf (env : Environment) (roots : List Name) : NameSet := Id.run do
  let mut seen : NameSet := {}
  let mut stack := roots
  while true do
    match stack with
    | [] => break
    | n :: rest =>
      stack := rest
      if seen.contains n then continue
      seen := seen.insert n
      match env.find? n with
      | none => pure ()
      | some ci =>
        let mut used := ci.type.getUsedConstants.toList
        match ci.value? (allowOpaque := true) with
        | some v => used := used ++ v.getUsedConstants.toList
        | none => pure ()
        match ci with
        | .inductInfo i => used := used ++ i.ctors ++ i.all
        | .ctorInfo c => used := used ++ [c.induct]
        | .recInfo r => used := used ++ r.all ++ (r.rules.map (·.rhs.getUsedConstants.toList)).flatten
        | _ => pure ()
        for u in used do
          if !seen.contains u then stack := u :: stack
  return seen

def modOf (env : Environment) (n : Name) : String :=
  match env.getModuleIdxFor? n with
  | some i => toString (env.header.moduleNames[i.toNat]!)
  | none => "<current>"

#eval show CoreM Unit from do
  let env ← getEnv
  -- suffixes of names: engine (Shea/Mamano), x-mutable's mutable-source engine, the C588M certificate checks, the seam;
  -- and Shea's ORIGINAL copy-carrying C588 objects (must NOT be on the path of the C588M theorems)
  let probes : List String := [
    "engine_of_ginv_401", "GXD.ginvOfZ", "GXD.ginvOf", "GXD.gscalOfZ", "GXD.glabOfZ", "GX.GLabZ.inv", "GLab.ValidE.glabZ",
    "GS.ScalOKZ", "GS.scalCheckZ", "GS.yCheckZ", "GS.eCheckZ", "GS.eCheckZ_sound", "GS.XRunZ.of_check", "GS.YRunZ.of_check",
    "GXD.C588M.ginv", "GXD.C588M.scalOKZ", "GXD.C588M.validE", "GXD.C588M.es_ok", "GXD.C588M.price", "GXD.C588M.sched",
    "GXD.C588M.scal0_ok", "GXD.C588M.scal3_ok", "GXD.C588M.ychk_ok", "GXD.C588M.echk0_ok", "GXD.C588M.echk2_ok", "GXD.C588M.sched_ok",
    "fold_C588Mx", "C588M.unit_cost_x", "SeamC588M.engine_C588Mx", "RAM.hillsC588M_program",
    "GXD.C588.ginv", "GXD.C588.scalOK", "RAM.hillsC588_program", "RAM.hills_program", "transform_mainC588", "fold_C588x"]
  for root in [`OAI.PowerSaving.transform_mainC588M, `OAI.PowerSaving.convolution_mainC588M] do
    let cl := closureOf env [root]
    IO.println s!"ROOT {root}: closure size {cl.size}"
    let names := cl.toList.map toString
    for p in probes do
      let hits := names.filter (fun n => n == p || n.endsWith ("." ++ p))
      let inEnv := env.constants.toList.any (fun (n, _) => let t := toString n; t == p || t.endsWith ("." ++ p))
      IO.println s!"  {if hits.isEmpty then "absent " else "REACHES"} {p} {hits.take 3}{if inEnv then "" else "  (no such constant in env)"}"
    let mut counts : Std.HashMap String Nat := {}
    for n in cl.toList do
      let m := modOf env n
      if (m.startsWith "Work.XMut") || (m.splitOn "C588M").length > 1 || (m.splitOn "C588").length > 1 then
        counts := counts.insert m (counts.getD m 0 + 1)
    for (m, c) in counts.toList.toArray.qsort (fun a b => a.1 < b.1) do
      IO.println s!"  module {m}: {c} constants reached"
