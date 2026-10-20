# afc-gate: can A-fc (the full cut on #2979's words, fill u1-2979sep) be built in Lean, and how?

Round wht45/afc, 20 Oct 2026. Agent afc-gate. Not public text. No Lean on the Mac; light Python only (each run under
2 s and 0.25 GB); GitHub not touched; nothing posted. Labels: [Lean-checked Cn] (compiled on a Colab, log path),
[computed] (command in section 7), [read] (file named), [argued].

**v2, 11:16 EDT (v1 10:25) — VERDICT: GO. The 48-line cut fact on #2979's words is now [Lean-checked C4].**
It is fc-engine's `bph_cutF`, which is stated and proved for any class type `T`, replica count `nR`, frames `base` and
head set `hb`, on modules that are byte-identical in design A's closure (L2979H) and design B's (L3945F). My test
module, compiled on a copy of the built design-A tree TA (L2979H + CoreF), proves it at `T = Fin 112`, `nR = 74`,
composes it with x-inst's heads-first certificate of the #2979 warp network through `MX.cut_heads_full`
(`xnetGeomW_cut_full`: target rank 25 on every head), and feeds A-fc's numbers through CoreF's `full_bundle` together
with afc-cut's generated rate fact (`bundle_Fc`: `E.BridgeHyp R typ 25 ∧ R.RateQ zRate2979Fc 1`). Everything compiled
on the first try, standard axioms only (section 1d). The full-cut machinery (CoreF, `cut_heads_full`, `splitWord`,
`Rule.twoTarget 73 25`, `Eng.engine_seam_final`) is generic, and with 16,576 heads = every rank-48 block,
`splitWord 24 24` splits nothing (`hcutF_Fc`, by `decide`). A-fc's data go through design A's generators with pin
changes and two small code fixes (section 2). The rate generator proves N = 17887045 and refuses 17887046
(section 3). No STOP written. What remains for the round is data and assembly (afc-data, afc-cut, integration), with
no open mathematical step.

## 1. The 48-line cut fact on #2979's words (the `hz` of `cut_heads_full`)

### 1a. Where design B got it [read]

`fullcut/fc-engine/lean/Work/Mixed/Bridge/HFTest3945FW.lean`:
```lean
def cutSF (l : H) : BLiveT (Orth (BC.LC (Fin 8) H l)) T 0 nR → Finset (BC.LC (Fin 8) H l)
  | Sum.inl (_, d) => a4Lines l (dDecD d).2.1          -- ALL 48 A4 lines of the head's side s
  | Sum.inr (_, b) => b.elim0
lemma cutSF_card (l) (hH : card H = 9) (r) : (cutSF l r).card = 48            -- a4Lines_card
lemma bph_cutF (base : T → OBase H) (l : H) (hH) (hb : Finset (Fin nR × Bool × T)) (r) (hr : r ∈ headsL (nb := 0) hb) (x) :
    MX.bph (zGeom base l hH r) x = MX.bph ((cutA base l r).fam (cutSF l r)) x
  -- proof: mem_headsL, then MX.bph_cast (a4Lines_card hH l s).symm ((a4Base base l g t).fam (a4Lines l s)) x
```
It is "Theorem 1G's hz" of fc-gate GATE.md section 2 ("`hz` is `HeadM1.bph_cast`" with `S r := a4Lines l s`).
The reason it is so short: p-hfirst's head directions ARE the frame lines, cast to `Fin 48`:
`zGeom r = a4Dir base l hH g (s,t) = fun i => (a4Base base l g t).fam (a4Lines l s) (Fin.cast _ i)`
(`Work/XHostT/HF/GeomA4.lean` line 67, `GeomInst.lean` line 36, `zGeom_head`).

### 1b. Where design A's half cut got its 24-line version [read]

`allheads/aA-data/lean/Work/Mixed/Bridge/HFTest2979HW.lean` (`engine2979H_geomW`) passes p-hfirst's
`bph_cut base l hH hb` (`Work/XHostT/HF/GeomInst.lean` line 157: the FIRST 24 directions of each head,
`MX.partA (a := 24) (b := 24) (zGeom ..)`, equal the 24 lines `cutLines l hH s` of the frame, by `MX.bph_partA_fam`)
and `HFTest2979H.cutS_card` into `engine2979H_of_cut`, i.e. into `MX.cut_heads_half`. The heads-first certificate it
feeds is x-inst's warp version `XHostTW.xnetGeomW_heads_first` (`Work/XHostT/HFW2979/GeomInstW.lean`), with head word
`MX.headWord Sum.inl (headsL hb) (fun _ => 48) (zGeom base l hH) _` — the SAME head directions `zGeom` as design B.
x-inst's `a4t_climb_eq` (`HFW2979/GeomA4W.lean`) is the word-side fact: on the warp network's frames the A4 climb of
class `t`, side `s`, gauge `g` is ONE block of the 48 lines `a4Lines l s` of `a4Base base l g t`, then a free shift.

### 1c. What must be proved for 48 lines on #2979's template [read; argued]

`MX.cut_heads_full he L hL rk z ind A S hz body hk` (HeadM1 line 503, byte-identical in both closures) needs:
1. `hk : GLiveKernel e (headWord e L (fun _ => 48) z ind ++ body)` — the heads-first certificate. On #2979's warp
   network it is `xnetGeomW_heads_first I Ic hIc base l hbase hH hOk S hfit hb hunh h48` — a theorem that is generic
   in the schedule `S` and already Lean-checked inside L2979H (design A) [read: aA-final NOTE section 3]. For A-fc it
   is the same theorem at `S = sched2979Fc` (needs only afc-data's `sched2979Fc_fits`; `hb = headBins S`,
   `hunh = headBins_unh S`, `h48 = headBins_48 S` are generic). The chronology of the 9 spare roles plays no role
   here: the certificate only asks that head bins are never hosted and have remainder 48.
2. `hz : ∀ l ∈ L, ∀ x, bph (z l) x = bph ((A l).fam (S l)) x` with `z = zGeom base l hH`, `A = cutA base l`,
   `S = cutSF l` — exactly `bph_cutF base l hH (headBins S)` at `T = Fin 112`, `nR = 74`. Nothing in its statement or
   proof mentions #3945 (no replica count, no class count, no provider): it uses only `zGeom`, `cutA`, `headsL`,
   `mem_headsL`, `a4Lines`, `a4Lines_card`, `a4Base`, `dDecD` and `MX.bph_cast`.
3. `hS : (cutSF l r).card = 48` — `cutSF_card`; the head target rank is `card α - 48 = 73 - 48 = 25`
   (`card_LC401 l hH = 73`).

So the cut fact is "already proved" up to restating `bph_cutF` in a module that does NOT import design B's files
(fc-engine's `HFTest3945FW` imports `HFTest3945HW`, which pulls in #3945's data). [computed] Every file common to
design A's release overlay (286 files) and design B's full-cut overlay (412 files) is byte-identical: 122 common
paths, 0 differ, among them `HF/GeomInst`, `HF/GeomA4`, `HFW2979/{GeomInstW,GeomA4W,Heads,...}`, `TGeomW`,
`TProvW`, `Mixed/HeadM1`, `Mixed/Bridge/Core`. So `bph_cutF`'s dependencies have the same text in both closures.

### 1d. The Lean try [Lean-checked C4]

`lean/Work/Gate/CutFact2979Fc.lean` and `CutFact2979FcRate.lean` (mine; namespace `OAI.PowerSaving.PB.Gate2979Fc`;
never imported by the chain), run by `colab/gate-lean.sh ta` on C4 in my own `cp -a` of the built TA
(`/content/work/afc-gate/T`, kept, 12 GB: TA + afc-cut's `Rate2979Fc.lean` + my two modules, all built), 15:10-15:14 UTC
(11:10-11:14 EDT). Inputs hash-checked (`colab/in/SHA256SUMS`); in-tree CoreF 7aa5f6c2, HeadM1 caac2d8c,
GeomInstW 715b89ef (= both overlays). Logs `logs/C4/` (gate.out, build1/2, ax1/2).

| name | statement (short) | axioms |
|---|---|---|
| `bph_cutF` | for any `T`, `nR`, `base`, `hb`, head role `r ∈ headsL hb`: `bph (zGeom base l hH r) x = bph ((cutA base l r).fam (cutSF l r)) x`, `cutSF` = all 48 lines `a4Lines l s` | propext, Classical.choice, Quot.sound |
| `bph_cutF2979` | the same at `T = Fin 112`, `nR = 74` (A-fc's network) | standard |
| `cutTgt_rank25` | every head role's target rank is `card LC - 48 = 25` | standard |
| `xnetGeomW_cut_full` | #2979's warp network: `xnetGeomW_heads_first` + `MX.cut_heads_full`: `∃ body`, proper, cost `card Γ * unitPriceTW Ic S - |headsL hb| * φ 48`, `GLiveKernelM Sum.inl (blockMat (cutTgt (headsL hb) (cutA base l) (cutSF l) r).2) (body ++ shiftWord ..)` for ANY schedule `S` with `hfit`, any head set with `hunh`, `h48` | standard |
| `hcutF_Fc` | CoreF's identity: `histFn cutFc r = [r ≠ 48] histFn unitFc r + 2 [r = 24] (histFn unitFc 48 - 16576)`, by `hcutF_of_check` + `decide` | standard |
| `unitFc_48`, `hcutF_Fc_fewer`, `unitFc_ne_cutFc` | 16,576 rank-48 blocks; the identity FAILS with 16,575 heads (one block would split); unit ≠ cut at 48 | propext |
| `bundle_Fc_hyp` | CoreF `full_bundle` at W0 66,304, Hd 16,576, c 48, a = b = 24, rmax 23, kc 25, the rate fact as hypothesis | standard |
| `cutFc_eq`, `bundle_Fc` | my `cutFc` is afc-cut's `MX.cut2979Fc` (`rfl`); `bundle_Fc` = `bundle_Fc_hyp` with `rate2979Fc_moment`: `∃ E R typ, E.e = e ∧ E.mx = 73 ∧ E.BridgeHyp R typ 25 ∧ R.RateQ zRate2979Fc 1` | none / standard |

Builds: 110 s (CutFact2979Fc) and 37 s (CutFact2979FcRate, incl. afc-cut's Rate2979Fc), at most 6.9 GB; no `sorry`
warning; source grep: only doc-comment hits. One deprecation warning (`if_pos`), harmless. So the chain's
`HFTest2979FcW` (afc-cut) can restate `bph_cutF` verbatim, exactly as planned in afc-cut's NEEDS section 3.
(A first try in an own copy of ours/T at 10:49 was stopped by me at 81 built modules: it duplicated the TA build
and C3 then had 7 heavy Lean processes; `logs/C3/gate-try1.out`.)

## 2. The 9 spare roles through design A's generators (dry runs on the Mac) [computed]

Inputs: A-fc's FILL (7fc766c6) and PROFILE (5df79411); ORDER, DISPATCH, PROVIDERS = #2979's files as in aA-data's
`in/alloc2979H` (A-fc's ORDER.json and DISPATCH.json are content-equal to them, `==` on the parsed JSON, but not
byte-equal; the converter sha-pins #2979's bytes, so use #2979's files). Outputs in `dry/` (small).

| generator (aA-data / x-tables) | result on A-fc | what assumes design A |
|---|---|---|
| `own_to_host2979H.py` (owner table) | **ALL OK** with only the FILL and PROFILE pins changed (`tools/own_to_host2979Fc_dry.py`): fill rules (hub window, our window, family once per call, compact prefix, capacity), all 271,728 slots hosted, rows injective, remainders = PROFILE.json minus the 666 calls, profile = PROFILE.json, rank total, dispatch counts 416 / 250; 0.8 s, 0.25 GB | only the two pins (FILL 764395fc, PROFILE 94a81974). It reads each role's bin from the FILL row: **no planned placement anywhere**; the spare roles are ordinary rows and pass the same window / same-call / capacity checks |
| `heads2979.py` (head list) | heads 16,576 (8,288 + 8,288), rem-48 bins = heads exactly, one head per array, FILL recount equal, template facts true; **but `nonheads_count` FALSE** with expected `8288,8288,0,0` | a code bug, not a data fact (afc-data found and fixed it too: its `tools/heads2979Fc.py`, 10:17 EDT): `byn == {12: E12, 13: E13}` compares an EMPTY Counter with `{12: 0, 13: 0}`, which is False on the Mac's Python 3.9 (True on Python >= 3.10). Fix in afc-data's copy: compare with the zero entries dropped (`dict(byn) == {k: v for k, v in {12: E12, 13: E13}.items() if v}`, same for `byc`) |
| `mk_parts.py` (host check parts) | wrote 27 parts, 666 kinds, records 271,728; `--perturb load` and `--perturb win` REFUSED (exit 3) | none (rems 0..63, any histogram) |
| `mk_disp2979.py` (dispatch) | Lean content identical to x-tables' `Data/D2979/Disp.lean`; only the doc comment differs (it prints the input files' sha256) | none; run it on #2979's DISPATCH/ORDER bytes and the module is byte-identical, so `D2979.*` is reused unchanged |
| host-table call fields vs design A and #2979 (`tools/fill_indep2979Fc.py`) | `disp, time, order_k, windows, res, local_profile, child_of` the same; only `hrows`, `offs`, `hrem` differ (271,672 of 271,728 slots hosted elsewhere than in design A) | none |
| `mk_inst2979H.py` (renaming generator) | **cannot be used as is** (not run: it hard-codes aA-data's paths; outcome read from the code) | it is design-A-specific by construction: (a) it REFUSES unless the remainder histogram equals x-tables' (A-fc's is {1: 2,269, 6: 5, 7: 3, 17: 8, 19: 1, 23: 25,905, 48: 16,576}, #2979's {1: 17,310, 3: 4,956, 5: 1,728, 48: 28,190}), so `sched2979Fc_hcost` and `unit2979Fc_arith` must be GENERATED from the MANIFEST's `rem_hist`, not renamed; (b) F1 = 0: `f1_cnt : cntK rt 48 6 74 = 0`, and the prose table indexes `nonheads_by_child['12']` (KeyError on `{}`); (c) it expects a half-cut rate JSON (`children['24'] = H`, `children['48'] = F1`); (d) its bridge files (`Inst2979`, `HeadsFirst2979`, `HFTest2979`, `HFTest2979W`) are half cut and must NOT be generated for 2979Fc (afc-cut's NEEDS section 1 agrees) |
| `mk_profile2979H.py` | not applicable (half cut) | replaced by `tools/mk_profile2979Fc.py` (section 3) |
| `mk_l2979H.py` (statement stage, i-skel's generator) | **REFUSES the suffix `L2979Fc`**: `re.fullmatch(r'L[A-Z0-9]+', S)` (lower-case `c`). With only that regex relaxed to `L[A-Za-z0-9]+` (`dry/gen/mk_l2979Fc_dry.py`; e-r45m1's `mk_fourier_mixed.py` already accepts `[A-Z][A-Za-z0-9]*`): final `L2979Fc` N 17887045 OK, 20 files, both challenges byte-identical to the derivation from OpenAI's statement (1 - 17887045/10^10, 1 - 17887046/10^10); 17887046 REFUSED ("the rate point 1 - 1788704602/10^12 lies above the mixed root"); stub `L2979FcS` (N 16256629) OK, 21 files | (a) the suffix regex (relax it, or use `L2979FC`); (b) the generated `Seam.lean` doc comment (line 15) still says "half-cut ... the 16,576 literal rank-48 heads cut by 24 of their 48 axes": wrong text for A-fc, needs fc-engine's `--word` option (comment only, no Lean change); (c) the `Stmt.lean` header says "generated by aA-data's tools/mk_l2979H.py, round wht45/allheads" (provenance text) |

Conclusion: no generator assumes "every role sits in its planned slot". The only places that assume design A are
pins, the heads generator's empty-Counter comparison, and mk_inst2979H's fixed histogram / F1 / half-cut checks.

## 3. Cut profile, targets, separation, N [computed]

`tools/mk_profile2979Fc.py PROFILE.json 16576 tools/profile2979Fc.json --sha 5df79411...` (refuses with 16,575 heads:
"pieces above 24 left after removing the heads"):
- unit profile = PROFILE.json = {1: 1,630,639, 2: 470,196, 3: 154,512, 4: 71,928, 5: 11,988, 6: 5, 7: 3, 8: 5,328,
  17: 8, 19: 1, 23: 25,905, 48: 16,576} = 666 x {1: 2,445, 2: 706, 3: 232, 4: 108, 5: 18, 8: 8} + remainders
  {1: 2,269, 6: 5, 7: 3, 17: 8, 19: 1, 23: 25,905, 48: 16,576} (the converter recomputes the remainders from the raw
  FILL with our lowering's residual widths and gets exactly these); rank mass 4,816,512 = 74 x 65,088;
- all 16,576 rank-48 pieces are the literal heads (heads2979.py: rem-48 bins = head bins); they are removed, their
  arrays get target 25;
- cut profile `cut2979Fc` = {1: 1,630,639, 2: 470,196, 3: 154,512, 4: 71,928, 5: 11,988, 6: 5, 7: 3, 8: 5,328, 17: 8,
  19: 1, 23: 25,905} (2,370,513 blocks, mass 4,020,864 = 4,816,512 - 48 x 16,576); nothing to split;
- targets {73: 49,728, 25: 16,576} per gauge copy, W0 = 66,304 = 74 x 896, Hd = 16,576;
- **largest paid rank 23 < 25**: separated. (The Lean rule may keep rmax = 24 as on design B: `23 <= 24 < 25`.)
- `mk_rate_mixed.py scan 2979Fc=tools/profile2979Fc.json`: root 0.0017887045643030791939300419693164749977017732938...,
  **N = 17887045 proved, 17887046 refused**, slack 1.5e-10 relative, root - N = 0.643 grid steps, rate point N + 0.02.
  (Equals hx-price's three routes and hx-skeptic's root4.)

## 4. The Lean plan (suffix 2979Fc)

Base tree: design A's (ours/T + `lever-l2979H-overlay.tgz`, sha256 e768d770..., 286 files) — it already contains
HeadM1 (`cut_heads_full`, `splitWord`), p-bridge's Core, e-r12/e-r3cost (Rule, `Eng.engine_seam_final`), x-words'
providers, D2979 dispatch, the HFW2979 warp heads-first certificate, x-inst's generic `HF/Inst2979`.

**Carried over unchanged (from L2979H's closure)**: everything generic or fill-independent: engine R1-R4, warps,
x-words' stage (`PWG2979.*`), dispatch `D2979`/`D2979W`/`D2979Net`/`D2979Geom`/`Disp2979W`, p-a3's warp network
(`TProvW`, `TDefW`, `TStageW`, `TCertW`, `TGeomW`), p-hfirst (`HF/GeomA4`, `HF/GeomInst`), x-inst's `HFW2979/*` and
`HF/Inst2979`, `HF/Count2979C`, `Mixed/HeadM1`, `Mixed/Bridge/Core`, the i-skel statement generator and scripts.

**Carried over from L3945F unchanged**: `Work/Mixed/Bridge/CoreF.lean` (fc-engine; generic, imports only Core).

**New, regenerated data (afc-data)**: `Data/H2979Fc/*` (mk_parts on the new host table, ~45 files incl. Calls, P0..P26,
Glob with the new `rem_hist`, Whole, Order), `Data/H2979Fc/{Heads,NegHeads,DispLink}` (heads_cnt = 16,576,
f1_cnt = 0), `Inst/S2979Fc` (`sched2979Fc`, `_fits`, `_hrem`), `Inst/S2979FcFacts` (`sched2979Fc_heads`,
`sched2979Fc_hcost` with RHS `2269 φ1 + 5 φ6 + 3 φ7 + 8 φ17 + 1 φ19 + 25905 φ23 + 16576 φ48`, `resK2979Fc_*`),
`Inst/S2979FcUnit` (`unit2979Fc_arith`: `416 c0 + 250 c1 + (that idle sum) = listCost φ unit2979Fc`), `S2979FcHB`.
NOT regenerated: anything half-cut (`Inst2979H`, `HeadsFirst2979H`, `HFTest2979H(W)`, `Rate2979H*`).

**New, the cut (afc-cut)** (fc-engine's 3945F files with #2979's numbers; signatures as fc-engine INTERFACE section 2-4):
1. `Mixed/Rate2979Fc` (+`Sharp`, `R5Joint2979Fc`): `mk_rate_mixed.py emit 2979Fc=profile 17887045`:
   `rate2979Fc_moment : moment 73 zRate2979Fc cut2979Fc ≤ (66304 - 16576) + T 73 zRate2979Fc 25 16576`,
   `zRate2979Fc = 1 - 1788704502/10^12`, `cut2979Fc_proper` (ranks in [1, 24]).
2. `Mixed/Bridge/Inst2979Fc`: `PB.unit2979Fc` (= PROFILE.json list), `hcutF2979Fc` by
   `hcutF_of_check unit2979Fc cut2979Fc 24 24 16576 [1,2,3,4,5,6,7,8,17,19,23,24,48]` (key list must contain 24 and
   48 and every rank of both lists), negatives (16,577 heads; unit ≠ cut at 48); `bundle2979Fc_of_cut` /
   `engine2979Fc_of_cut` = Inst3945F's with `80640 → 16576`, `322560 → 66304`, `cut3945F → cut2979Fc`,
   envelope `1 - 1788704501/10^12`, `rmax = 24` (or 23).
3. `Mixed/Bridge/HFTest2979FcW`: `cutSF`, `cutSF_card`, `bph_cutF` restated (no import of HFTest3945HW: use
   x-inst's `I2979.hLlen_card` and `I2979.card_live2979`, as HFTest2979HW does); `engine2979Fc_geomW(_auto,_unit)`
   with `Fin 74`, `Fin 112`, head count 16,576.
4. `L2979Fc/Stmt` (generated), `L2979Fc/EngineOfGeomW` (= L3945F's with `Fin 74`, 16576, `unit2979Fc`),
   `L2979Fc/Prov` (= aA-final's `L2979H/Prov` with the 2979Fc schedule, `hfit2979FcW`, `hpos2979FcW`, `hU2979FcW`),
   `L2979Fc/Engine` (one term, as `MX.L2979H.engine`).
5. Statement stage: `mk_l2979H.py`-style generator `final out/final-L2979Fc --suffix L2979Fc --num 17887045
   --profile 2979Fc=...` (+ fc-engine's `--word` comment change); hyp stage; stub dry run; refusal at 17887046.
6. Integration: aA-final's recipe (assemble, build, axioms, greps, Compare pos / neg at 17887046, official comparator,
   bigospeedrun merged with h-bos's tool — L2979H's closure was 581 modules, so this one is about 581 + ~10 [argued],
   over 500, merge needed — and the unmerged stage on another Colab).

Expected sizes [argued from fc-engine and aA-data's logs]: the cut layer ~600 lines (mostly copies), the data layer
generated; data build ~11 min (aA-data: 644 s), cut layer ~5 min, integration ~3-4 h (comparator ~1 h, bigospeedrun
~1-1.5 h per phase).

## 5. Theorem 1G and the Lean rule for this job model [read; argued]

Same as fc-gate section 3 with A-fc's numbers: targets {73, 25}, paid ranks {1..8, 17, 19, 23}, separation 23 < 25;
Lean `Rule.twoTarget 73 25 (card σ) |L| P 24 rk` meets `rk_le ≤ 24`, `24 < 25`, `2 ≤ 25`, `25 ≤ 73`, `H = G x 16,576
≥ 1`. The hub engine's halving guard (m^2 > 2 r_max^2) is not a hypothesis of Theorem 1G or of the Lean rule
(hx-skeptic item 5; design B is proved in Lean with the same 25 target).

## 6. Refinements

- 10:24 EDT: statement generator dry run (section 2, last row): suffix regex and the half-cut doc comment of Seam.lean.
- 10:27 EDT: the rate module regenerated independently from my profile (`mk_rate_mixed.py emit`, `dry/rate/`): its
  Lean text equals afc-cut's `Rate2979Fc.lean` except two doc lines (the generator says "half-cut", afc-cut wrote
  "full-cut"); `check 17887046` REFUSED.
- 10:40 EDT: closure size [argued from aA-final's MODULES.txt]: L2979H's closure has 66 design-A modules (39 host-table
  data, 16 stage, 11 others); L2979Fc replaces them one for one (data ~39, stage 17, `Rate2979Fc`, `CoreF`,
  `Inst2979Fc`, `HFTest2979FcW`, `S2979Fc`, `S2979FcFacts`, `L2979Fc/{Stmt,EngineOfGeomW,Prov,Engine}`; no
  `HeadsFirst2979H`/`HFTest2979H` counterparts), so about 580 modules: the bigospeedrun stage needs h-bos's merge
  (L2979H: 581 -> 447). No generic module of A's overlay imports a 2979H module [computed: grep of the overlay].
- 10:42 EDT: `L2979Fc/Prov.lean`'s `hU2979FcW` must regenerate the literal idle sum in its `show` line
  (`(17310 * φ 1 + 4956 * φ 3 + 1728 * φ 5 + 28190 * φ 48)` in L2979H) and import `Inst2979Fc` instead of `Inst2979H`;
  `D2979W.unitPriceTW_const` and `IcW_cost` are generic [read].
- 11:16 EDT: section 1d [Lean-checked C4]; verdict v2.

## 7. Commands

```
cd <work>/wht45/afc/afc-gate
D=<work>/wht45/fullcut/fc-reopt/r2/designs/u1-2979sep/design
A=<work>/wht45/allheads/aA-data
HA=<work>/wht45/lever/p-alt/dl/main2979/research/rank7-joint-endframes
V=<work>/wht45/lever/g-heads/dl/pr3945/owner/inputs
python3 -I -B tools/mk_profile2979Fc.py $D/PROFILE.json 16576 tools/profile2979Fc.json --sha 5df7941120150492c64e62dc3362095eff3dfac9facff041ebb5c25b2e78f75e
python3 -I -B tools/mk_profile2979Fc.py $D/PROFILE.json 16575 out/never.json            # REFUSED exit 2
python3 -I -B <work>/wht45/lever/e-r45m1/tools/mk_rate_mixed.py scan 2979Fc=$PWD/tools/profile2979Fc.json
# dry runs (dry/alloc = #2979's ORDER/DISPATCH/PROVIDERS + A-fc's FILL/PROFILE)
python3 -I -B tools/own_to_host2979Fc_dry.py dry/alloc $HA/data $V $A/in/residuals-2979.json dry/H2979Fc
python3 -I -B $A/tools/heads2979.py dry/H2979Fc-hosttab.json.gz dry/alloc/FILL.json.gz $V/TEMPLATE.json.gz - dry/heads2979Fc.json.gz dry/heads2979Fc-summary.json 8288,8288,0,0
python3 -I -B $A/tools/mk_parts.py dry/H2979Fc-hosttab.json.gz H2979Fc 25 dry/parts/H2979Fc     (+ --perturb load|win: REFUSED)
python3 -I -B $A/tools/mk_disp2979.py $D/DISPATCH.json $D/ORDER.json 74 416 250 dry/disp/Disp.lean
python3 -I -B tools/fill_indep2979Fc.py dry/H2979Fc-hosttab.json.gz $A/out/H2979H-hosttab.json.gz <work>/wht45/lever2979/x-tables/out/H2979-hosttab.json.gz
# statement stage, regex relaxed (dry/gen/mk_l2979Fc_dry.py = aA-data's mk_l2979H.py with L[A-Za-z0-9]+)
python3 -I -B dry/gen/mk_l2979Fc_dry.py final dry/final-L2979Fc --suffix L2979Fc --num 17887045 --profile 2979Fc=$PWD/tools/profile2979Fc.json
python3 -I -B dry/gen/mk_l2979Fc_dry.py final dry/never --suffix L2979Fc --num 17887046 --profile 2979Fc=...   # REFUSED exit 2
python3 -I -B dry/gen/mk_l2979Fc_dry.py stub dry/stub-L2979FcS --suffix L2979FcS --profile 2979Fc=...
python3 -I -B <work>/wht45/lever/e-r45m1/tools/mk_rate_mixed.py emit 2979Fc=$PWD/tools/profile2979Fc.json 17887045 dry/rate/lean
# the Lean try (C4): inputs colab/in (SHA256SUMS), script colab/gate-lean.sh
COPYFILE_DISABLE=1 tar --no-xattrs -czf - -C colab gate-lean.sh in | <ssh to the build machine> 'mkdir -p /content/work/afc-gate/scripts /content/work/logs/afc-gate && tar xzf - -C /content/work/afc-gate/scripts && mv /content/work/afc-gate/scripts/in /content/work/afc-gate/in'
<ssh to the build machine> 'tmux -L afc-gate new -d -s gate "bash /content/work/afc-gate/scripts/gate-lean.sh ta > /content/work/logs/afc-gate/gate.out 2>&1"'
```
Logs: `logs/` (Mac runs), `logs/C4/` (the Lean try), `logs/C3/gate-try1.out` (the stopped first try).
