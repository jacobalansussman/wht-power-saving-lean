# hx-admit: A-fc through every hub check of #2979's line (round wht45/hubx, 20 Oct 2026)

Agent hx-admit. Not public text. Hub code ran only on C5, in `/content/work/hx-admit/`, from a FRESH blobless clone
of CrocSwap/integer-mult-bounds (origin/main = `12fa8c5c844be21bfefb38a19503fbd867a11244`, our release's base) plus
the eight predecessor commits #2979's README names. #2979's package is a `git archive` of that commit; its tree
(`9981a8c1`) is the same as at 20ecf289, which ah-gate used. Venv `vA` has the authors' recorded versions (CPython
3.12, numpy 2.4.5, scipy 1.16.1, sympy 1.14.0, mpmath 1.3.0, gmpy2 2.3.2, psutil 7.2.1). No Lean. Nothing posted.
Labels: [hub code on C5] (unmodified), [our glue] (our code calling unmodified hub modules), [computed].

## Verdict (05:40 EDT): A-fc passes every fill check; it is rejected only where the hub asks for a higher hub b

- **Step 1** [computed]: A-fc's fill re-created from its class solution with fc-reopt's recipe on the fresh 12fa8c5c
  exports is the same fill. The decompressed content is identical. The `.gz` bytes differ only in the gzip header's
  4-byte time stamp, and with the published time stamp the sha256 is exactly 7fc766c6....
- **Step 2** [hub code on C5]: on the unmodified hub code, A-fc passes:
  - #2979's PR2839 owner engine (in `check_allocation.py`, 5 controls; standalone with its own 7 controls);
  - `verify.py` (86 controls);
  - the complete replay's fresh natives and the admission's owners, parent lift, rebind, discovery (7 steps) and
    outer ledger (63 positive checks, independent check);
  - the release's own `verify.py` (16 sub-stages, `PASS_T74_EXACT_FILL_COMPLETE_FINITE_ADMISSION`, 139 strict
    admission controls rejected, plus bit 12, controller 15, independent-owner 8, geometry 6) and its outer checks.
- **Where it is rejected.** Four places, every one a comparison of A-fc's hub b (1.6897030e-3) with a higher figure:
  1. `check_allocation.py`'s moment auditor: `ValueError: strict frontier improvement and lower endpoint`.
  2. `check_cap.py`: `AssertionError: recorded gap`. The literal root must lie within 1e-6 of the fixed-inventory cap
     1.7290358e-3; A-fc's is 3.93e-5 below it. Its four strict-gain asserts come after this one.
  3. The complete admission (`admit_extension.py` line 388): `RuntimeError: candidate must strictly improve full exact
     pinned frontier headline before full release`. The other half of that `need`, kappa = C - 1e-70 - 1e-46, holds.
  4. The release, run by our glue after (3): `ledger kappa does not exceed PR2269`. This is the first of its
     strict-gain list over earlier packages that lies above 1.6897e-3; the 14 lower entries before it pass.
- **Item 2 is not literally a strict-gain assert.** It is a near-optimality bound on the hub b (gap to the relaxation
  cap below 1e-6). It rejects A-fc for the same reason, its lower hub b, and not for any property of the placement.
  All of check_cap's other checks hold [our glue, rest_afc.py]: the cap certificate, rank mass, local profile paid,
  root at most the cap, and the finer-split control. The hub text should name it beside strict gain. **No STOP
  written**: no fill check failed.
- **Positive controls** (#2979 published, design A), same harness, same fresh clone:
  - every check above passes;
  - both complete replays end in `PASS_RANK7_JOINT_ENDFRAMES_COMPLETE_REPLAY`, κ exactly #2979's claim
    0.00172903398264651351771699..., `PASS_RELEASE`, with a fresh retained bootstrap;
  - `run_registered.py` passes on the published package. On design A it stops only at the stored receipt binding of
    the old fill hash, as in ah-gate.

## Step 1: A-fc's fill re-created from its class solution [computed, fc-reopt's code on C5]

- Recipe: fc-reopt's `r2/run2979.sh` exactly (same command, same environment variables), with the hub inputs taken
  from the fresh clone: #3815's construction `research/staged544-clock-ceiling` at 12fa8c5c (tree `e867b76f`, the same
  tree as #3815's head 77ffbf28) and #2979's package at 12fa8c5c. Class solution `SOL-MIP0.json` (96afbf48...).
  Script `colab/recreate.sh`, log `logs/recreate.out`, inputs `logs/recreate-inputs.sha`, output `out/recreate/`.
- Result: the tolerant placement needed the same 9 spare roles (8 of width 6 at call 553, 1 of width 4 at call 0) and
  dropped 8 planned width-6 slots, as FOUND.md says. The built-in engine check passed with its 7 controls.
- **The fill is the same.** The decompressed JSON is byte-identical (sha256 a79f48af...). The `.gz` file differs only
  in the 4-byte gzip header time stamp (Python's `gzip.compress` writes the current time). With the published file's
  time stamp put back, the sha256 is exactly **7fc766c6c3e1561cee60b258138a918020b96d85ea782a2b01ad0ee035fd082e**.
  PROFILE.json, ORDER.json and DISPATCH.json are byte-equal to the published A-fc files. A-fc's ORDER and DISPATCH are
  content-equal to #2979's (666 calls each); the package copies keep #2979's own ORDER/DISPATCH bytes.
- All checks below run on the published A-fc file (7fc766c6...), the one the hub PR would carry.

## Step 2: every hub check that applies to a #2979 fill

Package copies (`colab/setup.sh`, `colab/repin.py`): `pub` = #2979 unmodified (positive control); `A` = design A's fill
764395fc... (positive control); `Afc` = A-fc's fill and its PROFILE.json; `AfcC` = as `Afc`, plus CLAIM.json's
`profile` field set to A-fc's profile (the only CLAIM field `check_allocation.py` binds before its strict gain). In
each copy only these data files change, and ALLOCATION-PINS.json / MANIFEST.json are re-pinned in their own format.
No hub code is changed. Every copy's MANIFEST still matched after the checks ("manifest-equal True").

### 2a. Registered checks [hub code on C5] (`colab/reg.sh`; logs `logs/reg/`, summary `logs/reg/SUMMARY.txt`)

| check (script, hub code) | #2979 published | design A | **A-fc** |
|---|---|---|---|
| `verify.py` (three-stage local replay of the words, 86 controls) | PASS `PASS_RANK7_JOINT_ENDFRAMES_THREE_STAGE_REPLAY` | PASS | **PASS** (the fill does not enter it) |
| `check_allocation.py`: pins, words, PROVIDERS | pass | pass | **pass** |
| `check_allocation.py`: PR2839 owner engine (`allocation/vendor/independent_owners.py` `check`, claimed profile = PROFILE.json) | pass | pass | **pass** |
| `check_allocation.py`: physical 66,813 | pass | pass | **pass** |
| `check_allocation.py`: 5 rejection controls (missing_role, duplicate_helper, same_call_family_collision, capacity_overflow, changed_paid_profile) | 5/5 rejected | 5/5 | **5/5 rejected** (the run reached line 105, after them) |
| `check_allocation.py`: PR2917 exact moment auditor (`tools/audit_allocation_moment.py`), rank mass, root bracket, exact bisection | pass, root 0.00172903398264651351771... | pass | **pass, root 0.00168970301113743521378...** = 1.6897030e-3 |
| `check_allocation.py`: auditor's `need(target>prior and target<lower, ...)` | pass (gain +4.9988e-6 over PR2969) | pass | **REJECTED (strict gain)**: `ValueError: strict frontier improvement and lower endpoint` (gain -3.4332e-5; the second conjunct, target below the lower endpoint, holds) |
| `check_allocation.py`: claim binding, `strict gain over the PR2969 headline nomination`, input drift | pass | pass | not reached (after the strict gain) |
| `check_cap.py`: cap certificate, rank mass, local profile paid, literal root at most the cap | pass | pass | **pass** |
| `check_cap.py`: `recorded gap` (cap minus literal root < 1e-6) | pass (1.78e-9) | pass (1.78e-9) | **REJECTED**: `AssertionError: recorded gap` (gap 3.9333e-5: the same hub-b deficit, measured against the cap) |
| `check_cap.py`: finer-idle-split control; strict gain over PR2969, PR2966, PR2957, PR2942 | pass | pass | not reached; computed with the package's own modules [our glue, `colab/rest_afc.py`]: control holds, all four gains negative (-3.43e-5, -3.43e-5, -3.40e-5, -3.31e-5) |
| `run_registered.py` (all of the above, then the stored-receipt bindings) | **PASS** `PASS_RANK7_JOINT_ENDFRAMES_LOCAL_OWNER_MOMENT_CAP_AND_STORED_RELEASE` | verify, check_allocation, check_cap PASS; stops at the stored receipt binding `admit['fill_sha256'] == sha(FILL) == claim['fill_sha256']` (as in ah-gate; the complete replay replaces it) | verify PASS; stops at check_allocation's strict gain (same message) |
| `AfcC` (claim profile = A-fc's) | | | identical verdicts: the strict gain comes before the claim binding |

### 2b. PR2839's `independent_owners.py` with its own seven controls [hub code on C5, our driver `colab/own7.py`]

The engine's `main()` reads PR2802's provider words (`vendor/frontier`), which #2979 does not ship. So `own7.py` loads
the vendored engine unchanged, builds the providers exactly as `check_allocation.py` does, and runs the body of
`main()` (copied verbatim: `check` plus its seven `reject` controls, each with its required reason).

| | #2979 published | design A | **A-fc** |
|---|---|---|---|
| status | `PASS_INDEPENDENT_T_OWNER_CHRONOLOGY_MASKS` | same | **same** |
| roles / T | 271,728 / 74 | same | **271,728 / 74** |
| idle | {1: 17,310, 3: 4,956, 5: 1,728, 48: 28,190} | same | **{1: 2,269, 6: 5, 7: 3, 17: 8, 19: 1, 23: 25,905, 48: 16,576}** |
| controls rejected with their exact reasons | 7/7 | 7/7 | **7/7**: false_replica_count, missing_literal_helper_role, unbound_provider_at_actual_call, unavailable_source_window, same_call_physical_family_collision, unpaid_hole_in_compact_pending_prefix, changed_total_paid_profile |

Logs: `logs/reg/own7-{pub,A,Afc,AfcC}.{out,json}`.

### 2c. The complete replay, `reproduce.py` [hub code on C5] (`colab/full.sh`; logs `logs/full-*.out`, `logs/rp-*/`)

Each run: `reproduce.py --repo <fresh clone> --work rp-<v> --jobs 2`, nice 5, the three runs in parallel on C5
(04:44-05:38 EDT). Each run makes exact Git exports, a fresh retained bootstrap and fresh natives, then the complete
admission.

| stage (reproduce.py / admit_extension.py) | #2979 published | design A | **A-fc** |
|---|---|---|---|
| export-pr2894, export-pr2833 (exact exports, manifest-pinned) | PASS | PASS | **PASS** |
| fresh-retained-bootstrap | PASS (1505 s) | PASS (1695 s) | **PASS** (1418 s) |
| four-kernel, preload, parent-paired, paired-extension stages (PR2917 derived source 98d57c82) | PASS | PASS | **PASS** |
| rank7-overlay-stage (overlay 751b4636) | PASS | PASS | **PASS** |
| derive-closed-adapters | PASS | PASS | **PASS** |
| fresh-current-native (fresh five-stage natives, both providers) | PASS | PASS | **PASS** |
| admission: owners (`reconstruct.py`) | PASS | PASS | **PASS** `PASS_OWNER_VERIFY`, idle {1: 2,269, 6: 5, 7: 3, 17: 8, 19: 1, 23: 25,905, 48: 16,576} |
| admission: parent-lift | PASS | PASS | **PASS** `PASS_LITERAL_GENERAL_T_PARENT_PLAN_SPEC` |
| admission: rebind | PASS | PASS | **PASS** (new profile = A-fc's PROFILE.json) |
| admission: discovery (bind, emit, audit, invoice, assembly, freeze, independent) | PASS | PASS | **PASS, 7/7**: `PASS_GAUGE44_NATIVE_EXACT_GLOBAL_ASSEMBLY` (complex_lower 0.001689703011137435213783, bit_lower 0.001733427039589844234749); `PASS_INDEPENDENT_VARIABLE_T_LITERAL_WINDOWS_CHARTS_PREFIX_ROW` (271,728 roles) |
| admission: outer `check_ledger.py`, `independent_check.py` | PASS | PASS | **PASS**: FINAL KAPPA 0.00168970301113743521378299..., 63 positive checks; independent check passes |
| admission: `need(kappa == C-1e-70-1e-46 and kappa > previous)` | PASS | PASS | **REJECTED (strict gain)**: `RuntimeError: candidate must strictly improve full exact pinned frontier headline before full release` (kappa 1.6897030e-3 vs PR2969 1.7240352e-3; the first conjunct holds) |
| admission: release (`verify_release.py`) | PASS `PASS_RELEASE` (967 s) | PASS `PASS_RELEASE` (944 s) | not reached in the replay; run by our glue, see below |
| replay result | **`PASS_RANK7_JOINT_ENDFRAMES_COMPLETE_REPLAY`**, κ = claim | **same**, κ = claim, fill_sha256 764395fc... | rc 1 at complete-admission (as expected) |

**The release on A-fc, with only the strict-gain step skipped** (our glue `colab/release_afc.py`, 982 s, log
`logs/rp-Afc/candidate/release.log`, `logs/rp-Afc/HX-RELEASE.json`). The glue repeats what the adapter does between
the step and its release call: the OUTER-CLAIM.json update, the provenance copy and MANIFEST.json, with the same
values recomputed from the same files. It then runs the package's own `verify_release.py` unchanged, with the
adapter's arguments, environment and working directory. This is not a hub admission.

- The supplier's own `verify.py` passes all 16 sub-stages, ending in `PASS_T74_EXACT_FILL_COMPLETE_FINITE_ADMISSION`
  (strict-t74 `PASS_STRICT_MIXED_T74_EXACT_FILL`). Its controls: strict admission 139 rejected, 0 accepted
  (`PASS_EXACT_MIXED_048_ALLCOLD_NATIVE_ASSEMBLY_ADMISSION_NEGATIVE_CONTROLS`); bit 12; controller 15; independent
  owner 8; geometry 6 (`PASS_T74_GEOMETRY_AND_CHARGE_CONTROLS`). These are the same counts as the two controls.
  - The 16 sub-stages: fine-engine-oracle, retained-provider-and-parent-admission, rebuild-words-566,
    reconstruct-owners, lift-parent, rebuild-borrow-words, shared-kernel-independent-review,
    preload-move-independent-review, two-sided-helper408-review, paired-extension-review, audit-local-566-base,
    audit-local-566-extra, rebind-566, lemmag-labels, strict-t74, geometry-controls.
- The outer checks pass: `PASS_UNION_SCHEDULE`, `PASS_EXACT_REPLAY 19`, `PASS_OUTER`.
- The release's binding checks pass: the strict assembly, the endpoint bracket, the profile and bit endpoint, and
  ledger κ = OUTER-CLAIM.
- The release then stops at `ledger kappa does not exceed PR2269`. That is its list of strict gains over earlier
  hub packages: the 14 entries below 1.6897e-3 (PR1673 ... old union row) pass, and PR2269 (1.6962e-3) is the first
  one above. The supplier κ that `verify.py` prints, 0.001686852730998521, is the discovery's frozen witness. For
  design A the same field is 0.001726049584259595, so it is not the headline.

So the hub's whole chain judges A-fc's fill valid, and rejects it only because its hub b is lower.

## Files

- `STATUS.md` (milestones), `COMMANDS.txt` (every command), this NOTE.
- `colab/`: setup.sh, repin.py, full.sh, reg.sh, own7.py, rest_afc.py, release_afc.py, recreate.sh, `fcr/` (fc-reopt's
  tools, byte-identical to its C1 copies).
- `in/` + `in/SHA256SUMS`: the fills under test. `out/recreate/`: the re-created A-fc design.
- `logs/`: everything copied back from C5 (`/content/work/logs/hx-admit/`).
