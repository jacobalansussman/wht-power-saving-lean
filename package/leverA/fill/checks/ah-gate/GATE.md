# GATE for the all-heads fills (ah-gate, round wht45/allheads)

First verdict 19 Oct 2026, 14:58 EDT, after step (1). **Final verdict 15:58 EDT: all three steps done, both designs
PASS, no STOP written.** Agent ah-gate, independent of o-model (who made the fills).
Hub code ran only on Colabs, in my own folders `/content/work/ah-gate/`:
- C1: the hub's checks and the complete replays;
- C2 and C5: the operator replay.

The hub code came from fresh `git archive` exports of the hub's objects. I used my own uv venvs with the authors' pins.
On the Mac I ran only light Python: my code and our own tools. No Lean. Nothing posted.

Labels: [hub code on Colab], [computed] (our code), [read on hub], [argued], [what-if] (priced, not built).

## Verdict: B and A are both sound in every check the hub has (GO for both)

| | design B (#3945, T360) | design A (#2979, T74) |
|---|---|---|
| fill under test (sha256) | final 9c06677d...; #3815-stage 0e61fce2... | 764395fc... |
| (1) the hub's fill checkers, fresh copies, with controls and published positive controls | **PASS** (#3945 `owner/verify.py` plain and `-O`; PR3734 `owner_audit.py`; #3815 `verify.py` up to its stored-receipt binding) | **PASS** (`check_allocation.py`; `verify.py`, 86 controls; `check_cap.py`) |
| (1) heads: 5 counters (mine, recount.py, heads_fill.py, o-eval, our lowering G0) | **80,640** each (40,320 + 40,320) | **16,576** each (8,288 + 8,288) |
| (1) paid profile = published | yes, digit for digit | yes, digit for digit |
| (1) half cut, two routes, 24 digits, separated 48 < 49 [what-if] | **1.7581520e-3** (0.001758152016764631555031) | **1.7564708e-3** (0.001756470800681586217656) |
| (2) the hub's COMPLETE admission with the new fill (fresh natives, owners, parent lift, discovery, outer, release) | **PASS** `PASS_TERMINAL_WIDTH5_COMPLETE_CONDITIONAL_RELEASE`, 119/119 controls, κ = claim, `PASS_RELEASE` | **PASS** `PASS_RANK7_JOINT_ENDFRAMES_COMPLETE_REPLAY`, 139/139 controls, κ = claim, `PASS_RELEASE` |
| (3) a hub rule that forbids or prices the placement | **none found**; the operator replay (#4046's model) gives 0 differences on 2,029,945 characters, all 360 owners | **none found**; operator replay 0 differences on 594,132 characters, all 74 owners |

So design B (the headline) gives 1.7581520e-3 and design A gives 1.7564708e-3 [what-if, mixed price; base b unchanged
at 1.7306632e-3 and 1.7290340e-3]. Both are now backed by the hub's complete admission on the exact new fills. No kill
rule fires. The port notes are in section 3. The main one: B is the first #3945-line fill that hosts in children 12/13.

### Step (1) summary (first verdict, 14:58 EDT)

| check | design B (#3945, T360) | design A (#2979, T74) |
|---|---|---|
| the hub's own fill checker, fresh copy | **PASS** (#3945 `owner/verify.py`, plain and `-O`) | **PASS** (#2979 `check_allocation.py`) |
| its positive control (unmodified published files, same harness) | PASS | PASS |
| its built-in rejection controls | 4 of 4 rejected | 5 of 5 rejected |
| second hub checker | **PASS** (PR3734 `owner_audit.py` on the #3815-stage fill; 3 of 3 controls) | **PASS** (`verify.py` 86 controls, `check_cap.py`, inside `run_registered.py`) |
| my negative controls | old/new fill mismatch (both ways): rejected; PR3734's clock with the B fill: rejected | (built-in controls) |
| heads: my own recount (`tools/ahcount.py`) | **80,640** (40,320 + 40,320) | **16,576** (8,288 + 8,288) |
| heads: lever-reprice `recount.py` / `heads_fill.py` | 80,640 / 80,640 HEAD-literal (bridge-first control: 0) | 16,576 / 16,576 (control: 0) |
| heads: o-eval `evalfill.py` | 80,640, recount agrees | 16,576, recount agrees |
| heads in OUR lowering (G0, `g0_staged.py`, our gcert words, our hosted network) | 80,640, 0 violations, call blocks <= 8; controls: gate first gives 0 | 16,576, 0 violations; controls likewise |
| paid profile vs the published PROFILE.json | **equal, digit for digit** (3 tools; also the hub's verify.py) | **equal, digit for digit** (3 tools; also check_allocation.py) |
| uniform b (hub moment engine re-implementation) | 0.001730663236534549309146 = #3945's stated root | 0.001729033982646513517717 = #2979's stated root |
| **half cut, route 1** (o-eval evaluator) | **0.001758152016764631555031** | **0.001756470800681586217656** |
| **half cut, route 2** (release tools/lever: pricer.py enclosures, mixed.py) | **0.001758152016764631555031** (1e-24 grid certificate OK) | **0.001756470800681586217656** (certificate OK) |
| separated (48 < 49) | yes | yes |

On the published fills the same tools reproduce the known figures: #3945 has 57,312 heads and 1.7501109e-3; #2979 has
16,086 heads and 1.7556473e-3.

## 1. Details of step (1)

### Design B [hub code on C1]

Logs: `logs/regB/` (copies of `/content/work/logs/ah-gate/regB/`). Script: `colab/regB.sh`. Inputs: fresh exports
of pr3945 91f81d13, pr3815 77ffbf28 (with #3734's `staged544-producer-return`); file list in `logs/src-sha256.txt`.

- **#3945 `owner/verify.py`** on copies of `owner/`. In the design copy only `inputs/old-FILL.json.gz` and
  `FILL.json.gz` are replaced, with their pins in `SOURCES.json` and `MANIFEST.json` updated. `PROFILE.json` is kept
  as published, so the checker compares the new fill's histogram with #3945's own file.
  - Published (control): `PASS_LITERAL_ENLARGED_OWNER_FILL`, plain and `-O`.
  - Design B: `PASS_LITERAL_ENLARGED_OWNER_FILL`, plain and `-O`. Idle {3: 23,328, 5: 4,320, 48: 137,952}; changed
    bins 84,240; profile = #3945's. The four built-in controls are rejected (same_call_family, offset_hole,
    closed_owner_window, overflow).
  - Negatives: design B's fill against #3945's published old-FILL, and the reverse. Both are rejected with
    "owner assignment moved", which shows the checker is bound to the pair of fills.
- **PR3734 `owner_audit.py`**, byte-identical copy, in the layout #3815's `verify.py` uses: #3815's ORDER, FILL and
  PROFILE beside #3734's DISPATCH, NATIVE, PLAN and TEMPLATE.
  - #3734's own inputs, #3815's published fill and the design B #3815-stage fill all give
    `PASS_STANDARD_LIBRARY_INDEPENDENT_LITERAL_AUDIT`. For B: idle {1: 84,240, 3: 23,328, 5: 4,320, 48: 137,952}, the
    same as #3815's; controls same_call_family, duplicate_helper, unpaid_hole rejected.
  - Negative: PR3734's clock with the B fill is rejected ("strict lower and closed upper window"), as #3815's own
    control requires.
- **#3815's `verify.py`** on a package copy with the B #3815-stage fill. The positive control passes. The design copy
  passes steps 1-4a: clock, owner audit with controls, the block-norm alias premise (no role of call (s, r) in owner
  r's pair stage_pair[s]), exact price, ceiling certificate and ideal idle histogram. It stops at step 4b, "release
  input binding": the stored release receipt names #3815's fill hash. That is a stored-evidence binding, which the
  complete admission (step 2) replaces.
- **#3945's eight registered commands** on a package copy with the B fills and `owner/` re-pinned.
  - Pass: `independent_frames.py --controls`, `exact_root.py --controls`, `owner/verify.py`, `check_gain.py`.
  - `verify.py` and `test_verify.py` (plain and `-O`) first stop at the top-level `INPUTS.json` pin of
    `owner/inputs/old-FILL.json.gz`. With that one pin updated they pass: `PASS_TERMINAL_GAP_WORDS` and
    `PASS_6_RECIPE_CONTROLS` (`logs/regB2.txt`).
  - `admission/authenticate.py` stops at "complete admission frozen inputs": the stored admission evidence binds the
    old fill. That is what step 2 reruns.

### Design A [hub code on C1]

Logs: `logs/regA/`. Script: `colab/regA.sh`. Fresh export of main 20ecf2891. In the design copy only
`allocation/FILL.json.gz` is replaced, with its entries in `ALLOCATION-PINS.json` and `MANIFEST.json`.

- **`check_allocation.py`** (PR2839 owner engine, five rejection controls, PR2917 exact moment auditor):
  - published (control): `PASS_RANK7_LITERAL_OWNERS_AND_EXACT_MOMENT`;
  - design A: `PASS_RANK7_LITERAL_OWNERS_AND_EXACT_MOMENT`, with the same idle {1: 17,310, 3: 4,956, 5: 1,728,
    48: 28,190}, physical 66,813, κ target 0.00172903398264651351771699... and gain over #2969 4.998757898408639e-6;
  - five controls rejected: capacity_overflow, changed_paid_profile, duplicate_helper, missing_role,
    same_call_family_collision.
- **`run_registered.py`**:
  - Published (control): `PASS_RANK7_JOINT_ENDFRAMES_LOCAL_OWNER_MOMENT_CAP_AND_STORED_RELEASE`.
  - Design A passes `verify.py` (`PASS_RANK7_JOINT_ENDFRAMES_THREE_STAGE_REPLAY`, 86 controls), `check_allocation.py`
    and `check_cap.py` (gap 1.78e-9).
  - It then stops at the stored admission receipt binding, `admit['fill_sha256'] == sha(FILL) ==
    claim['fill_sha256']`. That is expected, and step 2 reruns the admission.

### Heads, profile, price [computed, Mac, light]

- **My recount** (`tools/ahcount.py`, my own code written from the data formats):
  - It rebuilds the owner windows from the template schedule and the widths from the words (dim helper_final -
    dim start_extra), then walks every role (identity order, window, one role per array per call, compact prefix,
    capacity).
  - It counts a bin as a literal head when its child's idle sweep is the first template touch of its array and no
    role sits on that array at a call at or before the sweep.
  - Controls in every run: one role re-bound into a head bin gives one head fewer; a bridge on pair 3 inserted before
    the first idle sweep gives 0 heads.
  - Outputs: `out/count-*.json`.

  | fill | heads (child 8 + 10) | empty rank-48 bins in children 12 / 13 | profile = published |
  |---|---|---|---|
  | #3945 published | 57,312 (28,672 + 28,640) | 40,320 / 40,320 | yes |
  | **B final** | **80,640 (40,320 + 40,320)** | 28,672 / 28,640 | **yes (#3945's)** |
  | #3815 published | 57,312 | 40,320 / 40,320 | yes |
  | **B #3815-stage** | **80,640** | 28,672 / 28,640 | **yes (#3815's)** |
  | #2979 published | 16,086 (8,064 + 8,022) | 6,512 / 5,592 | yes |
  | **A** | **16,576 (8,288 + 8,288)** | 6,342 / 5,272 | **yes (#2979's)** |

- **lever-reprice** `recount.py` and `heads_fill.py`:
  - They give the same numbers: B 80,640 HEAD-literal, A 16,576.
  - The bridge-first control gives 0 heads on every fill.
  - `heads_fill` reproduces the full profile of each design, and it equals the published one. Outputs:
    `out/hf-*.json`.
  - Side observation [computed]: on B, the HEAD-net (argued) class is 4,272 rank-5 blocks. It was 27,648 on #3945:
    the 23,328 rank-3 net heads have become literal heads.
- **Price.**
  - Route 1: o-eval `evalfill.py` on design specs that point at my fresh exports (`designs/`, `out/eval-*.json`).
  - Route 2: `tools/price2.py`, run on MY recount's profile and heads. It uses the release candidate's
    `tools/lever/pricer.py` for the uniform b (certified on the grid). It uses `lever_price.variants` (mixed.py,
    60 digits) for the half cut. It also certifies the half-cut root on the 1e-24 grid, with pricer.py's own rational
    enclosures applied to the mixed balance (`out/price2-*.json`).
  - Both routes agree to all 24 digits on all four fills.
  - Side note for the planner [what-if]: with 80,640 heads the full cut on B prices at 1.7806935e-3 (A: 1.7797095e-3).
    It is no longer at break-even as it was on #3945. Whether the Lean route takes target 25 is not my question.

## 2. Step (2): the hub's complete admission with the new fills: BOTH PASS [hub code on C1]

The runs used fresh copies of the packages in which only the fill changes, my own venvs with the authors' pins, and
`tmux -L ah-gate`. They ran 14:24-15:53 EDT on a shared C1. The scripts are `colab/fullB.sh` and `colab/fullA.sh`.
Copies of the logs and results are in `logs/frB/` and `logs/rpA/`.

| | B: #3945 `admission/full_replay.py` (as g-sound ran it for #3945) | A: #2979 `reproduce.py` (as x-gate ran it for #2979) |
|---|---|---|
| result | **`PASS_TERMINAL_WIDTH5_COMPLETE_CONDITIONAL_RELEASE`**, rc 0 | **`PASS_RANK7_JOINT_ENDFRAMES_COMPLETE_REPLAY`**, rc 0 |
| fill the admission ran on (its own `fill_sha256`) | 9c06677d... (design B) | 764395fc... (design A) |
| natives | fresh five-stage runs on both providers (`new_native_reuse` false/false); the cache is fresh | fresh retained bootstrap (`explicit_retained_reuse` false, 40 min) and fresh current natives |
| stages | prepare (the exact 666-file source closure); natives; owners `PASS_OWNER_VERIFY`; parent lift `PASS_LITERAL_GENERAL_T_PARENT_PLAN_SPEC`; rebind; discovery (bind, emit, audit, invoice, freeze, independent); outer `check_ledger` / `independent_check`; release | the exports; bootstrap; four-kernel; preload; parent-paired; paired-extension; rank-7 overlay; adapters; fresh natives; complete admission (owners, parent lift, rebind, discovery, outer, release) |
| release verify (inherited supplier) | all 21 sub-stages PASS, including `staged544-independent-word-review`, `terminal-width5-independent-word-review`, `strict-t74` and `geometry-controls`, then `PASS_T74_EXACT_FILL_COMPLETE_FINITE_ADMISSION` | the same supplier chain, `PASS_T74_EXACT_FILL_COMPLETE_FINITE_ADMISSION` |
| strict admission controls | **119 rejected, 0 accepted** (five widths, as #3945) | **139 rejected, 0 accepted** (six widths, as #2979) |
| outer ledger | 63 strict inequalities, 17 rejecting controls, 19 exact evidence replays | the same |
| κ | exactly 0.00173066323653454930914599... = #3945's claim | exactly 0.00172903398264651351771699... = #2979's claim |
| profile (admission) | = #3945's `owner/PROFILE.json` | = #2979's `allocation/PROFILE.json` |
| release receipt | `PASS_RELEASE` | `PASS_RELEASE` |
| wall time / peak memory | 1 h 29 min / 6.1 GB | 1 h 18 min / 1.4 GB |

- The positive controls for these runs are g-sound's fresh replay of the unmodified #3945 (09:47 EDT, the same
  script on C1) and x-gate's fresh replay of the unmodified #2979 (14:17 EDT).
- The κ is unchanged, as it must be: only the placement of the roles changed, and the profile did not.
- As with the published packages, the replay still reuses the hub's retained unchanged ancestry (the PR1297 provider
  and parent receipts) under hash checks. That is part of the scripts' design.
- Scope [read]: these are the hub's finite checks. They are conditional on its inherited uniform / all-size, hosting,
  precision, analytic and BIT obligations. They are not a Lean proof, and they do not check that the hub's native
  model is what our Lean engine executes. That link is the port's job.

## 3. Step (3): does any rule of the hub forbid or price the forced placement?

### Why the hub's solvers leave the choice free [read, o-model MODEL.md sections 0, 4, 6, 7; checked against the code]

- **The price sees only the profile.**
  - The hub's b is the root of Σ_r n_r r^(1-b) = 896T·73^(1-b): `exact_root.py`, `check_allocation.py`'s moment
    auditor and my pricer re-implementation.
  - At a fixed inventory, both construction models minimize a cost that is linear in the counts of bin patterns per
    capacity. These are #3815's `hall_lp.py` (Σ g_a(l_j)·y[k, j]) and #2979's `fill17.py` (a fixed face mix,
    feasibility only).
  - An empty cap-48 bin costs the same in children 8/10 (pair class 4) and in children 12/13 (pair class 6). Neither
    model has a variable for "head", so the LP vertex or CP-SAT solution chooses the split arbitrarily.
  - The controls show this. A fresh CP-SAT solve of #2979's own model gives 14,526 heads instead of 16,086, at the same
    price (MODEL.md, A2). Re-running #3815's construction on x86 gives the same 57,312 heads but different rows (B1,
    B1v2).
- **Forcing the split costs nothing on the hub's side.**
  - Both designs keep the profile exactly: my recount, heads_fill, o-eval, the hub's own checkers and our lowering
    (G0) all agree.
  - So the hub's b, κ assembly and ceiling certificates are unchanged. The two certificates are #3815's
    `ceiling_certificate.py`, which passed inside #3815's `verify.py` on B's #3815-stage fill, and #2979's
    `check_cap.py`, which passed on A.

### Hub rules checked against the forced placement

| rule (where) | B | A |
|---|---|---|
| owner windows lo < t <= hi, one role per array per call, compact prefix, capacity, full histogram (#3945 `owner/verify.py`, PR3734 `owner_audit.py`, PR2839 engine in `check_allocation.py`) | pass | pass |
| line-specific: exactly 84,240 changed bins, each +0/+1, old idle exact; i.e. every width-4 role of the 26 widened helpers keeps a free cell (#3945 `owner/verify.py`) | pass | n/a |
| clock: PR3734's clock must reject the fill, #3815's must accept it (#3815 `verify.py`, `owner_audit.py`) | pass | n/a |
| block-norm alias premise: no role of call (s, r) in owner r's pair stage_pair[s] (#3815 `verify.py` step 3) | pass | no separate check in #2979's package; the complete admission, with its block-norm checker, passes (step 2) |
| price / ceiling: profile-only (`exact_root.py`, `ceiling_certificate.py`, `check_cap.py`) | unchanged, pass | unchanged, pass |
| exact eight-stream operator replay, the PR3187 / PR3296 model with which #4046 refuted the #3934-#3959 chain (gist 75dd2bb8's method) | **pass**: 360 owners, 2,029,945 moved characters, 0 differences (section 3.1) | **pass**: 74 owners, 594,132 moved characters, 0 differences |
| hub rule documents (AGENTS.md, CONTRIBUTING.md, autoresearch.md, README.md at main 20ecf289) | nothing about which bins may stay empty. A changed construction needs its complete cold replay: done, PASS (step 2) | same |

**No rule of the hub forbids the forced placement, and none prices it** [read, computed]. The hub does not know about
heads: the lever is ours. The hub's price is placement-blind.

### 3.1 The operator replay (C2 for B, C5 for A; our driver `colab/opgate.py`, hub modules imported unchanged)

- **Model.** PR3187's `algebra.py` and `source.py` (Nick Petalas) and PR3296's `obstruction.py` (Matarener), from a
  fresh export of pr3945.
- **Glue.** From ChaloupeMusic's gist: the (stage, replica) dispatch key and the owner-parametrised replay.
- **Characters.** For every owner, every axis whose hosting event differs between the published fill and the design.
  The published replay must be the scalar operator (baseline). The design replay is compared with the baseline and
  with the scalar.
- **Sample of owners.** Every axis of every idle child plus PR3187's 147 base modes, against the scalar.
- **Controls.**
  - Identity: 0 moved.
  - Offset +1: rejected.
  - Each slice that the design hosts in a bin the published fill leaves empty is tested twice: moved to its child's
    window lo (wrong time), and moved to the partner stream (wrong stream). Both must change the operator.
- **Results (19 Oct, 15:37 EDT).**

  | | B (vs #3945 published) | A (vs #2979 published) |
  |---|---|---|
  | owners | all 360 | all 74 |
  | moved characters | 2,029,945 | 594,132 |
  | published replay not scalar (baseline) | 0 | 0 |
  | design replay not scalar | 0 | 0 |
  | design differs from published | 0 | 0 |
  | all-axes sample (every axis of every child + 147 base modes) | owners 10, 334, 356, 0: 98,112 characters, 0 not scalar | owners 0, 37, 73: 73,584, 0 not scalar |
  | controls (owner 10 / owner 11) | identity 0; offset +1 rejected; 1,120 newly hosted slices in children 12/13 are scalar, and wrong time changes 1,120/1,120, wrong stream 1,120/1,120 | the same, 1,120/1,120 |

  Logs: `logs/op/` (B on C2, about 30 min with 6 processes at nice 19; A on C5, about 10 min with 3).
  - The moved characters are axes of about 290,000 bins in B and 120,000 bins in A whose role sets differ from the
    published fill. The solvers re-placed most roles, not only the [9^5] bins.
  - The published fill is scalar at every one of these characters as well. For #3945 that is wider than g-sound's
    earlier replay, which covered #3945's own 115,473 moved characters.

### Notes for the port (not soundness) [argued]

- **Design B hosts 116,640 roles in children 12/13; the published #3945 fill hosts none there.** These are the 23,328
  [9^5] bins of pair 0, with window (τ2, τ8].
  - The #3945-line generators (p-own, p-hostck, p-a3 schedule) have never seen hosting in these children. They must
    accept that window.
  - The #2979 line already hosts there (22,366 roles in the published #2979 fill), so the mechanism exists in our
    #2979 port.
  - In our lowering these bins are TAIL-literal blocks of remainder 3 (G0), paid at rank 3 as before.
- **On B, the "net heads" of lead 2 fall from 27,648 to 4,272.** The 23,328 rank-3 net heads became literal heads.
  This is only relevant if lead 2 is ever priced on top.
