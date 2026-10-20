# fc-gate: is the FULL cut on design B a valid, correctly priced lever?

Round wht45/fullcut, 19 Oct 2026. Agent fc-gate (independent of bB-data, who priced it). Nothing posted; no Lean on
the Mac; GitHub not touched. Labels: [computed] (command in section 6), [read] (file named), [argued].

**v2, 22:14 EDT (v1 22:09) — VERDICT: GO.** The full cut is a valid program transformation in our lowering, Theorem 1G (both the
job-model proof and its Lean form, e-r12's `Rule.twoTarget`) covers it with separation 24 < 25, and it prices at
**1.7806935e-3** (N = 17806935; N + 1 refused) by four routes (A-D below; B and C share no code with bB-data's price_ours.py,
D only the hub engine's log/exp enclosures, which price_ours.py used for the uniform price alone), 1.28 % above the half cut's 1.7581520e-3.
No STOP written. Refinements below are marked with their time when added.

## 1. Price [computed]

Profile of the full cut, re-derived from #3945's PROFILE (= design B's, `DS5H.rem_hist`, G0) and design B's head
classes (G0 `out/g0-3945H.json`, Lean `DS5H.heads_cnt = 80640`, `DS5H.f1_cnt = 57312`):
- unit paid profile {1: 8,041,680, 2: 2,306,880, 3: 729,648, 4: 353,160, 5: 69,120, 8: 25,920, 48: 137,952};
- the 80,640 literal heads (hub children 8/10, kind A4) are removed: their arrays get target 73 - 48 = **25**;
- the other 57,312 rank-48 blocks (hub children 12/13, kind F1, remainder 48) are each split 24 + 24;
- `cut3945F` = {1: 8,041,680, 2: 2,306,880, 3: 729,648, 4: 353,160, 5: 69,120, 8: 25,920, **24: 114,624**};
  targets {73: 241,920, 25: 80,640}, W = 322,560, H = 80,640. Rank mass 19,560,960 = 23,431,680 - 48 x 80,640;
  11,641,032 paid blocks. Max paid rank 24 < 25: separated. Profile file: `out/profile3945F.json` (the generator's
  JSON format; identical numbers to fc-rate's `tools/profile3945F.json`).

| route | full cut b | N (largest grid numerator) | N + 1 |
|---|---|---|---|
| A. lever-reprice `lever_price.py` + `mixed.py` (Decimal prec 60, bisection) | 1.780693540576e-3 | - | - |
| B. my own exact rational interval arithmetic (`tools/route2.py`: atanh logs, Taylor exps, outward 2^-330 rounding; no shared code) | [1.780693540576287216656517872997541608626131335188768004413535325e-3, ...5398e-3] (64 digits) | **17806935** proved at z1 = 1 - (100N+2)/10^12 | **refused** (certified false at z1(N+1) and at 1 - (N+1)/10^10) |
| C. the Lean data generator e-r45m1 `mk_rate_mixed.py scan` (x-pack's enclosures, what Lean re-checks) | root 1.78069354057628721665651787299754160862613133518876800441354028e-3 | **17806935** (slack 8.9e-11 relative; root - N = 0.406 grid steps) | refused (rc 2) |
| D. (v2) the hub moment engine's own enclosures (release `tools/lever/pricer.py`: `log_bounds`, `exp_bounds`, 224-bit outward rounding), rate fact divided by W in the engine's normal form `Σ (r n_r)/(mW) (m/r)^b ≤ (W-H)/W + (H kc)/(mW) (m/kc)^b` (`tools/route_moment.py`) | 1e-24 grid: **0.001780693540576287216656** (certified at b and b + 1e-24, as the hub's `certify`) | **17806935** | false at z1(N+1) and at (N+1)/10^10 |

Controls (all reproduce):
- uniform #3945: route B gives 1.730663236534549309146766...e-3, the hub's stated root [0.001730663236534549309146, ...147] digit for digit;
- design B half cut: 1.7581520168e-3 (A), 1.758152016764631555...e-3 (B), 0.001758152016764631555031 (D), N = 17581520, N + 1 refused (B, D) = bB-data's Lean figure;
- #3945's published fill (57,312 heads, 80,640 others) full cut: 1.7307679642e-3 (A, B, D; N = 17307679) = 1.7307680e-3;
  break-even 57,262 heads (A) = lever-reprice's;
- uniform #3945 by the hub engine `price_b` (D's code): 0.001730663236534549309146, the hub's root;
- design B's best separated d (lever_price's sweep d = 1..48) is d = 48, the full cut itself.

Full cut vs half cut: +2.2541e-5 (+1.28 %); vs #3945's b: +5.0030e-5 (+2.89 %).

## 2. Validity in OUR lowering [read / argued]

**What the 57,312 blocks are.** G0 on design B (`../../allheads/bB-data/out/g0-3945H.json`, classes): hub children 12
and 13 = our kind **F1** (climb step `fin`, the final sweep, comps 1/2 of our network = hub pair 0), bin rank 48,
**remainder 48** (28,672 in child 12, 28,640 in child 13). Remainder 48 = the whole bin is idle:
no helper role is hosted in these bins (the F1 bins that host roles have remainder 3: 23,328 of them). They are not
heads (stage0/stage1 locals and gateA2 touch the array first) and not tails (gateC2 reads/writes after): class
"neither". In the price they are ordinary paid blocks of rank 48 (`sched3945H_hcost`: 137,952 φ 48 per gauge copy,
of which 80,640 are the heads).

**The split is an exact identity, independent of what the blocks are.** `Work.Mixed.HeadM1` (compiled, same file in
the leverB package, md5 86343055):
- `split_block : blockMat z = blockMat (first a) * blockMat (last b)` for any block of rank a + b;
- `splitWord a b w` replaces EVERY block letter of rank a + b by two block letters on the SAME role (first a, last b
  directions); `gwalk_splitWord : gwalk f (splitWord a b w) x = gwalk f w x` at every column count f, every input;
  `GLiveKernelM.splitWord` keeps every mixed certificate; `hist_splitWord` (1 ≤ a, 1 ≤ b) moves rank a + b to a and b;
  `proper_splitWord` (1 ≤ a, b < u). **No hypotheses on the block** (hosted or not, head or not, phase or frame).
So no data condition is needed for the split. Which blocks get split is forced by the histogram: after the heads are
removed, the body's rank-48 letters number exactly G x (137,952 - 80,640) = G x 57,312 (from the port's price
identity `hcost`, as `hist_headsFirst` reads it); every other paid rank is ≤ 8, 3 or 5. So `splitWord 24 24` splits
exactly the 57,312 F1 bins per gauge copy and nothing else.

**The heads: `cut_heads_full`** (`HeadM1` line 503): hypotheses `he` (e injective), `hL` (L nodup),
`hz : ∀ l ∈ L, ∀ x, bph (z l) x = bph ((A l).fam (S l)) x` (the WHOLE head, not just its first 24 directions) and
`hk : GLiveKernel e (headWord e L rk z ind ++ body)` (rk any function). Output:
`GLiveKernelM e (fun l => blockMat (cutTgt L A S l).2) (body ++ shiftWord e L (cutShift A S))`, target rank
`card α - |S l|` (`cutTgt_rank`).
Do the data hold them? Yes, with NO new port data:
- `hk` is exactly what the half cut already uses: p-hfirst's `xnetGeom_heads_first` (heads of rank 48 first, on
  `headsL (headBins S)`), unchanged.
- `hz`: p-hfirst's head directions are `zGeom r = a4Dir base l hH g (s,t) = fun i => (a4Base base l g t).fam (a4Lines l s)
  (Fin.cast (a4Lines_card hH l s).symm i)` (GeomA4.lean line 67). So with `A := cutA base l` (unchanged) and the new
  `S r := a4Lines l s` (the head's 48 lines), `hz` is `HeadM1.bph_cast` (compiled, line 115). Rank:
  `card_LC401 = 73`, `a4Lines_card = 48`, so target 25.
- (In the half cut, `hz` was on `partA z l` with 24 lines `cutLines`; the full cut is simpler.)

## 3. Theorem 1G for this job model [read / argued]

PROOF.md's Theorem 1G: targets W_k at ranks K = {73, 25} (W_73 = 241,920 G, W_25 = 80,640 G), paid rounds in word
order with ranks r_j ∈ {1, 2, 3, 4, 5, 8, 24}, r_max = 24, k_min = 25.
- (S) separation r_max < k_min: **24 < 25** holds.
- Each unit must be a correct word at every column count f: `GLiveKernelM` is quantified over every f (cut_heads_full
  and splitWord are both "at every column count").
- The heads must be first on their own array, alone: same as the half cut (p-hfirst / G0: 80,640 literal heads, one
  per array, unhosted, first touch). The split blocks need no such property: splitting a job a + b into jobs a and b
  is allowed in the job model (PROOF.md section 0 preamble), and in Lean it is `gwalk_splitWord`.
- theta = ln 24 / ln 25 = 0.98732 < 1 (Lemma D's exits still shrink); eta_73 = 0.99063 at tau = 1 - 1.7806935e-3;
  kappa = ceil(ln(2/eta)/ln(1/theta)) = 56 [computed], finite. For comparison the Lean-checked half cut has
  theta = ln 48 / ln 49 = 0.99470 and kappa = 133, so the full cut is the better-conditioned of the two: the log factor only changes the constant's exponent of ln n, absorbed by any exponent above
  tau (the statement's exponent is above alpha, which is above z1).
- The Lean form (e-r12 `Rule`, e-r3cost `Eng.bridge_twoTarget`, `Eng.engine_seam_final`): its hypotheses are
  `rk_pos`, `rk_le : rk j ≤ rmax`, `sep : rmax < kmin`, `two_le_kmin`, `kmin ≤ kr ≤ m`, `kmin_att` (H ≥ 1) and the
  rate fact `RateQ z1 1` with 0 < z1 < z ≤ 1. Nothing quantitative in rmax/kmin (e-r12 uses
  `rate_mono` with q = (rmax/kmin)^(z - z1) < 1; for 24/25 and z - z1 = 1e-12, q = 1 - 4.1e-14 < 1, only the
  constant grows). `Rule.twoTarget 73 25 W H P 24 rk ..` meets all of them: 24 < 25, 2 ≤ 25, 25 ≤ 73, H = G x 80,640 ≥ 1.
- Eng's own fields (Core.lean `cutEng` pattern): `sep : every paid rank < tr l for all l` (24 < 25 ≤ 73), `tr_pos`,
  `tr_le`, `three : 3 ≤ card α`, `hP` proper (`proper_splitWord 73 24 24` + body proper + shifts proper).
- The hub engine's halving constraint m^2 > 2 r_max^2: 5329 > 1152.

## 4. The Lean plan (for fc-rate and fc-engine)

Everything the port proves is reused; only the bridge between the port certificate and the mixed engine changes.

**Unchanged, reuse as built (bB-data's / bB-final's trees):** ENG (385 files: engine R1-R4, warps, words, dispatch,
network, generic bridge Core/Inst/HeadsFirst/HFTest, p-hfirst HF, HeadM1 with `cut_heads_full`, `splitWord`),
bB-data's 176 files (S5H data, `DS5H.*`, `Inst.S5H` / `sched3945H` / `S5HFacts` / `S5HHB` / DispLink,
`Work.L3945H.Prov`'s fill facts `hfit3945W`, `hpos3945W`, `hU3945W` with `PB.unit3945H` — these do not mention the cut).

**New (suffix 3945F / L3945F):**
1. R5 (fc-rate, running): `Work.Mixed.Rate3945F` (+`Sharp`, `R5Joint3945F`) from `mk_rate_mixed.py emit
   3945F=profile 17806935`: `rate3945F_moment : moment 73 z1 cut3945F ≤ 241920 + T 73 z1 25 80640` at
   z1 = 1 - 1780693502/10^12; `Rule.twoTarget 73 25 322560 80640 P 24 rk`.
2. Bridge, generic (fc-engine), e.g. `Work.Mixed.Bridge.Full`: a full-cut analogue of HFNet with `z : σ → Fin c → Space α`,
   `S`, `hz : ∀ l ∈ L, bph (z l) = bph ((A l).fam (S l))`, `hS : (S l).card = c`, body, hP, hk, G, U, hcost;
   `fullWord := splitWord 24 24 (body ++ shiftWord e L (cutShift A S))` (shifts have no blocks, `hist_shiftWord`;
   splitting before or after appending the shifts is the same);
   `fullEng` = `cutEng` with `wd := fullWord`, `hw := GLiveKernelM.splitWord 24 24 (cut_heads_full he L hL (fun _ => c) z ind A S hz body hk)`,
   `hP := proper_splitWord 73 24 24 ..`; histogram: `hist_splitWord` + `hist_headsFirst` give
   `fullWord.hist r = G * histFn cut r` from a profile identity such as
   `hcutF : ∀ r, histFn cut r = (if r = 48 then 0 else histFn U r) + (if r = 24 then 2 * (histFn U 48 - Hd) else 0)`
   with `Hd ≤ histFn U 48` (or prove `fullWord.hist` per key by `decide` like `hcut_of_check`).
   Then **reuse `bundle_of_eng` / `bridge_of_eng` unchanged** with a = 48 (`ha`, `haα : 48 < 73`,
   `hsep : 24 < 73 - 48`, `h2 : 2 ≤ 25`), and `E.engine_seam_final` as in `cut_engine`.
3. Port helper (tiny): `fullS l hH r := a4Lines l s` on head roles and `bph_full : bph (zGeom base l hH r) x =
   bph ((cutA base l r).fam (fullS l hH r)) x` by `bph_cast`; `fullS_card = 48` by `a4Lines_card`.
4. Instance: `Inst3945F` (`hcut3945F` for `unit3945H` and `cut3945F`, negatives: one head more, #3945's 57,312),
   `HeadsFirst3945F` / `HFTest3945FW.engine3945F_geomW_unit` = copies of the 3945H files with `cut3945F`,
   kc 25, rmax 24, envelope `1 - 1780693501/10^12`, the full-cut bridge; same hypothesis list (I, Ic, base, S, hfit,
   hpos, hheads = 80,640, hU).
5. Statement stage L3945F: `mk_l3945.py final out/final-L3945F --suffix L3945F --num 17806935 --profile 3945F=..`
   (alpha `1 - 1780693501/10^12`, target `1 - 17806935/10^10`, negative `1 - 17806936/10^10`); `Work.L3945F.EngineOfGeomW`
   and `Work.L3945F.Engine` = copies of L3945H's calling the 3945F engine; `Prov` can be imported from L3945H as is.
6. Integration as bB-final: Compare (pos PASS / neg FAIL at 17806936), official comparator (~1 h), merged bigospeedrun
   (~2.5 h, start by 03:30 EDT), axioms, greps, challenge derivation.

## 5. Refinements (v2, 22:14 EDT)

- **Recount from the raw design-B FILL** (`tools/recount_full.py`, own code, no chronology; FILL sha256 9c06677d,
  1,321,920 rows) [computed, `out/recount_full.json`]: every bin of children 8 and 10 is fully idle (40,320 + 40,320
  = 360 x 2 x 112: the heads); children 12 / 13: 28,672 / 28,640 fully idle bins and 11,648 / 11,680 bins hosting
  exactly 5 roles each (5 x 9 = 45 of 48: the remainder-3 bins). Idle 48s = 80,640 + 57,312 = 137,952 = PROFILE[48].
  The derived full-cut children equal `out/profile3945F.json` exactly. Agrees with G0 and with the Lean facts
  `DS5H.heads_cnt = 80640`, `DS5H.f1_cnt = 57312`, `sched3945H_hcost` (137,952 φ 48) [read, bB-data logs/C3/build].
- **Lean shapes checked against the compiled sources** [read]: `Eng` (Work/Mixed/Prog/Desc.lean) needs only
  `tr_le`, `tr_pos`, `hP`, `sep : ∀ paid r, ∀ l, r < tr l`, `W_pos`, `three`, `hw` (= `GLiveKernelM` unfolded);
  `hist_headsFirst` gives `body.hist r + [r = 48] |L| = G histFn U r`; p-bridge's `bundle_of_eng` / `bridge_of_eng`
  are generic in the cut size `a` (here 48) and need `1 ≤ a < card α`, `rmax < card α - a`, `2 ≤ card α - a`;
  `Eng.bridge_twoTarget` and `Rule.twoTarget` need `rmax < kc`, `kc ≤ mx`, `2 ≤ kc`, `1 ≤ H`. All hold for 48 / 24 / 25.
- **Gate review of fc-engine's drafts** (`../fc-engine/lean/Work/Mixed/Bridge/CoreF.lean`, `Inst3945F.lean`,
  `HFTest3945FW.lean`, `Work/L3945F/Engine.lean`, as of 22:12): they implement exactly section 4 — `fullWord =
  splitWord a b (body ++ shiftWord ..)`, `hw := GLiveKernelM.splitWord a b (cut_heads_full ..)`, `hcutF` in the
  form `[r ≠ 48] U r + [r = 24] (U 48 - Hd) + [r = 24] (U 48 - Hd)` (gives 114,624 at 24, 0 at 48: correct, since
  a = b both fire), `bph_cutF` by `bph_cast` with `cutSF = a4Lines`, `rate3945F_moment` at kc 25, and the L3945F engine
  reusing `L3945H.Prov`'s fill facts unchanged. I found nothing wrong in them. (fc-engine's STATUS reports CoreF compiled on
  C3 at 22:20.)
- Head-rank vs split-rank: CoreF's `hcost` subtracts the heads at rank `a + b`, i.e. it assumes the heads have the
  same rank as the split blocks (48 = 24 + 24). True for design B; noted because the generic theorem is not more
  general than that.

## 5b. What would change the verdict (none observed)

- A paid rank above 24 in the full-cut word: impossible, the histogram is forced (`fullWord_hist_eq`), and the unit
  profile's only rank above 8 is 48.
- A head that is not first/alone on its array: that is the half cut's hypothesis too (p-hfirst `xnetGeom_heads_first`,
  Lean-checked in L3945H); the full cut adds nothing to it.
- A Lean-side quantitative condition on rmax/kmin: none in e-r12/e-r3cost (only `rmax < kmin`, `2 ≤ kmin`).

## 6. Commands

```
cd <work>/wht45/fullcut/fc-gate
python3 -I -B tools/route1.py > out/route1.json                       # route A + controls (< 1 s)
python3 -I -B tools/route2.py out/route2.json                         # route B, 5.5 min, < 100 MB
python3 -I -B <work>/wht45/allheads/bB-data/tools/mk_rate_mixed.py scan 3945F=out/profile3945F.json   # route C
python3 -I -B tools/route_moment.py out/route_moment.json             # route D, 13 s
heavy.sh python3 -I -B tools/recount_full.py <work>/wht45/allheads/bB-data/in/owner/FILL.json.gz \
    <work>/wht45/allheads/bB-data/in/owner/PROFILE.json 360 > out/recount_full.json   # FILL recount
```
