> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# an-build: two-sided slack tried exactly — it does not survive

Anatomy round, build track, 15 October 2026, 09:22 to 11:35 EDT (Mac clock). Task: test an-synth's top proposal
("two-sided slack": a second unpaid coordinate per family, predicted ×2.03), after first reproducing the source unit's
saving with our own code. All computations are our own code on our own data files, on the Mac (each under 10 s, under
200 MB). Nothing was posted. No Colab was used. No third-party code was run. No skeptic has read this note.

Labels: [computed] = a program in `calc/` and its output; [read] = taken from a named file; [argued] = my reasoning.

## In short

1. **Calibration passes** [computed].
   - The replica frame model built from our Lean layout gives D = 320 per replica (every data array pays 72 of 73)
     and the exact-fill price 1.6647327e-3.
   - Our hosting tool reproduces its own same-run run: W 7544, 376 banks, 1.6225278e-3.
2. **Two-sided slack is well defined, and it is forced to have one shape.**
   - The 112 ports span F₂⁹, and the scatter's reader graph is connected.
   - So the Y line must be one extra coordinate s, shared by every Y array, at m = 74. The X arrays never hold s.
   - The gain is real on paper: each data array saves 2 of 74 instead of 1 of 73, i.e. +896 moves per replica.
     [computed, argued]
3. **But the network must also move Y data into X.** Our nine stages alternate the forward local word (y += x) with
   the backward one (x −= y). Both run on one invocation with one frame word [read: our Lean].
   - In a backward call, every channel of y information has to leave s.
   - The y-to-x map has rank 112, so at least 112 crossings are needed per backward call, at ≥ 1 move each.
     [computed, argued]
4. **On our word (U8X) at its fixed gate order, the crossings eat the whole gain** [computed, exact]. The Lean word has
   two "dirt" gates, one at the start frames and one at the end frames [read: our `InvThm.lean`].
   - Both need s-crossings: rank 112 at the start and rank 112 at the end.
   - That is 224 per backward call, and 4 backward calls per replica make 896. That equals the gain exactly.
   - Net: D stays 320 per replica, now at m = 74. Priced with our hosting ILP, the price falls to 1.29e-3 to 1.33e-3,
     ×0.79 to ×0.82 of the one-sided price.
   - **Verdict: kill.**
5. **The other lowerings are illegal outright.** On U8X0 (the word of our Lean-checked line) and on both of our T74
   lowerings (the hub's record word), "mutable sources" put 113 to 124 direct x–y gates in phase B. Under two-sided
   frames x avoids s and y holds s, so these gates can never share a frame [computed].
6. **Ceiling for any redesigned local word in this template.**
   - D ≤ 768 per replica (×2.4), by the rank bound in item 3.
   - That would price at about 3.1e-3 to 3.2e-3 with our hosting model.
   - Reaching it needs a new local word that routes all of a backward call's y information through a single
     112-channel crossing, and whose forward calls pay nothing for theirs.
   - No such word is known. This is an upper bound, not a construction [computed price; argued bound].
7. **Proposal 2 (a port design beyond E8): BW16 fails the screen** [computed].
   - Junk per port is 1,304 (mark ≤ 120).
   - An explicit junk clique of 30 lines forces a window of at least 30 (mark ≤ 12).
   - The Leech family fails the same way, through the clique {4e₁ + 4e_j} of 23 lines. Coxeter–Todd was not screened.

## 1. Calibration

- `calc/calib.py` builds the frame flag of one family from our Lean layout and checks it.
  - Sources: `HDef.lean` gluing x01..y78, host chains `hH0_*` (X side on the line, Y side at 0), and the
    idle-children table of `host_price.py`.
  - Result: every pair's X and Y pay 72 of 73 (in-call plus idle children equal to the table), and D = 320 per
    replica.
- Pricing uses an-ours' fold, the Lean fold hypothesis with s = 40:
  - U8 exact fill: 1.6647327e-3;
  - the U8X0HQ profile: 1.6256631e-3 (Lean grid point 1.6256629e-3).
- `calc/hosty.py` imports our same-run hosting tool (`fsp-price/tools/host_price.py`) unchanged. It reproduces that
  tool's recorded run exactly: W 7544, 376 banks, root 1.62252782e-3; our fold gives 1.6225278e-3.

## 2. What two-sided slack has to be

- In our complex engine the X line is free: x starts on its port line, and Y ends on its hyperplane.
- Two-sided slack adds a free line for Y: y starts on a line ⟨s⟩, and x ends in s^⊥.
- `calc/twosided_local.py` shows that s is forced to be a single new coordinate shared by all y's:
  - The ports span all of F₂⁹, so s cannot lie in the window.
  - Each total serves 56 y's at the scatter, its descent must pass through their frames, and the reader graph is
    connected. So all y's share one s.
- Then m = 74. The global chain is the old one with ⟨s⟩ added on the Y side, so every gluing equation still holds.
- Identity: D per word = 2·7168 − 72·64 − (crossings) = 9728 − (crossings), against 2560 now. [computed, argued]

## 3. Where it breaks: the backward stages

I first missed this, and my first price (5.2e-3) is withdrawn.

- **The stages** [read]:
  - The template's scalar chain applies GF, GB, GF, … (`Gates.lean`: `liftQ0 GF`, `liftQ1 GB`, …).
  - GF is y += x and GB is x −= y (`Reframe/SNet.lean`).
  - Both come from one invocation with the same frames (`Bridge/Inv.lean`, `istage_fwd` and `istage_bwd`).
- **The backward word** [read]:
  - It is the transposed inverse of each gate, applied in the same time order (`GCert/Chain/Phased.lean`).
  - It also has two dirt gates (`Carrier/Live.lean`):
    - g1 = 1 − PY·M·PS ("y −= K s") at the start frames;
    - g8 = PX + PY + PS·Mi (the slots set from the inverse) at the end frames.
  - These are real joins: `InvThm.lean` builds them as `GRoute.gate` at `start_frames` and `end_old`.
- **The data** [computed]:
  - `calc/dirtjoin.py` replays our U8X word exactly in fractions. The Lean facts hx, hy and hid hold on the replay.
  - The start dirt gate joins 328 slots with the y's. 112 of those slots meet an x in phase A.
  - With that gate, the backward reading delivers y information to all 112 x's; without it, to none
    (`calc/bwdpaths.py`).
  - `calc/endstep.py`: only 82 of 616 slots are restored by phases A and B, so the end gate joins the 490
    x-dependent slots with the x's. The y gates of phase B alone force 354 of them to hold s.
  - `calc/rank0.py`: the two maps that must cross s have rank exactly 112 and 112.
- **The cost** [computed]:
  - A crossing costs at least one move: an erased copy stepping ⟨s⟩ → 0. Any other register pays 2.
  - At this gate order: 224 per backward call, and 4 backward calls per replica (stages 1, 3, 5, 7) make 896.
  - The gain is +896. So the net gain is zero, at a larger block.
  - `calc/hosty.py 100 224` prices it at 1.290e-3 (one switch per bank) or 1.332e-3 (one per bank and run), against
    1.6225e-3 one-sided at the same hosting.
  - This price is optimistic: it takes the forward calls' upward crossings as free.
- Before this miss was found: `calc/mintrans.py` and `calc/mintrans2.py` (exact LPs, integral) showed the frames
  legal with 112 crossing copies per backward call. `calc/price_twosided2.py` and `hosty.py` then priced D = 768 at
  3.1e-3 to 3.4e-3. Those figures stand only as the ceiling of item 6.

## 4. Why it cancels so exactly [argued]

- Every Y family's data must end in an X array: ρ swaps them.
- Its line saves 2. Its channel must also leave s once, and in this template it leaves twice: once through the start
  dirt gate and once through the end dirt gate.
- So the second line returns exactly what the two dirt gates take: 2 × 112 per backward call. This is the deficit
  identity again (paid = W·m − 2·Σ source dims + 2·dropped dims), with the drops now forced on the Y side.
- A redesign that crosses once would keep half: D ≤ 768.

## 5. Predictions (STATUS.txt, written 09:31) against outcomes

| prediction | outcome |
|---|---|
| B1: the calibration reproduces D = 320 and 1.6647327e-3 (80 %) | **held** (and the hosting tool's 7544 / 376 / 1.6225e-3) |
| B2: not legal at zero extra; the exchange forces a descent (65 %) | **held**: the backward stages force 224 crossings per call at this order |
| B3: extra ≥ 64 per call; kill 45 %, between 30 %, pass 25 % | **kill**: 896 per replica, 99.6 per call on average, 224 per backward call; net gain exactly 0 |
| B4: nothing above 6.5e-3, nothing near 1e-1 (97 %) | **held**: the ceiling for any redesign is about 3.2e-3 |

## 6. What this leaves

- **For the hub's coverage-law card** ("explain the single unpaid coordinate per family; stop if the endpoint
  contract forces it"): the endpoint contract does not forbid a second line. The exchange does.
  - The Y data must come back across the new coordinate.
  - In this template it crosses twice per backward call (the two dirt gates), and that cancels the gain exactly.
  - The single unpaid coordinate per family is what survives a backward stage for free.
- **Open, as an upper bound only:** a local word whose backward stage sends all y information through one
  112-channel barrier would give D ≤ 768 per replica, about 3.2e-3. The dirt gates of every word we have stand in the
  way.
- **Proposal 3** (a small-block shared-window network) was not tried. The backward-stage cost found here applies to
  it too: each family's second line is worth 2 and costs at least 1 to cross back. So u_eff ≤ 1.5 there, which is
  well short of what the 1e-1 bar needs [argued].

## Files

All in `checks/wht43/anatomy/an-build/`:

- `calc/calib.py` → `calib.out`, `calib.json` (calibration).
- `calc/twosided_local.py` → `twosided_local.out` and `.json` (s is forced; the first-pass frames).
- `calc/hostclass.py` → `hostclass.out`; `calc/viol_probe.py` (direct x–y gates of U8X0 and T74).
- `calc/mintrans.py`, `calc/mintrans2.py` → `mintrans_U8X*.json`, `mintrans2_U8X.json` (exact LPs).
- `calc/bwdpaths.py`, `calc/dirtjoin.py` → `dirtjoin_U8X-parsed.json`; `calc/endstep.py`; `calc/rank0.py` →
  `rank0.json` (backward paths, dirt gates, exact ranks).
- `calc/price_twosided.py` (withdrawn first price), `calc/price_twosided2.py` → `price_twosided2.out` (ceiling,
  idealised).
- `calc/hosty.py` → `hosty.out`, `hosty.json` (ceiling at same-run hosting); `hosty_fixedorder.out`,
  `hosty_fixedorder.json` (the verdict price).
- `calc/bw16.py` → `bw16.out`, `bw16.json` (proposal 2 screen).
- Sources read:
  - our Lean in this repository: `community/u8xmhh/overlay/Work/XHost/HDef.lean`,
    `HStage.lean`; `Work/GCert/Chain/{Win,NetDef,Phased,InvThm}.lean`; `Work/Reframe/{SNet,Shared}.lean`;
    `Work/Bridge/Inv.lean`; `Work/Carrier/Live.lean`;
  - the archived template `Gates.lean` (`wht43/bound-check/h-L3-design/gen/template/` in the private archive of
    the agents' working notes);
  - pr717 README (P1, P4, P5, S1 to S4, F2; read as text);
  - the hub's `03-motifs.tex` and `04-swap.tex` (read as text);
  - our lowerings `k-lower/out/U8X-parsed.json.gz`, `fsp-cc/out/U8X0.gcert.json.gz`,
    `lg-a1/out/T74E.gcert.json.gz` and `T74G.gcert.json.gz`.
