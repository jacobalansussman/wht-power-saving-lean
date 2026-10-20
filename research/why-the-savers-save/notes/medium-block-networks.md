> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# cx-medium: a crossing unit at a medium block (13 to 25), built and priced

Anatomy round 3, track cx-medium, 15 October 2026, 12:25 to about 13:35 EDT (Mac clock).

**Task** (BRIEF3.md):
- build a crossing-slack network at blocks 13 to 25, with E7 as the port design and P112, D6 and A5 for comparison;
- start from the hub's template, shrunk to the smaller block;
- price it for real, with junk cleaning paid;
- say whether anything beats the hub's 1.68e-3 or reaches 1e-2;
- check any candidate exactly before pricing.

**Labels:**
- **[computed]**: a program of ours in `calc/` and its output.
- **[proved]**: a complete short argument written here. It has been checked by me only; no skeptic has read this note.
- **[read]**: taken from a named file.
- **[argued]**: reasoning that no program has checked.

All programs are our own. Small runs were on the Mac. The integer programmes, the E7 call and the large subspace
enumerations ran on Colab machines through a private remote-run helper script (not included). No third-party code was run, and nothing was posted.

## In short

1. **Nothing in blocks 13 to 25 beats the hub, and nothing built comes near 1e-2.**
   - The best real price in the range is **7.39e-4**, for the 112-port design at block 25 with the hub's own record
     word (T74 gauge44). With U8X0 (the word of our Lean-checked line) it is **7.31e-4**. Both are ×0.44 of the hub's
     1.6797550e-3.
   - Both figures are upper bounds, because the banks are packed ideally. First-fit packing gives 6.88e-4 to 7.09e-4
     [computed].
2. **Shrunk to blocks 13 to 25, the hub's template has five points, and only one of them saves.**
   - With P bank pairs the template has 2P + 1 stages and block m = 1 + σ(2P + 1), where σ = k − 1 for a window of
     k bits.
   - The points in range: E7 at 19 (one pair); D6 and A5 at 13 (one pair) and at 21 (two pairs); the 112-port design
     and E8 at 25 (one pair).
   - Per replica, D = 2Pv − (2P + 1)c, where v is the number of ports and c the centre charge per call. With one pair,
     the centre must cost less than 2/3 of a level per port per call.
3. **New floors on the centre charge** [proved; computed by exhaustive enumeration].
   - They hold for *every* decoder whose totals gather, descend once to frame 0 and scatter (the hub's accounting).
   - **D5, D6, A4, A5: at least 0.8 per port per call.** That is exactly the two-pair break-even, so these four designs
     can never save with one or two pairs (blocks 9 to 21). This closes an item that round 2 left open.
   - **E7: at least 35 per call.**
     - Every decoder of the kind used so far (totals Σ_b (f·r_b) x_b) costs at least 42, the one-pair break-even.
     - The one known decoder costs 49, so **D = −21** at block 19.
     - The only 7-total decoders near the known one are rescalings of it [computed: the tangent space has exactly the
       dimension of the obvious deformations].
     - A different decoder costing 35 to 41 is not excluded.
   - **112-port design: at least 52.4 per call** against the hub's 64. E8: at least 54.5.
4. **What sinks the one design that saves is junk cleaning without idle levels.**
   - With one pair, every data array is busy in all three calls, so there is no idle climb to host a helper on.
   - At block 25, U8X0's helpers need 453 banks per 224 data arrays. They carry **69 % of the narrow mass**
     [computed].
   - With junk cleaning free, the same network would price at **2.33e-3**, ×1.39 over the hub [computed].
   - To beat the hub, the helpers' narrow mass must fall to **17.9 % of U8X0's** (5.6 times less), with no idle level
     anywhere [computed requirement; argued that no such circuit exists].
   - When side information travels only on copies that never climb, free meetings (a source and a target at the same
     frame) can clean at most 57.7 % of the 112-port junk and 43.8 % of E7's [proved]. Data arrays that carry other
     sources' content on their own climbs (the hub's mutable sources) are not covered by that bound.
5. **1e-2.**
   - Excluded at block 25 with the hub's centre: even with junk cleaning free and the best piece sizes the stage gates
     allow, the ceiling is **8.26e-3** [computed; the piece bound is argued].
   - Not excluded outright: a 112-port centre of at most 61 per call (the proved floor is 53) together with nearly free
     junk cleaning would pass 1e-2. On E7, a decoder of at most 36 per call with free junk cleaning would pass it.
   - Neither ingredient is known, and the best word known still needs 504 helpers per call.
6. **Shrinking the hub's template only loses.** With the U8 word and generous hosting:
   - block 25: 7.31e-4; block 41: 1.467e-3; block 57: 1.652e-3;
   - block 73: 1.6647327e-3, the calibration, exact; block 89: 1.649e-3 [computed].

   Fewer stages leave fewer idle climbs to host on, and that loss outweighs the smaller block.
7. **Exact checks** [computed].
   - The candidate's local word (U8X0): our own checker was rerun, and every exact check passed. The regenerated
     certificate is identical in content to the stored one.
   - New crossing calls, built here and checked on raw complex arrays with real frame operators: **A5 exact at
     block 5**, **E7 exact at block 7**.
   - Not checked on raw arrays: a whole three-stage run with its gauges between stages.

## Words used here

- **Template, stage, call.**
  - The hub's network runs a local word in stages. In each stage one bank pair is active (a *call*), and the other
    pairs climb idle.
  - Calls alternate between GF (y += x) and GB (x −= y).
  - With one pair, the calls GF, GB, GF swap the pair: X ends with −y and Y with x.
- **Window, σ.** The local frames of a call live in a k-bit window. A data array climbs σ = k − 1 levels in each stage,
  active or idle. So m = 1 + σ × (number of stages).
- **Centre, total, decoder, break-even.**
  - A *total* is a sum of sources. Its copy, a fresh array, descends from the total's frame to frame 0. It is added
    into the targets there and then thrown away. It costs the dimension of the frame it starts from.
  - The *decoder* is the choice of totals.
  - One pair gains 2v free directions per run against 3c of centre, so c must be below 2v/3.
- **Exact pairs.**
  - A pair of ports must be delivered exactly when their labels are non-orthogonal (orthogonal roots).
  - Otherwise it may carry ±1/2 junk.
  - For E7 the 63 ports are the non-zero vectors of a 6-dimensional symplectic space over F₂ (label u = 1 + e).
    Exact pairs are the commuting ones.
- **Junk, free meeting, carrier, bank.**
  - A target's junk term from source b can be removed for free if the target, on its way up, passes through a frame at
    which a copy of x_b was taken. This is a *free meeting*.
  - Otherwise some array must carry x_b up to the target's frame. It is either a *carrier*, a fresh array whose climbs
    are paid and lost, or a helper hosted on a *bank*, an extra array that pays for its own transform and saves
    nothing.
- **Narrow mass N.** Σ r ln(m/r) over all paid pieces. To first order, the saving is D/N.

## 1. The hub's template, shrunk

**What changes** [read: the template `Gates.lean`, an-build `calib.py`, an-ours `stages.py`; argued]:
- *Stages.* A swap needs three calls on some pair (an2-smallnet's bank lemma). The hub gets the other pairs through
  with two calls each, by free twin gates, so P pairs take 2P + 1 stages. One pair is the smallest template: three
  stages.
- *Windows.* Each stage costs every data array σ fresh levels, so m = 1 + 3σ with one pair. This is also
  an2-smallnet's argued floor for a swap: 19 for E7, 25 for the 112-port design and E8, 13 for D6 and A5.
- *Bridges.*
  - Between stages, each family's frames are glued by a free address map, as in the hub (Φ_k(1) = Φ_{k+1}(t) for X,
    Φ_k(t⊥) = Φ_{k+1}(0) for Y).
  - X climbs 1 → 1 + σ → 1 + 2σ → 1 + 3σ, and Y climbs 0 → σ → 2σ → 3σ.
  - With one pair there are no twin gates.
- *Idle levels.* None. Every data array is active in all three calls, so all its levels are spent inside calls. This
  is the change that matters.

E7 at block 13 and the 112-port design at block 17 (an2-smallnet's proved floor 1 + 2σ) are not template points. A swap
there would need a backward call and the forward call after it to share a window, which an2-smallnet proved
impossible. Points between template blocks only add idle coordinates. That could help hosting, but it cannot help a
design whose D is not positive.

**The template points in 13 to 25, before any junk cleaning** [computed, `calc/items.out`, `calc/cliquefloor.json`,
`calc/cliquefloor9_*.json`]:

| design | pairs | block | data arrays | best centre known (per call) | proved floor (per call) | break-even | D per replica |
|---|---:|---:|---:|---:|---:|---:|---:|
| E7 | 1 | 19 | 126 | 49 | 35 | 42 | **−21** |
| D6 | 1 | 13 | 60 | 30 | 24 | 20 | **−30** (≤ −12 at the floor) |
| A5 | 1 | 13 | 30 | 30 | 12 | 10 | **−60** (≤ −6 at the floor) |
| D6 | 2 | 21 | 120 | 30 | 24 | 24 | −30 (≤ 0 at the floor) |
| A5 | 2 | 21 | 60 | 30 | 12 | 12 | −90 (≤ 0 at the floor) |
| E8, 120 ports | 1 | 25 | 240 | 72 | 54.5 | 80 | +24 |
| 112-port | 1 | 25 | 224 | **64** (the hub's, Lean-checked) | 52.4 | 74.7 | **+32** |

## 2. The centre

**Coordinate decoders** [computed, `calc/natural.out`]. These are the totals T_f = Σ_b (f·r_b) x_b over an orthonormal
basis f of the root space. Each is exact on every pair that must be exact, and leaves ±1/2 on the 60° pairs. A total
costs the dimension of the span of the labels it gathers.

| design | totals | cost per total | per call | per port |
|---|---:|---:|---:|---:|
| 112-port | 8 | 9 | 72 (the hub's own decoder: 64) | 0.643 (0.571) |
| E8, 120 ports | 8 | 9 | 72 | 0.600 |
| E7 | 7 | 7 | 49 | 0.778 |
| D6 | 6 | 5 | 30 | 1.000 |
| A5 | 6 | 5 | 30 | 2.000 |

The "hub-type" charges t(k − 1) used in rounds 1 and 2 (E7 42, D6 24, A5 20) are a formula, not decoders anyone has.

**The clique floor** [proved here; `calc/cliquefloor.py`, `calc/cliquefloor9.py`].
- *Setting.* Every total gathers, descends once to frame 0 and scatters (the hub's and our accounting). It then costs
  at least c(S) = dim span{u_b : b ∈ S}, where S is its gather set. Frames only grow along a route without a descent,
  so the total's frame contains every gathered label.
- *Exact cliques.* Take a clique C of ω ports whose labels are pairwise non-orthogonal. Inside C nothing may be left
  over, because no side information passes between such ports (an2-smallnet 2.1). So the delivered map restricted to
  C × C is the identity I_ω.
- *Counting.* Each total adds a piece of rank at most 1 there, and only if its gather set meets C. So at least ω
  totals meet each maximum clique. Summed over all maximum cliques:
  Σ_k T(S_k) ≥ ω · (number of cliques), where T(S) is the number of maximum cliques that S meets.
- *Ratio.* Each total satisfies c(S_k)/T(S_k) ≥ dim W / T(labels in W), with W the span of its labels. So the centre
  charge per call is at least ω · #cliques · min_W dim(W)/T(W).
- *Computation.* The minimum was found by enumerating every subspace W:
  - every subspace of F₂⁵ and F₂⁷;
  - for the 9-bit designs, every subspace of each dimension 1 to 8 (up to 3.3 million per dimension).

| design | ω | maximum cliques | best W (dim; cliques met) | floor per call | per port |
|---|---:|---:|---|---:|---:|
| A4 | 2 | 15 | 4; all 15 | 8 | **0.800** |
| A5 | 3 | 15 | 4; all 15 | 12 | **0.800** |
| D5 | 4 | 15 | 4; all 15 | 16 | **0.800** |
| D6 | 6 | 15 | 4; all 15 | 24 | **0.800** |
| E7 | 7 | 135 | 5; all 135 | 35 | 0.556 |
| 112-port | 8 | 1,260 | 6; 1,155 of 1,260 | 52.36 | 0.468 |
| E8, 120 ports | 8 | 2,025 | 6; 1,785 of 2,025 | 54.45 | 0.454 |

- **D5, D6, A4 and A5 meet 0.8 per port exactly.** That is the two-pair break-even (4/5), so with one or two pairs they
  can never save, whatever the decoder (in this accounting).
  - This is the floor of about 0.8 that an2-smallnet asked for (its section 9, item 3). The touch-count relaxation
    could not give it; counting cliques does.
  - D4 and A3 were already excluded.
- **For the 112-port design** the floor is 52.4 against the hub's 64. So pr717's floor of 64 (proved there for
  hyperplane-flat decoders) is not proved for all decoders by this argument. Its gap is at most 11 per call.

**E7 in particular.**
- *Symmetric frame decoders* [proved]. Totals Σ_b (f·r_b) x_b, the class of every decoder known:
  - each gathers every root not orthogonal to f;
  - the most lines that one hyperplane can hold is the 36 of an E6 (root subsystems of rank 6; confirmed by sampling,
    `calc/e7hyper.py`), so each total gathers at least 27 ports;
  - 27 odd labels cannot lie in a 5-dim space (it holds at most 16), so each total costs at least 6;
  - at least 7 totals are needed (the exact cliques have 7 ports);
  - so these decoders cost at least 42, the one-pair break-even.
- *Rigidity* [computed, `calc/e7rigid.py`].
  - Take the 7-total decoders M = Ψᵀ Φ with M exact on the 1,953 constrained entries.
  - At the root decoder their tangent space has dimension exactly 111 = 49 (changes of basis) + 62 (diagonal
    rescaling). So every 7-total decoder near the known one is a rescaled copy of it, and costs at least 42.
  - Other components (and decoders with 8 or more totals) are not covered.
- *The window that is left.* A decoder costing 35 to 41 per call. If one with 7 totals of cost 5 existed, each total
  would gather the 15 ports of 1 + V for a 4-dim subspace V. The local rigidity says it cannot be close to the known
  decoder.

**The integer programme** [computed, `calc/decoders.py`; Colab runs in `cxm-dec-*`]. This was a MILP over symmetric
frame decoders built from a pool of small integer directions.
- It found nothing better than the coordinate decoders.
- Its bounds were weak (E7: best 63, bound 14; D6: best 30, bound 10).
- Only for A5 did it close the question within the pool (bound 20 per call, against break-even 10 and 12).
- Its calibration failed. On the 112-port design it did not even find the coordinate decoder in 900 s (best 225,
  bound 15). So it is not evidence for anything beyond A5.
- The clique floor and the E7 argument above supersede it.

## 3. Junk cleaning without idle levels

There are three ways to remove junk, and with one pair all three are bad.

1. **Helpers on banks** (the hub's words, minus the hosting). Banks are D-neutral, but each costs a whole array of
   narrow mass. This is the price in section 5.
2. **Carriers.** Every carrier level is subtracted from D [computed, `calc/e7call.py`]:
   - A5 call: 76 of the 120 junk entries need carriers, which cost 181 levels per call. A5 gains 10 per call.
   - E7 call: 1,780 of the 2,016 junk entries need carriers, 6,415 levels per call. E7 gains 42 per call.
3. **Free meetings.** These cost nothing. But a frame of dimension i holds at most 2^(i−1) odd labels, and at most
   2^(k−i−1) targets can stand in it (their labels lie in its complement). So, when side information travels only on
   copies that never climb, at most v · Σ_i min(2^(i−1), 2^(k−i−1)) junk entries can be met for free [proved, in that
   scope]:

   | design | junk entries per call | free-meeting bound | share | greedy flags built here |
   |---|---:|---:|---:|---:|
   | 112-port | 5,824 | 3,360 | 57.7 % | not built |
   | E7 | 2,016 | 882 | 43.8 % | 229 to 243 |
   | A5 | 120 | 90 | 75 % | 42 to 44 |

   Data arrays can also carry other sources' content on their own climbs, which they pay for anyway (the hub's
   mutable sources and direct delivery). That spreads a source's content over several frames of one dimension, so the
   bound does not cover it. U8X0 uses exactly that and still needs 504 helpers and 3,778 helper levels per call.

So in the shrunk template, most of the junk must go through something that pays, unless a new way of carrying it is
found. The hub avoids the payment by hosting, and hosting needs exactly the idle levels that one pair removes.

## 4. The calls built and checked

Each call is an explicit GF call:
- sources on their lines, targets ending on their hyperplanes;
- the coordinate centre: totals gathered at the full frame, copies descending to 0, scattered, then thrown away;
- free meetings along greedy flags, and paid carriers for the rest.

The check runs on raw complex arrays at the window's block, with frame operators Φ_U = J diag(i^{w(proj_U x)}) J
(every frame used has an orthonormal basis). Every gate asserts one common frame and every step asserts nesting. The
targets must end with Φ_{u⊥}(y + Φ_u⁻¹ x_in) and the sources with Φ_full Φ_u⁻¹ x_in.

| call | block | outcome | paid levels: sources / targets / copies / carriers | total |
|---|---:|---|---|---:|
| A5 (`calc/call_A5.json`) | 5 | **exact** | 60 / 60 / 30 / 181 | 331 |
| E7 (`cxm-call-e7/call_E7.json`) | 7 | **exact** | 378 / 378 / 49 / 6,415 | 7,220 |

**The 112-port word, U8X0** [computed, `calc/csx-rerun/`]. Our checker `bound-check/fsp-cc/tools/csx_check.py` was
rerun. It works in exact rationals, on both the label side and the scalar side; it is the Python mirror of the Lean
checks of the U8X0HQ line.
- Every check passed: hid, hy, exit rows and columns, random-dirt word tests.
- The per-call histogram is {1: 2548, 2: 338, 3: 446, 4: 162, 5: 72, 8: 8}, mass 5,634.
- The regenerated certificate is identical in content to `fsp-cc/out/U8X0.gcert.json.gz` (sha256 87d6080e…).

This is the forward identity per call of the word priced in section 5. The same word is Lean-checked in U8X0HQ.

## 5. Real prices

**Model.**
- One pair, three calls per replica, with the local word's per-call pieces.
- Every helper piece goes on a bank: at most one piece per call per bank, 25 levels per bank, the rest a pad.
- Priced by our fold (`an-ours/anatomy_ours.price`, the Lean fold hypothesis with s = 40) at m = 25 [computed,
  `calc/price_p112.out`, `calc/price_t74_p1.out`].

| word, packing | replicas | W | banks | D | price | against the hub |
|---|---:|---:|---:|---:|---:|---:|
| T74 gauge44 (the hub's record word), ideal | 25 | 16,802 | 11,202 | 800 | **7.3920e-4** | ×0.440 |
| U8X0 (our Lean line's word), ideal | 25 | 16,934 | 11,334 | 800 | 7.3057e-4 | ×0.435 |
| U8X0, first-fit, banks shared across replicas | 25 | 17,484 | 11,884 | 800 | 7.0854e-4 | ×0.422 |
| U8X0, first-fit, banks inside one replica | 25 | 19,450 | 13,850 | 800 | 6.8779e-4 | ×0.409 |
| U8X0, helpers as thrown-away scratch | | | | −11,302 per replica | no saving | |

The banks ignore one constraint: a bank's frame dimension must match the stage it serves in. So these are upper bounds.

**Every cost, per replica (U8X0, ideal packing)** [computed, `calc/items.out`]:

| item | arrays | paid levels | narrow mass at m = 25 | share |
|---|---:|---:|---:|---:|
| sources: 8 levels × 3 calls each (the free line saves one) | 112 | 2,688 | 6,798 | 15.5 % |
| targets: 8 levels × 3 calls each (the free hyperplane end saves one) | 112 | 2,688 | 6,652 | 15.2 % |
| centre: 8 copies of rank 8 per call | (scratch) | 192 | 219 | 0.5 % |
| helpers: 504 per call, 3,778 levels per call | 453.4 banks | 11,334 | 30,089 | 68.8 % |
| D = 25 W − paid = 224 − 192 | | | | **32** |

First order D/N = 7.313e-4, against the fold's 7.306e-4.

**E7, D6 and A5**: D ≤ 0 before any junk cleaning at every decoder known, and for D6 and A5 at every decoder in the
hub's accounting. There is no price to give.

## 6. What shrinking the template costs

The model is the U8 word with generous hosting: every piece that fits by capacity is hosted, the idle remainder goes in
pieces of the largest idle child, and the rest goes on ideal banks [computed, `calc/shrink.out`]. At P = 4 it
reproduces an-ours' exact fill, 1.6647327e-3, exactly.

| pairs | stages | block | D per replica | idle cover per call | pieces hosted per call (of 3,778) | price | against the hub |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 3 | 25 | 32 | 0 | 0 | 7.306e-4 | ×0.435 |
| 2 | 5 | 41 | 128 | 1,792 | 1,749 | 1.467e-3 | ×0.873 |
| 3 | 7 | 57 | 224 | 3,584 | 3,443 | 1.652e-3 | ×0.984 |
| 4 | 9 | 73 | 320 | 5,376 | 3,778 | 1.6647e-3 | ×0.991 |
| 5 | 11 | 89 | 416 | 7,168 | 3,778 | 1.649e-3 | ×0.982 |

an-ours' first-order guess for five stages (×1.8, about 2.9e-3) assumed that every helper piece fits into idle climbs.
Only 46 % do. With the rest priced on banks, five stages give ×0.87. Nine stages is the best layout for the U8 word.

## 7. What it would take

- **To beat the hub at block 25** with the hub's centre [computed, `calc/items.out`]:
  - The helpers' narrow mass per replica must fall from 30,089 to at most 5,381 (17.9 % of U8X0's), keeping U8X0's
    source, target and copy pieces.
  - Per call that is 1,794 against U8X0's 10,030, which spends 1.72 per junk entry.
  - With side information on non-climbing copies only, free meetings take at most 3,360 of the 5,824 entries. The rest
    would have to cost at most 0.73 each: 2.4 times U8X0's efficiency, with no idle level.
  - With mutable data arrays carrying content, no bound is known, but U8X0 already uses them.
  - The hub's whole record chain moved helper mass by about 1 % per step [read: an-hub]. I do not think a 5.6-times
    cheaper helper circuit exists [argued].
- **The most generous ceilings** [computed, `calc/ceilings.out`]. These assume junk cleaning free, no helpers, and
  every data array's levels in the two largest pieces its gates allow, {σ, 2σ}. The piece bound is argued: a data array
  has a gate in each of its three calls, so one cut falls in the middle call's range.

  | design, block | centre per call | ceiling |
  |---|---:|---:|
  | 112-port, 25 | 64 (the hub's) | **8.26e-3** |
  | 112-port, 25 | 61 | 1.06e-2 |
  | 112-port, 25 | 53 (at the proved floor) | 1.68e-2 |
  | E8, 25 | 72 (coordinate) | 5.81e-3 |
  | E7, 19 | 42 / 41 / 38 / 36 | 0 / 1.76e-3 / 7.05e-3 / 1.06e-2 |
  | 112-port, 25, U8X0's own data pieces | 64 | 2.33e-3 |

- So 1e-2 at these blocks needs two things at once. The first is a centre below anything known: at most 61 per call
  for the 112-port design against the hub's 64, or at most 36 for E7 against 49. The second is junk cleaning that is
  almost free.
- **The one stack of assumptions that would beat the hub at block 19** [argued, first-order estimate]. Suppose:
  - an E7 decoder at the floor of 35 existed (D = 21 per replica);
  - free meetings reached their bound;
  - the residual junk were cleaned at U8X0's efficiency on banks.

  Then E7 would come to about 2e-3. Each of the three is unproven, and the first is the crux: it needs a decoder that
  local rigidity puts away from the only one known.

## 8. Predictions (STATUS.txt, written 12:35) against outcomes

| prediction | outcome |
|---|---|
| C1 (90 %): nothing in 13 to 25 beats 1.68e-3 at its real price; nothing reaches 1e-2 | **held**: best 7.39e-4. Nothing built reaches 1e-2; at block 25 it is excluded with the hub's centre (8.26e-3 even with free junk cleaning) |
| C2 (80 %): E7 at 19 has D = 0 at the natural charge (42); nothing below 42 found (60 %) | **held, and worse than predicted**: the real decoder costs 49, so D = −21. 42 was a formula, and it is now the proved floor for every decoder of the known kind. The general floor is 35, and no decoder below 49 was found |
| C3 (75 %): the 112-port design at 25 with U8 on banks prices between 4e-4 and 9e-4; D = 32 | **held**: 6.88e-4 to 7.31e-4 (7.39e-4 with T74's word); D = 32 |
| C4 (85 %): D6, A5 at 13 and 21 have D ≤ 0 with every decoder I can construct; a cheaper one found: 25 % | **held, more strongly**: D ≤ 0 is now proved for every decoder in the hub's accounting (floor 0.8 per port, break-even 2/3 and 4/5) |
| C5 (85 %): an explicit call passes an exact raw-array check; free meetings ≤ 882 of E7's 2,016 | **held**: A5 and E7 exact. The bound holds, but only when side information rides on non-climbing copies, a scope I had not stated; the greedy flags reached 229 to 243 |
| C6 (70 %): the best real price is the 112- or 120-port design at 25, below 1e-3 | **held**: 112-port, 7.39e-4 |

## 9. What this leaves

1. **Medium blocks are closed for the record inside the hub's template, as far as anyone's words go.** The one design
   that saves there (112 ports at block 25) saves 0.14 of a direction per data array. With no idle level to host
   helpers, banks eat it down to ×0.44 of the hub. Fewer stages lose for the same reason (section 6).
2. **For 1e-2:** D5, D6, A4 and A5 are excluded at blocks 9 to 21 for every decoder (in the hub's accounting). E7 at 19
   and the 112-port design at 25 are excluded unless both of these hold: a new centre (E7: 35 to 36 per call;
   112-port: 53 to 61), and junk cleaning far cheaper than any word has.
3. **Cheap next steps, if anyone wants them:**
   - an exact search for an E7 decoder with 8 or more totals of cost at most 5 (it would need a decoder far from the
     root decoder);
   - the same for the 112-port design between 53 and 63, which would also lift the hub's own record. Each unit of
     centre is worth about 2.8 % at block 73 (an-ours);
   - a multi-layer version of the clique floor;
   - a bound on free junk cleaning that covers mutable data arrays.
4. **The new tools** run on any port design in seconds:
   - the clique floor;
   - the free-meeting bound (for non-climbing copies).
5. **A skeptic should check:** the clique-floor argument (section 2), the E7 symmetric-decoder bound and its rigidity
   computation, and the piece bound behind the 8.26e-3 ceiling.

## Files

All in `checks/wht43/anatomy/cx-medium/`:
- `STATUS.txt` (dated predictions, progress, outcomes), `RUNLOG.txt`, this `NOTE.md`.
- `calc/designs_cx.py`: designs with explicit F₂ labels.
- `calc/natural.py` → `natural.out`, `natural.json`: coordinate decoders.
- `calc/cliquefloor.py` → `cliquefloor.json` (k ≤ 7); `calc/cliquefloor9.py` → `cliquefloor9_*.json` (k = 9, by
  dimension; d = 4, 5 and 6 in `cxm-clq-a/`, `cxm-clq-b/`).
- `calc/e7pool.py`, `calc/e7hyper.py`, `calc/e7rigid.py` → `e7rigid.json`; `calc/e7sparse.py` → `e7sparse.json`. That
  last search used generic values, which always fail the column condition, so it is uninformative and is not cited as
  evidence.
- `calc/decoders.py` → `decoders_A5_D6.json`; Colab runs in `cxm-dec-e7/`, `cxm-dec-d6p112/`, `cxm-dec-d6pair/`,
  `cxm-dec-p112/`.
- `calc/price_p112.py` → `price_p112.out`, `.json`; `calc/price_t74_p1.py` → `price_t74_p1.out`, `.json`.
- `calc/shrink.py` → `shrink.out`, `.json`; `calc/items.py` → `items.out`, `.json`; `calc/ceilings.py` →
  `ceilings.out`, `.json`.
- `calc/e7call.py` → `calc/call_A5.json`; E7 in `cxm-call-e7/` (`call_E7.json`, `e7call.out`).
- `calc/csx-rerun/`: our checker's rerun on U8X0 (`csx_checks.json`, the regenerated certificate).

## Sources

- Briefs: `../BRIEF3.md` (section cx-medium), `../BRIEF.md`, `../BRIEF2.md`; `../NOTE.md` (rounds 1 and 2).
- Round 2:
  - `../an2-smallnet/NOTE.md` and `calc/` (the frame operator in frames.py, ceiling.py, designs.py,
    classfloor.json), `../v-an2-smallnet/NOTE.md`;
  - `../an2-word/NOTE.md`, `../v-an2-word/NOTE.md`.
- Round 1:
  - `../an-hub/NOTE.md` (T74 profiles, idle-child layout, caps);
  - `../an-ours/NOTE.md`, `stages.py`, `u8word.py`, `anatomy_ours.py` (the fold);
  - `../an-build/NOTE.md`, `calc/calib.py`, `calc/hosty.py`.
- Our data:
  - `bound-check/fsp-price/data/words.json` (U8 E1 profile and leftover widths);
  - `bound-check/fsp-cc/out/` (U8X0 certificate) and `bound-check/fsp-cc/tools/csx_check.py`;
  - `bound-check/fsp-price/tools/host_price.py` (header: the idle-child table).
- Our Lean template, read as text: `wht43/bound-check/h-L3-design/gen/template/Gates.lean` in the private
  archive of the agents' working notes (the GF/GB order, the twin gates, the run's result).
- This repository:
  - `tools/gx/gxcore.py` (gcert/1 semantics: y at frame 0 at the scatter; a copy costs its total's frame dimension);
  - README sections 5 to 7; `notes/scratch-copies.md`.
- The community note on general frames, `general-clifford-frames.tex` (local copy under `wht16/rs-crocswap/sources/`),
  read as data only.
- The 63 E7 lines as the points of a symplectic F₂⁶ (commuting = orthogonal roots) are recomputed here from the
  labels, which are all odd vectors of F₂⁷ except the all-ones vector.
