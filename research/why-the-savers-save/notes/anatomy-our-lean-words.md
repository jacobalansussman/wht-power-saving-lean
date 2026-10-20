> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# an-ours: why our Lean-checked words save, and what caps them

Anatomy track, 15 October 2026, 08:12 to 08:40 EDT by the Mac clock. Source: our own Lean-checked words C588MHH
(1.5699305e-3), U8XMHH (1.6197705e-3) and U8X0HQ (1.6256629e-3), with the rest of their ladder (Shea's C588,
the C588 repack, C588M, U8X, U8XM, U8X0, U8X0HH) and, for comparison, our banked lowering of the hub's T74 word.
Nothing was posted. No third-party code was run. The arithmetic is in `anatomy_ours.py`, `identity.py`,
`knobs.py`, `stages.py` and `u8word.py` in this folder: our own code, small, on the Mac, under 1 s each. Labels: [computed] = a script here
and its output; [read] = taken from a named file; [argued] = my reasoning, not checked by anyone else.

## In short

1. **Our words do not contain a saving object of our own.** In every one of them the saving comes from the
   community's network and its helper circuit (the #401 nine-stage network of Dugongue and Mamano, with Kothe's
   atm3 circuit or the U8 circuit). What we added (exit frames, mutable sources, no centre copies, same-run
   hosting, bank packing, a sharper rate lemma) removes cost that our Lean realisation charged and the hub's
   pricing never charged. Prediction P1 held.
2. **The saving is a fixed deficit divided by a variable "narrow mass".** Every word on both of our lines, and
   the hub's T74, has the same deficit per data array, **5/14 of a move** (2560 moves on 7168 data arrays). The
   exponent is close to D / N, where N = sum over paid blocks of r ln(73/r). All of our 31.8 % gain over Shea,
   and all of the hub's gain from atm3 to U8 to T74, is a reduction of N. Not one step changed D. [computed]
3. **Where D comes from, move by move:** in each run, every data array pays 72 of its 73 directions and gets
   the 73rd free. That free direction is the pair exchange: two arrays at the same frame swap contents through
   free additions. The eight central totals of each call give back 64 moves (eight blocks of rank 8). Per run:
   +896 - 576 = 320. Once their leftover pieces are hosted, the helpers cost no moves net. Each call's helpers
   pay exactly as many in-call moves as their leftover pieces later remove from the idle climbs (5634 - 3778 =
   1856 = 224 x 8 + 64 on U8; 5916 - 4060 = 1856 on atm3). [computed from our own tables]
4. **What caps our words:** the hub's word fixes D and about 95 % of N. Our layer can only shrink the rest. On
   U8 that rest is 85,217 of N's 1,570,083 in U8X0HQ, against 48,238 at an exact fill. The U8 line's ceiling in
   this network is therefore **1.6647327e-3** (+2.40 % over U8X0HQ), and only cross-run hosting moves towards
   it (1.6438659e-3 priced). Inside same-run hosting, U8X0HQ is as good as our tools can make it: the host ILP was
   solved to optimality for its first-order objective (HiGHS, gap 0), and the bank table has one pad, which
   x-pack's concavity argument makes optimal. Prediction P3 held.
5. **What caps the mechanism itself:** one free direction per data array, out of m = 73. Even with free central
   totals (D = W) and the U8 helper circuit, the price is only **4.67e-3** [computed what-if]. To reach 1e-1, the
   same "one free direction per array" has to happen at block m <= 6 with rank-1 moves (s = ln(m/(m-1))/ln m:
   0.139 at m = 5, 0.102 at m = 6), or arrays must get more than one free direction each. Within this family
   the block size is set by the helper circuit itself: its leftover pieces must fit into the idle climbs of the
   other pairs, and for U8 that first happens at nine stages (m = 73). Five stages (m = 41) would be worth about
   a factor 1.8, if the helper leftovers were half as wide [argued, first order].

## 1. What the saving object is

- **Block and arrays.** m = 73 directions per array per pass. W = 7168 data arrays: 8 runs (replicas) x 2 sides
  (X, Y) x 112 ports x 4 pairs, which is W = 2PTv with P = 4, T = 8, v = 112. Our Lean words add banks to W:
  Shea 5040, C588M 4032, U8XM 3771, C588MHH 111, U8XMHH 310, U8X0HQ 323. [read: Rate*.lean, pr717 P4,
  ATTEMPT.md]
- **Calls.** 9 stages per run, 72 calls in all. Each call runs the local helper circuit on the active pair's
  224 data arrays (PAIR = 1,1,1,2,2,3,3,4,4) with 504 helpers on both atm3 and U8 (616 with the source copies)
  and 8 central totals. The
  other pairs sit idle and climb in fixed "children" of rank 16, 8, 32, 8, 48, 8, 48, 32, 16 (our x-hosting
  model of the Lean layout). [read: `bound-check/x-hosting/tools/host_price.py` header; fsp-price `words.json`]
- **Moves.** A block of rank r is one recursive call on r of the 73 directions. The price is the root a of
  sum_r (r n_r / 73W) (73/r)^a = 1, with a 2^-40 fill term (the Lean fold `fold_*`). The per-call profile of U8
  is {1: 2548, 2: 338, 3: 446, 4: 162, 5: 72, 8: 8}; the helpers' leftover pieces have widths {4: 80, 5: 32,
  7: 54, 8: 122, 9: 216}. [read]
- **Mixed targets?** No mixed targets in the sense of the moonshots. All 224 active arrays are targets of the
  same kind. The sharing is between targets of one call: the 8 central totals are read by every target, and
  target rows carry 36 to 97 terms (pr717 C2 on the #626 word). [read]
- **Game and free steps.** The engine game of our generalised engine: free additions between arrays at equal
  frames (any subspace of F_2^m, the frame lemma), address permutations, multiplication by an address-dependent
  power of i, entrance gauges (a helper enters part of the way up its climb), exit frames (a helper ends below
  the full frame), and the undo step that restores the helpers. The circuit's scalars are rational (register
  scales absorbed into rational units, ATTEMPT.md) [read: README section 5; ATTEMPT.md].
- **Halving: yes.** In our lowering of U8 (`k-lower/out/U8X-parsed.json.gz`), 2042 of the 5919 single additions
  (in 7500 gates) have a denominator: 1501 have 2 (1494 of them the coefficient 1/2), 401 have 4, 112 have 8 and
  28 have 16. Pairs such as y -= x; x += y/2 form a + b and a - b in place. The scatter gives each of the 112 targets 4 of the 8 totals, with coefficients +-1/2 and +-1
  [computed, `u8word.out`].

## 2. Where the saving comes from

**The identity** [computed, `identity.py`]. Per run, with 896 data arrays:

| item | moves per run |
|---|---|
| every data array: 8 levels per active stage + its idle children | 72 each, 64,512 in all (pair 1: 24 + 48; pairs 2-4: 16 + 56) |
| 8 central totals per call, one rank-8 block each, 9 calls | 576 |
| helpers: in-call blocks minus the idle rank their hosted leftover pieces relieve | 0 (9 x (3778 - 3778)) |
| total paid | 65,088 = 73 x 896 - 320 |

The local circuit's rank minus its leftover widths is 1856 = 224 x 8 + 64 on atm3, on U8 and on U8 with five
centre copies alike. That is why D = 2560 on every word: the helper circuit, where all of the community's
ingenuity goes, cannot change D in this network. It changes only N. The centre copies of U8XM are also neutral
in rank (+45 moves and +45 of leftover width per call); they cost only narrow mass.

This also holds helper by helper, as distributions [computed, `u8word.out`]. In our lowering of U8, the span of
each helper in the call (the dimension of its final frame minus that of its start frame) is distributed
{4: 80, 5: 32, 7: 54, 8: 122, 9: 216}, plus 112 copy helpers of span 9. These are exactly the leftover widths
of the U8 word. A helper pays its span once inside the call, in small blocks. Its leftover piece of the same
width is then climbed inside a host's idle climb, which shrinks by that width. So a helper adds no rank. All it
does is replace part of a big idle block by small in-call blocks.

**The deficit, move by move.** +1 per data array per run is the free direction, the "+1" in m = 1 + 8(2P+1).
It is done by the pair exchange through free additions. Inside one call this shows in the local word's frames
[computed, `u8word.py`]:
- each source x starts on the line of its own port label and climbs 8 levels to the full local frame;
- each target y starts at 0 and ends on the hyperplane perpendicular to its port, also 8 levels.

So within a call each data role pays 8 levels on a local frame of dimension 9 (h = 9). How the local frames are
glued from stage to stage is the network's part (the frame equations x01, y01, ... of `GNetDef`). Over a run the
count is 72 paid of 73 per array, as in the table above. The local price identity closes exactly: rank =
R h + 2v(h-1) + cst - entrance levels - exit levels = 616 x 9 + 1792 + 64 - 408 - 350 = 6642 (U8X, with copies;
126 helpers end at exit frames).

The -64 per call is the eight central totals, each climbing one rank-8 block. The hub's own formula agrees:
D = T[2Pv - (2P+1)cst] = 8 x (896 - 576) (pr717 P4, CHECKED there as a model). Our own E8 unit has the same accounting: "of the 480 moves that the 120 exchanges
could save per unit, the scratch copies of the 8 totals take back 360" (README section 6). There the deficit per
data array is 1/4; on the paired cube it is 0.583; on the P112 line it is 5/14.

**The price, move by move** [computed, `run1.out`]. N for U8X0HQ is 1,570,083. Of this, 1,484,866 (94.6 %) is
72 x the hub's U8 local profile, and 53 % of that is rank-1 blocks at ln 73 = 4.29 each. The other 85,217 is
idle children left unhosted, the small remainders where a leftover piece did not fill a child (2320 rank-1,
864 rank-2, 256 rank-3 and 896 rank-4 pieces cost 29,000), and the 323 banks' unused children. Of the local
word's narrow mass, the helpers carry 77.0 % in our lowering of U8, the targets 13.8 %, the sources 8.6 % and the
centres 0.6 % [computed, `u8word.out`; U8X, with the 112 copy helpers]. The helpers' blocks are
{1: 2308, 2: 282, 3: 414, 4: 98, 5: 56}. pr717 C1 gives 68 % on the #656 gauge44 word, U8's ancestor [read].

**How it escapes our near misses' failure modes** [argued, from the sources above]:
- *Zero-divisor dilemma.* No leftover is killed by a factor such as 1 + u. The exchange leaves junk on each
  target: in our E8 unit, halves of the 56 sources it "knows" (README section 6); on P112 the pattern is the
  E8 inner products (pr717's decoder model). Helpers holding exact signed partial sums subtract that junk with
  free additions. The undo step then restores the helpers by running the exact inverse of the forward run, so
  the readout of each target stays invertible. The cancellation is subtraction plus uncomputation, never a
  zero divisor.
- *Pattern identity* (every junk row is touched only by moves that also reach the target). Paid blocks here act
  on whole arrays and are never asked to separate junk from target. The separation is done by free additions
  from separate helper arrays.
- *One-route trap.* Each target receives its value over many routes (36 to 97 terms). The central totals are
  shared by all 224 targets, and one helper set serves every stage.
- **The price of the escape** is the helpers' climbs. In rank they are free: the identity above makes them
  neutral, because the multi-stage layout has idle climbs of other pairs that the helpers' leftover pieces
  replace. In narrow mass they are not free: they turn cheap big idle blocks (rank 8 to 48, ln(73/48) = 0.42)
  into expensive small in-call blocks (rank 1 to 5). A small unit at block 3 to 6 has no idle arrays to absorb
  its junk removal. I take that to be the structural reason why the network saves and the small units do not.

## 3. What caps the saving

| level | figure | what caps it | status |
|---|---|---|---|
| U8X0HQ (our best) | 1.6256629e-3 | same-run hosting: 4032 of 36,288 leftover pieces have no idle child in their own run, so 323 banks remain; the host ILP is optimal for its first-order objective and the table has one pad | Lean-checked; ILP gap 0 [read: A_h869.padone.json] |
| cross-run hosting, no copies | 1.6438659e-3 (1.6459856e-3 deferred finals) | needs a cross-replica hosted network in Lean (1,500-3,000 lines, unwritten) | Python price [read: fsp-price] |
| U8 at an exact fill | **1.6647327e-3** (T = 8; 1.6647298e-3 at T = 37) | D fixed, the U8 local word's N fixed, the idle remainder all at rank 48: the floor for this word in this network | [computed; matches bh-plan] |
| the hub's T74 word | 1.6797550e-3 | a better helper circuit: rank-1 blocks per call 2548 -> 2340, N per call 20,623 -> 20,421 (-1.0 %) | hub price [read] |
| free central totals (D = W) on the U8 word | 4.67e-3 | one free direction per array at m = 73 | what-if [computed] |
| one free direction per array, rank-1 moves | ln(73/72)/ln 73 = 3.2e-3 | the same, at the plainest price | [computed] |

**Tight certificates that exist** (all from pr717 and CHECKED by its author, not by us): the centre decoder
floor cst >= 64 is attained by the coordinate decoder, on 8-term hyperplane-flat decoders (D1). The 112 ports
are a cap in F_2^9, with every hyperplane missing at least 48 of them (C3). gauge44's 44 cold ends equal the
gauge bound (G2). For a fixed gate order the frames are optimal (#698 exactly, P5 sampled). None of these bounds
the exponent itself. Each bounds one ingredient at fixed structure.

So our words are capped in three nested ways:
- **by pricing:** our layer can only move from 1.6257e-3 to 1.6647e-3, and every rung of that is overhead the hub
  already leaves out;
- **by the hub's helper circuit:** about 95 % of N is the hub's local word;
- **by the mechanism:** b is about (D/W) / (m l), where l = N / (paid rank) is the mean of ln(73/r) and is 2.89
  on U8X0HQ. D/W can never exceed 1 here, because there is one free direction per array.

**Knobs, priced on U8X0HQ or taken from the hub's checked knob table:**
- one more unit of deficit per call (for example a centre climbing 7 levels instead of 8): +2.82 % [computed];
- one centre fewer: +22.6 %; one port more: +2.15 % [read: pr717 P1, CHECKED there];
- every local rank-1 block paired into a rank-2 block: 1.8154e-3, +9.1 % over the exact fill (hypothetical,
  for scale) [computed].

## 4. How it was found, and what each step changed

The word is the hub's: Kothe's atm3 (#588); his #799 and then the U8 circuit (#869), with Boukhalfa's frame
changes (#741), icekylinx's gauge44 (#645) and the triangle-hybrid of NewBArch and ColdWitness/icekylinx; the
whole line rests on the #401 network of Dugongue and Mamano. It was found by the community's search and
structure work, which an-hub covers. Our layer was found by Lean engineering and by integer programmes: x-pack's
concave bank packing, and hosting ILPs laid out after DaysSky's hosting_core.

| step | figure | change | what it removed (D = 2560 throughout) [computed] |
|---|---|---|---|
| Shea C588 (start, 13 Oct) | 1.2331685e-3 | - | N 2,071,832; 5040 banks; 112 source copies per call |
| bank repack, one pad | 1.2383291e-3 | +0.42 % | 2448 rank-1 pad blocks -> 724 + one rank-37 block (N -8,645) |
| C588M: mutable sources | 1.3552606e-3 | +9.44 % | the copy helpers' rank-1 entrance and the sources' rank-8 climb, 112 per call (N -178,435; W -977) |
| C588MHH: same-run hosting | 1.5699305e-3 | +15.84 % | 3921 banks' unused idle children (N -258,736; rest 17 % -> 6 % of N) |
| word atm3 -> U8, on the hosted rung | 1.6197705e-3 (U8XMHH) | +3.17 % | the hub's better helper circuit (N -50,207) |
| U8X -> U8XM -> U8XMHH | 1.2704023e-3 -> 1.3868856e-3 -> 1.6197705e-3 | +9.17 %, +16.79 % | the same two steps on U8 |
| U8X0HH: no centre copies | 1.6256506e-3 | +0.36 % | 5 copies per call (N -5,718) |
| U8X0HQ: one-pad table, term_le4 | 1.6256629e-3 | +0.0008 % | a rate-lemma grid loss of 123 points; one rank-1 piece |

Of our factor 1.318 over Shea, a factor 1.032 is the hub's change of word (about 11 % of the gain on a log scale).
The remaining factor 1.277 is overhead that our Lean route had added and our layer removed. Prediction P2 held, with one correction: the steps were not "less than one unit's
worth" each. Hosting alone is +16 %, but it is still pure overhead removal.

## 5. Knobs that would scale it, and what stops each

- **Finish our layer** (cross-run hosting, then an exact fill): at most +2.4 % on U8. It is stopped by the Lean
  network, which is unwritten, and in the end by the ceiling 1.6647e-3. Porting T74 only ties the hub
  (1.6797548e-3). The recursion lever (grouped heads, Theorem 1G) is +0.90 % over the hub as a what-if; it is
  unbuilt and needs a new engine (HANDOFF, 14 Oct). It is also a narrow-mass lever: it reprices 11,010 rank-48
  heads, and D does not move.
- **Cheaper central totals (cst below 64).** Worth +2.8 % per unit per call, and +180 % if the totals were
  free. It is stopped by the decoder floor: cst >= 64 for every decoder whose costly terms are hyperplane-flat
  supports (pr717 D1, CHECKED there; decoders that mix terms of beta <= 7 are open). These centres are the
  "scratch copies of the totals" of our own notes (`notes/scratch-copies.md`), and the same lever moved the
  lineage before. Copy rank per invocation went from 241 (our first circuit, h = 16) and 440 (the star rule,
  h = 22) to 72 (our E8 unit: 8 totals at dimension 9) and then 64 (the hub's P112: the coordinate decoder, 8 star terms of rank 8,
  on the E8 root lines with the 8 frame lines removed). Here they still take back 576 of every 896 free directions (64 %).
- **More ports per run (v).** +2.15 % per port. It is stopped by the label geometry: the 112 ports are the E8
  root lines minus the 8 frame lines, a cap in F_2^9 (pr717 C3, D1 model).
- **Cheaper helper circuit (lower l).** This is the knob the whole hub record chain turns: atm3 -> U8 -> T74
  gains about 1 % per step, and the price is dominated by rank-1 helper blocks. It is stopped by the search;
  no certificate is known. P5 in pr717 shows that frame annealing at a fixed gate order gains 0.0 %.
- **Smaller m with the same one free direction per array.** This is the only knob with a large effect,
  because b is about (D/W) / (m l). It is stopped by the layout. With P pairs the run has 2P + 1 stages and
  m = 1 + 8(2P + 1). The helpers' leftover pieces must fit into the idle climbs of the other pairs, and those
  come to 16 v (P - 1) per call. So the helper circuit's leftover width per call must satisfy Rw <= 16 v (P - 1).
  The formula reproduces pr717 P4's cover of 3584 per call at 7 stages [computed, `stages.py`]:

  | P | stages | m | D per data array | (D/W)/m | idle cover per call | U8's 3778 fits? |
  |---|---|---|---|---|---|---|
  | 1 | 3 | 25 | 0.143 | 0.0057 | 0 | no |
  | 2 | 5 | 41 | 0.286 | **0.0070** | 1792 | no (0.47) |
  | 3 | 7 | 57 | 0.333 | 0.0059 | 3584 | no (0.95) |
  | 4 | 9 | 73 | 0.357 | 0.0049 | 5376 | **yes (1.42)** |
  | 5 | 11 | 89 | 0.371 | 0.0042 | 7168 | yes |

  So nine stages is the smallest layout in which U8's helper leftovers fit. A helper circuit with half of U8's
  leftover width (at most 1792 per call) would allow the five-stage layout. To first order, with
  l(41) = l(73) - ln(73/41) ≈ 2.31, that is worth about a factor 1.8, roughly 2.9e-3 [argued, first order; no such
  circuit is known]. The hub's model at 7 stages gains only 0.1 to 0.2 % (pr717 P4), because of this overflow.
  Small blocks have no idle cover at all. Even the P -> infinity limit of D per data array, 1 - cst/v = 0.43,
  keeps the family far from 1e-1.
- **More than one free direction per array** (nested exchanges): no construction is known. The per-array
  deficit of every network in our lineage is at most 0.583, and 5/14 on this line.

**For an-synth, the plain lesson** [argued]: the real saver is a layout fact (one free direction per data
array, through pair exchanges), paid for by junk removal whose rank cost is absorbed by idle climbs that the
layout already has. A 1e-1 construction needs that same exchange at block 5 or 6, where there are no idle
climbs, or exchanges that give more than one free direction per array.

## Files

- `anatomy_ours.py`, `run1.out`, `anatomy_ours.json`: prices of all 11 words with our own code (they match the
  Lean kernel numerators exactly for the term_le4 words, 12383292 and 16256630, and lie 12 to 33 grid points
  above the term_le words, as expected); N decompositions; ceilings.
- `identity.py`, `identity.out`, `identity.json`: the per-run rank identity (72 + 1 per array; local rank minus
  leftover width = 1856 on every word).
- `knobs.py`, `knobs.out`, `knobs.json`: what-ifs (exact fill, free totals, paired rank-1 blocks, small-block
  formula, per-array deficits).
- `stages.py`, `stages.out`, `stages.json`: the stage-count knob (D per data array, m, idle cover per call).
- `u8word.py`, `u8word.out`, `u8word.json`: class-by-class replay of our lowering of the U8 local word
  (`k-lower/out/U8X-parsed.json.gz`): blocks and narrow mass of sources, targets, helpers and centres, the local
  price identity, the coefficient census.
- Sources read: the Rate*.lean files listed in `anatomy_ours.py`; `bound-check/fsp-build/out/RESULT.json`;
  `bound-check/fsp-price/out/A_h869.padone.json` and `data/words.json`; `bound-check/x-hosting/tools/host_price.py`
  (header); `bound-check/d-lean-route/shea-src/research/complex-transfer-588-596/ATTEMPT.md`;
  `bound-check/g-p112-decoder/pr717/README.md` (P1, P4, P5, C1-C3, D1); `bound-check/bh-plan/PLAN.md`;
  the r6 candidate README (sections 5 and 6); `meta/MAP.md` 1.3; `moonshot/NOTE-wave4.md`; `meta2/NOTE.md`.
