# Idea rounds: 10 October evening and 16 October

Most of the research since the fourth revision was construction: searches for a small program that saves work. This
section tells the two rounds in which I asked for thinking instead. The first ran on the evening of 10 October 2026,
right after the fourth revision (saving 8.7625e-4): I asked where the effort should go next, and two follow-ups ran
the same night, a floor theorem for programs of linear size and the first bounds on very small units. The second ran
on 16 October, after the moonshot waves had found no unit: I asked for one serious attempt at a new mechanism, with no
big coding projects, generation first and criticism last. The ten idea waves of the night of 10 to 11 October, which
produced the 5/6 circuits, are public in
[wht-fewer-operations](https://github.com/jacobalansussman/wht-fewer-operations) (its `waves/` folder) and are not
repeated here.

I did this work with a team of AI agents (Claude) that I directed: I set the goals, chose which ideas to pursue and
decided at each stage what to do next; the agents read the papers, wrote the proofs and programs, and tried to refute
one another's work. **No person has read any of the proofs in this section, and nothing in it is checked in Lean.**
Each finding below carries its status in words: *proved in Lean* (none here), *checked by computer*, *argued by hand,
not checked*, *refuted* (and by what), or *open*.

## In short

- **Neither round found a unit, a new mechanism or a route to linear size.** What came out is a list of closed doors,
  a few structural facts, and leads that are still open, each with a first experiment.
- **Closed on 16 October** (each argued by hand, not checked, and none of them examined by a skeptic):
  - the second mechanism inside OpenAI's existence proof: its 25 % margin lives only in a "price world" that real
    programs do not obey, the win it forces is astronomically small, and even if realising its laws cost nothing, a
    fixed-matrix version could not certify more than 1.97e-3 (1.58e-3 with OpenAI's own count of the extra pairs in
    its chain); it also cannot be made right for every number of digit columns;
  - streaming (feedback) units: a streaming cell is an ordinary unit with unboundedly many arrays, and feedback
    through free steps is worth nothing;
  - "busy" cells in which every register moves in every layer: they need at least m moves per input, so any saver must
    leave some register idle;
  - phases that are not powers of i (eighth roots, any complex scaling): the hook certificates are kept by every
    diagonal, so the bounds 2, 14/5 and 41/13 moves per array at blocks 2 to 4 stand;
  - the Schur-Weyl basis of the transform: too expensive to reach.
- **New on 16 October:**
  - no exactly tight unit exists at block 4 in the engine game: 0 of the 30 positions one move from the start and 0 of
    the 280 two moves away pass the termwise test that exact tightness requires. This closes the grouped 2 in 7
    (blocks of rank 2, 2, 2, 1; it would have saved 1.675e-1) and 13 in 41 (1.71e-1). The computation is checked by
    computer (one program); the argument around it is argued by hand, not checked, and has no skeptic;
  - the best saving of Kronecker powers is the same for every kernel width (argued by hand, not checked);
  - an exactly tight 5 in 14 at block 3 (6.28e-2) lives in a small sub-game whose group is the symmetric group S8, and
    it can carry very little symmetry (checked by computer; the structure was repeated by a skeptic);
  - E8 is the best possible centre for the record units among root systems, by Haemers' bound from index coding
    (argued by hand, not checked);
  - the triangular kernel [[1,0],[1,−1]] is the one binary kernel that the block-3 hook does not pin down (2.1625
    moves per array instead of 2.8), but the same certificate is too weak at block 2, so that room is uncalibrated.
- **The two skeptics of 16 October:** the lead "tight block-3 family" survives as a question, but the plan to test it
  does not; the lead "targets on a subspace, block-4 cores first" is **refuted**.
- **10 October, floor theorem:** if some pair of digit columns keeps its diagonal under every cheap step, a unit needs
  more blocks than live arrays. So inside the Lean engine neither linear size nor one block per array is reachable.
  Argued by hand, not checked (three agents by three routes, three skeptics re-derived it); weak on the exponent
  itself.
- **10 to 11 October, tiny units:** in the sub-game without address maps, nothing saves at block 3, and block 4 needs
  at least 3.058 moves per array (both argued by hand, not checked). Later work absorbed or superseded most of it:
  the block-3 certificate became half of Theorem C3 (15 October), the old 1e-1 target 3 in 8 at block 3 is excluded
  in the whole engine game by the 14W/5 bound (proved in Lean, 13 to 14 October), and 2 in 6 at block 4 by 41/13
  (checked by computer, upheld by skeptics). The night's exclusion of 2 in 7 at block 4 still stands, in the sub-game
  only.
- **Open, each with a first experiment** (table near the end): the tight block-3 family (start with one 6-move
  problem), the triangular kernel (first find a lower bound that is right at block 2), a defect budget for near-tight
  block-4 schedules, the 7-array problem at block 5, an idle-register carousel at block 6, and the escape hatches of
  the floor theorem.

## Words used here

The glossary of [../README.md](../README.md) applies; the words this section leans on most are these.

- **Unit, "W in M", saving.** A unit is a small fixed program that transforms W arrays with M *paid moves*. At *block
  size m* each array is owed m directions, so the plain cost is m·W. A paid move applies one direction's 2-point
  transform to one array (inside the recursion it is one recursive call); everything else is a *free step*
  (relabelling addresses, phases, adding multiples of arrays). Used recursively, a unit gives the transform in about N
  (log N)^(1−s) operations with saving s = ln(mW/M)/ln m.
- **Games.** The *engine game* has the free steps of this repository's Lean engine: linear address maps, shifts,
  Clifford phase patterns (powers of i) and rational mixing of arrays. The *sub-game without address maps* of 10 to 11
  October keeps everything except address maps. The *no-phase game* keeps everything except phase letters.
- **Digit columns.** An address is f digit columns of m bits; a free step acts the same way on every column, and a
  unit must be right for every f.
- **Frame, depth.** Every term of an array's contents sits at a point of a finite graph, its *frame* (135 frames at
  block 3, 2,295 at block 4). Free steps leave a term at its frame; one paid move takes it to a neighbouring frame. A
  term's *depth* is its distance from the start frame, where the raw input sits.
- **Certificate, cap, demand, hook.** A certificate is a subspace that every free step keeps. One paid move can push
  at most *cap* dimensions out of it, and each array's target pushes out its *demand*, so M·cap ≥ W·demand (this is
  "Theorem R"). The *hook* is the best family of certificates known; it gives 14/5 moves per array at block 3 and
  41/13 at block 4.
- **Exactly tight.** A shape is exactly tight when M·cap = W·demand: every move must then gain the full cap and waste
  nothing. 5 in 14 at block 3 is exactly tight (14 × 80 = 5 × 224).
- **Skeptic.** An independent AI agent that tries to refute a claim with its own reasoning and its own code.
- **The bar.** A saving of 1e-1 is the standing target for a new route; since 15 October 1e-2 counts as a
  stepping-stone tier.

---

## Part 1. The evening of 10 October

### 1.1 Where to put the effort

After the fourth revision I asked four questions: what would it take to reach 1e-2, 1e-1, linear size, and a faster
FFT on a real machine? Five agents each looked through one lens (very small units, levers on the present mechanism,
imports from other fields, programs outside the engine, real machines); five skeptics checked them; a sixth agent
wrote the plan. Its answer, in its own words, was "one candidate, no mechanism": a tiny unit in the game the engine
really accepts would clear the bar, nothing excluded one, and nobody knew how one would work. It judged that the
present mechanism (one free direction per data array, through an E8 centre) could not get there, that linear size had
no route inside the engine, and that on a real machine the E8 unit makes 167 to 852 whole-array additions per saved
move where a real machine allows about 3 to 10. The odds it attached to each avenue were the agents' guesses, not
measurements, and I do not repeat them. What happened to each avenue:

| Avenue of the plan | What it hoped for | What happened |
|---|---|---|
| Tiny units, bounds first: 3 in 8 at block 3, 2 in 6 and 2 in 7 at block 4 | savings 1.07e-1, 2.075e-1 and 9.63e-2 | run that night (section 1.3). Later: 3 in 8 at block 3 is excluded in the whole engine game by 14W/5 (proved in Lean); 2 in 6 at block 4 by 41/13 (checked by computer, upheld by skeptics); 2 in 7 at block 4 is still open ([../moonshot-block4-block5/](../moonshot-block4-block5/)) |
| A 16-point identity: the transform on 16 entries below the textbook count | a kernel that could be timed | not run as a search. The ten idea waves that started at 21:00 that night, looking for linear size, found small circuits with better constants on the way (the 5/6 result, public in wht-fewer-operations, 16-point circuits among them); a search for a 16-point circuit in OpenAI's count has still not been made (that repository's `waves/LEADS.md`, lead 1.1) |
| The two network codes posted on 7 October (Braverman and He; Zhang, Li and Li) | a design principle from network coding | read that evening: neither code can be turned into a unit (argued by hand, not checked; one skeptic). Told in [../other-fields/](../other-fields/) |
| Fat labels; labels on frames that are not subspaces; tiny units with at most two additions | 1.7e-2 to 2.5e-2; 1e-2 to 2e-2 (a guess); a real-time gain of about 5.6 % per cache tile | not run in these rounds (a screen earlier on 10 October had found no fat-label family that gains) |
| Linear size | | the floor theorem below (section 1.2), and the ten idea waves (public) |

### 1.2 A floor theorem for linear size

**The question.** Could a program of the engine's kind compute the transform in a number of operations that grows only
like N? The engine charges one recursive call per paid block, so this is a question about how many blocks a unit
needs.

**The statement.** Pick two digit columns. Their *diagonal* is the data sitting where the two digits are equal. Call a
step *cheap* if it keeps "every array lies in the diagonal" true: sums of whole arrays with any coefficients,
relabelling and scaling done alike on both columns, swapping the two columns, anything on other columns. Then, if
every cheap step of a program keeps that diagonal, the paid steps satisfy

  sum over blocks of (1 − 2^(−rank)) ≥ W · (1 − 2^(−m)),

where a block of rank k transforms k of the m bits. Since a block has rank at most m − 1, a unit needs at least
W·(2^m−1)/(2^m−2) blocks: **more blocks than live arrays**. One size with two columns is enough, and coefficients,
scratch arrays and phases may be anything.

**Why it holds** (one paragraph). Feed every array data on the diagonal. The potential is the number of dimensions by
which the arrays' contents, all together, stick out of the diagonal. It starts at 0; a cheap step cannot raise it; a
paid block of rank k raises it by at most 2^m − 2^(m−k); and at the end each array holds the transform of its input,
which sticks out by 2^m − 1. Per bit, the one-bit transform keeps one of the two diagonal vectors inside and sends the
other outside, which is where the halving comes from.

**What it means.** In the Lean engine every letter acts on all columns alike, so its hypothesis holds; the engine
would certify linear size only with fewer blocks than live arrays, and that cannot happen. One block per array (about
N log log N) is out of reach too, at every fixed block size. **Status: argued by hand, not checked.** Three agents
reached the inequality by three routes; three skeptics re-derived every step; exact tests (606 random exact programs
at block sizes 1 to 3, 2,880 random letters, 7,020 letters in a widened class where a planted fault was caught) did
not break it. A Lean statement of the core lemma was written as text and never compiled.

**What it does not do.** It is weak on the exponent: at block size 110 it asks for only 1 + 7.7e-34 blocks per array,
where the published unit uses 24.6, and at block 3 it caps the saving at 6.2e-1, far above the cap of 6.28e-2 that the
14W/5 bound now gives (proved in Lean). It counts blocks, not work.

| Statement from that round | Status |
|---|---|
| The diagonal floor above | argued by hand, not checked (three agents, three skeptics) |
| Inside the Lean engine linear size and one block per array are out of reach | argued by hand, not checked; the Lean statement was read, not compiled |
| Child calls cover each array at least once, when the calls act on parts chosen by other columns | argued by hand, not checked; never run on a whole program |
| The same, when the parts are chosen by the pair's own digits ("known zeros") | **open** |
| Cost at least about N log log N in OpenAI's count, for a wide class | **open**: it needs "a sum costs N, a block a full call" for partly zero arrays |
| "Each array gets its own call" | **refuted** in the same round by a counterexample (which saves no call) |
| "Treating two columns alike is enough" | **refuted**: the right condition is "keeps the diagonal" |
| The plainest half-touch identity, transform = (multiple of the identity) + cheap · (half-size transform) · cheap, does not exist | argued by hand, not checked (one agent); its support bound was checked by computer at lengths 4 and 16 |

**The escape hatches** are the open part. A program of linear size must, for every pair of columns, contain a cheap
step that moves that pair's diagonal: a mixing step on one column only, a relabelling that mixes columns, fan-out and
fold, blocks run only where the pair's digits agree, calls to a different transform, batched calls, or a unit that
changes with N. Every such program the agents found was the textbook method in disguise. No mechanism for linear size
is in sight. The one small experiment the round proposed, a "half-touch" identity (one half-length call plus a sparse
remainder), was run later that night in the second idea wave, public in wht-fewer-operations: none of its identities
beat the textbook count (open table).

### 1.3 Tiny units without address maps (night of 10 to 11 October)

**The question.** The plan's first avenue: can 3 arrays share 8 moves at block 3 (saving 1.07e-1), or 2 arrays share
6 or 7 moves at block 4 (2.075e-1, 9.63e-2)? The night attacked it from both sides: searches for a unit, and certificates
that would exclude one.

**The sub-game.** Most of the night's results hold in the *sub-game without address maps*: additions with rational
coefficients, scratch arrays, shifts and every phase letter are free, but no step may rename the addresses of a block.
At two digit columns an array at block 3 is an 8 × 8 matrix, written on 64 "labels", and the free steps of the
sub-game cannot move a label out of its group (15 groups in all). That turns the sub-game at two columns into a
routing problem, and gives a "joint test": for every set of label groups, the moves must be able to carry into it at
least what the W targets require.

What came out, in order of the night:

- **The joint test excludes 3 in 8 at block 3 in the sub-game** (argued by hand, not checked; checked by computer by
  five separate programs). Every one of the 5,764,801 ordered patterns of 8 directions has a cut certificate, and a
  ten-line hand argument with six fixed inequalities covers all of them; the same test gives at least 13W/5 moves for
  W arrays and excludes a saving of 1e-1 for 1 to 5, 7 and 13 arrays. In the full game, with address maps, the method
  gives only 7W/3 and excludes nothing.
- **Three digit columns close block 3 in the sub-game: W arrays need at least 3W moves**, which is the plain cost
  (argued by hand, not checked; every step has a hand proof, read line by line by a second agent; the ranks were also
  computed exactly by three programs). The certificate is the span of the antisymmetric combinations of the 24 "face
  triples" (three addresses on one face of the cube): every direction pushes at most 8 dimensions out of it, and the
  transform pushes out all 24. It also holds when the address maps all permute one common basis, so a saving unit at
  block 3 needs, for example, an address map of order 4 or 7. **Later:** this same "face certificate" became half of
  Theorem C3 (15 October; told in [../reed-muller-and-one-move/](../reed-muller-and-one-move/)): without phase
  letters, whatever the address maps, nothing saves at block 3 either. So a block-3 saving needs address maps and
  phase letters together.
- **Block 4 in the sub-game needs at least 3728/1219 = 3.058 moves per array** (argued by hand, not checked; the
  certificate was recomputed exactly by a second agent). This excludes 2 in 6, 3 in 8 and 3 in 9 there. **Superseded**
  by the hook bound 41/13 = 3.154 (13 October), which holds in the whole engine game.
- **2 in 7 at block 4 is excluded in the sub-game** (11 October morning; argued by hand, not checked). Of the 116,280
  multisets of 7 directions only 120 pass every order-free test, and every ordering of them fails an order-aware form
  of the same lemma, derived twice; a second agent recomputed the certificates in whole numbers. This does not touch
  the full game: 2 in 7 with address maps is still open (see
  [../moonshot-block4-block5/](../moonshot-block4-block5/)).
- **At every block size, units whose directions all lie in one basis need m·W moves** (argued by hand, not checked).
- **Searches in the full game found nothing**, but each covered a thin slice: one family of 8-move schemes for 3 in 8
  (two of its stages complete, a third stopped at about a quarter), the same with one relative address map taken from
  the six powers of one map of order 7 (complete: about 645 million nodes), and a staged search for 2 in 7 at block 4
  (of its nine stages, one complete, three at 9 of 12 shares and one at 3 of 12; the rest were not reached, and the
  cloud machine, a Linux machine of Google's Colab service, that ran the other shares expired). Checked by computer, one program each, for those slices only. For 3 in
  8 they are now moot: 14W/5 excludes it.
- **Symmetries of the sub-game:** only the 6 coordinate permutations of the 168 address maps keep the transform, so
  the 7^8 patterns fall by a factor of at most 12, not 168; running a word backwards with transposed additions is a
  symmetry. Checked by computer (exact replays), and re-run by a checking agent.

Two cautions that apply to all of it. Every lower bound of that night rests on one lemma of the Theorem R kind (three
written proofs by agents), and that lemma was never tested on a real unit that shares work between arrays, because
none is known at block 3 or 4; random programs stood in for it. And a statement about the sub-game says nothing about
units that use an address map, which is where every known saver lives.

---

## Part 2. 16 October: one serious attempt at a new mechanism

### 2.1 The set-up

By 16 October the small shapes had been searched and bounded from many sides and no unit had appeared. I asked the
agents to stop stepping through small searches and to think hard instead: no big coding projects, small exact
experiments only to test an idea on its smallest case, and every method calibrated on a known fact before it speaks.
The coordinating agent wrote six seed observations:

1. the best saving should be the same for every kernel width, so a unit found in any "playground" (another kernel,
   another set of free steps) would transfer;
2. OpenAI's existence proof ("Finite tensor savings and exact Fourier circuits", 25 September 2026) contains a second
   mechanism, never built as a unit, with a local margin of 25 %;
3. without phases the game becomes ordered products of reflections, and the block-3 frame group is the Weyl group of
   E7 up to sign;
4. OpenAI's model allows free steps that none of our games use (any precomputed diagonal, eighth-root phases);
5. streaming units, with state carried from one period to the next, might collapse the unknowns;
6. the Schur-Weyl picture of the transform (computed by hand: negative).

Six tracks took these up, one more agent quantified OpenAI's existence proof, a recombiner put the tracks together,
and two skeptics attacked the recombiner's two best leads.

### 2.2 What each track found

**OpenAI's second mechanism is a price-world artefact** (two tracks; argued by hand, not checked; no skeptic).
- The whole 25 % margin comes from one object, the "dense pair" B(u)⊗B(u), priced at 6 instead of 8; every other law
  in the chain (feedback, specialisation, pivots) only converts a cost from one description to another. Priced at 2
  per pair, the chain is a tie. The feedback law itself implies price ≤ number of nonzeros for every invertible
  matrix, which caps the margin at 1 − π/2 + O(1/d) for dense-pair price π. Checked by computer: every identity of the
  reflection argument at d = 1 to 4.
- The price world is not reality: the proof assumes no win exists, and in reality wins exist. The win it forces is
  astronomically small (below 1e-23 even with the most generous reading of the Lean packing; a tower of exponentials
  with the laws' own run lengths). Even if realising the laws cost nothing, a fixed-matrix version could certify at
  most 1.58e-3 (with OpenAI's own count of the chain's extra pairs) to 1.97e-3 (with each of those extra pairs priced
  like a dense pair). These are prices
  computed by a program under an idealisation, not proved; they are of the same size as the records already built by
  explicit units (the best Lean-checked saving that day was 1.6256629e-3), and the real object is far smaller.
- It cannot be made right for every number of digit columns f: feeding back per digit column needs state that grows
  with f (the McMillan degree of (I + zJ)^(⊗f) is exactly f·4^f, computed for f = 1 to 4), and closing one loop around
  the whole power gives the wrong matrix (checked exactly at f = 2).
- Our engine has no analogue of OpenAI's "sparse but costly" object: on one array a Clifford operator costs exactly
  log2 of its nonzeros per row (the Bruhat decomposition; checked by computer on 6,400 random letters), so sparse
  means exponentially cheap here.

**Streaming cells add nothing** (argued by hand, not checked; no skeptic). A cell that carries state from one period
to the next, run for T periods, is a block unit with p·T arrays, so every certificate binds on its rate. If no paid
move sits on a feedback loop, evaluating the delay at one generic number turns the cell into a block program with the
same counts. The smallest cells are "stage programs" in which arrays are owed partial transforms that chain across
periods; that points to one new small target, (2,2,4) in 7 at block 4 (9.63e-2 if streamed). In the commutative slice
one stream never saves (checked by computer on all 37,059 direction sets below m moves at blocks 3 to 5).

**The Coxeter carousel, and the busy-register wall** (the step is checked by computer; the wall is argued by hand, not
checked; no skeptic). Borrowed from coding theory (the accumulator of repeat-accumulate codes) and Lie theory (the
Coxeter element): "one paid move along direction 1, then relabel the directions cyclically" is a step c with c^m equal
to the whole transform, exactly, for m = 2 to 8. A carousel of r carrier arrays then advances m inputs in flight with
one move each round, the most packing any unit can have, and computes the transform exactly at the plain rate. The
wall: if every register gets a paid move in every layer, the trivial representation turns the cell into a scalar
system whose Hankel matrix must have rank p(δ + 1) yet factors through the r registers, so r ≥ p(δ + 1) and the rate
is at least m moves per input. So packing alone buys nothing; a saver must leave some register idle, as every known
saver does.

**Root systems, and why E8 is the best centre** (no skeptic; the status of each item is given with it).
- The block-3 frame game is the Weyl group of E7 up to sign: the 135 frames are the 135 sets of seven orthogonal E7
  root lines, and each paid move is a 4 × 4 block "J − 2I" in frame coordinates, which is OpenAI's dense pair (checked
  by computer, one program).
- A crossing centre of the record units is a linear index code. Its copy ratio is (decoder rank) × (window) / (ports):
  E8 0.600, E7 0.778, E6 1.167 (checked by computer, one program). E8 and E7 are perfect codes (Lovász theta =
  independence number = rank), so by Haemers' bound no decoder of any kind, with these labels, brings E6, D_n or A_n
  below E8 (argued by hand, not checked). Used as a whole unit it is still crossing slack, capped near one free direction per array.
- Without scratch, programs whose gates are reflections obey short lemmas: the rank-one matrices z·zᵀ of the gate
  directions must add up, modulo 2, to the sum of the targets (so no "plain minus one move" shape exists without
  scratch), a form may never switch on exactly one gate that is not a sign flip, and two layers of commuting gates
  never save (argued by hand, not checked). The track then reduced the block-5 "7-array problem" to a block-4 saving
  core; **that reduction was refuted** by the second skeptic (section 2.4).

**A tight 5 in 14 lives inside S8** (checked by computer; the structure repeated by a skeptic with its own code; the
symmetry screen by one program). An earlier lemma, which no skeptic has read, says every term of an exactly tight
block-3 unit sits on 30 of the 135 frames. These 30 are exactly the maximal totally singular subspaces of the quadric
x·z = 0: the points and planes of the projective space PG(3,2), by the Klein correspondence. Their symmetry group is
O6+(2), which is the symmetric group S8; in S8 the paid move is a transposition, the target is a product of three
disjoint transpositions, and the free steps that respect the quadric form AGL(3,2). This gives a closed sub-game 36
times smaller than the frame game whose programs are genuine units. The symmetry answer is negative: a tight 5 in 14
can carry a symmetry group of order at most 24; no group of order 24 or 12 survives the screen; among groups of order
6 or more, only one shape with D8 symmetry (which needs S-gate sectors outside the pure sub-game) and six shapes with
S3 symmetry are left. If it exists it is nearly asymmetric. One more fact: the one-piece relaxation (the
56-dimensional irreducible) behind part of the numerical evidence for 5 in 14 is already satisfied by a scheme whose
moves act on raw inputs only, which no real unit can do, so that relaxation is blind to depth.

**Other kernels and other free steps** (the "playgrounds"; argued by hand, not checked; no skeptic).
- **The saving is the same for every kernel width.** OpenAI's generation lemma gives all nonmonomial q × q kernels one
  saving; the Kronecker power of H is a 2^b × 2^b kernel; adding a 1 × 1 block never lowers it. So s(q) = s(2) for
  every q. A unit found for any legal kernel transfers to the Walsh-Hadamard exponent, and the Fourier exponent is at
  least that along OpenAI's lengths. That the Fourier saving *equals* it does not follow from these steps (see section
  2.5).
- **Block 2 is dead in every playground**: Theorem 2W and the two-layer bound W + m hold for every kernel once the
  free steps act the same way on every column (in their group-algebra form, which does not need the move to be
  "spread").
- **The hook applies to every kernel.** Every dense 2 × 2 kernel gets exactly 14/5 moves per array at block 3 (the cap
  80 argued by hand for all; the demand 224 computed for six cross-ratios). The single exception is the **triangular
  class** (the shear [[1,0],[1,1]], the subset-sum kernel, the binomial involution [[1,0],[1,−1]]): 173/80 = 2.1625,
  exact. That would leave room for 2 in 5 (1.66e-1), 3 in 7 and 4 in 9 at block 3.
- **But that room is uncalibrated.** At block 2 the same hook gives the triangular class only 7/4, where 2 holds by the argument above,
  and every certificate of this kind that was tried stays at or below 7/4 (all pattern subspaces at four columns:
  1.673). So the 2.1625 is a fact about the certificate, not evidence for a unit. With only diagonals free, a
  triangular unit never saves at any block (argued by hand, not checked), so a triangular saver would need
  permutations that break the subset order.
- **Non-Clifford phases open nothing at blocks 2 to 4.** Every hook is kept by every diagonal and every permutation of
  column states, so 2, 14/5 and 41/13 hold with eighth-root phases, controlled phases of any angle and any scaling.
- The qutrit kernel is uncharted at block 3: the hooks within reach are weak, and its finite frame game (1,120 frames)
  has never been searched.

### 2.3 What only appeared when the tracks were put together

**No exactly tight unit at block 4 in the engine game.** An exactly tight unit at one f is tight at every f, and in
the engine game, where the group is finite, it is then tight in every irreducible piece; in the regular representation
that means every single term of every move must pass a termwise filter on its own. An earlier block-4 test had found
that sampled terms all fail, but left an escape through terms that only pass in combination; the termwise statement
closes that escape. The recombiner enumerated every engine position at block 4 by depth (reproducing the distance
counts 1, 30, 280 of the frame graph as a calibration) and ran the earlier test's own filter, exactly modulo the prime
33,554,393: **0 of 30 positions at depth 1 pass, and 0 of 280 at depth 2**, while the controls (the identity, the
target, the target followed by free letters) pass. So the first paid step of any program already breaks exact
tightness. This closes every schedule exactly on the 41/13 bound: 13 in 41 (1.71e-1), the grouped 2 in 7 with blocks
2,2,2,1 (1.675e-1), and the others in the table below. *Status:* the computation is checked by computer (one program,
not repeated by anyone); the argument that turns it into an exclusion is argued by hand, not checked; the lemma it
rests on has no skeptic. The schedules with rank-2 blocks are excluded only with high probability, because their
tightness uses a rank-2 cap computed at two primes but not certified. Near-tight schedules, such as the plain 2 in 7
with rank-1 moves, are not covered.

| Exactly tight schedule at block 4 (W arrays; moves) | Would have saved |
|---|---|
| 13 in 41 (and 26 in 82, ...), single-direction moves only | 1.71e-1 |
| 5 arrays: 14 single moves and 1 rank-2 block | 1.71e-1 |
| 7 arrays: 15 single moves and 4 rank-2 blocks | 1.70e-1 |
| 2 arrays with blocks 2,2,2,1 | 1.675e-1 |
| 4 arrays: 2 single moves and 6 rank-2 blocks; 6 arrays: 3 and 9 | 1.67e-1 |

**Smaller results of the recombination:**
- the stage variants of 5 in 14 that the streaming track tabulated are all excluded by an earlier theorem (no tight
  block-3 unit gives any array a two-row target);
- (2,2,4) in 7 at block 4 is impossible without scratch in the multi-target slice (argued by hand, not checked; the
  second skeptic re-derived it with a shorter proof: its three targets sum to zero, so its 7 gates form a 3-flat, and
  a target of rank 4 cannot live there);
- the S8 sub-game of the tight 5 in 14 is the A7 reflection subgroup of the E7 picture (the counts agree, and the
  quadric is forced by the start and target frames);
- four tracks agree that OpenAI's second mechanism and streaming are closed.

The recombiner ranked three leads: **D1**, the tight block-3 family (5 in 14, and a 7-fold symmetric 15 in 42) inside
the S8 sub-game, worth 6.28e-2; **D2**, "targets on a subspace" with block-4 cores first, worth 4.66e-2 and then
7.54e-2; **D3**, the triangular kernel, worth up to 1.66e-1 but with no cheap decisive test.

### 2.4 The two skeptics

**D1 survives, narrowed** (skeptic 1). The question "is there an exactly hook-tight block-3 engine unit?" is open;
every structural fact D1 states checks out, recomputed from scratch with the skeptic's own code (the 30 frames, their
7-regular bipartite step graph, the group of order 40,320, the 8 sectors). One small new fact: every step between the
30 frames changes the depth by exactly 1, so a tight unit has no sideways steps. But the plan as written does not
survive:
- its first test was to calibrate the only numerical evidence (a descent that "reached 2.2e-8 and was still falling"
  on 5 in 14; its log shows 5.4e-9 at iteration 1500, when it was stopped). One control had already been run and
  misfired: its "plateau" was still falling when the job was stopped. The exact half of the other control comes out the uninformative way: two moves already span the whole
  relaxation algebra;
- the control that matters is missing, and everything predicts that it fails: shapes that an earlier theorem excludes
  for a junk reason (4+1 in 13, 3+2 in 12, 2+3 in 11) pass every certificate tried on these two pieces, just as
  5 in 14 does;
- the relaxation is underdetermined (3.5 unknowns per equation) where the real problem is overdetermined (0.34 in the
  sub-game, 0.076 in the frame game), and the plain 1 in 3, which exists, shows the same slow tail. "Still falling" is
  that optimiser's signature on anything that exists;
- the proposed decider is a search, which can find a unit but never refuse one; and "36 times smaller" holds only for
  the pure sub-game, one of the 8 sectors a tight unit may use.

If a tight unit exists it saves exactly 6.28e-2, with no hidden dilution. The skeptic's guess for its chance is "well
under 1 in 10", and it advised planning with zero until an exact sub-instance is decided. The smallest instance that
can be refused is the trivial-character 6-move problem of the 7-fold 15 in 42 (open table).

**D2 is refuted** (skeptic 2). The step that fails is "without scratch and without interleaving, the 7-array problem
needs a block-4 core". The core must deliver a total target *rank* of at least 19, not a total *length*; every core of
rank 19 or more is already excluded by Theorem R with the block-4 hook values (414,688 needed against 411,840
available); and the cores D2 listed have total rank 16 (18 for the mixed ones), short of 19, so they are not cores
of the 7-array problem at all. The same
count shows that the outside gates of such a 7 in 31 must span all five directions, so the problem is genuinely at
block 5. The listed block-4 cores cannot exist without scratch either: their four targets would have to sum to zero,
and none of the 1,612 surviving classes does; all 27,256 zero-sum target sets fail the earlier exact tests. Three of
D2's claimed consequences do not follow. *Status:* refuted, by the skeptic's argument and exact computations (one
skeptic, not repeated). The two new negative facts (outside gates must span F_2^5; a 4 in 15 on the 15 directions of
F_2^4 needs scratch and a nonzero target sum) are argued by hand and checked by computer, with one skeptic. The
7-array problem itself stays open; it is told in [../reed-muller-and-one-move/](../reed-muller-and-one-move/).

Neither skeptic examined D3 or the block-4 exclusion of section 2.3; those have no skeptic.

### 2.5 The Fourier transform at every length

The playground track showed that the Fourier saving is at least the Kronecker-power saving along OpenAI's lengths, and
that equality does not follow from its steps. At the end of the round the coordinating agent gave a short argument
that the Fourier saving at *every* length equals the Kronecker-power saving, through Rader's algorithm over finite
fields F_(p^k). The argument is not written out in these notes, and no skeptic has read it. *Status:* argued by hand,
not checked; I list writing it down as an open item.

---

## What is still open, and the first experiment for each

None of these has been started.

| Lead | Worth if it exists | First experiment |
|---|---|---|
| A hook-tight block-3 unit (5 in 14, or the 7-fold 15 in 42) | 6.28e-2 | Decide the trivial-character 6-move problem of the 7-fold 15 in 42 exactly: deliver 2 full arrays and 1 trivial piece in 6 moves, every move gaining exactly 80, every term on the 30 singular frames. Start from the case where every move touches one array, and show that the two zero-divisor transfers cannot both keep +80. A refusal kills that instance; a solution is a unit. If machine time goes to the old descent anyway, run the one control that could rescue it, 4+1 in 13 on the same two pieces |
| The triangular kernel at block 3 (2 in 5, 3 in 7, 4 in 9) | up to 1.66e-1 | First find a lower-bound method that gives 2 moves per array at block 2 for [[1,0],[1,−1]] (every subspace certificate tried gives 7/4); the group-algebra rank arguments are the only ones known to. Then apply it at block 3: at least 14/5 closes the playground. Only then search for 2 in 5 with AGL(3,2) maps, bit flips and continuous diagonals |
| Near-tight block-4 schedules (the plain 2 in 7; 6 in 19 at 1.69e-1) | 9.63e-2 to 1.69e-1 | A "termwise defect budget": charge each failing term's defect against the schedule's slack. Every depth-1 term was measured at a defect of at least 3,920, above the 2,112 slack of the 1.6e-1-class schedules; decide whether defects add up per lineage |
| The 7-array problem at block 5 (7 in 31) | 7.54e-2 | No block-4 shortcut exists (section 2.4). With scratch, run the Reed-Muller round's exact route (899 unknowns per time order), starting with the local Gröbner test on the move sets of the weight-8 forms; see [../reed-muller-and-one-move/](../reed-muller-and-one-move/) |
| An idle-register carousel at block 6 | 1.02e-1 | In the no-phase game with the cyclic address group of order 6: five busy carriers with the Coxeter step plus one idle delay register, 2^21 forms (about 350,000 orbits), lags up to about 12. Calibrate on the rate-6 carousel (must pass) and a deleted-move control (must fail) first. A job for a large cloud machine, not a laptop |
| (2,2,4) in 7 at block 4, with scratch or address maps in the middle | 9.63e-2 streamed | Calibrate first (find the plain (2,2,4) in 8; refuse (2,2) in 3 at block 4 and (2,3) in 4 at block 3), then extend the exact six-move frame solver of the block-4 rounds to three inputs with targets (K, B, C) |
| The escape hatches of the floor theorem | linear size | The half-touch experiment that the floor round proposed (one half-length call plus sparse steps) was run the same night in the second idea wave, public in wht-fewer-operations: exact one-call identities at 4 to 64 points, none below the textbook count. What it left open is one half-length plus one quarter-length call, where the known obstruction stops; the only construction of that shape found so far has rows that grow like N^(1/4). First experiment, the one that wave's note names: for that construction, compute the input side B by exact linear algebra at 16 and 64 points and count its nonzeros, which decides whether that construction gives such an identity at all (that repository's `waves/LEADS.md`, lead 2.10) |
| The diagonal floor lemma in Lean | (a proof) | Compile the two-theorem statement written in the floor round (a rise lemma by dimension counting, and the rank-k value 2^m − 2^(m−k)), then a bridge from the engine's letters |
| The Fourier exponent at every length | (a theorem) | Write out the Rader argument of section 2.5 and give it to a skeptic |
| The qutrit kernel at block 3 | unknown | Port the exact frame-game tools to the 1,120 Lagrangians of F_3^6 and compute frame certificates there; hooks are too weak for q = 3 |

## Method lessons

- **Calibrate before a method speaks.** Two numerical searches of 10 October failed their own controls (one found a
  known 9-move solution once in 11 tries; one sweep found its planted answer once in 8 starts) and were stopped or
  discounted. On 16 October the descent offered as evidence for 5 in 14 turned out to behave the same way on schemes
  that exist and on the control it was meant to beat.
- **A failed search refuses nothing.** A decider must be able to say no; a plan whose "decisive test" can only find a
  unit is not decisive.
- **Read the scope of every exclusion.** Most of the night of 10 October holds only without address maps, and every
  known saver uses them. Several sentences in the agents' notes said "excluded" without that scope, and the checking
  agents had to restore it.
- **The later verdict wins.** The reduction of the 7-array problem to block-4 cores was stated as proved and then
  refuted the same afternoon (rank, not length); the tiny-unit bounds of 10 October were overtaken within days by
  bounds in the whole engine game.

## Programs and data

The programs are in [programs/](programs/), grouped by finding; their saved outputs are in [data/](data/) under the
same folder names, and each program folder's README says what each file shows and how it was called. Most need only
Python 3.9 or later with numpy (some also scipy); the two C programs need a C compiler. Run each command from inside
its program folder; a program that writes its own output file writes it there, to compare with the copy in `data/`.

| Folder | What it shows | Re-run for this revision |
|---|---|---|
| [openai-mechanism-price-world](programs/openai-mechanism-price-world/) | the dense-pair accounting, p ≤ nnz, the McMillan degrees, the Clifford cost check, the zero-overhead ceilings | all six re-run: same output |
| [streaming-cells](programs/streaming-cells/) | the streaming evaluator with plain cells and controls; the stage shapes and their hook room | both re-run: same output |
| [coxeter-carousel](programs/coxeter-carousel/) | the Coxeter step, the carousel, the border check, the busy-register wall on scalar shadows, the carry identities | all six re-run: same output (`busy_wall.py` took about 10 minutes at background priority) |
| [root-systems-and-index-codes](programs/root-systems-and-index-codes/) | the E7 dictionary, the index-code copy ratios, the reflection-form probes | `s1`, `s2`, and `s4` with `pos 1` and `neg 1` re-run: same output; not re-run: `s3_subspace.py` (6.5 minutes), the `s4` probe seeds, `s5`, `s6` |
| [tight-5-in-14-subgame](programs/tight-5-in-14-subgame/) | the hook calibration, the 30 singular frames and S8, the symmetry screens, the depth-1 staircase | all 19 programs re-run in order: same output for the 18 with saved outputs (`c11_staircase.py` differs only in the source lines that numpy's overflow warnings print) |
| [kernel-playgrounds](programs/kernel-playgrounds/) | the hook for other kernels (14/5 for dense kernels, 173/80 for the triangular class), the block-2 calibration | `exact_hook.py`, `hookk.py 3 4`, `allinv.py`, `patterns.py 2 3` and `patterns.py 2 4` re-run: same output (`patterns.py 2 4` took 16 minutes at background priority); the other scans not re-run |
| [no-tight-unit-at-block-4](programs/no-tight-unit-at-block-4/) | the block-4 positions by depth and the termwise filter (0 of 30, 0 of 280); the smallest dependency sets; how to rebuild the filter's data file | `b4_admissible.py` at depth 1 and in both depth-2 halves, `enum_test.py 3` and `dep_supports.py` re-run: same results (0 of 30, 0 of 140 and 0 of 140; the depth-1 output has one extra progress line that the saved run, made with an earlier copy of the program, did not print); the data file E rebuilt from `make-E/`, byte-identical |
| [skeptic-tight-block3-family](programs/skeptic-tight-block3-family/) | the first skeptic's depth cuts, labelled certificates, singular frames and type graphs | all four programs re-run (`s1_depth_cuts.py` in both modes): same output |
| [skeptic-targets-on-a-subspace](programs/skeptic-targets-on-a-subspace/) | the second skeptic's zero-sum filters, harness calibration and the Theorem R table for cores | all five re-run: same output |
| [floor-theorem](programs/floor-theorem/) | the floor round's exact tests by the three authors (`f-*`) and the three skeptics (`g-*`) | 22 runs re-run with the same output; `span_clo.py 3 2 65537` was stopped after 5 minutes before its first result, as the original run was; not re-run: `span.py`, `g2.py` |
| [tiny-units-without-address-maps](programs/tiny-units-without-address-maps/) | the hand transcription of the engine game (`regmachine.py`), the joint test and its certificates, the face certificate and the block-4 certificate, the symmetry tests | `model_cert.py`, `hand.py`, `static.py 3 8`, `staticlp.py 9`, `check_cert.py`, `sk_six.py`, `facet.py` for six (m, f) pairs, `sk_general.py` in its three modes, `sk_cert.py`, `symtest.py` and `gaugetest.py` re-run: same output (`facet.py 2 3` stops with a division by zero, because its subspace is empty, exactly as the saved run did); `certfind.c` and `certcheck.c` not re-run |

For this revision I re-ran the programs marked above from a clean copy of this folder: 98 runs. Every run with a
saved output in `data/` gave that output again, except one long run that was stopped before its first result, as the
original run had been; three runs have no saved output. The list, with times, is
[data/rerun-19-october.txt](data/rerun-19-october.txt). I also rebuilt the block-4 filter's data file, byte for byte.

Not published: the search programs of the night of 10 to 11 October and their outputs (about 560 MB, mostly search
states and logs), the 160 MB of binary cut certificates (rebuilt by `certfind.c`, see its folder), and the 53 MB data
file of the block-4 filter (rebuilt by the programs in `no-tight-unit-at-block-4/make-E/`).

## Credits

I did this work with a team of AI agents (Claude) that I directed: I set the goals, chose which ideas to pursue and
decided at each stage what to do next; the agents read the papers, wrote the proofs and programs, and tried to refute
one another's work.

Outside work used or cited:
- OpenAI, "Finite tensor savings and exact Fourier circuits" (preprint, 25 September 2026) and its Lean development:
  the generation, amplification and transfer lemmas, the price laws, the dense-pair lemma and the reflection argument
  that section 2.2 analyses. Nothing of OpenAI's is copied here.
- Braverman and He (arXiv:2610.10108) and Zhang, B. Li and Z. Li (arXiv:2610.09367), the two network-coding preprints
  of 7 October that the plan of section 1.1 proposed to read.
- Standard mathematics, cited from memory by the agents: Haemers' minrank bound and the Lovász theta function (index
  coding); the classification of line systems at 60 and 90 degrees (Cameron, Goethals, Seidel and Shult); the Klein
  correspondence; the Bruhat decomposition of the symplectic group; Coxeter elements (Bourbaki); repeat-accumulate
  codes (Divsalar, Jin and McEliece, 1998); Hamming-weight phasing (Gidney, 2018); Kalman's realisation theory
  (McMillan degree and Hankel rank); Rader's algorithm (C. M. Rader, 1968); Young's orthogonal form for the symmetric
  group; Agler's transfer-function realisations; and, as comparisons, the rigidity work of Alman with Rao and with
  Williams and the multiplication algorithm of Harvey and van der Hoeven.
- The engine, the games, the hook certificates and the earlier lemmas this section builds on are those of this
  repository and its earlier revisions and of the other sections of research/. The matrices of the small irreducible
  representations of Sp(6,2) used by the first skeptic were computed in the theory round of 14 October
  ([../theory-of-theories-and-logic/](../theory-of-theories-and-logic/)); the second skeptic uses the Reed-Muller
  round's programs, which are published in [../reed-muller-and-one-move/](../reed-muller-and-one-move/).

## Sources

All of this section comes from my working notes, kept in a private archive; paths below are relative to its `checks/`
folder. Where a track's note and a later skeptic disagree, this section follows the skeptic: the reduction of the
7-array problem to a block-4 core (track e7-rm, section 4) gives way to skeptic-2's verdict, the recombiner's figure
"2.2e-8" for the 5-in-14 descent gives way to skeptic-1's 5.4e-9 at iteration 1500, and the 10-11 October bounds give
way to the later theorems named in section 1.3 (archive entries L7, L9, L22, L24).

**10 October evening.**
- Plan: `wht37/avenues/NOTE.md` (the plan, by agent v-plan), `FIRST-EXPERIMENTS.md` (three briefs), the five lens
  reports `v-imports.md`, `v-levers.md`, `v-outside.md`, `v-realworld.md`, `v-tiny.md` and their skeptics
  `k-imports.md`, `k-levers.md`, `k-outside.md`, `k-realworld.md`, `k-tiny.md`, and `v-plan.md`.
- Network codes: `wht38/codes/c-read.md`, `c-skeptic.md` (told in other-fields).
- Floor theorem: `wht38/floor/NOTE.md`; authors `f-engine.md`, `f-model.md`, `f-break.md`, `f-note.md`; skeptics
  `g-engine.md`, `g-model.md`, `g-break.md`. Program folders here: `floor-theorem/f-*` and `floor-theorem/g-*` keep
  those labels.
- Tiny units: `wht38/tiny/STATUS.md` (dated summaries from 20:40 on 10 October to 09:42 on 11 October),
  `h-skeptic.md`, `SYMMETRY.md`, `SYMMETRY-STATUS.txt`, `y-skeptic.txt`; the routing filter `filter/NOTE.md`,
  `filter/argue.md`; the joint test `filter/joint/NOTE.md` with `cert.md`, `proof.md`, `full3.md`, `block4.md`,
  `ratio.md`, `skeptic.md`; three columns `next/NOTE.md` with `cutf3.md`, `netcode.md`, `amap.md`, `skeptic.md`; the
  sub-game closure `close/NOTE.md` with `general.md`, `wedge4.md`, `skeptic.md`; 2 in 7 at block 4 `b4w2k7/NOTE.md`
  with `patterns.md`, `skeptic.md`, `full.md` and the search logs in `b4w2k7/full/`. The engine transcription
  `regmachine.py` comes from `wht35/wider/w-engine/`.
- The superseding results: archive `RESULTS.md` entries L7 (14W/5, Lean), L8, L9 (41/13), L22 (Theorem C3), L24 (block
  4 with maps and phases at the ends). The fat-label screen of earlier on 10 October: entry R30.
- The half-touch experiment and the 16-point lead: the public repository wht-fewer-operations, `waves/wave-02/e-7.md`,
  `waves/wave-02/b-2.md`, `waves/wave-02/NOTE.md` (idea D) and `waves/LEADS.md` (leads 1.1 and 2.10).

**16 October.** `wht44/ideas/`: `BRIEF.md` (the six seeds), the track notes and the recombination `RECOMBINE.md`, and
the two skeptics. Track labels and the folders here:

| Track label | What | Folder here |
|---|---|---|
| unwind, quantify | OpenAI's second mechanism; the existence proof as a number | `openai-mechanism-price-world/` |
| stream | streaming cells | `streaming-cells/` |
| wildcard | the Coxeter carousel and the busy-register wall | `coxeter-carousel/` |
| e7-rm | root systems, index codes, reflection forms | `root-systems-and-index-codes/` |
| tight514 | the S8 sub-game of 5 in 14 | `tight-5-in-14-subgame/` |
| playgrounds | other kernels and free steps | `kernel-playgrounds/` |
| recombine | the block-4 check, the dependency sets | `no-tight-unit-at-block-4/` (with `make-E/` from `wht43/wave16/tb4-tight/`) |
| skeptic-1 | verdict on D1 | `skeptic-tight-block3-family/` (input `irreps_small.npz` from `wht43/meta2/a2-gct/a2gct-mixed/`) |
| skeptic-2 | verdict on D2 | `skeptic-targets-on-a-subspace/` (modules from `wht43/rmfam/rm2-phase/`) |

Earlier notes the 16 October round leans on: `wht43/px/px-t514/PROOF.md` (Lemmas A to D),
`wht43/wave16/tb4-note/NOTE.md` and `tb4-tight/STATUS.txt`, `wht43/tier2/t2-block3/PROOF.md`,
`wht43/meta2/a2-gct/STATUS.txt` and its logs, `wht43/rmfam/rm2-phase/PROOF.md`, `wht43/anatomy/an-openai/NOTE.md`;
archive entries L1, L4, L5, L7, L9, L11, L12, L21 to L26, M1, M18, M23, M24, M27, N1, N3, N14, N17, R31 to R37, X7,
X8, X36 to X44. The verdicts of the day are summarised in the state block of 16 October, 14:05 EDT, of
`wht33/HANDOFF.md`; the Rader argument of section 2.5 is recorded only there.

The programs were copied from the folders above with private paths replaced by relative ones; a few path lines were
changed so that each program finds its sibling files and its data in this repository (`E_DIR`, `IRREPS_NPZ`,
`RM_PHASE_DIR`, `RM_PHASE_DATA` and `CERT_JSON` override the defaults). The re-runs were made on 19 October 2026 on an Apple M4 Mac
mini, at background priority under a memory cap.
