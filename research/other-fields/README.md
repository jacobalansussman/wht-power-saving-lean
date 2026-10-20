# Famous problems, network coding and quantum computing

This section covers four research rounds, from 10 to 15 October 2026, that asked how my problem relates to other
fields. Is it a famous hard problem in disguise? What does it look like to someone who works on quantum circuits, or
on network coding? Can the tools of those fields prove anything here, and can anything here help them? I wanted the
answers for two reasons: so that I would not spend effort on something mathematicians already know to be hopeless,
and because parallels between fields are useful to both communities.

I did this work with a team of AI agents (Claude) that I directed. I set the goals, chose which ideas to pursue and
decided at each stage what to do next. The agents ran the searches, wrote the proofs and programs, and tried to
refute one another's work. Nothing found in these rounds is checked in Lean; the one Lean-proved fact used here, the
block-3 bound at the end of section 10, comes from another section. No person has read any of these arguments.
Where I say "a skeptic", I mean one independent AI agent that tried to break a result with its own reasoning and its
own code; it is not a human referee. Most results here had one skeptic; a few had two or three; some had none, and I
say so.

**In one paragraph.** My problem is not another famous problem in disguise: it is a famous problem in its own right,
the old question whether the fast Fourier transform's N log N operations can be beaten, which OpenAI's preprints of
September 2026 answered with yes. What is left, how large the saving can be, has no famous name that the surveys
found. Some reshapings are exact, but all of them are known: the Fourier saving is the saving for multiplying
polynomials, and the Walsh-Hadamard saving is the saving for XOR convolution and for every other Kronecker-power
transform. One half of the problem is famous-hard: proving that the saving can never reach 1 would solve Valiant's
old problem of a superlinear lower bound for an explicit linear map. The half I work on, building units and proving bounds in restricted games,
meets no famous barrier; each open unit question is a finite, decidable algebra problem. In quantum language an
array's coordinate frame is a stabiliser state, the free steps are the Hadamard-free Clifford group and a paid move
costs one Hadamard; every fact in that dictionary is textbook, and my earlier belief that the link was new was wrong. The
dictionary still did work: it closed one route (labels on frames outside the subspace frames) and gave a normal form
and exact counting tools. In network-coding language a unit whose additions only join arrays at the same frame is a
linear network code on the frame graph; entropy inequalities, at every size small enough to solve, prove no more than
the plain cut bound there, far below my certificates, though cuts alone exclude most three-array layouts. Finally,
neither of the two network codes posted on 7 October gives a saving unit: one never divides by 2, which every saving
unit must; the other lives over the field with 9 elements and loses almost all its free symbols when its rule is
read over the Gaussian integers (one kind of placement, at block 3, is not ruled out, but it would first need a
rational version of the code, which nobody has). No saving moved because of any of this.

## Contents

- [The question](#the-question)
- [Words used here](#words-used-here)
- [What I tried](#what-i-tried)
- [What I found](#what-i-found)
  1. [It is a famous problem, not another one in disguise](#1-it-is-a-famous-problem-not-another-one-in-disguise)
  2. [Exact reshapings, all of them known](#2-exact-reshapings-all-of-them-known)
  3. [Which half is famous-hard](#3-which-half-is-famous-hard)
  4. [The network coding conjecture and the saving](#4-the-network-coding-conjecture-and-the-saving)
  5. [What only looks the same](#5-what-only-looks-the-same)
  6. [The quantum dictionary](#6-the-quantum-dictionary)
  7. [What the dictionary did for the problem](#7-what-the-dictionary-did-for-the-problem)
  8. [Common-frame units are network codes, and entropy bounds stop at the cut bound](#8-common-frame-units-are-network-codes-and-entropy-bounds-stop-at-the-cut-bound)
  9. [Few arrays: when cuts alone exclude a saver](#9-few-arrays-when-cuts-alone-exclude-a-saver)
  10. [The two network codes of 7 October, read against the game](#10-the-two-network-codes-of-7-october-read-against-the-game)
  11. [Claims corrected along the way](#11-claims-corrected-along-the-way)
- [What is still open, with a first experiment for each](#what-is-still-open-with-a-first-experiment-for-each)
- [Programs, data and proofs](#programs-data-and-proofs)
- [Credits](#credits)
- [Sources](#sources)

## The question

A program has *saving* s if it computes the transform on N points in about N·(log N)^(1−s) operations. The plain fast
transform has s = 0. In October 2026 the best savings anyone had proved were around 1.6e-3 to 1.8e-3 (the record
figures of this revision are in the [top-level README](../../README.md)). I build savings from *units*: a fixed
program that transforms W arrays at block size m with M paid moves, fewer than the plain W·m, gives the saving
ln(W·m/M)/ln m.

Four questions:

1. **Famous problems (14 October).** Is the whole problem equivalent, or isomorphic, to another really hard problem
   that mathematicians have been trying to solve? If it is, which half is hard, and is any moonshot target blocked by
   a famous barrier?
2. **Quantum computing (13 October).** I had noticed that the coordinate frames of the arrays look like stabiliser
   states. Is that exact? Is it new? Does it give anything, in particular for one route I had kept open: labels on
   frames other than the subspace frames that all certificates use?
3. **Entropy bounds (15 October).** Common-frame units are network codes. Can the entropy inequalities of network
   coding (Shannon, Ingleton, Dougherty-Freiling-Zeger) prove lower bounds stronger than my certificates?
4. **The two preprints of 7 October (10 October).** Two network-coding preprints posted that day give codes that beat
   routing on undirected networks. Can either code be read as a saving unit?

## Words used here

- **Array, paid move, free step.** An array is one vector being transformed. A paid move applies one direction's
  2-point transform to one array; it is what is counted. Free steps cost nothing in the count: relabelling addresses,
  phases, bit flips and sign flips, and mixing arrays with rational coefficients.
- **Unit "W in M" at block m.** A fixed program that applies the m-direction transform to W arrays with M paid moves.
- **Games.** The rules of what is free. The *engine game* has the letters of my Lean engine. The *equal-frame game*
  is the narrower game in which additions only join arrays at the same frame. Wider games add more free steps.
- **Digit columns f.** A unit must be correct for every number f of digit columns. This makes an array's content a
  formal sum of terms, each one input passed through one fixed product of steps.
- **Frame, frame graph Γ_m.** The coordinate system an array is written in, up to free steps. Each frame is labelled
  by a Lagrangian subspace of F_2^(2m) (a subspace of dimension m on which the symplectic form vanishes). There are 3,
  15, 135 and 2,295 frames at blocks 1 to 4. Two frames are adjacent when one paid move takes one to the other; the
  distance between frames A and B is m − dim(A ∩ B). Every input starts at the start frame, and its transform sits at
  the frame opposite it, at distance m. The frame graph is the *dual polar graph* of Sp(2m, 2).
- **Subspace frames.** The frames built from a subspace U of the address bits: 2, 5, 16 and 67 at blocks 1 to 4. All
  certificates so far, mine and the community's, use only these.
- **Common-frame (coherent) unit.** A unit with a choice of coordinates in which every paid move acts on an array
  whose terms all sit at one frame. OpenAI's unit and the community's units are of this kind as designed.
- **Session, routing, coding gain.** In network coding a session is a message that must go from a source node to a
  sink node. Routing sends each message along a path, without ever adding two messages. The coding gain is the cost
  of the best routing divided by the cost of the best code.
- **Cut bound.** The lower bound on a code's cost that says: whatever enters a set of nodes must carry the messages
  whose sinks are inside and whose sources are outside. Entropy inequalities can in principle prove more than cuts.
- **Certificates (Theorem R).** My lower bounds on moves per array: 2, 2.8 and 41/13 = 3.154 at blocks 2, 3 and 4.
  They hold for every unit in the engine game and cap the saving at 0, 0.0628 and 0.1714. The bounds at blocks 2 and
  3 are proved in Lean for the engine game; the bound at block 4 is checked by computer and upheld by skeptics, not
  checked in Lean. They are described in [../block-bounds-and-new-rules/](../block-bounds-and-new-rules/).
- **Copy wall.** In the community's record units, scratch arrays called *totals* collect terms and deliver them to
  their targets; the moves they spend in between are the copy cost. The community's complex side pays 64 copy moves
  per call (8 totals of 8). See [../why-the-savers-save/](../why-the-savers-save/).
- **Skeptic.** An independent AI agent that tried to refute a result with its own reasoning and its own code.

## What I tried

- **10 October, the two preprints.** One agent read the two network-coding preprints of 7 October against the rules
  of my game: it rebuilt one preprint's code from its published data, computed over two fields, and tested whether
  its graph could be laid on a frame graph. A skeptic checked the reading the same evening.
- **13 October, the quantum round.** Two tracks. One searched every frame and every symplectic twist at blocks 2 and
  3 (all 720 and all 1,451,520) for cancellations that the community's rules miss, then tried whole label families
  at E8 scale with an exact SAT encoding. The other built the dictionary to the stabiliser formalism by exact
  computation, read the quantum literature first-hand, and asked what transfers each way. Four skeptics: two on the
  mathematics and the model, one recomputing the dictionary, one on novelty.
- **14 October, the famous-problems round.** Five surveys, one per area (circuit lower bounds and rigidity; matrix
  multiplication and Strassen's asymptotic spectrum; network coding, integer multiplication and convolution; quantum
  gate counts and simulation; combinatorial objects such as Hadamard matrices, Reed-Muller codes and E8), each with a
  skeptic who re-derived every relation it claimed and corrected its type; then a completeness critic that looked for
  famous problems all five had missed. Each relation got a type: *equivalent* (maps both ways that keep the
  exponent), *one-way*, *same shape* (the same structure, no reduction known), *shared tools*, or *none*.
- **15 October, the entropy round.** One agent wrote the entropy linear programs for codes on Γ_2, Γ_3 and Γ_4,
  calibrated them on networks with known answers, proved why no solvable LP can beat the cut bound there, and then
  asked when cuts alone settle a layout with few arrays. One skeptic rebuilt the graphs and LPs with its own code and
  solved the hard cases the author had left.

The computations are exact where they can be: integers, fractions, Gaussian rationals, finite fields, exact group
elements, exact SAT. The LPs are floating point with rigorous dual bounds or exact rational certificates. Most ran on
a Mac in seconds; the large LPs ran on Linux machines of Google's Colab service.

## What I found

Each finding carries its status in words. "Argued by hand, not checked" means a written argument that no person and
no proof checker has checked; I say whether a skeptic re-derived it. "Checked by computer" means an exact
computation; I say whether a skeptic repeated it with its own code. A relation to a published theorem is only as good
as the reading of that theorem, which no person has checked here either.

### 1. It is a famous problem, not another one in disguise

- **The question is the classic one.** Whether N log N is optimal for the Fourier transform is an old open question.
  Morgenstern proved in 1973 that it is optimal when every constant has modulus at most 1; without that restriction it
  stayed open. Ailon calls a superlinear lower bound for the Fourier transform in linear circuits "a long standing
  open problem", and Alman and Rao call N log N "popularly conjectured" to be optimal. OpenAI's finite-tensor preprint of
  September 2026 says its result "refutes the Omega(n log n) lower bound in the unrestricted complex linear-circuit
  model"; this repository's Lean proofs give explicit positive savings for both the Walsh-Hadamard transform and the DFT.
- **The link runs one way.** A positive saving refutes the N log N conjecture, but a best saving of 0 would not have
  proved it: a cost of N log N / log log N also has saving 0.
- **What is left has no famous name.** The open question now is the exact value of the best saving s*. As far as the
  five surveys and the critic found, nobody has named or studied that refined question as a problem of its own.

*Status: argued by hand, not checked. A literature reading by five survey agents, each read by one skeptic, plus a
completeness critic; no person has checked it. The skeptic of the combinatorics survey added the main point (that
the problem is itself the famous one) when that survey had missed it.*

### 2. Exact reshapings, all of them known

These are my problem in other words, not other problems. Each transfers the exponent defined over all sizes; an
exponent that holds only on a sparse set of sizes does not transfer automatically.

- **The Fourier saving is the saving for multiplying polynomials (cyclic convolution).** One direction is the
  convolution theorem (three DFTs and pointwise products); the other is Bluestein's trick, which writes a DFT of
  length n as scalings plus one product with a fixed polynomial, and a circuit for that product can be linearised in
  its fixed factor at a cost of at most 3 times. *Status: argued by hand, not checked (a standard argument; OpenAI's
  DFT preprint proves its convolution corollary the same way); the identities were checked by computer for every n
  from 1 to 40 by a skeptic.* Polynomial multiplication has its own N log N question, solved only with bounded
  constants (Bürgisser and Lotz 2004); it is the same question, not a second one.
- **The Walsh-Hadamard saving is the saving for XOR convolution, and for every Kronecker-power transform.** For bits,
  (−1)^(x·y) = i^|x| · i^|y| · i^(−|x XOR y|), where |x| counts ones. So the WHT is one XOR convolution with a fixed
  vector between two diagonal scalings. Every invertible 2 × 2 kernel that is not a permutation times a diagonal gives
  the same exponent for its Kronecker powers: H = U·diag(1, −2)·P·U·P, and three calls to H give U up to diagonals,
  where U is the zeta kernel and P a swap. So Yates's zeta transform, Alman and Li's disjointness matrix and OpenAI's
  kernel C share the Walsh-Hadamard exponent. *Status: argued by hand, not checked; the identity was checked by
  computer exactly up to N = 128 and over finite fields, by two agents independently.* A side consequence, subset
  convolution over C in n^(2−s*+o(1))·2^n operations, has not been compared with the subset-convolution literature.
- **The saving for H is the saving for one quadratic polynomial, within a factor 15.** A circuit for H gives one for
  q(x, y) = y^T·H·x with 2N more operations. Conversely, reverse-mode differentiation of any circuit for q gives all of
  H·x at 5 times its size (Baur and Strassen); setting y = 0 and linearising costs another factor 3. So the lower-bound
  half of the problem sits at degree 2, where algebraic complexity has no method. *Status: argued by hand, not
  checked; no skeptic has read it; one exact test (H on up to 32 points, and random matrices) found the identity
  with constant 12 on those instances.*
- **Applying a layer of Hadamard gates to a state vector on a classical computer is the WHT**, and the quantum
  Fourier transform circuit, with its controlled phases merged into twiddle factors, is the radix-2 FFT. So the
  problem is exactly the cost of simulating a Hadamard layer or the QFT on a classical computer. The units are not the
  first simulations below n·2^n operations: Alman and Rao's 23/24 and my 5/6 (in
  [wht-fewer-operations](https://github.com/jacobalansussman/wht-fewer-operations)) already lowered the constant.
  They are the first whose saving grows with n. For simulation the effect is a curiosity: at 40 qubits, 40^(1.63e-3)
  is about 1.006, and even a saving of 0.1 would give only about 1.45 times. *Status: argued by hand, not checked (a
  literal restatement); a skeptic upheld it and corrected the "first below n·2^n" claim.*
- **There are two numbers, not one.** The best saving for the DFT (all lengths) and the best saving for the WHT are
  not known to be equal. OpenAI's DFT preprint turns unit-type constructions for Kronecker powers into DFTs of every
  length with the same exponent (at the cost of a (log log N)^4 factor), so WHT-type units count for the DFT; no
  reduction is known for arbitrary circuits, and none from the DFT back to the WHT. On 14 October my best
  Lean-checked figures were 1.6257e-3 for the DFT and 1.3553e-3 for a WHT corollary; on each kernel where both were
  carried through, the two figures agree, so that gap reflects which kernels have been checked, not a known
  difference between the two numbers.
  A later argument (16 October) that the all-length DFT exponent equals the Kronecker exponent through Rader's trick
  over finite fields is in [../idea-rounds/](../idea-rounds/); it is argued by hand and not checked. *Status: argued
  by hand, not checked; the statement that the two numbers are not known to be equal is a skeptic's correction of an
  earlier summary.*

### 3. Which half is famous-hard

- **Proving the saving can never reach 1 is famous-hard.** If s* < 1, the cost divided by N is unbounded: a
  superlinear lower bound for the WHT or the DFT. Nobody has proved a superlinear lower bound for any explicit matrix
  in linear circuits (Valiant's problem of 1977), not even at logarithmic depth. A skeptic added a nuance: s* < 1
  gives "above c·N infinitely often", where the textbook problem asks for "above c·N for all large N"; either would be
  far beyond anything known.
- **It would also give a lower bound for matrix multiplication.** The WHT on n² points is X ↦ H·X·H^T, two products
  with a fixed matrix, so WHT(n²) ≤ 6·MM(n) (Alman 2021, Lemma 8.3, with linearisation); a lower bound N·(log N)^ε for
  the transform gives about n²·(log n)^ε for matrix multiplication, where no superlinear bound is known either (Raz
  2002 has n² log n only with bounded constants). This follows from the previous item, so the two problems are nested,
  not independent. The other direction: matrix multiplication in n²·(log n)^c operations would give a saving of at
  least 1 − c; ω = 2 alone says nothing. *Status: argued by hand, not checked; a skeptic built Strassen's algorithm as
  a circuit (n = 2 and 4), fixed one factor to H and linearised it gate by gate: exact, within 6·C(n) gates.*
- **The usual tools are blocked.** H and F are not rigid (Alman and Williams for H; Dvir and Liu for F),
  linear-size superconcentrators exist so wiring arguments give only linear bounds, and determinant and entropy arguments need bounded or well-conditioned
  constants. Cut and flow arguments, network coding included, can never prove s* < 1.
- **The half I work on is not famous-hard.** Building units and proving bounds inside restricted games matched no
  famous open problem. Each open unit question is finite: fix W, M, m, the game and a bound on scratch arrays, and
  there are finitely many program shapes, with correctness a system of polynomial equations in the mixing
  coefficients, which is decidable over C. These questions are like "what is the rank of 3 × 3 matrix
  multiplication?" (between 19 and 23): hard by size, not by depth. My certificates are restricted-model results like
  Morgenstern's and Ailon's; they would become famous-hard only if someone proved that units capture every circuit.

*Status: argued by hand, not checked; the skeptics of three surveys upheld the link to Valiant's problem and two
upheld the link to matrix multiplication, with the nuances above.*

### 4. The network coding conjecture and the saving

- **The published implication.** The undirected multiple-unicast network coding conjecture (Li and Li 2004;
  independently Harvey, Kleinberg and Lehman) says that coding never beats routing in undirected networks. Afshani,
  Freksen, Kamma and Larsen (2019) showed how a fast circuit for the transform gives a network: H·diag(H·e_s/N)·H is
  the XOR shift by s, so two copies of a circuit with a diagonal between them carry every shift; reduce the constants
  to a finite field and bring every node to degree at most 3; routing all the shifts then costs about N log N. So the
  conjecture forces circuits of size about N log N, which means a best saving of 0.
- **What I checked.** Only the implication: one agent re-derived the map step by step, and a skeptic checked its
  finite-field step by computer (the shift identity for the cyclic DFT over the field with 97 elements, and the
  WHT-to-C identity modulo 10009).
- **What I do not claim.** Braverman and He, and the public repository `chris1234555/ncc-certificate`, state what
  this implication means for the conjecture now that positive savings are proved. I make no claim about the
  conjecture myself, I have not checked their statements, and I have not checked the outside claim that OpenAI's
  September theorem already settles it.
- **The networks are astronomically large.** For a saving of 1.68e-3 to show up as an advantage of coding over
  routing in the networks the map produces, log2 N must be at least 2^595 even in the most favourable accounting, and about 2^1786 with the constant loss that is actually
  guaranteed for degree-3 networks (about 8; the skeptic's computed ratios are 11.9, 9.6, 8.7, 8.4 and 8.0 at
  log2 N = 8, 16, 32, 64 and 1024). *Status: checked by computer (arithmetic of the map's constants), one skeptic;
  the guaranteed factor of about 8 is the skeptic's correction of a first estimate of about 2.*
- **The other direction.** The map gives s* ≤ γ, where the undirected coding gap grows like (log k)^γ on degree-3
  graphs; read backwards, any theorem that the gap is at most (log k)^(1−δ) would give a superlinear lower bound.

### 5. What only looks the same

None of these is a reduction. Each either has the same mathematical shape or shares tools with my problem. The
types are the skeptics' corrected ones. *Status of every row: argued by hand, not checked; one skeptic per row unless the row says otherwise.*

| famous problem or object | relation | what it means here |
|---|---|---|
| Matrix rigidity (Valiant's route; H and F are not rigid) | shared tools for constant factors; none for the exponent | rigidity gave Alman and Rao's 23/24; it cannot decide whether units exist, and a unit must hold for every f, which makes it an identity in a semigroup ring where a low-rank-plus-sparse split of one matrix has no meaning |
| Strassen's asymptotic spectrum and ω | same shape, different order | without its order, my composition semiring is the same free object as the one generated by matrix-multiplication tensors; capping s* is like proving ω > 2, and neither theory has the monotone that would do it |
| Direct-sum problems (Strassen's additivity conjecture, refuted by Shitov) | an analogy only | one array alone needs exactly m moves, so a unit is a failure of direct-sum additivity; no shared structure was verified |
| Alman 2021 (depth d), Alman and Li 2025 (depth 2) | same shape | the closest theory and the cleanest precedent; Alman and Li's depth-2 duality theorem is the twin of the theorem that the best saving is that of the best single unit, and their rebalancing lemma always leaves about (e ln 2)·c·N log2 N, so the units are not that lemma in disguise |
| Morgenstern (bounded constants); Ailon (well-conditioned steps) | necessary conditions that run from their theorems to mine | a saving needs large or ill-conditioned constants; vacuous at every size anyone will run |
| Integer multiplication below n log n (OpenAI's Turing-machine result) | shared tools: the same unit recursion in a different model | the labels disagree: the skeptics say shared tools or same shape, the critic says one-way; I treat it as unresolved (open) |
| Index coding and minrank over different fields | shared tools (no skeptic) | where OpenAI's network and the community's label decoder come from; a place to import mechanisms |
| Quantum Hadamard count of a Clifford circuit | equivalent only for one array with no scratch and no mixing, where the answer is m | with two or more arrays, quantum copies are tensor copies, where the Hadamard count adds exactly (checked), so quantum amortisation is trivial; my copies are direct sums with free copying and mixing, where savings exist |
| Stabiliser rank, magic monotones, the Clifford hierarchy, coherence theory, ZX | tools or vocabulary | |
| E8 mod 2, Reed-Muller cosets, dual polar spaces | same shape or shared tools; finite objects where everything is classified | see the next list |
| The group FFT conjecture (Clausen, Baum, Rockmore) | same shape (no skeptic) | it asks my question for every finite group; the abelian case is beaten, the nonabelian case untouched |
| Hadamard conjecture, orthogonal designs, asymptotic rank conjecture, Cohn-Umans, VP vs VNP, P vs NP, log-rank, sunflower, Kakeya, Hilbert's 10th and 13th problems | none | they concern first-order or polynomial-scale growth, high degree, or another object |

Three facts about shared finite objects came out of the combinatorics survey:

- **The community's E8 labels and the block-3 frames are the two point orbits of one quadric.** The 120 labels and
  the 135 block-3 frames are the nonsingular and singular points of the quadric O+(8, 2), with explicit
  isomorphisms computed by a skeptic. The packing and kissing properties of E8 play no part. *Status: checked by
  computer, one skeptic computed the isomorphisms with its own code.*
- **Reed-Muller cosets.** For relabellings that form a 2-group, the "frustration" classes of two-array separations
  are exactly the nonzero cosets of RM(1, k), also with X letters allowed. An order-3 letter breaks this picture, since
  two pieces then suffice (used in [../moonshot-small-blocks/](../moonshot-small-blocks/) and
  [../moonshot-block4-block5/](../moonshot-block4-block5/)). *Status: argued by hand, not checked; checked by computer
  on small cases, one skeptic.*
- **Hirzebruch's inequality** rules out the perfect 27-line-type certificate at F_16 with multiplicity 1 (the
  arithmetic is in the data folder). That is a settled theorem doing work for me, not a link to an open problem. The
  survey's first claim, that the Sylvester-Gallai theorem itself bears on the problem, was downgraded by its skeptic to
  "none". *Status: argued by hand, not checked; one skeptic.*

What my results would be news for, according to the surveys: exact polynomial multiplication and the whole
Kronecker family (it is the same number, so any gain in s* is a gain there); two sentences in Alman's papers that
predate OpenAI's preprints (Alman 2021, section 1.2: a circuit of size o(N log N) for H_n "is currently beyond our
techniques", now answered in the general linear-circuit model, though not checked for his synchronous circuits; Alman
and Li 2025, section 1.1: "All the smallest-known (unbounded-depth) circuits for Kronecker powers can be achieved in
this way", which no longer holds); group-FFT researchers; and, only together with a better bit-side unit, integer
multiplication. The qualitative news is OpenAI's; what my work adds is the size of the saving.

### 6. The quantum dictionary

I had read the frames as stabiliser states on 13 October and asked a round to make that exact. It is exact on the
complex side. Each row below was computed exactly at small m by one agent and recomputed by a second agent with its
own code, written from the definitions of my Lean engine. *Status of the table: checked by computer at the sizes
named; one skeptic repeated it with its own code.*

| in the transform game | in quantum language | checked at |
|---|---|---|
| the m address bits of an array | m qubits (the data, not a physical state) | notation |
| shift; sign pattern (−1)^(a·x); label (x\|z) | Pauli X; Pauli Z; Pauli X^x Z^z, with the symplectic form | m ≤ 3 |
| free steps on one array | the Hadamard-free Clifford group ⟨X, CNOT, CZ, S, i⟩; on labels, the stabiliser of the start frame. Exact up to scalars: e^(iπ/4) is Clifford but not free | all elements at m = 1, 2 (32 and 3,072); label images at m ≤ 3 (2, 48, 10,752) |
| one paid direction | √X on one qubit after a CNOT circuit; one Hadamard up to free S gates and the scalar e^(iπ/4) | m ≤ 3 |
| the kernel; the WHT | (√X)^⊗m; the WHT is H^⊗m; they differ by free S layers, up to scalars | m ≤ 3 |
| a frame | a Lagrangian: the stabiliser group, signs dropped, of a stabiliser state; equivalently an orthonormal stabiliser basis | 3, 15, 135, 2,295 frames and 6, 60, 1,080, 36,720 states at m ≤ 4 |
| the subspace frames (2, 5, 16, 67) | CSS states twisted by S on every qubit; as sets, the Lagrangians spanned by Y and Z, which S^⊗m maps onto the CSS ones | m ≤ 4 |
| the other frames | the other quadratic forms in Dehaene and De Moor's description (a subspace frame is the dot product; CSS is 0) | m ≤ 4 |
| the frame graph | the dual polar graph of Sp(2m, 2): nearest-neighbour stabiliser states, counted up to sign | m ≤ 4 |
| frame distance m − dim(A ∩ B) | the least Hadamard count of the transition Clifford (the rank of its lower-left block); also the overlap \|⟨ψ_A\|ψ_B⟩\|² ∈ {0, 2^(−d)}. It is not the cost of steering one state to another: S takes \|+⟩ to \|+i⟩ for free, yet their distance is 1 | all of Sp(2m, 2) at m ≤ 3 |
| distance between subspace frames | the Koetter-Kschischang subspace distance of random network coding | m ≤ 4 |
| which ports of the community's E8 unit can serve each other (120 root lines; orthogonal = can serve) | the commutation graph of the 120 four-qubit Paulis with an odd number of Y's; E7's 63 lines are the 63 three-qubit Paulis | all pairs |
| E8 octads (2,025); E7 heptads (135) | four-qubit Lagrangians that are not totally singular; the three-qubit Lagrangians | exhaustive |
| DFT twiddle factors | non-Clifford from order 4 on (controlled-S and up) | n ≤ 4 |

**Where it stops being exact.** Scalars: move counts ignore them, exact correctness does not. Mixing across arrays:
the engine copies, adds and erases arrays with rational coefficients, and a sum of two Clifford matrices is in
general not a multiple of a unitary (I + C is a computed example), so sums of Cliffords are an analogy only. The bit
side of the community's units has labels over Q, not F_2, and no stabiliser reading. The dictionary covers the
kernels, not the DFT's outer twiddle layers. And a name clash: the community's "Clifford rank" means Hadamard count,
while "stabiliser rank" in quantum information means something else.

**It is not new.** Every single fact in the table is standard in quantum information: stabiliser counts and overlaps
(Aaronson and Gottesman), the Hadamard-free group and its normal form (Bravyi and Maslov), Hadamard-optimal normal
forms with the count given by a block rank (Maslov and Roetteler; the same rank argument in Vandaele, Martiel, Perdrix
and Vuillot), nearest-neighbour stabiliser states (Garcia, Markov and Cross), the quadratic-form description (Dehaene
and De Moor; Pllaha, Tirkkonen and Calderbank), H^⊗n|A⟩ = |A^⊥⟩ (Aaronson and Christiano), three qubits and E7
(Cerchiai and van Geemen), the 135 heptads (Lévay, Planat and Saniga), the 120 E8 rays with their 2,025 orthogonal
bases in Kochen-Specker work (Lisoněk, Raussendorf and Singh), and the subspace metric (Koetter and Kschischang). The
Clifford reading of frames was also already in the community: icekylinx's general Clifford frame lemma (written with
GPT-6 Astra) has the frame label, the transition cost as a block rank and the subspace frames with their distances;
the community's exact-Clifford scripts and its Lean file "Exact Clifford circuits on n qubits" use Pauli tableaux;
eumemic's pull request speaks of a Clifford model, stabiliser-state modules and the dual polar graph. And it was in my
own notes of 9 October (least Hadamards = rank of the lower-left block, computed on Sp(4, 2) and Sp(6, 2)). My
update of 13 October said the link "appears new as a connection"; that was wrong (*refuted by the novelty skeptic*).
OpenAI's preprints use no Clifford or stabiliser language; the idea of a frame is theirs.

What may be new is narrow and is not new mathematics: the explicit match with the quantum literature, checked by
computer; the E8-line to odd-Y-Pauli map written out (classical O+(8, 2) geometry, for which no source was found);
and the consequences for my routes, in the next section. A page written for quantum readers, with the dictionary and
four questions, is [notes/for-quantum-readers.md](notes/for-quantum-readers.md).

### 7. What the dictionary did for the problem

**The interval lemma.** Take two frames X and Y at distance d. The frames on their shortest paths are exactly those
Z with Z = (Z ∩ X) + (Z ∩ Y), and they form a copy of the lattice of subspaces of F_2^d. In quantum words: a
stabiliser group lies on a shortest path from |0…0⟩ to |+…+⟩ exactly when it is CSS. So the shortest paths from the
start to the end pass only through subspace frames, and an array that visits any other frame has paid at least one
extra move. *Status: argued by hand, not checked; two skeptics re-derived it; checked by computer on every pair at
m = 2 and 3 and on 94,095 pairs at m = 4. Its use for the copy wall rests on the community's assumption that ports
run without detours; whether a detour through another frame could pay for itself when several paths share it is
open.*

**Labels on the other frames are dead in their simple form.** This was the route I had kept open ("avenue 4"). The
start frame has 2(2^h − 1) neighbours, half of each kind: each subspace neighbour is h − 1 moves from the end, each
other neighbour h moves. So a port labelled by a non-subspace neighbour costs h + 1 instead of h, and no such frame
can be a thin label. *Status: argued by hand, not checked; checked by computer for m ≤ 4 and h ≤ 9; one skeptic. At
h = 9 the neighbours split 511 + 511; the route's earlier description, "1,022 against 256", compared all neighbours
with the odd-label subspace ones only.*

**A new kind of cancellation, in the wider game only.** Two terms riding in one array keep a fixed twist between
them. With the right twists, a target can remove a source's content using the source's own register even when their
labels are not orthogonal, which the community's rules forbid. At m = 3, labels 001 and 111 share one total at copy
cost 2. The twisted relation covers exactly twice as many ordered pairs as the orthogonal one at odd sizes (12
against 6 at m = 3; 65,280 against 32,640 at h = 9) and less than twice at even sizes (56 against 32 at m = 4). It has
nothing to do with the non-subspace frames (one route through the example uses only subspace frames), and it is legal
only in my engine's wider game, not in the community's format, which requires a common frame at each addition.
*Status: checked by computer (exhaustive over all 720 and 1,451,520 twists at m = 2, 3); both skeptics rebuilt the
example with their own code and corrected the "exactly twice" claim at even sizes.*

**But a pair saves nothing.** A total's copy cost is at least the F_2-rank of a twisted cross-Gram matrix with one
entry for each pair (source read, target served): its diagonal is all ones, and each pair cancelled directly is a
zero. So a total whose interference is all cancelled directly costs at least one move per port; for a pair that is
2, the same as two separate totals. A skeptic extended the bound to totals that walk through several frames.
*Status: argued by hand, not checked; two skeptics; checked by computer on the searched shapes.*

**Whole families: no gain at E8 scale.** On the community's E8 frame decoder with a different twist on each total, an
exact SAT encoding of a generous relaxation finds copy cost 63 impossible and reaches 64 with the community's own
labels; with one twist for the whole family, no total drops from 8 to 7; small families (h = 3 and 4 with the exact
frame condition, and every 4-port pattern at any h with payload F_3 or F_5) all cost one copy move per port. With
payload F_2 the copy cost is at least one per port for every family and every twist (a short proof). *Status: checked
by computer, not repeated by anyone; the F_2 lemma argued by hand, not checked; nobody has re-checked the SAT
encoding. In a looser version of the encoding, with one modelling condition dropped (that the pairing values inside one
twist class agree), 62 becomes satisfiable, so the result rests on that condition.*

What it would have been worth, priced with this project's copy of the community's pricing program (its control
matches the community's record of the time, 1.6644e-3, to the last digit): one copy move fewer per call would give
1.7112e-3 (+2.81 %); one whole total fewer per call 2.0385e-3 (+22.47 %, where the bit side then binds at 1.7304e-3).
*These are prices computed by a program for constructions that do not exist, not results.*

**Tools that transfer from quantum information.**
- The description of frames by a subspace U and a symmetric form β gives exact enumeration by distance and by
  *excess* (the extra moves any path through a frame pays). The excess histogram at m = 3 is {0: 16, 1: 35, 2: 56, 3:
  28}. *Status: checked by computer; two independent recounts confirm the histogram; the formula excess =
  rank(β − dot product) was checked by one program for m ≤ 4 and not repeated.*
- Bravyi and Maslov's normal form: a run of letters on one array equals free · one block of rank r · free, so a search
  needs at most one block per array between additions that touch it. *Status: checked by computer at m ≤ 2 in the
  engine's own matrices; not repeated by anyone.*
- Vandaele and coauthors' rank-potential argument is the quantum twin of my potential bounds; the method transfers,
  not the numbers.
- A modelling check: at m = 1, all 96 group elements stay linearly independent once up to 6 digit columns are tested,
  so the "every f" rule makes cancellation exact there (*checked by computer, one program*). Larger m is untested.

**What does not transfer.** Single-copy identities (superposition lemmas, Pauli expansions, linear combinations of
unitaries): the every-f rule makes my game formal, and at f = 1 the kernel would even be a free combination
(computed at m = 1 by one program). Local-complementation classes use local Hadamards, which I pay for, and leave out CNOT, which is
free for me (on the graph states of 3 qubits the two classifications differ: local complementation puts the path
and the triangle together, the distance to the end frame puts the empty graph and the path together; computed by one
program). Quantum network coding: no-cloning
is not my copy wall. Hadamard-optimal gate counts are numbers from a different cost model.

**The network-coding parallel, firmed up.** The distance between subspace frames is exactly the Koetter-Kschischang
subspace distance, and the full frame distance (the dual polar distance) extends it to all frames (checked by
computer for m ≤ 4). The copy floor is a linear index-coding (minrank) problem, whose side-information graph at E8 is
the commutation graph of the 120 odd-Y four-qubit Paulis; the largest set of them that pairwise anticommute has 8
members, which bounds that minrank below by 8 over every field (checked by computer, one program).

### 8. Common-frame units are network codes, and entropy bounds stop at the cut bound

**The reading.** Take a unit with a choice of coordinates in which no paid move acts on an array whose terms sit on
two frames. Then it is a linear network code on Γ_m: every frame is a node that keeps everything it receives and
combines it for free, every paid move sends one symbol along one edge, and input a enters at its source frame and
must be read out at the opposite frame. The symbols may be matrix rows, so it is a vector-linear code over C; routing
costs exactly m moves per input, and a saver is a code that beats routing. This reading is not mine first: a pull
request by eumemic in the community's repository states it, and the fourth revision's
[note for network-coding readers](../../notes/for-network-coding-readers.md) builds on it. *Status: argued by hand,
not checked; several AI readings in October found no gap.*

**The worst layout is the uniform one** (Lemma 1). The symmetry group Sp(2m, 2) moves any pair of opposite frames to
any other, so copying a code once for every group element gives a code for the layout with every opposite pair
equally often, at the same cost per input. So one number per block caps every common-frame unit at that block.
*Status: argued by hand, not checked; the skeptic checked the symmetry facts by computer (transitive on all 120 and
8,640 ordered opposite pairs at blocks 2 and 3, and on arcs).*

**Cuts against certificates.** On the uniform layout the cut bound is 1.5, 1.75 and 1.875 moves per array at blocks 2,
3 and 4. As ceilings on the saving of common-frame units that is 0.415, 0.491 and 0.547; my certificates allow 0,
0.0628 and 0.1714. (*Checked by computer; the slack round's skeptic recomputed all three cut constants, and the
entropy round's skeptic the block-2 value again.*)

**No entropy LP small enough to solve proves more than the cut bound on Γ_2, Γ_3, Γ_4.** The full entropy LP has one
variable per set of sessions and arcs: 2^210 sets at block 2, before symmetry. Any LP small enough to solve must lump
frames together, and lumping caps the result: route every session along a shortest path between groups, sending
nothing inside a group; that routing satisfies all the LP's rows and every inequality real linear codes satisfy,
because it is itself a linear code on the lumped network. So whatever inequalities the LP uses (Shannon, Ingleton,
Dougherty-Freiling-Zeger, or any other valid for linear codes), it can never certify more than plain routing on the
lumped network costs (Lemma 3); and lumping into four groups or fewer never beats the cut bound, because every metric
on four points is a sum of cut metrics (Lemma 4). *Status: argued by hand, not checked; the skeptic reproduced the
values but wrote no separate check of these two proofs.*

How much lumping costs, from searches (*checked by computer for the lumpings found; that no smaller lumping does
better is a search result, not a proof*):

| block | cut bound | certificate | smallest lumping that beats cuts | what reaching the certificate needs |
|---|---|---|---|---|
| 2 | 1.5 | 2 | five groups (one frame, the three lines through it, the eight frames opposite it; a K_{2,3}) reach 1.6, but their LP has 26 variables (2^26 sets), already beyond reach | all 15 frames kept apart |
| 3 | 1.75 | 2.8 | none below 14 groups found | about 124 of the 135 frames kept apart |
| 4 | 1.875 | 41/13 | lumpings into up to 300 groups stay below 1.875 | about 2,027 of the 2,295 frames kept apart (2,026 cells still fall short) |

**The method itself works where it should, and already stops short on slightly larger networks.** On the
Okamura-Seymour network (K_{2,3} with four one-way sessions) cuts allow 6, routing needs 8, and the entropy LP proves 8
(rigorous bound 7.9999986; costs are whole numbers), matching Harvey, Kleinberg and Lehman (cited from memory); over all rates of its four
session types it proves no coding gain. On the directed butterfly it finds the coding gain 8/7. But the same K_{2,3}
with sessions in both directions (the core of the five-group lumping of Γ_2) gives only 70.8 against routing 80. The
skeptic's own Shannon LP found the butterfly's gain, a second known gain (the combination network, 9 against 10) and
reproduced the three-group lumping of Γ_2 (80), and confirmed that the Ingleton and DFZ rows are encoded correctly
(the Ingleton row gives +1 on the Vámos matroid and the DFZ row +1 on the Fano matroid over GF(2), as they must; both
are never positive on random linear rank functions).
*Status: checked by computer, one skeptic repeated it with its own code.*

So network-coding entropy bounds cannot sharpen the packing ceilings of
[../frame-packing-two-arrays/](../frame-packing-two-arrays/) at blocks 2 to 4. A common-frame unit with saving
0.1 is still not excluded at block 4: it would need a coding gain of about 15 % on its layout, at least three arrays,
and a layout that is neither single-source nor single-sink. Everything in this section and the next is about
common-frame units; it reaches every unit only through Proposition 6 (no unit costs less than the cheapest
common-frame design with the same number of arrays), which is open, see
[../frame-packing-two-arrays/](../frame-packing-two-arrays/).

### 9. Few arrays: when cuts alone exclude a saver

The open targets have three to five arrays, where the uniform layout is far too pessimistic. For one fixed layout of
W opposite pairs the cut rows alone give an exact LP that a computer solves in about a minute at block 3.

- **Lemma W (two-class stretch).** Orient each pair. A set of pairs is *stretchable* if every start in it is opposite
  every end in it. If the pairs of a layout split into two stretchable sets, every common-frame unit with that layout
  needs W·m moves: no saving. The proof uses one distance function per set and the fact that the plane with the
  max-norm is the plane with the sum-norm turned by 45 degrees. *Status: argued by hand, not checked; the skeptic
  re-derived it and checked its packings exactly on every two-class stretchable layout among 300 random ones for
  each W = 3 to 6 at each block 2 to 4, with no failure.*
- **Lemma W+ (three arrays).** If a three-pair layout falls short of compatibility by e steps, every common-frame unit
  with it needs at least 3m − e moves. The proof moves one end of a pair one step at a time, using the fact that on
  every line of Γ_m the point nearest a given frame is unique. *Status: argued by hand, not checked; the skeptic
  re-derived it and checked the near-polygon fact by brute force (all 840 and 613,440 triples at blocks 2 and 3,
  381,300 at block 4).*
- **What that covers.** e = 0 for 85 %, 78 % and 74.5 % of random three-pair layouts at blocks 2, 3 and 4, and e ≤ 1
  for about 99.45 % at block 4. So a common-frame 3 in 10 at block 4 is excluded by proof on every layout with e ≤ 1.
- **An exact screen for the rest.** Map each frame to its distances from one end of each pair; cut bounds on that
  small image (at most (m + 1)^W points) are cut bounds on Γ_m, and a rational dual certificate is checked exactly.
  At block 4, of 709 random three-pair layouts, 531 are covered by the lemma and the other 178 were all certified (177
  by a small image, one by the full LP). Of 46 layouts hill-climbed toward the worst case, the images certify 22; the
  author solved 8 of the others on all of Γ_4 and the skeptic solved the remaining 16. In all, 37 three-pair layouts
  at block 4 have been solved on all of Γ_4, most of them hill-climbed toward the worst case, and every one has cut
  value exactly 12. *Status: checked by computer; the skeptic repeated the values with its own LPs and solved the
  layouts the author had left. One correction: the author's full Γ_4 solves use the directed LP, which holds for
  the sampled orientation; the skeptic's undirected solves cover every orientation.*
- **Four and five arrays: cuts are not enough.** With four arrays, cut gaps exist: at block 3 a triangle of mutually
  opposite frames plus one more pair has cut value 11 against 12 (harmless there, since the certificate already
  demands 12), and at block 4 the skeptic found two such layouts with cut values 14 and 15 against 16. So cuts cannot
  exclude a common-frame 4 in 14 or 4 in 15 at block 4; 4 in 13 is not settled either way. At block 3 with five
  arrays, the eight random layouts tried all have cut value above 14, but layouts built on three mutually opposite frames reach
  exactly 14 and one reaches 13, so cuts cannot exclude a common-frame 5 in 14, and my certificate (14 moves for five
  arrays) is the stronger bound. Calibration at block 2, where Theorem 2W proves that routing is optimal, shows the
  same: cut values fall short from five arrays on (9 against 10, 11 against 12), so beyond four arrays a low cut value
  says nothing about the truth. *Status: checked by computer; one skeptic.*

### 10. The two network codes of 7 October, read against the game

On 7 October 2026 two preprints stated codes that beat routing in undirected multiple-unicast networks: one by
Xindan Zhang, Baochun Li and Zongpeng Li (arXiv:2610.09367) and one by Mark Braverman and Zhongtian He
(arXiv:2610.10108). I asked whether either code can be read as a saving unit. **Nothing here says whether either
preprint is right, and nothing here confirms or refutes anything in network coding.** Every statement about a
preprint is what the preprint states; my computations run on my own rebuild of one preprint's published data and
are not a check of its theorem.

**How a code would become a unit.** A network code has symbols at vertices, so it can only be read into the
equal-frame part of my game: a vertex is a frame, a session is an array owed the transform (its ends m apart), one
symbol over one edge is one paid move, a vertex combining what it holds is a free addition, and zeros are scratch
arrays. The engine's additions have rational coefficients and nothing else (a first reading thought Gaussian
coefficients were available through the phase letter; the skeptic showed that letter multiplies an array by i^f, not
by i, so the lift a code needs is to the rationals).

**Zhang, Li and Li: no placement can work.** Their code has integer coefficients, works over every finite field, and
never divides. Every unit the engine accepts that saves anything must divide by 2 in a free step (statement B of the
fourth revision's note for network-coding readers). So no placement of their code is a saving unit. *Status: argued
by hand, not checked; it rests on an argument read only by AI agents; one skeptic.*

**Braverman and He: the gain comes from the field with 9 elements.** Their code lives on the 182-vertex incidence
graph of the projective plane over the field with 9 elements: 157 sessions, 735 remaining edges, every session 5
edges long,
so routing needs 785 edge uses against the code's 735. Its feature, for a reader outside the field: a symbol on every
edge, with the rule that around each vertex the ten symbols fit a pattern in which any five fix the other five;
counting says only the zero assignment should survive (910 symbols, 910 conditions), but in this graph the
conditions repeat each other so heavily that 166 symbols stay free, and a receiver fills in its missing symbols
without the message ever travelling the five steps.

- My rebuild from the preprint's data agrees with every number it states that I tested: 175 deleted edges, 157
  sessions, 5 in and 5 out at every vertex, all sessions exactly 5 apart, and 166 free symbols over the field with 9
  elements (131 without the cubing the rule applies at lines, below the 152 a gain needs). That is agreement of my
  reading with their data, nothing more.
- The same integer rule read in the Gaussian integers modulo a large prime (10007 and 1000003, two digit
  conventions) keeps only 25 free symbols. That is an upper bound for this lift over the Gaussian rationals, a more
  generous ring than the rationals the engine needs, so the obstruction only gets stronger; a gain needs at least
  152. No rational code with these sessions is stated or known. A reason to expect no good lift: the symmetries of the plane that fix a vertex move its 10 neighbours as the full group
  of 2 × 2 matrices, and the degree-4 pattern is invariant only in characteristic 3 (a recalled textbook fact, not
  looked up).
  *Status: checked by computer, one program, re-run for this revision with the same result; the skeptic did not
  repeat the rank computations.*
- **Placement at block 5 with every edge one move: impossible.** Sessions 5 apart force block size 5; a route of 5
  edges must go to a shortest path, which keeps the distance of every pair on it; then 1,524 forced cases need 8 or 9
  different first steps between vertices 3 apart, where every frame graph has exactly 7. The cause is plain: a line
  of the plane over 9 elements has 10 points, while a frame graph over 2 elements branches 3, 7, 15 ways. *Status:
  argued by hand, not checked, with a computer count: the step from routes to shortest paths is an argument, the 7
  is computed at blocks 3 and 4 and is a standard fact, recalled, at block 5, and the 1,524 cases are counted by
  computer; the skeptic reproduced the count with separate loops and found the argument sound for this case.*
- **Blocks 3 and 4.** The first reading dismissed them because 735 moves exceed 157 × 4 = 628 and 157 × 3 = 471;
  the skeptic showed that step was invalid, since two vertices placed on one frame exchange symbols for free, so the
  move count is 735 minus the collapsed edges (*the dismissal refuted*). The skeptic then closed block 4 by another
  argument: every two edges at a vertex lie on one session route, so the collapsed edges form a matching, at most 91,
  leaving at least 644 moves against 628 (*argued by hand, not checked, with a computer count; no second skeptic*).
  Block 3 is not closed: a saving there needs at least 265 collapsed edges and 0.1 needs 313, and the counts decide
  nothing (*open*; no such placement is known, and every one would need the missing rational code first).
- **Even repaired, it would not reach 0.1.** With 157 arrays at block 5, 735 moves in blocks of rank 1 give a saving of
  4.09e-2 (4.23e-2 by the whole-block bound). That is right for the placement described and below 0.1, but it is not
  a ceiling over all placements: one that stretches some edges and collapses others is not bounded by it. At block 5,
  0.1 needs 14.9 % saved (4.26 moves per array); a code of this design would need more than 19.0 % of its edges to be
  free symbols, against the preprint's 17.6 % (causal) and 18.2 % (certificate). *Prices computed by a program for a
  construction that does not exist.*

**The other direction.** In my engine's wider game an array holds a formal sum of symbols sitting at different
vertices, and one paid step moves every one of them along its own edge for the price of one. In network words,
transmissions are charged per bundle and generator, not per symbol and edge, and a vertex cut cannot even be stated,
since one array can lie on both sides. So flow and cut bounds are not arguments there. The matching question is a
direct-sum question in a group ring: does serving W messages cost W times one? *Status: argued by hand, not
checked.* The fourth revision's network-coding note listed three messages at block 3 in 8 paid moves as open; that
case is now closed: the block-3 bound of 14/5 moves per array, proved in Lean on 14 October, gives at least 8.4, so 9
moves for three arrays (see [../block-bounds-and-new-rules/](../block-bounds-and-new-rules/)).

### 11. Claims corrected along the way

Failures are as useful as successes, so here are the claims of these rounds that turned out wrong, with what broke
each. In every case the later verdict is the one used above.

| earlier claim | what is true | broken by |
|---|---|---|
| the quantum link "appears new as a connection"; "no prior work linking operation counts to Clifford theory" | every fact is textbook, and the community already had the Clifford reading | the novelty skeptic, 13 October |
| "a paid move is one Hadamard" | as a cost, yes; the gate is √X | the quantum round's dictionary check |
| the twisted relation "exactly doubles" the orthogonal one | only at odd sizes; 56 against 32 at m = 4 | both sandwich skeptics |
| avenue 4 as I first described it: "1,022 non-subspace neighbours against 256" | 511 + 511; the comparison mixed two classes | the quantum round's own count |
| "with timing respected, even the old rule cannot share one total" | true only for the searched shape; a longer walk shares one total at cost 3 | a sandwich skeptic's explicit schedule |
| the units are "the first exact Hadamard-layer simulations below n·2^n" | Alman and Rao's 23/24 and my 5/6 came first; the units are the first whose saving grows with n | the quantum-survey skeptic |
| "with one array the game is a solved quantum problem" | only with no scratch arrays and no mixing, which is not one of my games | the quantum-survey skeptic |
| the DFT and WHT savings are one number | two numbers, not known to be equal | a skeptic of the lower-bound survey |
| the shift networks are "about 2 times" below the plain FFT | that holds only on the plain FFT's own graph; the guarantee for any degree-3 network is about 8 | the network survey's skeptic |
| Gaussian coefficients are free through the phase letter | the letter multiplies by i^f; the lift must be to the rationals | the preprint-reading skeptic |
| blocks 3 and 4 dismissed for the preprint code because 735 > 628, 471 | collapsed edges are free; block 4 closed by another argument, block 3 open | the preprint-reading skeptic |
| "cut gaps at W = 4 show cuts cannot settle 4 in 13 to 4 in 15" (block 4) | the gap was only known at block 3; now 4 in 14 and 4 in 15 are shown (cut values 14, 15), 4 in 13 is not | the entropy skeptic |
| the Sylvester-Gallai theorem, design matrices and cap sets as relations | none, or tools at most | the combinatorics skeptic |
| Strassen's additivity conjecture as "same shape" | an analogy only | the tensor skeptic |

## What is still open, with a first experiment for each

1. **Three arrays at block 4: every layout.** "Every three-pair layout on Γ_4 has cut value above 11" is proved for
   layouts with e = 0 (Lemma W, about 74.5 % of random layouts); for e = 1 Lemma W+ gives only 11, which excludes
   3 in 10 but not 3 in 11. It is checked on every layout anyone reached (37 solved on all of Γ_4, all exactly 12,
   and every other layout screened was certified above 11). A proof for the rest would make "every common-frame unit with at
   most three arrays needs W·m" a theorem at every block, extending the two-array theorem to three. Solving every
   layout is out of reach: dividing the number of three-pair layouts by the order of the symmetry group Sp(8, 2)
   leaves at least about 5.7 × 10^6 classes (4.6 × 10^7 if each pair's direction counts), and still about
   1.5 × 10^6 among the layouts with e ≥ 1, at 20 to 60 minutes of LP each. *First experiment:* take the optimal cut
   packings (the LP duals) of the 37 solved layouts, see which families of cuts they use (`w3probe.py` does this at
   blocks 2 and 3), and look for one family, like the level sets of distance functions behind Lemma W+, that gives
   3m − 0.5 or more on every layout with e ≥ 1; test a candidate on the hill-climbed layouts before trying to prove
   it.
2. **Four arrays at block 4: 4 in 13.** Cuts cannot exclude 4 in 14 or 4 in 15, and 4 in 13 is unsettled. *First
   experiment:* run `layout_search.py` at block 4 with W = 4 seeded by the triangle-plus-pendant family, looking for
   a layout with cut value at most 13; if none appears, try Lemma W+ style moves on the four-pair layouts.
3. **An entropy LP that does not lump.** Keep every frame separate and restrict instead which sets of arcs and
   sessions get variables (a partial polymatroid on the full ground set). Lemma 3 does not cap it; nobody knows how
   strong such LPs are. *First experiment:* at block 2 (90 arcs, 120 sessions in the uniform layout), give variables
   only to the sets of arcs at one frame together with the sessions starting or ending there, solve, and see whether
   it passes the cut bound 1.5.
4. **The E8 twisted decoder, properly checked.** Nobody has re-checked the SAT encoding behind "63 impossible" or the
   F_2 lemma, and one case is untested: decoders other than the frame decoder with a different twist per total.
   *First experiment:* a skeptic rewrites the encoding from the model's definition and re-runs `min 7`;
   then extend the community's label-decoder SAT with twist groups. Low odds.
5. **Detours through other frames.** The interval lemma shows a lone port gains nothing from a non-subspace frame; it
   is open whether a detour shared by several ports could pay for itself. *First experiment:* the exhaustive twist
   search of `twist.py` at m = 3, extended to totals with one detour step, on families of three and four ports.
6. **Minrank of the Pauli commutation graph.** The copy floor of the community's E8 unit is a minrank problem on the
   commutation graph of the 120 odd-Y four-qubit Paulis, with lower bound 8 over every field. *First experiment:*
   compute the Lovász theta function and Haemers-type rank bounds of that 120-vertex graph (seconds), then a SAT search
   for small-rank fitting matrices over F_2 and F_3; the related E8 centre bounds are in
   [../idea-rounds/](../idea-rounds/) and [../why-the-savers-save/](../why-the-savers-save/).
7. **A rational code of Braverman-He type at block 3.** The only route to 0.1 from that preprint would need a rational
   code on its sessions and a block-3 placement with at least 313 collapsed edges; nobody has either, and my odds are
   under 1 %. The underlying idea, local rules with far more solutions than counting allows, can be tried directly:
   *first experiment:* local rules on the block-3 frame graph (945 edges, 4,320 opposite pairs in one orbit) that are
   invariant under two or three subgroups of its symmetry group, with the dimension of the solution space by a rank
   computation; stop if no rule keeps more than a third of the edges free. Odds under 2 %.
8. **The quadratic polynomial.** The factor-15 equivalence between H and y^T·H·x has had no skeptic. *First
   experiment:* a skeptic re-derives it and extends `k1_baur_strassen.py` to random straight-line circuits for q with
   cancelling junk, checking the reverse-mode construction exactly.
9. **Integer multiplication.** Whether a better complex unit says anything about OpenAI's integer-multiplication
   exponent is unresolved (the labels disagree). *First experiment:* write out which of the two procedures that cap
   that exponent (the bit-side address interchange and the complex Hadamard layer) a complex unit would replace, using
   the preprint's own parameters.
10. **Synchronous circuits.** Whether the units give o(N log N) circuits in Alman's synchronous model, where size
    counts wires, is unchecked. *First experiment:* count the wires of one level of the recursion for the smallest
    published unit, with the scratch arrays laid out as Alman's model requires.
11. **The coding-gap exponents.** Whether a saving of 1.68e-3 is above or below the coding-gap exponents that follow
    from the two preprints' constructions is open; it needs one parameter of their amplification that nobody here
    computed. *First experiment:* compute that parameter from the published amplification theorem; this is a
    comparison of exponents, not a check of either preprint.
12. **Subset convolution.** The side consequence of section 2 has not been compared with the subset-convolution
    literature. *First experiment:* a literature search for the best known operation count of subset convolution
    over C.

## Programs, data and proofs

Every program here was written by the AI agents of these rounds, or by their skeptics. Each folder has a README with
the command, the time and memory, and whether it was re-run for this revision.

| folder | what it holds | re-run for this revision |
|---|---|---|
| [programs/famous-problems/](programs/famous-problems/) | 20 small exact checks of the reshapings and reductions of sections 1 to 5 (lower bounds, matrix multiplication and kernels, the shift networks, quantum simulation, finite geometry, the quadratic polynomial) | all 20, same outputs |
| [programs/quantum-dictionary/](programs/quantum-dictionary/) | the dictionary checks and their independent re-check, the literature checks, the interval lemma and twist search with both skeptics' programs, the E8 family searches, the transfer checks; two helper modules in `lib/` | all light ones, same outputs; see the folder README for the long ones |
| [programs/entropy-bounds/](programs/entropy-bounds/) | the entropy LP, the lumping searches, the cut LP of a fixed layout, Lemma W and W+, the exact screen, the Colab drivers; the skeptic's own graphs and LPs | the light calibrations and the skeptic's checks, same outputs; the large Γ_4 LPs were not re-run |
| [programs/preprint-reading/](programs/preprint-reading/) | the rebuild of Braverman and He's code from its published lists (which you extract from your own copy of the preprint), the field computations, the placement tests, the skeptic's collapsed-edge test | all five, same outputs |
| [data/](data/) | the outputs of the original runs, in folders named like the program folders | |
| [proofs/](proofs/) | the entropy round's proof file, its report and its skeptic's record, each with a status header | |
| [notes/for-quantum-readers.md](notes/for-quantum-readers.md) | the dictionary and four questions, written for people who know the stabiliser formalism | |
| [results.json](results.json) | every finding with its status, figure and evidence file, the open questions and the programs, machine-readable | |

Requirements: Python 3.9 or later; numpy for some programs, scipy (with its HiGHS solver) for the LPs, `highspy` for
the large first-order LPs, `python-sat` for one SAT program, `pdftotext` (poppler) to make a text copy of the
preprint. The re-runs used an Apple Mac with Python 3.9.6, numpy 2.0.2 and scipy 1.13.1.

## Credits

- **OpenAI** for the preprints "Finite tensor savings and exact Fourier circuits", "An explicit power saving for the
  exact discrete Fourier transform" and "Integer multiplication below n log n" (September 2026), the frame idea, the
  kernel C, and the Lean development `openai/math` on which this repository builds.
- **The community of `CrocSwap/integer-mult-bounds`**: icekylinx (working with GPT-6 Astra) for the general Clifford
  frame lemma (pull request #130); eumemic for pull request #288, the common-frame reading as network codes, the
  Clifford model and the dual polar graph; the community's exact-Clifford scripts and Lean file; the E8 frame decoder
  and labels of the community's record units, and the records #799 (1.6644e-3) and #942 (1.6647e-3) used as controls
  for the pricing what-ifs.
- **The two preprints of 7 October**: Xindan Zhang, Baochun Li and Zongpeng Li; Mark Braverman and Zhongtian He.
- **Network coding**: Li and Li, and Harvey, Kleinberg and Lehman, for the conjecture; Afshani, Freksen, Kamma and
  Larsen for the map from circuits to networks; Okamura and Seymour, and Harvey, Kleinberg and Lehman, for the
  calibration network; Dougherty, Freiling and Zeger for the inequality used; Deza and Laurent for the fact that on
  four points every metric is a sum of cut metrics (Lemma 4); Hu, and Papernov and Lomonosov, for the two-commodity and
  small-demand cut theorems; Koetter and Kschischang for the subspace distance.
- **Circuit complexity and algebra**: Morgenstern, Valiant, Ailon, Alman, Alman and Li, Alman and Rao, Alman and
  Williams, Raz, Bürgisser and Lotz, Baur and Strassen, Strassen, Shitov, Fritz, Hirzebruch, Kühne and Yashfe.
- **Quantum information**: Aaronson and Gottesman; Bravyi and Maslov; Maslov and Roetteler; Vandaele, Martiel,
  Perdrix and Vuillot; Garcia, Markov and Cross; Dehaene and De Moor; Pllaha, Tirkkonen and Calderbank; Aaronson and
  Christiano; Cerchiai and van Geemen; Lévay, Planat and Saniga; Lisoněk, Raussendorf and Singh.
- The work was done with AI agents (Claude) that I directed, as described at the top.

## Sources

**My rounds and their working names.** The famous-problems round of 14 October (`iso`: surveys `iso-lower`,
`iso-tensor`, `iso-network`, `iso-quantum`, `iso-combi`, one skeptic each `v-iso-*`, and the critic `iso-critic`); the
Lagrangian and quantum round of 13 October (`lagrangian`: `t-sandwich`, `t-family`, `x-dict`, `x-lit`, `x-transfer`,
skeptics `v-sandwich-math`, `v-sandwich-model`, `v-dict-compute`, `v-dict-novelty`); the entropy round of 15 October
(`nc-gain`, skeptic `v-nc-gain`, part of the packing round); the reading of 10 October (`c-read`, skeptic
`c-skeptic`). Folder map: `programs/famous-problems/` ← iso checks (`skeptic_*` from the `v-iso-*` tracks);
`programs/quantum-dictionary/` ← lagrangian (`dictionary-check` = x-dict, `dictionary-skeptic` = v-dict-compute,
`literature-checks` = x-lit, `sandwich-and-twists` = t-sandwich, `sandwich-skeptic-math` = v-sandwich-math,
`sandwich-skeptic-model` = v-sandwich-model, `twisted-families-e8` = t-family, `transfer` = x-transfer, `lib/frames.py`
from the tiny-game round of 10 October, `lib/copyfloor.py` from the copy-wall round of 13 October);
`programs/entropy-bounds/` ← nc-gain (`lp-and-layouts`) and v-nc-gain (`skeptic`); `programs/preprint-reading/` ←
c-read and c-skeptic. Private catalogue entries: M2, M3, M4, M5, M6, M7, M8, M28, L2, L3. The verdicts used are the
skeptics' of 10, 13, 14 and 15 October; none was overturned later.

**Earlier public material.** The fourth revision's
[notes/for-network-coding-readers.md](../../notes/for-network-coding-readers.md) and
[notes/network-coding/](../../notes/network-coding/) (statements A, B, C and the open case of three messages at block
3); [wht-fewer-operations](https://github.com/jacobalansussman/wht-fewer-operations) for the 5/6 result.

**Community.** `CrocSwap/integer-mult-bounds`: research/native-complex-moment552/notes/general-clifford-frames.tex
(pull request #130, icekylinx); research/hyperplane-end-lean/HyperplaneEnd/Clifford.lean; docs/research/
quadratic-affine-vectors.md; pull request #288 (eumemic); records #799 and #942. `chris1234555/ncc-certificate` (its
README says it has not been peer reviewed).

**Preprints and papers.** OpenAI, "Finite tensor savings and exact Fourier circuits", "An explicit power saving for
the exact discrete Fourier transform", "Integer multiplication below n log n" (September 2026). X. Zhang, B. Li,
Z. Li, "Network Coding Can Beat Routing in Undirected Multiple-Unicast Networks", arXiv:2610.09367. M. Braverman,
Z. He, "The Multiple Unicast Conjecture is False", arXiv:2610.10108 (version 1, 7 October 2026). Afshani, Freksen,
Kamma, Larsen, arXiv:1902.10935. Braverman, Garg, Schvartzman, arXiv:1608.06545. Alman, arXiv:2102.11992. Alman and
Li, arXiv:2509.14489. Alman and Rao, arXiv:2211.06459. Ailon, arXiv:1305.4745 and arXiv:1403.1307. Alman and
Williams, arXiv:1611.05558. Dvir and Liu, arXiv:1902.07334. Raz, ECCC TR02-012. Bürgisser and Lotz, arXiv:cs/0301016. Fritz, arXiv:2112.05949. Zuiddam,
arXiv:1504.05597. Kühne and Yashfe, arXiv:1912.06123. Björklund and Kaski, arXiv:2310.11926. Huang, Newman and Szegedy,
arXiv:1804.10368. Amy, Glaudell and Ross, arXiv:1908.06076. Dougherty, Freiling and Zeger, arXiv:1311.4601. Harvey,
Kleinberg and Lehman, IEEE Trans. Inf. Theory 52 (2006). Okamura and Seymour, J. Combin. Theory B 31 (1981). Hu,
Operations Research 11 (1963). Deza and Laurent, Geometry of Cuts and Metrics (1997). Aaronson and Gottesman, quant-ph/0406196. Bravyi and Maslov,
arXiv:2003.09412. Maslov and Roetteler, arXiv:1705.09176. Vandaele, Martiel, Perdrix and Vuillot, arXiv:2302.07040.
Dehaene and De Moor, quant-ph/0304125. Garcia, Markov and Cross, arXiv:1210.6646 and arXiv:1711.07848. Pllaha,
Tirkkonen and Calderbank, arXiv:2102.12384. Aaronson and Christiano, arXiv:1203.4740. Cerchiai and van Geemen,
arXiv:1003.4255. Lévay, Planat and Saniga, arXiv:1305.5689. Lisoněk, Raussendorf and Singh, arXiv:1401.3035. Koetter and
Kschischang, arXiv:cs/0703061. Bailey and Spiga, arXiv:1711.04012. Sources read only as abstracts, as the rounds
recorded them: Alman and Williams, Dvir and Liu, Raz, Bürgisser and Lotz, Shitov 2019. Sources cited from memory:
Valiant 1977, Morgenstern 1973 (through Ailon), the Strassen and Baur-Strassen theorems, Hirzebruch's inequality, the
rank of 3 × 3 matrix multiplication; and, in the entropy round, Hu, Harvey-Kleinberg-Lehman, Okamura-Seymour,
Deza-Laurent, Papernov-Lomonosov and the DFZ inequality (whose encoding a skeptic tested on known cases).
