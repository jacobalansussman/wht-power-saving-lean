# Moonshot waves at blocks 2 and 3: the everything-free game and mixed units

This section covers the parts of the moonshot construction waves of 14 and 15 October 2026 that worked at block
sizes 2 and 3. I ran eight waves of moonshot rounds in those two days; the work at blocks 2 and 3 came in waves 1,
2, 3, 4 and 7, and the work at blocks 4 and 5 is in [../moonshot-block4-block5/](../moonshot-block4-block5/). I did
this work with a team of AI agents (Claude) that I directed: I set the goals, chose which ideas to pursue and decided
at each stage what to do next; the agents ran the searches, wrote the proofs and programs, and tried to refute one
another's work. No person has read any of the proofs in this section, and nothing in it is checked in Lean.

**In one paragraph.** I looked for small programs ("units") that transform two or more arrays together with fewer
paid moves than doing them one at a time, in the games where the usual lower bounds did not reach: the
*everything-free game*, where Toffoli gates, controlled-S phases and the running XOR across digit columns cost
nothing, and *mixed units* at block 3, where some arrays get the full transform and others only a two-row part of
it. No unit was found anywhere. What came out instead are bounds: block 2 can never save, in any game whose free
steps touch only a few entries per row (Theorem 2W, argued by hand and re-derived by one skeptic); every unit with
two or more arrays needs at least W + m moves; mixed units that are tight on the block-3 counting bound cannot exist;
and an approximate unit can be worth anything at block 2 only in one narrow form that is still open. Several claims
made on the way were refuted, and they are listed as plainly as the results.

## Contents

1. [The question](#1-the-question)
2. [Words used here](#2-words-used-here)
3. [What I tried](#3-what-i-tried)
4. [What came out, question by question](#4-what-came-out-question-by-question)
5. [What failed or was refuted](#5-what-failed-or-was-refuted)
6. [What is still open, with a first experiment for each](#6-what-is-still-open-with-a-first-experiment-for-each)
7. [The programs and how to run them](#7-the-programs-and-how-to-run-them)
8. [Sources](#sources)

## 1. The question

The Walsh-Hadamard transform of an array of length 2^n can be computed block by block: split the address bits
into groups of m (a *block*), and transform each group in turn. A *unit* is a short fixed program that applies one
block's transform to W arrays at once with M *paid moves*, where the plain way needs mW. If such a unit exists and
is right for every array length, using it inside the recursion saves the exponent

    saving = ln(mW / M) / ln m,

so "2 in 3" at block 2 (two arrays in three moves) would save ln(4/3)/ln 2 = 0.415. My standing bar for a new route
was a saving of 0.1; since 15 October a saving of 0.01 also counts, as a stepping stone. For comparison, the record
savings in this repository are below 2e-3, so any unit above the bar would be a very large step.

On the morning of 14 October the lower bounds already closed most small shapes in the *engine game* (the free steps
of the Lean engine in this repository): block 2 cannot save at all (proved in Lean), and block 3 needs at least 14/5
moves per array (that morning proved in Lean from one numeric fact taken as a hypothesis; later that day proved in
Lean with no hypothesis left), so block 3 can save at most 6.28e-2. Those proofs live in
[../block-bounds-and-new-rules/](../block-bounds-and-new-rules/). Two places were left where large savings were not
excluded:

- **The everything-free game.** If the running XOR across digit columns is free, the counting certificates behind
  the old bounds stop working, because that step keeps no subspace except 0 and the whole space. There, 2 in 3 at
  block 2 (0.415), and 2 in 5 (0.166), 3 in 7 (0.229) and 3 in 8 (0.107) at block 3, were open.
- **Mixed units at block 3.** The block-3 bound counts each array as needing the full three-direction transform. A
  unit that gives some arrays only the two-row part of it has more room under the same count: up to 0.1255.

So the waves asked four questions:

1. Can arrays share moves at block 2 in any game?
2. Is there a block-3 unit when everything is free?
3. Can mixed units beat the block-3 counting bound?
4. Are approximate units (families whose error only tends to zero) worth anything here?

**Short answers.**

1. **No.** At block 2, W arrays need 2W paid moves in every game whose free steps touch a bounded number of entries
   per row, the everything-free game included (Theorem 2W; argued by hand, not checked; one AI skeptic re-derived
   it and found no error).
2. **None was found.** An exact test found no way for a shared move to cancel the other array's work, within the
   menus of free words it tried, except through readouts built from unboundedly long free words. Every unit with
   W ≥ 2 arrays needs at least max(W + 3, 2W) moves, so 2 in 5 and 3 in 6 sit exactly on the bound and 3 in 7 and
   3 in 8 above it. These four are open, and any of them would need a collapse (argued, not checked).
3. **No unit, and the room is shrinking.** Mixed units tight on the block-3 hook are impossible (Theorem T). Of the
   five targets with a little room that the waves attacked, two are excluded by Theorem 2W′; 1+5 in 12 (0.102),
   widened 2+9 in 22 (0.110) and engine 10+36 in 93 (0.114) are open, pinned into narrow shapes. One conjecture (FT),
   together with versions of it for partly touched arrays and for moves that touch several arrays at once, would
   close all of them. A few more targets with room were listed but never attacked.
4. **Almost nothing.** An approximate unit can beat the two-layer bound W + m by at most one move, and only with
   exploding coefficients. At block 2 the only form that could still save is a "kept" approximate unit with between
   W + 1 and 2W − 1 moves, such as an approximate 2 in 3. That case is open; every search for it found nothing.

## 2. Words used here

- **Array, block m, digit columns, f.** An address of the array is split into f digit columns of m bits each. A
  unit must be right for every f; agreement at a few small f proves nothing.
- **Paid move, free step.** A paid move applies the two-point transform along one direction of the block (one bit
  of every digit column) and costs one recursive call. A free step costs a constant per entry and is not counted:
  relabelling addresses, multiplying by phases, adding arrays together. Which steps are free depends on the game.
- **Unit, "W in M".** One fixed program that transforms W arrays with M paid moves. The plain cost is mW.
- **Games.** In the *engine game* the free steps are those of the Lean engine of this repository (rational mixing
  of arrays, shifts, and the engine's Clifford address maps and phases). The *widened game* also makes every address
  permutation (Toffoli included) and every pattern of the phases 1, i, −1, −i free. The *everything-free game* also
  makes free the controlled-S phase, Toffoli on any column, the *chain maps* that XOR each digit column into the next,
  and their inverses, the *running XORs* (each column XORed with all the earlier ones). These steps mix digit
  columns, and as f grows they generate infinite groups. All these games have *bounded free steps*: each free step
  touches a bounded number of entries per row and column, however large f is.
- **Kernel array, two-row array, mixed unit.** At block 3 a *kernel array* gets the full transform K (3 moves alone);
  a *two-row array* gets R, the transform on two of the three bits (2 moves alone). "1+5 in 12" means one kernel
  array and five two-row arrays in 12 moves, against 13 plain. Its saving is 1 − τ, where W_K·3^τ + W_R·2^τ = M
  (W_K kernel arrays, W_R two-row arrays).
- **Certificate, cap, demand, room, the hook.** A certificate is a counting bound: one paid move can raise a
  certain dimension by at most the *cap*, and each array needs its *demand* by the end. The *room* is
  cap × M minus the total demand; a unit needs room at least 0. The *hook* is the best family of these; at block 3
  it has cap 80 per move and demand 224 per kernel array and 144 per two-row array. A unit is *tight* when its room
  is 0: then no move may waste anything.
- **Collapse.** A product of paid moves and free steps that equals something with fewer paid moves. The trivial
  one: the same direction applied twice is a free flip (D·D = X).
- **Readout.** The free combination of moves that forms an output.
- **Approximate (border) unit.** A family of programs depending on a small number ε, whose error tends to 0 with ε
  while some coefficients grow like 1/ε^h (h is the *pole order*). It is the same idea as border rank in fast
  matrix multiplication (Bini, Capovani, Romani and Lotti, 1979). To be used in the recursion it must be made exact
  (*exactified*), which usually costs moves.
- **Kept identity, kept approximate unit.** Inside the recursion each paid move is computed by a smaller call, which
  need not be exact. A *kept* identity holds whatever the smaller calls return. A *kept* approximate unit blows up
  only along kept identities; such a unit can be run at its own move count (an exactification theorem proved in a
  separate round on approximate units; see [../theory-of-theories-and-logic/](../theory-of-theories-and-logic/)).
- **Skeptic.** An independent AI agent that writes its own code and tries to refute a claim.

**Status words.** Every finding carries one: *proved in Lean*; *checked by computer* (an exact or exhaustive
computation, with whether anyone repeated it); *argued by hand, not checked* (a written proof, with whether a
skeptic re-derived it; no person and no proof checker has checked any of them); *refuted* (and by what); *open*.
Numerical searches that found nothing are evidence, not proof, and are reported as such. When two notes disagree,
the later skeptic verdict wins, and I say where that happened.

## 3. What I tried

| Wave (date) | Block-2 and block-3 tracks | What each did |
|---|---|---|
| 1 (14 Oct, morning) | everything-free block 3; mixed units; mixed units with packing | least squares over 574 pairs of program shape and direction pattern; an exact test of cross-talk between arrays that share a move; exhaustive tables of first steps into two-row arrays; boundary searches along the hook; separation costs of packed two-row arrays |
| 2 (14 Oct, midday) | the two-layer bound; symmetric mixed units; the everything-free game at block 2; the two-layer bound in the everything-free game; a theory guide and an analyst | proofs of W + m in all games with bounded free steps; symmetric programs split by symmetry type; a 2-adic lattice bound; search rules from necessary conditions |
| 3 (14 Oct, evening) | approximate units in the everything-free game | exact twins for well-behaved approximate units; a search for an approximate 2 in 3 |
| 4 (14 Oct, night) | approximate units in the everything-free game, again | which identities survive the recursion; where an approximate unit must blow up |
| 7 (15 Oct, early morning) | approximate units at block 2; mixed units at block 3 | Theorem 2W; rank budgets of approximate units; a search for a "triangular" 2 in 3; prices of the first touch of each array; f = 2 screens of 12-move programs |

Most tracks had a skeptic afterwards (the exceptions are named with their results), and the standing rules were
the same throughout: predictions
written before computing, calibration on known cases (every search had to find a planted answer first), only exact
checks count, and every numerical fit is re-tested on fresh inputs.

## 4. What came out, question by question

### 4.1 Block 2: no saving in any game with bounded free steps

**Theorem 2W: at block 2, W arrays need at least 2W paid moves, in every game with bounded free steps.**
*Argued by hand, not checked; one AI skeptic re-derived it line by line and found no error, and checked its key
numerical fact with its own code.* The plain cost is 2W, so block 2 has no saving at all, exactly, in the engine,
widened and everything-free games alike. With the two-layer bound below, every unit with W ≥ 2 arrays at block m
needs at least **max(W + m, 2W)** moves.

The idea of the proof, in words. Let f grow. A paid move spreads every entry over many positions (each entry of the
block-2 move has size 2^(−f/2)), while a free step touches only a few. So, in the limit, a paid move is orthogonal to
every free word, and the target is orthogonal to everything with fewer than two paid moves. Split each path of the
program into its part with exactly one paid move and the rest. The one-move parts form a triangular matrix with the
paid moves themselves on the diagonal, so it is invertible. The outputs' one-move part must vanish (the target has
none), while their top part must equal the target for every array. A rank count in the von Neumann algebra of the
limit (Murray and von Neumann's direct finiteness, Sylvester's rank inequality) then gives:
(rank of the readouts) + (rank of the inputs) ≤ M, and each of the two is at least W. So M ≥ 2W. Nothing beyond
bounded support is used: no finiteness of any group, no certificate, no 2-adic structure, so the infinite groups of
the everything-free game play no role. The coefficients may even be infinite combinations of free words.

- It closes, in the everything-free game: 3 in 5 (0.263), 4 in 6 (0.415) and 5 in 7 (0.515) at block 2; 4 in 7
  (0.491), 5 in 8 and 5 in 9 at block 3; 5 in 9, 6 in 10 and 6 in 11 at block 4. At larger blocks it adds to W + m
  only when W > m, so 2 in 5, 3 in 6, 3 in 7 and 3 in 8 at block 3 are untouched.
- It also covers approximate units with bounded coefficients, formally polynomial ones, and (given an earlier
  theorem judged plausible but not checked) those whose cancellations use a finite group of free steps. It does not
  cover kept approximate units that use the infinite groups of the everything-free game (section 4.4).
- The skeptic's own orthogonality check used a wider menu than the author's: chain maps and their inverses,
  Toffoli across columns, CZ, S, parity phases and non-monomial one-column steps. The overlap still falls, like
  2^(−(f−1)/2) instead of 2^(−f/2) (*checked by computer, re-run for this revision*; the program is
  [programs/block2-two-moves-per-array/check_orth.py](programs/block2-two-moves-per-array/)).

This supersedes earlier statements. On the morning of 14 October the note
[notes/open-directions-3.md](../../notes/open-directions-3.md) §8 listed the
plain chain map with its running-XOR inverse as open at block 2 (a certificate built on the 27 lines of a cubic
surface argued it, with two lemmas unproved); Theorem 2W covers that game. It also supersedes, at block 2, the
2-adic lattice bound of 5/3 moves per array below.

**Theorem 2W′: targets that each apply at least two directions need two moves each.** *Same proof, same status
(one skeptic, not refuted).* So two-direction ("rank-2") targets never save in these games, and a mixed unit needs
at least 2(W_K + W_R) moves. This excludes 7 in 13 and 9 in 17 at block 3 with a spare row, 5 in 9 at block 4 with
two-row targets (0.152), and the widened mixed units 6+35 in 80 (0.123; it would need 82) and 9+52 in 119 (0.124;
it would need 122). It leaves 1+5 in 12 and 2+9 in 22 exactly on the bound and 10+36 in 93 with one move to spare.

**The two-layer bound: every unit with W ≥ 2 arrays needs at least W + m moves.** *Argued by hand, not checked;
one skeptic upheld it and reproduced its numbers with its own code.* It was first proved for the engine and widened
games (the wave-2 track on blocks 5 and 6; that version never had a skeptic of its own), then for every game with
bounded free steps by a trace-limit argument (Theorem T′). The everything-free version excluded, before Theorem 2W,
the largest open target of that game: **2 in 3 at block 2 (0.415)**, with or without division by 2, and also 2 in
4 and 3 in 5 at block 3 and 2 in 5 and 3 in 6 at block 4. At blocks 5 and 6 it excludes 2 in 6 at block 5 (0.317),
the case that had earlier "missed exclusion by 0.008": that near miss was a weak bound, not a hidden unit.

- **The easy proof route fails, the bound holds anyway.** The engine-game proof projects identities onto "layers"
  (how many paid moves a word has). In the everything-free game identities can mix layers: conjugating the
  parity-phase identity P = (1+i)/2 + (1−i)/2·Z by D2·CS·D1 equates a combination of one-move words with a free word
  (*checked by computer*: 72 such cases with four paid letters, none with three, in the program `percol.py`). The
  trace-limit argument avoids layers altogether.
- **Most earlier searches in this game could not have held a unit.** The first half of the proof says a unit needs
  at least W *deep* moves (at the end of a chain of m moves) and W *co-deep* moves that read inputs (*argued by hand,
  part of the proof the skeptic upheld*). The shared-top families searched earlier had one deep move for two arrays,
  so the recurring signature "the shared move serves one array" was this rank deficit, not bad luck.

**The parity lattice: units that never divide by 2 need 5/3 moves per array at block 2.** *Argued by hand, not
checked; no skeptic read it; its two key numbers were computed four independent ways.* Every free step of the
everything-free game keeps the lattice of vectors whose entries all agree modulo 1 − i and sum to an even number. At
f = 2 one move pushes out 6 units of it and each array needs 10, so a 2 in 3 would have to divide by 2. Theorem 2W
supersedes this at block 2 for every unit. What remains useful is the method: lattices built from Reed-Muller codes
give, in the "chi game" at block 3 (Clifford steps plus the running XOR), 2.327 moves per array at f = 3 and values
rising toward 2.5 at larger odd f, which would exclude a 3 in 7 that never divides by 2 there (*argued; the deciding
value at f = 5 comes from a formula fitted to 161 exact values and has not been checked exactly*).

### 4.2 Block 3 with everything free

**No shared move found a non-trivial collapse.** *Checked by computer; one skeptic re-ran the tests with the terms
the author had left out, and a second skeptic confirmed every part it finished.* When two arrays share a paid move,
each output receives a copy of the other array's work, and that copy must cancel. The cheap way to cancel it would be
a collapse. About 136 exact runs tested every combination of free words around 1 to 4 paid moves (menus of 9 to 115
free words, with Toffoli, controlled-S, chain maps and parity phases, for shared last, first and middle moves and for
two shared moves in a row), at f = 3 and 4 together and re-checked at f = 5. "Exact" here means linear algebra in
floating point with a gap of 10 to 14 orders of magnitude between the singular values kept and those dropped, not
modular arithmetic. The only collapses were the trivial ones
(a direction applied twice). The skeptic's re-run against the complete lower span flagged 34 collapses at f = 3 and 4
(42 counting one configuration twice, since one of its frames commutes with the move). Each was extracted and
re-tested at f = 5 and 6 with fresh inputs: 25 were coincidences of small f (they miss by 0.54 to 0.76 at f = 5),
and 9 were the known trivial identities in disguise (a repeated direction combined with the parity-phase identity,
proved and checked at f = 3 to 7). None is genuine. The re-classification of these flags is re-run for this
revision from the saved outputs.

**What a shared move would need.** *Argued by hand, within the menus tested.* With no non-trivial collapse, a shared
move can serve two arrays only if its readout is singular and the mix that cancels the readout is built from
unboundedly long free words. Such words exist: as permutations, the two-sided chain word A·B (a chain map on every
row, then its mirror image) has order 3, 7, 7, 31, 63 at f = 2 to 6 (*checked by computer*; the skeptic's own program
extends this to f = 18). So this door is narrow but open.

**Least squares found nothing.** Over 574 pairs of program shape and direction pattern (every shape first recovered a
planted target), the best fit missed the transform by 43 % (relative residual 0.427, a 3 in 8), and the best 2 in 5 by
0.474; a unit needs 0. This is evidence, not proof; and by the deep-move count above, the shared-top shapes among them
with a single shared last move could not have held a unit at all.

**What is left.** With 2W and W + m, every unit at block 3 in this game needs at least max(W + 3, 2W) moves, that is 2
moves per array for W ≥ 3, so the largest saving still allowed there is ln(3/2)/ln 3 = 0.369. Open: 2 in 5 (0.166)
and 3 in 6 (0.369) exactly on the bound, and 3 in 7 (0.229) and 3 in 8 (0.107) above it. Any saving there must use a
collapse (*argued, not checked*). The other bounds known in this game are in
[../block-bounds-and-new-rules/](../block-bounds-and-new-rules/).

### 4.3 Mixed units at block 3

**Theorem T: no exact unit tight on the block-3 hook can give any array a two-row target.** *Argued by hand, not
checked; one skeptic re-derived it line by line. Its finite input, a table of all first steps, is checked by
computer: the skeptic rebuilt it with its own code in three arithmetics, derived the 1,680 classes by hand and
enumerated them all.* It holds in the engine and widened games. In words: tightness at one f
forces tightness at every f, because the certificate is kept
when a column is added; a program that is the same at every f must then meet the tightness condition term by term
(a Vandermonde argument); and no single first step into a two-row array meets it (0 of 14 engine classes, 0 of 1,680
widened classes; the smallest miss is √8). This excludes 1+4 in 10 (0.119), 7+38 in 88 (0.125), 2+8 in 20, 3+7 in 21,
4+21 in 49 and 10+35 in 91. Scope: exact units, moves counted by rank, free steps that act the same way on every
column. In the engine game a separate certificate from the group Sp(6,2) (Theorem M, from a separate round on
certificates; see [../theory-of-theories-and-logic/](../theory-of-theories-and-logic/)) excludes four of them by
another route.

**Units with a little room.** Theorem T says nothing about units with room left on the hook. Room 16 means two of the
cheapest single-step defects (8 each) fit. The waves attacked five such targets above the bar: 1+5 in 12 (0.102; open
in both games, the engine version the more constrained), widened 2+9 in 22 (0.110), widened 6+35 in 80 (0.123),
widened 9+52 in 119 (0.124) and engine 10+36 in 93 (0.114). In the engine game 2+9 in 22 was excluded by the symmetric
track below, and 6+35 in 80 and 9+52 in 119 fail Theorem M's count; Theorem 2W′ (section 4.1) later excluded the
widened 6+35 in 80 and 9+52 in 119. The wave-1 skeptic also found targets the first list had missed: engine 12+44 in
113 (0.115) and 11+40 in 103 with room 16, and units with room 32 or 48 such as engine 14+53 in 135 (0.114), 13+49 in
125, 12+45 in 115, 11+41 in 105 and 14+54 in 137, and widened 9+53 in 121 (0.122), 10+57 in 131, 5+32 in 72 and 2+10 in
24. Against Theorem 2W′ the engine ones each have one move to spare, widened 2+10 in 24 sits exactly on
the bound, and widened 9+53 in 121, 10+57 in 131 and 5+32 in 72 fall below it: they would need 124, 134 and 74 moves.
Two of the engine ones, 12+44 in 113 and 14+53 in 135, are also exactly tight on Theorem M's certificate (20 per move
against 49 per kernel array and 38 per two-row array: 2,260 of 2,260 and 2,700 of 2,700), so the term-by-term argument
that excluded engine 2+9 in 22 (the symmetric track below, argued by hand with no skeptic) would exclude them too.
These last steps are my arithmetic for this write-up, not checked by a skeptic.

**The separation toll, in a strict model only.** The wave-1 packing track found that separating packed two-row arrays
always cost exactly 8 per array of the hook's room (at f = 5, f = 6 and in the large-f limit), and claimed a theorem
(S8). *Refuted as written* by its skeptic: the group of free steps that pass through the two-row target (the *Levi
letters*) has 73,728 elements, not the 8,192 the proof used, so its key step (a 2 × 2 grid of failure kernels; the
full group has 9) and its "exhaustive" check were wrong, and only 64 of its 128 permutations were engine steps. The
skeptic proved a narrower statement with a short new argument: in a strict two-level model with any finite group of
Levi letters, separation costs at least 8 per two-row array in the limit, so that model excludes the room-16 units
(*argued by hand, not checked; no one has read the repair*). The toll in general (H8) is open; if it held, the room for
block-3 mixed units would fall from 0.1255 to 0.074.

**Symmetric programs split into smaller ones.** *Argued by hand, not checked; no skeptic read it.* A program with a
symmetry is a direct sum of smaller programs, one per symmetry type, each paying the full depth of the transform, so
symmetry only shrinks the search. With Theorem M's certificate this excludes engine 2+9 in 22 (0.110), 10+35 in 91,
4+18 in 44 and 6+7 in 28: they are exactly tight on that certificate, and no single first step into a kernel array
passes its tightness test (defect 6 or 7, never 0). That defect table was later rebuilt over all 10,752 relevant
engine steps by the wave-7 track and, from an independent model, by its skeptic (*checked by computer, repeated*). One
loophole survives: for an order-4 twist, eigen-pieces of two-row arrays admit steps of defect 0. The wave-2 analyst
found that the symmetric screen overcounts "absorbed" moves, so its count of surviving symmetry types for 1+5 in 12 is
an undercount (*argued, not checked*; no exclusion is lost).

**The price of the first touch.** *Argued by hand, not checked; one skeptic read every lemma line by line. The tables
they rest on are checked by computer and were rebuilt by the skeptic from an independent model.* Group the arrays by
the move that first touches them. The waste these first touches cause adds up over the moves and can never exceed the
room (Lemma Q). A first touch by one free step has a fixed price; a touch by m steps costs at least its least stable
rank divided by m (Lemmas Q1 to Q3); a basis-free version also charges the later separation of packed arrays (Q″).
With exhaustive tables (Theorem M's certificate over all 10,752 engine steps: defect 6 or 7 on a kernel array, 2, 4 or
18 on a two-row array; hook stable ranks over all 1,680 widened classes: 8 for the cheapest two-row class, at least 20
for every other) they force:

- engine 1+5 in 12: no array may be first touched by a single step; the kernel array's first touch needs at least 5
  steps, each two-row array's at least 2, and all first touches together may waste at most 1 (Theorem M's
  certificate leaves 1 unit of room: the demand is 239 of a capacity of 240);
- engine 10+36 in 93: every kernel array's first touch needs at least 3 steps;
- widened room-16 targets: the two-row arrays first touched by single steps fall into at most two packed groups.

None of this excludes a target on its own: a two-step first touch can cost 4 instead of 8, and in Theorem M's
certificate a two-step touch can cost exactly 1, the whole room of engine 1+5 in 12 (the skeptic confirmed that such
steps exist).

**Conjecture FT: touching a two-row array completely costs at least 8, however the touch is split.** *Open.* A cheap
first touch is always a zero divisor in the samples (it leaves part of the array untouched), and every sequence that
touched an array completely paid at least 8 in the end: splits 4 + 4, 4 + 2 + 2, 2 + 4 + 2, 2 + 6, 4 + 3 + 1, 1 + 7
and 8, in letter groups of order 4 and 8, abelian and not; two packed arrays always paid at least 16. The skeptic
proved one special case: one two-row array, every touch drawn from an abelian group of Levi letters (*argued by hand,
not checked*; fractional subadditivity of rank over the group's characters). With versions for partly touched arrays
and for moves that touch several arrays at once, FT would give the toll H8 outside the strict model, and any toll
above 3.2 per two-row array excludes all the room-16 targets.

**The f = 2 screen: 12-move programs stop at 0.25, and 0.25 is not a near miss.** *Checked by computer; the skeptic
reproduced every floor with its own evaluator.* At f = 2 the data splits into three pieces that every widened free
step keeps (8, 28 and 28 coordinates), and "the program reaches the targets" becomes a row-space test per piece. The
plain 13-move programs pass and their families find them (7 of 7 fits). Every 12-move family for 1+5 in 12 stopped
at 0.25, 0.5, 0.75 or more; the best (a "kernel relay", where the kernel array's moves also carry the two-row arrays)
passes all five two-row arrays and misses exactly two directions of the kernel array (the Fourier characters 001 and
101). The smallest relay, 1+1 in 4, which the hook already forbids, stops at the same 0.25. So 0.25 is the signature
of one move shared between the kernel array's middle and a two-row array's top, and a trap for the optimiser.
Widened 2+9 in 22 stopped at 0.5 where the 24-move control reached 6.5e-5.

### 4.4 Approximate units

**At most one move of room, and only with exploding coefficients (Theorems B and B′).** *Argued by hand, not checked;
one skeptic upheld them.* The top half of the two-layer argument is a rank count, and rank counts survive limits. So
an approximate unit with bounded coefficients obeys the full bound W + m, and one whose coefficients blow up can reach
W + m − 1 at best. A shape with t < W deep moves has large-f error at least W − t whatever its coefficients: for two
arrays and one deep move, a relative residual of at least √(1/2) = 0.707, measured at 0.699 and 0.711 (*checked
numerically*; the skeptic reproduced the floor out to f = 8).

**Making approximate units exact.** Two earlier theorems, FE and FE′, give every formally polynomial approximate unit
an exact twin with the same moves, also when its cancellations use free steps that commute with the paid move and
form a finite group (*argued by hand, not checked; two skeptics upheld them*). They come from a round on approximate
units that ran beside the moonshot waves on 14 October, together with Lemma W: making a unit exact copy by copy gives
an exact unit of the same game. The wave-3 track restated both: Theorem F restates FE (its skeptic found it correct), and
Theorem F′ says the same as FE′ (judged plausible; its proof rests on a theorem cited from memory).

**"The one-move window holds no useful approximate unit in the everything-free game": refuted.** The wave-3 track
claimed that approximate units in the window W + m − 1 are worthless there. Its skeptic refuted it: the relative form
of the exactification theorem needs no finite group, only that every smaller call respects the commuting steps, and
the chain map on one bit row (and so its inverse, the running XOR) commutes exactly with the paid move on the other
row (checked at f = 1 to 8; its order 1, 2, 4, 4, 8, 8, 8, 8 grows without bound). A toy family whose cancellation uses that
commutation has error exactly
0 at f = 2 to 5 under the relative recursion (*checked by computer*). So 2 in 3 at block 2 stayed open for this kind
of approximate unit.

**What the recursion keeps at block 2.** *Argued by hand and checked by computer at f = 1 to 5; one skeptic rebuilt it
with its own code.* Whatever the smaller calls return, the recursion keeps exactly two kinds of identity: D1 D2 =
D2 D1 (the two paid moves act on different bit rows), and "D1 commutes with every free step that touches only the
other row". Several identities that hold for the exact transform are not kept (D1 commuting with a flip of its own
row, D1² = X1, the braid relation). A full-value approximate 2 in 3 must blow up along a kept identity that three moves
can compute (Lemmas H and K), and such identities exist without special design (R8, in directions 1, 2, 1, and its
mirror). One claimed narrowing was *refuted*: "a unit without an exact twin must use the running XOR with lamps (single
bit flips moved along by it)". The skeptic found a second commuting step of infinite order, the mirror running XOR;
together with the running XOR it generates a large group with no lamps (all 168 invertible 3 × 3 bit matrices at f =
3, 2,520 elements at f = 4, at least 3.08 million at f = 5). It also built a family on R8 with a different halving
that escapes the "exact program in disguise" proposition (*refuted*: that proposition holds only for the halving it
used). What stands: a unit without an exact twin must use some commuting step of infinite order, or D1 D2 = D2 D1.

**At block 2, only one narrow case is left.** *Argued by hand, not checked; it rests on Lemma W, which two skeptics
upheld in the earlier round; the skeptic of this track did not list this step among those it read line by line.*
Making an approximate unit exact copy by copy gives an exact program in the same game, so at least 2W moves by Theorem 2W: the
plain cost. Making it exact with separate scales counts at least 2(W + 1). So an approximate unit can save at block 2
only if it is *kept* and has between W + 1 and 2W − 1 moves, and then only if every limit of kept programs fails to be
a kept program (a failure of "closure"). For W = 2 this is an approximate 2 in 3 (0.415); for W = 3, 3 in 4 and 3 in 5.

- **Rank budgets.** A first-order approximate unit at block 2 needs at least 1.5W moves (Lemma Λ1: the ranks b of
  its exploding readouts and g of its base inputs satisfy W ≤ b + g and b, g ≤ M − W). The skeptic checked the
  smallest case exhaustively over the two-element field: 4,032 solvable configurations, every one with b = g = 1, and
  none in the exact case (*checked by computer, re-run for this revision*). With pole order h, M ≥ W(h + 2)/(h + 1)
  (Proposition Λh: proved when the coefficients are numbers or matrices; argued in the limit algebra), so an
  approximate W in W + 1 needs pole order at least W − 1.
- **The two layer equations alone allow an approximate 2 in 3** (Proposition Λ0, by explicit example), so the rank
  argument cannot finish the job by itself.
- **The simplest first-order 2 in 3 is an exact "triangular" 2 in 3** (Proposition Λ2): an exact 3-move program with
  outputs K x + T y and K y, for some cross-term T. It is worth 0.415 exactly when "array 2's output does not depend on
  x" holds for every smaller call. Without collapses such a program is impossible (*argued, not checked*).

**The searches found only a trade-off.** *Evidence, not proof; the skeptic re-evaluated them on new inputs.* The first
triangular runs fell to about 0.04, but the cross-term T grew without bound (17 times the signal, then 37), and the
track proved that with free steps that commute with both paid moves the residual can only fall as T grows like
1/residual: read back as an approximate unit, the two errors multiply to a constant. With T penalised, all 36 runs
sat on one curve, residual × |T| ≈ 0.5 to 0.8 (0.82 to 0.89 at f = 6). The ordinary 2 in 3, which is excluded, stayed
at 0.52 to 0.75; the plain 2 in 4 and a solvable 3-move target converged. The earlier wave-4 searches with the running
XOR in the menu ended, on fresh inputs at f = 5, at errors 0.80 to 0.85 with no 1/ε part left (they had converged to
exact programs), except one run that diverged (2.04).

### 4.5 Tools and by-products

- **Search rules from necessary conditions** (*argued by hand; not checked, except that one skeptic re-derived one
  corollary of the top rule*): W arrays need W moves reaching the target by independent routes, and, by transposing
  the program, W independent routes out of the inputs; a unit with exactly W + m − 1 moves is a single chain; the
  division by 2 must come after a move and again before a move in the reversed program; and a transposed unit is a
  unit (the transposition principle of Tellegen, Bordewijk and Fiduccia).
- **Two pieces suffice to separate two arrays with an order-3 address map** (*checked by computer, upheld*): from
  x1 + x2 and x1 − σx2, with σ the cyclic shift of three address bits, the transform of x1 is recovered exactly,
  because 1 + σ is invertible while 1 + u is a zero divisor for every order-2 step u. It is used in
  [../moonshot-block4-block5/](../moonshot-block4-block5/).

## 5. What failed or was refuted

Failures are as useful as the results, so here they are in one list.

- **"The one-move window holds no useful approximate unit in the everything-free game"** (wave 3): refuted by its
  skeptic, as in section 4.4; the relative exactification keeps the running XOR's commutation.
- **"A unit without an exact twin must use the running XOR with lamps"** and **"families built on R8 are exact programs
  in disguise"** (wave 4): refuted by the skeptic's mirror running XOR and its second halving.
- **Theorem S8 as written** (wave 1): refuted (the Levi group is nine times larger than used); a narrow version
  survives.
- **"Free words of the everything-free game have exponential order"** (wave 1): refuted as stated; the orders grow
  without bound but erratically (A·B has order 15 at f = 7, 63 at f = 10, 31 at f = 15 and 87,381 at f = 18). The
  argument it supported only needs unbounded order, so it survives.
- **"In the shared-move fits one array is served and the other is not"** (wave 1): false as a general statement; the
  best fits are balanced. Only the exact cross-talk test points at cancellation as the failure.
- **A "border family" for 2 in 3 at block 2** (found on 14 October while re-checking an earlier search of the
  everything-free game): its error kept falling while its coefficients grew, but it had been fitted to one random
  input vector. On fresh inputs one array missed by 3.2 to 5.7 times the target. It was memorising its input. Since
  then every fit in this section is re-tested on fresh inputs.
- **"Descent to 0.034"** in the wave-3 approximate search: an artefact of rescaling ε; errors are now measured at
  fixed ε.
- **The false descent of the triangular search** (wave 7): the residual fell only because the cross-term exploded
  (section 4.4).
- **"Every complete packed two-array sequence cost 16"** (wave 7): false as written; one block of the output came from
  another letter group, where complete sequences cost 80 and 88. Still at least 16, which is all FT needs.
- **"FT plus less than 60 % untouched would be enough"** (wave 7): does not follow; FT as stated prices only arrays
  touched completely, and partial touches are not priced in proportion to what they remove.
- Smaller counting errors the skeptics corrected (none changes a conclusion): 136 cross-talk runs, not 142; 59 layered
  12-move fits, not "about 110"; two search patterns of wave 7 never ran; the wave-1 open list of mixed units was too
  short.

## 6. What is still open, with a first experiment for each

1. **A kept approximate 2 in 3 at block 2 in the everything-free game** (0.415; also 3 in 4 and 3 in 5). It needs
   closure to fail among kept programs, a first-order version must have the triangular shape, and that shape needs a
   collapse. *First experiment:* run the triangular search
   ([programs/approximate-units-everything-free/](programs/approximate-units-everything-free/), `tri_als.py`) with its
   32-word `wide` menu, which holds the mirror running XOR, parity phases and boundary steps (that batch was lost to a
   machine outage before it reported), then with Toffoli-type steps across columns added; penalise the cross-term
   (`MU=0.001`), and count as a hit only a run with fresh residual below 0.05 at f = 5 and a cross-term under twice the
   signal. In parallel, prove Proposition Λh in the
   limit algebra; that would close every approximate W in W + 1 of low pole order.
2. **Block 3 with everything free: 2 in 5 (0.166) and 3 in 6 (0.369) on the bound, 3 in 7 (0.229) and 3 in 8 (0.107)
   above it.** *First experiment:* the one door the collapse test left, readouts of the form 1 − g with g a free word
   of large order (start with A·B), with contents built from powers of g and orbit sums truncated at length L; watch
   whether a solution appears as L grows ([programs/block3-collapse-test/](programs/block3-collapse-test/) has the
   exact test). On the proof side, the third-layer step (W + m + 1) is set up in the two-layer proof for W = 2,
   M = m + 2 and stops at a case split; finishing it would close 2 in 5 here and 2 in 7 at block 5 (3 in 6 would need
   the same step for W = 3).
3. **Conjecture FT** (with versions for partly touched arrays and for moves that touch several arrays at once, it
   would close all three remaining mixed targets). *First experiment:* extend the skeptic's proof from abelian Levi
   groups to the non-abelian letter group of order 8 already tested by `fulltouch4.py`, and test letter groups of order 16 and more and three packed arrays with
   [programs/first-touch-prices/](programs/first-touch-prices/); then formulate the version for partly touched
   arrays (a price of at least 8·g(x) for removing a fraction x, with g(0.4) > 0.4).
4. **Engine 1+5 in 12 (0.102).** It is exactly on the 2W′ bound, its first touches must use several steps (the kernel
   array's at least 5) and together may waste at most 1 (of Theorem M's capacity of 240), and its input and readout
   directions in the one-move layer must be exact complements. *First experiment:* build candidates from multi-step
   first touches that meet these conditions and price them with the 105-dimensional certificate tools before any search.
5. **Engine 10+36 in 93 (0.114) and widened 2+9 in 22 (0.110).** 10+36 in 93 has one move to spare on 2W′ and was
   never searched; 2+9 in 22 is exactly on the bound. *First experiment:* the f = 2 screen
   ([programs/first-touch-prices/fam.py](programs/first-touch-prices/)) on 10+36 in 93 with every kernel array first
   touched by at least 3 steps, and for 2+9 in 22 a search built around the order-4 twist loophole.
6. **The 0.25 floor as a theorem.** *First experiment:* prove exactly, starting from the 4-move relay 1+1 in 4, that a
   kernel array whose middle move is shared with a two-row array's top can never supply the characters 001 and 101.
   That would exclude every program of that shape for every choice of steps.
7. **Mixed units with more room** (engine 12+44 in 113 at 0.115 and 14+53 in 135 at 0.114, among others, each with one
   move to spare on 2W′). *First experiment:* have a skeptic check the term-by-term argument on Theorem M's certificate
   (section 4.3), which on my arithmetic already excludes engine 12+44 in 113 and 14+53 in 135, both exactly tight
   there; then apply the first-touch lemmas to the rest as was done for 10+36 in 93.

## 7. The programs and how to run them

Every program is Python 3.9 or later; most need only numpy, a few need scipy. Run each from its own folder with
`python3 file.py` (each program imports its sibling files from that folder). Times are on an Apple M4 Mac mini;
"a Colab machine" means a Linux machine of Google's Colab service
(8 processors, about 50 GB). Outputs from the original runs are in [data/](data/), under the same folder names. "Re-run
for this revision" means I ran the program again on 19 October 2026 and its output agreed with the saved one, except
where the notes below the table say otherwise; the re-run outputs are in `data/<folder>/rerun-2026-10-19/`.

| Folder | What it backs | Re-run for this revision | Not re-run (heavy, or ran on Colab machines) |
|---|---|---|---|
| [block2-two-moves-per-array/](programs/block2-two-moves-per-array/) | Theorems 2W and 2W′ (the orthogonality facts), Lemma Λ1 | `check_orth.py` (blocks 2 and 3), `lemD.py`, `check_lambda1.py` | none |
| [two-layer-bound-everything-free/](programs/two-layer-bound-everything-free/) | Theorem T′ in the everything-free game, the layer-mixing relation, Lemma O, the Theorem B floor | `sanity2.py`, `percol.py`, `lemO.py` | `orth.py`, `relfind.py`, `floor.py` |
| [block3-collapse-test/](programs/block3-collapse-test/) | the collapse test at block 3, the orders of free words, the follow-up classification | `sanity3.py`, `orders.py`, `orders_mine.py`, `ctl_annih.py`, `annih4.py` and `annih6.py` (one and two configurations), `frame_check.py`, `identity_check.py`, `classify.py`, `terms.py`, `second-skeptic/vident.py` | the full campaigns (about 136 runs), `krigid.py`, `xtalk_ext.py`, `path_ext.py`, `fu.py`, `second-skeptic/vchk.py`, the least-squares campaigns (`als3.py`, `fam.py`) |
| [mixed-units-tight-hook/](programs/mixed-units-tight-hook/) | Theorem T's table of first steps | `v1_entry.py`, `v2_contphase.py`, `c10_entry_exhaustive.py` (first prime), `c24_defects.py`, `c23_open_mixed.py` | none |
| [separation-toll-strict-model/](programs/separation-toll-strict-model/) | the refutation of Theorem S8 as written | none | `v1_struct.py` (about 7 minutes on a Colab machine) |
| [first-touch-prices/](programs/first-touch-prices/) | the first-touch tables, Conjecture FT's evidence, the f = 2 floors | `srank_classes.py`, `twostep_quick.py`, `fulltouch.py` (one batch), `small11c.py`, `skeptic/v105.py`, `skeptic/vplain.py`, `skeptic/vod.py`, `skeptic/vhook.py` | `fulltouch2.py` to `fulltouch4.py`, the 86 fits of `fam.py`, `skeptic/vlift.py` |
| [approximate-units-everything-free/](programs/approximate-units-everything-free/) | the kept identities, the refutations of waves 3 and 4, the triangular search | `check_chain.py`, `check_relX.py`, `keptrel.py`, `r8check.py`, `lamplighter.py`, `braid.py`, `propS_check.py`, `dual_toy.py`, `vcheck.py`, `vcheck7.py`, `rescore.py` (the five saved solutions) | the search campaigns of `tri_als.py` |
| [parity-lattice-block2/](programs/parity-lattice-block2/) | the parity lattice | `verify_lattice.py`, `montecarlo.py` | none |

Each folder has its own README with what each file shows, the commands, the time and memory, and the expected
output. Notes on the re-runs. The block-3 part of `check_orth.py` samples free words at random: the saved output
was made with 12 samples per value of f and is reproduced exactly with them, and a run with 40 samples gives larger
maxima (0.5, 0.354, 0.177 for the rank-2 target at f = 2, 3, 4, against 0.064, 0.088, 0.044 with 12) that still stay
within the bound 2^(−f/2), so both support the theorem. The current `annih4.py` and `annih6.py` also print a count
of non-trivial collapses that the earliest runs did not report; it is 1 for the configuration 1;2;4 (where a mix
from the menu turns one move into the direction of another) and 4 when both top moves are direction 4. Both are the
repeated-direction collapses of section 4.2; the table of runs (`out_crosstab.md`) records the same 4 for a later
run of that configuration. Three re-runs differ from the saved files: `fulltouch.py` and `small11c.py` draw random
letters or trials (every complete touch still costs at least 8; every trial still stops at 0.25 or above), and
`vcheck7.py` now also prints lines its saved file lacks (the lines in both agree). Memory:
`lemO.py` needs about 2.6 GB and `second-skeptic/vident.py` about 1.7 GB; every other re-run stayed under 1.6 GB.


## Sources

The working notes behind this section are the agents' notes of the moonshot rounds of 14 and 15 October 2026, kept in
my private research archive under `checks/wht43/moonshot/`. The wave summaries are `NOTE.md` (wave 1), `NOTE-wave2.md`,
`NOTE-wave3.md`, `NOTE-wave4.md` and `NOTE-wave7.md`; the catalogue of results is the archive's `RESULTS.md`, entries
L4, L5, L6, L13, L14, L15, X7, X8, X9, X10, X11, X34, N5, N10, T15 and C7. Where a wave summary and a later skeptic
disagree, this section follows the later skeptic: for example, the wave-7 summary drops the hedge "given F′" on
Theorem 2W's scope for approximate units, which its skeptic restored; and the collapse flags are counted as 34
distinct ones, as the follow-up's table shows (42 counts one configuration twice).

Track labels (the names the agents used, which also appear in comments inside the programs), with the folder here
that holds their programs and data:

| Track (author) | Skeptic | Wave | What | Folder here |
|---|---|---|---|---|
| ms-xw3 | v-ms-xw3, its follow-up xw3-extract (`v-ms-xw3/followup/`, `FOLLOWUP.md`) and a second skeptic of the follow-up (`v-ms-xw3/v-followup/`) | 1 | everything-free game at block 3: least squares, the exact cross-talk test | `block3-collapse-test/` |
| ms-mixed | v-ms-mixed | 1 | mixed units at block 3; Theorem T | `mixed-units-tight-hook/` |
| ms-round2 | v-ms-round2 | 1 | packed two-row arrays; Theorem S8 and its refutation | `separation-toll-strict-model/` |
| mw-big56 | none of its own | 2 | the two-layer bound in the engine and widened games | (proof only) |
| mw2-xwlayers (also `xwlayers`) | v-mw2-xwlayers | 2 | the two-layer bound in every bounded-support game; Theorems B and B′ | `two-layer-bound-everything-free/` |
| mw-xw2tight | none | 2 | the parity lattice at block 2 | `parity-lattice-block2/` |
| mw-sym3mix (orchestrator notes in `sym3mix/`, `mixed/`) | none; read by the analyst `mw-analyst` | 2 | symmetric mixed units; engine 2+9 in 22 | (its tables are superseded by `first-touch-prices/skeptic/`) |
| theory-guide (`RULES.md`) | none | 2 | search rules from necessary conditions | (rules only) |
| w3-border | v-w3-border | 3 | approximate units in the everything-free game; Theorems F, F′ | `approximate-units-everything-free/` |
| xw-border | v-xw-border | 4 | the kept identities at block 2; R8 | `approximate-units-everything-free/` |
| w7-xwborder | v-w7-xwborder | 7 | Theorems 2W and 2W′; Lemma Λ1, Propositions Λ0, Λh, Λ2; the triangular search | `block2-two-moves-per-array/`, `approximate-units-everything-free/` |
| w7-mixed3 | v-w7-mixed3 | 7 | first-touch prices, Conjecture FT, the f = 2 screens | `first-touch-prices/` |
| px-border | v-px-border-logic, v-px-border-compute | (14 Oct, beside the waves) | Lemma W and Theorems FE, FE′; the fresh-input test that refuted the "border family" for 2 in 3 | (told in sections 4.4 and 5 here; no programs shipped) |
| v-ms-fano | | 1 | the order-3 pair that separates two arrays | (in [../moonshot-block4-block5/](../moonshot-block4-block5/)) |

The main written proofs, by file in the archive: Theorem 2W and 2W′, `w7-xwborder/PROOF.md` §1 (the skeptic's verdict
`v-w7-xwborder/NOTE.md`); Theorem T′ and Theorems B, B′, `mw2-xwlayers/PROOF.md` §3 and §5 (verdict recorded in the
wave-3 brief as "upheld"); the engine-and-widened two-layer bound, `mw-big56/PROOF.md`; Theorem T, `ms-mixed/STATUS.txt`
(11:53 entries; the skeptic's verdict is in the wave-1 workflow record, summarised in `NOTE.md`); Theorem S8 and its
refutation, `ms-round2/PROOF.md` and `v-ms-round2/STATUS.txt`; the symmetric splitting, `mw-sym3mix/PROOF.md`; the
first-touch lemmas and Conjecture FT, `w7-mixed3/PROOF.md` and `v-w7-mixed3/NOTE.md`; the kept identities,
`xw-border/PROOF.md` and `v-xw-border/RUNLOG.txt`; the wave-3 approximate-unit theorems and their refutation,
`w3-border/PROOF.md` and `v-w3-border/RUNLOG.txt`; the parity lattice, `mw-xw2tight/PROOF.md`; the search rules,
`theory-guide/RULES.md`; the collapse test, `ms-xw3/NOTE.md` and `v-ms-xw3/FOLLOWUP.md`; Lemma W and Theorems FE and FE′,
`checks/wht43/px/px-border/PROOF.md`, summarised in `checks/wht43/px/NOTE.md` §2 (two skeptics); the morning state of
14 October, `notes/open-directions-3.md` §8 in this repository, and the full Lean proof of 14/5 at block 3,
`checks/wht43/px/px-lean/`.

Outside work used or cited: border rank (D. Bini, M. Capovani, F. Romani and G. Lotti, 1979); direct finiteness of
finite von Neumann algebras (F. J. Murray and J. von Neumann, *On rings of operators*, 1936) and the regular ring of
a finite algebra (S. K. Berberian), cited from memory by the agents; Sylvester's rank inequality; the transposition
principle (B. D. H. Tellegen; J. L. Bordewijk, 1956; C. M. Fiduccia, 1973); the lamplighter group and Reed-Muller
codes, both standard. Nothing in this section uses a pull request of the community's hub or OpenAI's code.

The programs were copied from the track folders named above with private paths removed; a few import lines were
changed so that each program finds its sibling files in its own folder, and two read-back programs (`classify.py`,
`terms.py`) take the data folder as an argument. In the saved outputs, file paths were shortened to file names, and
one machine name inside an output file name in `out_crosstab.md` was replaced by `m2`. The re-runs were made on
19 October 2026 on an Apple M4 Mac mini.
