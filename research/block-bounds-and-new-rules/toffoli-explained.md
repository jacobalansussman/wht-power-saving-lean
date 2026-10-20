# The Toffoli idea, explained from the start

Written on 13 October 2026, late afternoon, when the idea was new; revised for this release with the figures as
they now stand and a last section on what happened next. The rest of this folder has the details and the programs.

## 1. What is being computed

The Walsh-Hadamard transform (and the Fourier transform built on it) takes a list of N numbers and produces another
list of N numbers. The textbook fast method costs about N × log2(N) steps. The game is to prove that it can be done
in fewer: N × (log N)^(1 - s), where s is the **saving**.

- The best savings with a Lean proof are tiny. On the afternoon this note was first written, mine was 1.3553e-3;
  by the next day it was 1.6257e-3. The current record of this repository is at the top of its main README.
- The community's best figures were of the same size, about 1.7e-3.
- The Toffoli idea was about a route to savings around 0.1 to 0.25, a hundred times bigger. It was a long shot. Here
  is why it was even possible.

## 2. The game everyone plays

Think of the data as several **arrays** (working copies). At the smallest scale, each array has 8 entries, addressed
by 3 bits: 000, 001, 010, ..., 111. That is "block size 3".

The goal: every array must end up holding the transform of its own data. There are two kinds of steps.

- **Paid moves.** Apply one "direction" of the transform to one whole array. A direction mixes pairs of entries
  whose addresses differ in one chosen bit (000 with 100, 001 with 101, and so on). An array with 3-bit addresses
  needs 3 directions, and each costs 1.
- **Free steps.** Things that cost only one quick pass over an array, and are ignored in the count:
  - adding arrays together (with any fractions as coefficients),
  - multiplying entries by i (the square root of -1),
  - shifting or relabelling addresses: moving entries to new positions by a fixed rule.

The plain way transforms every array on its own: 3 arrays × 3 directions = 9 paid moves.

A **unit** is a trick that serves several arrays with fewer paid moves than the plain way, by mixing the arrays with
free additions in between. The engine of this repository turns any unit into a permanent saving in the exponent.

- Today's units save about one move per array, but they need hundreds of arrays and helpers to do it. That is why
  the savings are tiny.
- A tiny unit (2 or 3 arrays at block size 3) would be a completely different regime:
  - 3 arrays in 8 moves instead of 9: saving about 0.107;
  - 3 arrays in 7 moves: saving about 0.229;
  - 2 arrays in 5 moves instead of 6: saving about 0.166.

## 3. Why everyone thought tiny units were impossible

On 10 October I had shown that 2 arrays in 5 moves is impossible at block size 3 (argued by hand and upheld by AI
skeptics, not in Lean), and other proofs showed that tiny units need very special ingredients. But every one of
those proofs quietly assumed something about the free relabellings: **the only relabellings allowed are linear
ones.**

A linear relabelling moves address bits around with XORs. Examples on 3-bit addresses (a, b, c):
- swap two bits: (a, b, c) → (b, a, c);
- XOR one bit into another: (a, b, c) → (a, a XOR b, c);
- flip a bit (that is a shift): (a, b, c) → (a, b, NOT c).

These play nicely with the Walsh-Hadamard transform: transforming and then relabelling linearly is the same as
relabelling (in a matching way) and then transforming. Because of that, a linear relabelling can never change the
"shape" of the work still left to do. The old proofs track that shape (the "frames") and show that it can never shrink
fast enough. That is where the impossibility came from.

## 4. What a Toffoli relabelling is

A Toffoli is the simplest relabelling that is not linear:

> (a, b, c) → (a, b, c XOR (a AND b))

In plain words: flip bit c, but only when a and b are both 1. On 8 addresses it does exactly one thing: it swaps entry
110 with entry 111 and leaves the other six alone.

It is still just a permutation of the entries, a reshuffle. But it is not linear (the AND makes it non-linear), and it
does not play nicely with the transform: transforming a Toffoli-shuffled array does not give a shuffled transform. So
the shape bookkeeping of the old impossibility proofs does not apply to it. In the language of quantum computing,
linear relabellings are "Clifford" operations and Toffoli is the classic non-Clifford gate; the whole frame calculus
is built on Clifford structure, which is why everyone had used only linear ones.

## 5. Is a Toffoli relabelling really free?

That was the key question, and I checked OpenAI's own rule book (their Lean file that defines how cost is counted).
The answer looks like yes, by the same standard as the linear ones.

- In OpenAI's model every primitive operation costs 1: reading an array cell, writing one, doing a little integer
  arithmetic on an index (add, multiply, divide, remainder, ...).
- The linear relabellings everyone treats as free are not literally free either. Each costs one pass over the array,
  about N operations, which is negligible next to the transform. That is why the exponent ignores them.
- A Toffoli relabelling is also one pass over the array. Working out each entry's new address takes a fixed small
  number of bit operations. Same order of cost.
- Nothing in OpenAI's rules treats linear permutations differently from non-linear ones.

So the game the whole community had been analysing was probably smaller than the game OpenAI's model actually allows.
The impossibility proofs were correct, but for the smaller game.

What it would take in this repository: the engine computes every new address bit by looking it up in an XOR table.
A Toffoli also needs an AND, and an AND table can be built the same way at the same cost. That is about 400 to 450
lines of Lean, not written, because no unit appeared. This is argued by hand from the Lean files and OpenAI's model;
no one has checked the argument.

## 6. What it did and did not mean

It did not mean that a tiny unit exists. It only meant that the old proofs no longer ruled one out.

The first computation, with every reshuffle allowed, gave a lower bound of 7/3 = 2.33 moves per array at block size 3.
The plain way is 3, and a unit is worth something below about 2.69. So a valuable unit was no longer excluded. Not
excluded is not the same as found.

The odds were low. A Toffoli shuffle may just make the remaining work messier, not smaller. But nobody had searched
this wider game, and the prize was large.

## 7. What happened next (13-14 October)

- **The same evening, block size 3 closed again.** A better certificate (the "hook") showed that every unit at block
  size 3 needs at least 14/5 = 2.8 moves per array, even with every reshuffle and every phase pattern free. That rules
  out 3 in 8, 3 in 7 and 2 in 5, and caps any block-3 saving at 6.28e-2. In the engine's own game this bound is now
  proved in Lean. With every reshuffle free it is argued by hand and was upheld by three AI skeptics; its numbers are
  checked by computer.
- **Block size 2 never saves**, with or without Toffoli (proved in Lean in the engine's game; argued by hand for
  every game with bounded free steps, one AI skeptic).
- **Block size 4 is still open.** Every unit there needs at least 41/13, about 3.154 moves per array, with or without
  Toffoli (its numbers are checked by computer). 2 arrays in 7 moves (0.096) and 3 in 10 (0.132) are allowed by
  every bound. With Toffoli free one array has 736,785 shortest routes instead of 315, but no search found a way to
  share moves between arrays.
- **End-column tricks do not rescue it.** The bounds of blocks 2, 3 and 4 still hold when Toffoli and every other
  column-uniform step are free together with any fixed operators on a few end columns (argued by hand, two AI skeptics
  for each of two proofs).
- **The same question for other free steps.** Maps that mix neighbouring digit columns are also free in OpenAI's
  model; they are told in section 4.1 of the [README](README.md).

So the Toffoli long shot did not come in at block sizes 2 and 3. One door at block size 3 is still open through a
different free step: with the running XOR (a map that mixes neighbouring digit columns) free, 3 in 8 is not excluded
(section 4.1 of the [README](README.md)). Otherwise the long shot can only come in at block size 4 or above, where
nothing found so far excludes a saving of 0.1.

## 8. If it had worked

A single tiny unit would give a saving of roughly 0.1 to 0.25 for the Walsh-Hadamard and Fourier transforms, against
about 1.7e-3 today. It would change the picture for the whole field, not just a leaderboard. And because a tiny unit
needs no huge helper network, it is also the only kind of result here that could ever matter for real FFT code.

Even the failure is useful: an argument that tiny units are impossible with all shuffles allowed closes that question,
and that is what happened at block sizes 2 and 3 (in Lean for the engine's own game; by hand, upheld by AI skeptics,
with every shuffle free).
