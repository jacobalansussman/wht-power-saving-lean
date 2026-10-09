# The scratch copies: what is known

A note for anyone who takes up the open direction of the README, section 7. Nothing in this note is a Lean
theorem. Each statement is marked: (computed) exact arithmetic by the project's calculator; (paper) a derivation
on paper by the AI agents of this project, with computer checks at small sizes; (data) read from a certificate.
Savings are in units of 10^-8, whole-block accounting unless said otherwise, truncated.

## 1. The quantity

One *invocation* is one run of the helper circuit on v pairs of data arrays, with R helper arrays. At one moment
of the run (the *scatter*) h totals are added into all the outputs. Each total is held by a scratch copy, and a
copy makes as many moves as the dimension of the frame at which its total stands. cst is the number of copy
moves of one invocation.

- first result, h = 16: 15 totals at frames of dimension 15 and the grand total at dimension 16, so
  cst = (h-1)^2 + h = 241 (data);
- second result (the community's star rule), h = 22: 22 totals at frames of dimension h - 2 = 20, which are
  degenerate subspaces, so cst = h(h-2) = 440 (data).

In the bridged word with s = 5 stages a unit has W = 4v + R arrays and m = 5h directions. The plain schedule
needs W m moves; the network needs W m - D with

    D = 4v - 5 cst.

4v is what the exchanges save. 5 cst is what the copies of the five invocations take back: 2200 of 5280 in the
second result (41.7 percent), 1205 of 2240 in the first (53.8 percent).

## 2. Exchange rates on the first network (computed)

- One unit of copy rank per invocation changes D by 5, which is 5/1035 = 0.48 percent of the saving, in
  whole-block and in per-rank accounting alike: 261 units whole-block, 152 per-rank.
- One helper array fewer is worth 6.35 units on average (4.67 to 8.83, by the length of its chain). So one unit
  of copy rank weighs as much as 41 helper arrays.

## 3. The what-if (computed; no design behind it)

The block list of the same unit with the copy blocks deleted and everything else unchanged:

| network | as it is | without copies | factor |
|---|---|---|---|
| second result (h = 22, R = 9412) | 74745 | 128147 | 1.7145 |
| first result (h = 16, R = 7122) | 53992 | 116867 | 2.1645 |

`python3 tools/whatif_copies.py <certificate.json>` recomputes a row from the `blocks` field of a certificate.

## 4. Small steps that exist

- First circuit: retaining, in place of the grand total T, the total N_q over the triples through one point q
  (T = E_q + N_q, and N_q stands at a frame of dimension 15) gives 16 copies of rank 15: cst 241 to 240, D 1035
  to 1040, with every scatter coefficient still a half-integer. A certificate with this change was built and
  accepted by the Python reference checker; the 21 extra helper arrays it needs leave a net gain of 151 units
  (54143 against 53992). Not in Lean. (computed)
- With 16 retained totals that span the same scatter space, 240 is the floor: the only values of that space
  with a frame of dimension 15 are the multiples of E_q and N_q. (paper)
- An older variant retains E_0 .. E_(h-1) and no grand total: cst = h(h-1), scatter coefficients with the
  denominator 2(h-3). Worth a factor 1.008 at h = 16 in the three-stage word. (computed; tested on arrays at
  h = 8, 10)

## 5. What is excluded (paper)

Class: the copies move along coordinate directions and are read only by the outputs.

- Every identity over triples in this class has cst >= h(h-3) per invocation. One step of the argument quotes
  an outside theorem, the minimum rank of the complement of the Kneser graph K(h-1, 2).
- The identity of the first circuit has cst >= h(h-1) in this class; it has h(h-1) + 1.
- Even the floor h(h-3), which no construction reaches, would give at most a factor 1.249 at h = 16, 1.156 at
  h = 18 and 1.108 at h = 20 in the three-stage word. So halving the loss is excluded in this class.

## 6. What has not been examined

- Copies at arbitrary subspaces. The class of section 5 assumes coordinate directions. Since this revision the
  engine allows every subspace as a frame (README, section 5), and the star rule of the second circuit already
  puts its totals at degenerate frames. Whether other frames or other totals give copies of much lower rank is
  open.
- Other routes (fewer totals, a hand-over without copies, copies shared between the stages as the helper
  arrays are) do not appear in the working notes behind this file.

A circuit in the certificate format of this repository (`tools/gx/`, `Work/GCert`) is a data file for the
checker. A change of the scatter step itself would probably also need a change of the invocation theorem in
`Work/GCert/Chain`.
