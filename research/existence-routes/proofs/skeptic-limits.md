<!-- header: existence-routes -->
# The skeptic's record: limits

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** This is an AI skeptic's record of its check, as it wrote it. Its verdict is the status this revision uses for the track it checked, unless a later verdict is named in the section README.

---

# v-ex2-limits: skeptic's report on ex2-limits

14 October 2026, Mac clock 20:47 to 21:03 (first instance) and 23:01 onwards (second instance, this report). Nothing
posted. Labels: [checked] = I went through every step myself; [computed] = program and output named below;
[argued]; [gap] = a step that is missing or wrong but does not change the conclusion; [error] = a false statement.

## Verdict in five lines

1. **Not refuted.** The theorems that carry the obstruction are right as stated. I checked Theorems 1, 7(a), 8, 8′ and 9,
   Propositions 3(a)–(c) and 4, Lemmas 5 and 6, and §5.1 line by line, and found no break.
2. **The numbers reproduce.** An independent code base reproduced the E1 and E1b ranks and smallest eigenvalues to every
   digit printed (Gram matrices are the matrices of overlaps between the words). One new screen of the open case H1
   (CS combined with two paid letters, f = 2..4) shows no decaying eigenvalue.
3. **The calibration holds in substance, but two of its stated reasons are wrong.** In the engine game the letter group
   contains the scalar i·I. So the limit along one ultrafilter sees only the irreducibles of one central character,
   that is, one residue of f mod 4. It does not see "every irreducible at once". Hence:
   - the px-t514 explanation is wrong (f = 5 and f = 6 lie in different residue classes);
   - "the f-limit is the formal problem" is wrong. It is the formal problem twisted by one central character, which
     is what Theorem 1(d) itself says.
   Every refusal still holds, for reasons given below.
4. The other problems are small: one misdescribed relation in E1, an "iff" in Proposition 3(d) that is too strong, a
   missing global phase in Theorem 9, and a definition in Theorem 7(b) that the calibration example does not satisfy
   literally.
5. **By-product (trace pricing).** At odd f its guess holds for every Theorem R bound, by a short argument (§4 below).
   At even f no column-uniform free letter swaps the two parity sectors, so that case is still open.

## 1. Calibration: did it really certify and refuse, by a valid argument?

**Known units.**
- The plain schedule, the E8 unit and the 8-move control word are exact at every f, so they are limit units
  trivially (ρ = 0). The m-limit statement for the plain schedule (Lemma 5, a strongly convergent infinite product) is
  right [checked].
- For U^(⊗∞) (Theorem 7(c)) the restriction to the first n rows cannot be taken letter by letter. A diagonal phase on
  a row beyond n does not map functions of the first n rows into themselves, so it is not in the algebra of Lemma 6, and
  U^(⊗∞) is not "locally finite" in the sense of 7(b). The claim survives if one restricts block by block. Each block's
  sub-program is an exact identity whose composite is K on that block, and that composite acts as the identity on the
  first n rows [gap, presentational].
- As the author says, none of this is existence power. The limits recognise programs they are given.

**Known impossibles.**
- *Engine 2 in 5 at block 3*: refused along every ultrafilter. The reason is that the 14/5 hook certificate (cap 80,
  demand 224) lives at f = 5 and, tensored with V^(⊗(f−5)), bounds programs that are exact at any single f ≥ 5
  (px-t514 Lemma A). The reason the author gives ("L1 = the formal problem, where 14W/5 is proved") is not enough by
  itself, because a limit unit along one ultrafilter only has to be exact on one residue class mod 4.
- *px-t514's 14-move program*: refused, but for a trivial reason. It is not a unit at f = 5 either: px-t514 STATUS
  (19:16) records that "input 0 never passes depth 1 and only input 1 has a K-term". What it "passes" at f = 5 is
  certificate tightness. Its slot-0 residual has coefficient −1 on the class of K at every f, so the limit residual is
  at least 1 along every ultrafilter. The stated mechanism, "the limit contains every irreducible, including the one
  that first appears at f = 6", is false [error]: f = 6 has central character μ², f = 5 has μ¹, and no single ultralimit
  contains both.
- *mw2-absorb4's near-miss limit*: its 2/3 is a residual in the untwisted group algebra. Along one ultrafilter the twisted
  quotient is what matters. absorb4's own proof (§4 there: a unit-versus-zero-divisor argument in a finite-dimensional
  algebra) carries over to the twisted quotient: frames become cosets of F·Z instead of F, step 2 gives X ∈ F·Z, and a
  scalar does not change the frame distance used in step 3. So no fixed-coefficient point of that shape is exact on
  any residue class, and the refusal holds along every ultrafilter. ex2-limits does not give this argument; it is mine
  [argued].
- *Non-exactifying border families*: refused member by member. Each member has fixed coefficients and a nonzero formal
  residual, and the twisted caveat applies to them in the same way.

So calibration_holds = true, with the corrections above. The route is still blind for existence, for the reason the
author gives.

## 2. The theorems, line by line

- **Theorem 1** (column-uniform games) [checked].
  - (a)–(e) are right.
  - One simplification: every letter of E and Wd has entries in Z[i, 1/2], and an element of Z[i, 1/2] of modulus 1 is
    ±1 or ±i (2 = −i(1+i)², so a/2^k with |a|² = 4^k forces a = unit·(1+i)^(2k)). So the scalars are always of order
    dividing 4, and Wd needs no "bounded order" assumption.
  - The twist is real in E: (S·dir·S)² = i·I exactly [computed, `local/twist_check.py`]. So the word (S.D.S)² has
    family i^f·I, and 1 − (S.D.S)² vanishes exactly at f ≡ 0 mod 4.
  - The paragraph "What the limit object is" must read: the regular trace of the twisted algebra, which contains only
    the irreducibles with central character μ ↦ μ^(f_0) [error in the paragraph and in the report, not in the theorem].
- **Corollary 2** [checked].
- **Proposition 3** (a)–(c) [checked].
  - (d) is too strong as an "iff". If some modulus-1 root β of the recurrence is not a root of unity, the peripheral
    part Σ c_k β_k^f can tend to 0 along an ultrafilter without vanishing (e.g. |1 − β^f|²). Then the limit is softer,
    with no "decaying-only" relation.
  - Theorems 8, 8′ and 9 prove the direct statement (lim_U ρ = 0 ⇒ x_f = 0 for most f along U), so nothing proved is
    affected.
  - But the open case H1 should also be screened for quasi-periodic near-relations, which a "decaying eigenvalue" test
    up to f = 8 cannot see [gap].
- **Proposition 4** [checked].
- **Theorem 8** [checked]. Fat and thin pairs are fixed along the ultrafilter. Transitivity holds at those f where every
  thin codimension exceeds the largest sum of two fat codimensions, and that set is in the ultrafilter. The cut by
  separating hyperplanes, E′ ⊆ C_T, the Reed–Muller count on E′ and the final union-of-zero-blocks step are all right.
  This is a real proof.
- **Theorem 8′** [checked].
  - Layers separate (the largest entry is 2^(−f/2)).
  - F_ij is affine in (v, t), and the partner map is affine.
  - Codimensions are symmetric between i and j (the partner map is a bijection F_ij → F_ji).
  - Transitivity holds by taking preimages.
  - The bound ρ ≥ c·μ_Ω(bad) is right; the factor 1/n is even unnecessary.
- **Theorem 9** [checked, with one gap]. A word's nonzero entries are 2^(−d/2) times a power of i times a **global
  phase** that the proof does not mention. Because the entries lie in Z[i, 1/2], that phase is an 8th root of unity, so
  it takes finitely many values and is constant along the ultrafilter; the minimum-value step survives [gap, fixable].
  Every letter of XQ and Xchi is Clifford (ut-direct §7.1 defines Q as Clifford families; chi and P = i^parity are
  Clifford).
- **Lemmas 5 and 6, Theorem 7(a)** [checked], given Theorem D′ (mw2-xwlayers, upheld per moonshot/BRIEF3).
- **§5.1** (Noetherian compactness) [checked].
- **§5.2** (Lefschetz) needs "right at every f" to be one first-order sentence. That holds because a weighted
  automaton with D states is zero iff it is zero on all words shorter than D, over any field. With that line added it is
  right [gap].

## 3. Numbers

The first instance wrote and launched the rebuild code at 20:58 to 21:02. I collected its outputs at 23:16.

- **Independent rebuild** [computed]: my own block-2 letters (`vex2lim_rb/myxw.py`, written from the xw2lib docstring),
  my own exact-integer Gram code and my own exact paid-letter columns.
  - Control: 0 mismatches against xw2lib on 54 letters and words at f = 2..5.
  - My merge of words that agree up to a global phase gives the same 514, 1,178 and 722 families as E1.
  - Main menu: rank 152, 391, 400, 402, 402, 402 at f = 2..7; on the complement of the 111 fixed relations the smallest
    eigenvalues are 0.0286, 0.0345, 0.0372 (f = 5..7), the same as E1.
  - chain3: rank 890, 975, 994, 996, 994 at f = 3..7, 182 fixed relations, complement 0.0115, 0.0125, 0.0130.
  - One paid letter: rank 256, 711, 712, 712 at f = 2..5, 10 relations, complement 0.0115, 0.0135.
  - The relation the author attributes to X1.Pr0 is misdescribed. My exact check: X1.Pr0 and Pr0.X1 are equal at even f,
    and at odd f they differ by ±i **depending on v**, so they are not proportional (not "a sign (−1)^f"). The count of
    parity-only relations still matches, but the odd-f one is some other identity [error, harmless].
- **H1 screen** (CS combined with two paid letters, 594 families, `vex2lim_rb2`): rank 584 at f = 3 and 4; on the
  complement of the 10 fixed relations the smallest eigenvalues are 0.0146 (f = 3) and 0.0384 (f = 4). They rise, so
  there is no decaying eigenvalue (weak: f ≤ 4).
- **H1 to f = 5** (`vex2lim_rb3`, menu h1small: 1.D_z.M.D_y.1 and CS.D_z.M.D_y.1 with M among CS, A.CS, CS.A, lCS,
  CS.Pr0, Pr0.CS, CS.CS, A, Pr0, 1; 117 words, 102 families) [computed]: rank 101 at f = 2..5, 1 fixed relation; on
  its complement the smallest eigenvalues are 0.0067, 0.0630, 0.0949, 0.1063 at f = 2..5. They rise, with no decaying
  eigenvalue up to f = 5 (a small menu, and it cannot see quasi-periodic near-relations).
- **Theorem 9 screen** (Clifford letters only, two paid letters; menu cliff2: L.D_z.M.D_y.R, L in {1, A, Pr0, S1},
  12 Clifford middles, R in {1, Ai, Pr1}; 1,296 words, 1,141 families) [computed]: rank 256, 1,098, 1,114 at f = 2, 3, 4,
  still growing at f = 4 (coincidences at small f). 23 relations are exact at f = 3 and 4. f = 4 is too small to read a
  trend, so a reduced menu (cliff2s, 216 words) runs to f = 5 (§5).

## 4. The trace-pricing by-product (H4): partial check [argued, mine]

Take the setting of Theorem R: a certificate Y kept by every free letter, and a move D_z whose content is p·c, with
p = (1 + Z_r^(⊗f))/2 and z_r = 0.
- p is in the free algebra, so pY ⊆ Y and Y = pY ⊕ (1 − p)Y. D_z commutes with p.
- So the potential gain of that move is at most f_p(D) := dim(pY + D pY) − dim pY, and f(D) = f_p(D) + f_(1−p)(D).
- At odd f the free letter X_r^(⊗f) commutes with D_z and swaps p and 1 − p. Hence f_p = f_(1−p) = f(D)/2: the gain is at
  most τ(p)·f(D).
- The same works for products of such sector projections. A unit is right at infinitely many odd f, and the weights do
  not depend on f.
- So **every Theorem R bound survives trace-weighted move costs** for Pauli-sector supports.
- At even f, a column-uniform Clifford maps Z_r only to ±Z_r, and (±1)^f = 1, so no free letter swaps the sectors. That
  case, and Theorems D and T′ (which are support-growth bounds), are still open.

## 5. Screens run this evening

- **cliff2s** (`vex2lim_rb4`, Colab #5) [computed]: Clifford letters only, two paid letters,
  L.D_z.M.D_y.1 with L in {1, A}, 12 Clifford middles, z, y in {1, 2, 3}.
  - 216 words make 198 families.
  - Rank 145, 196, 196, 196 at f = 2..5, with 2 relations exact at f = 4 and 5.
  - On their complement the smallest eigenvalues are 0.0541, 0.1012, 0.1284 at f = 3, 4, 5: rising, with no decaying
    eigenvalue, as Theorem 9 requires.
- **h1small** (§3): rising to f = 5.
- **cliff2** (§3): too small an f to read a trend.

Every job of mine has finished (EXIT 0 on Colabs #4, #1 and #5). Nothing is left running.

## Sources

ex2-limits/PROOF.md, NOTE.md, STATUS.txt; px/px-t514/STATUS.txt (lines at 15:20, 19:16, 21:16 to 27:16:45);
moonshot/mw2-absorb4/NOTE.md §2–§4; moonshot/mw2-xwlayers/PROOF.md §0–§2; meta/reps/TABLE.md (game letters);
wave17/ut-direct/PROOF.md §7.1 (class Q); proofs/NOTE.md (14W/5, hook cap 80, demand 224); moonshot/BRIEF3.md line 10
(mw2-xwlayers upheld). Code: vex2lim_rb/ (ctl, gram, paid, analyze), vex2lim_rb2/ (chain3, H1 screen), vex2lim_rb3/
(Theorem 9 screen, H1 to f = 5), local/twist_check.py. Outputs: */remote-out/ssh4/ and vex2lim_rb3/remote-out/ssh/.
