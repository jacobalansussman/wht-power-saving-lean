> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# The open window for a cheaper 112-port centre

ctr-window, 15 Oct 2026. **Version 4, 16:50 EDT.** Earlier versions: v1 16:24, v2 16:31, v3 16:46. Later versions refine this file in place,
and the change log is at the bottom.
- Everything about the hub was read first-hand from its GitHub: PR bodies, READMEs at the PR heads, one doc on main
  and one Discussion.
- Copies are in `centre/ctr-window/hub-read/`. They are untrusted data, and none of the hub's code was run.
- Our own computations are in `centre/ctr-window/calc/`.

Labels: **[read]** from a named hub source; **[computed]** by our code; **[argued]** reasoning only.

## In short

1. **The hub has already proved "≥ 64" for every centre that keeps today's kind of decoder** [read]. That covers:
   - any number of terms whose factors lie in the column space of a rank-8 realisation (#2338 Theorem P, "case (a)");
   - any decoder with at most 8 paid terms (#2338);
   - 9 terms with one factor outside that space ("case (b)", **#2724**, opened 18:27 UTC today).

   Every label table of rank 9 is the dense one (#2338 Theorem R), so this covers the hub's whole frontier format for
   those shapes. Every one of these PRs except #2228 is open, not merged.
2. **Four single-layer shapes are open.** Nobody has closed them or exhibited a decoder under 64 in them [read]:
   - **O1:** a rank-8 realisation (so the junk is unchanged) with **at least 10 terms** and some factors outside its
     column space.
   - **O2:** 9 terms that realise a **rank-9** matrix (case (c)). The junk changes, and it contains all the co-star
     decoders.
   - **O3:** at least 10 terms with realisation rank 9 or more.
   - **O4:** label tables of rank 10 or more, which means a wider window. The Gram case of O4 is already closed.

   **O5**, decoders with more than one layer, lies outside every floor anyone has, ours included.
3. **Accounting, now settled from the hub's own text.**
   - The hub's *shipped* format charges a centre copy **the dimension of its source frame, with the copy descending to
     frame 0** [read: #2338 §3.7, "COPY rank equals the dimension of the source frame, temp at frame 0"]. That is
     cx-medium's one-sided accounting, so **its floor of 52.36 holds there**.
     - So, in the shipped format, **the window is 53 to 63 per call** in shapes O2 and O3, and 60 to 63 in O1 (next
       point).
   - The hub states its floor theorems in the cheaper **pairing rank** (`rank_F2 B[R,A]`), so they hold in the shipped
     format a fortiori.
   - A *projected* copy, charged only its pairing rank, is SPEC-legal but forbidden by the checker convention.
     - In that accounting the clique count gives **no floor above 16.9** [computed]. One explicit term costs 1 and meets
       598 of the 1,260 frames.
     - The hub priced one projected copy, on the bit side, at −3.57 % net, because readers pay extra climbs [read].
   - **So ctr-search should refuse below 52.36 only for one-sided costs.**
   - **New: O1 has a higher floor, 60** [computed, exact integer re-check]. A real rank-8 realisation is
     `D(2G)D⁻¹` (hub #1212, Theorem 1, merged). So every 8 lines that form a basis of R^8, not only the 1,260 frames,
     need 8 terms whose gather sets meet them.
     - A weighted count over the frames and 9,300 bases, checked against every one of the 8,283,457 subspaces of
       F2^9, gives **cost ≥ 59.07**, so at least 60 per call (costs are whole numbers).
     - This holds for every single-layer decoder of a rank-8 realisation in the one-sided accounting, with any number
       of terms and any factors.
     - So in the shipped format **O1's window is 60 to 63**. O2 and O3 keep 53 to 63, because only the frames are
       forced there (section 4).
4. **For ctr-search:**
   - **O1 is the cheapest to price, but now the narrowest.** Its junk is today's ±1/2 pattern up to diagonal gauge,
     so the U8/T74 helper circuits still apply, and so does the price rule of +2.8 % per unit. Its window is 60 to 63
     in the shipped accounting: at most 4 units, about +11 %. **Refuse any O1 candidate whose one-sided cost is below
     60.**
   - **O2 is wider (53 to 63) but changes the junk.** A hit there needs a helper re-price before it can be valued.
   - At 10 terms, a total of 63 or less needs an average cost below 6.3.
   - A known O1 decoder exists: #1633's exact two-sided ten-term P112 witness, with a copy bill of 82. Use it as a
     calibration point.
   - Stay out of the fixed-reader neighbourhood of the coordinate decoder, which is closed piece by piece (C7 to C9
     below).

## 1. The model, in the hub's words [read]

Sources: #2338 §1, #2228 §1, #2724 §2, and the hub doc `docs/research/complex-decoder-scope-map.md` on main.

- **Unit.** P112 is the 120 E8 root lines minus the frame F0 = {e_{2a} ± e_{2a+1}}. These are the hub's 112 ports and
  ours.
- **Dense pairing.** `B_ts = (<r_t,r_s> + 1) mod 2`. It is 1 on the diagonal and on orthogonal pairs, and 0 on the 60°
  pairs.
  - Its F2-rank is 9.
  - It equals `Ω(u_t, u_s)` for labels `u_t ∈ F2^9`. The rank of a block is that of Ω on the two label spans:
    `rank B[R,A] = dim U_R − dim(U_R ∩ U_A^⊥)`.
- **Realisation.** A realisation is a matrix N with `N_tt ≠ 0` and `N_ts = 0` on orthogonal pairs. It is free on the
  60° pairs: those entries are the junk, which the helpers clean. On each of the 1,260 orthogonal 8-frames, N is
  diagonal and invertible.
- **Decoder.** A decoder is `N = Σ_j p_j q_jᵀ`. It is single-layer: each term is one total, gathered from
  `R_j = supp p_j` and scattered to `A_j = supp q_j`.
- **Costs.**
  - *Pairing rank* (the hub's theorems): `Σ_j rank_F2 B[R_j, A_j]`.
  - *One-sided* (the shipped checker, cx-medium): `Σ_j dim(source frame of total j)`, which is at least
    `Σ_j dim U_{R_j}`.
  - One-sided is never below pairing rank.
- **Today's decoder.** The removed-frame coordinate decoder has 8 terms, gathering `N(f)`: the 56 live lines not
  orthogonal to a removed line f.
  - Each term costs **8 in both accountings**, 64 per call [computed, `calc/a2clique.py`].
  - It is rigid up to gauge (#717 D1(d), #1212) [read].
- **Rank and junk.**
  - Rank-8 realisations are gauge-rigid: `D1·Gram·D2` [read: #1212, #1815; #1389 to #1474 for punctures]. Their junk
    is today's pattern, rescaled.
  - A nine-term decoder of a rank-9 realisation adds a rank-one Schur defect on the junk [read: #2338's frame-module
    theorem]. Its junk is new, so its helpers must be re-priced.

## 2. What the hub has closed at ≥ 64 (P112, dense B; PR status at about 20:30 UTC today) [read]

| # | decoder shape | statement | status in the hub | where |
|---|---|---|---|---|
| C1 | any label table of F2-rank ≤ 9 | it equals the dense B, so everything below holds for every such B | CHECKED (16 DRAT-verified UNSATs) | #2338 Theorem R (open); P111 row is #620 / #1097 |
| C2 | factors in the live root code, **any number of terms** | ≥ 64 | PROVEN | #2228 Thm 4 (**merged**); #1180 (merged) |
| C3 | **any rank-8 realisation, factors in col M / row M, any number of terms** (case (a)) | ≥ 64; every nonzero term costs ≥ 8 | PROVEN from a CHECKED census: P112 has only 8 deficient label classes, none containing a plane | #2338 Theorem P (open); #2277 had 56 here |
| C4 | **≤ 8 paid terms, factors anywhere, any realisation** | square, so C3 gives ≥ 64 | PROVEN | #2338 §3.3 1.3 (open); #717 D1(b) for 8 hyperplane-flat terms |
| C5 | **9 terms, rank-8 realisation, exactly one factor outside col M** (case (b)) | ≥ 64 (P104: ≥ 56) | PROVEN from exact enumeration, not DRAT; d0 = 5, 6 re-implemented independently | **#2724** (open) |
| C6 | free terms (`B[R,A] = 0`); proportional normals | removed or merged without raising the cost | PROVEN | #2338 1.1-1.2; #1224 Lemma D (merged) |
| C7 | the coordinate decoder's 8 readers kept, plus one arbitrary extra total (real) | ≥ 64 | proof and finite checks | #845 (merged) |
| C8 | fixed readers plus two extra totals, source union up to 48 | ≥ 64 | scoped | #1787 (43), #2382 (47), #2581 (48) |
| C9 | #2021's fixed two-column rank-9 completion: 9 terms, any coefficient, any number of terms in its factor spaces, its right-space extension | ≥ 64, sharp | scoped proofs | #2321 (merged), #2385, #2398, #2515 |
| C10 | k independent calls on one table, factors in col M_c | ≥ 128 (k = 2); ≥ 64k (k ≥ 3, symmetric charts) | PROVEN | #2338 Thm B / Bk; #1714, #2294, #2342, #2395 |
| C11 | label table of rank 10 or 11, **Gram realisation, factors in the root code**, any m, any fitting B | rank 10: every term ≥ 8, so ≥ 64; rank 11: every term ≥ 7, so ≥ 56, above the tie price of 55 | RECORDED (rank 10: 687 of 688 DRATs; rank 11: kissat only) | #2338 §5.1, #2639 H1/H2 (open) |
| C12 | co-star nine-term decoders with 3 to 7 distinct normals | **no exact decoder exists**, whatever its cost | exact certificates in a fork, not a PR | Discussion 1267 (philippeVV, 14-15 Oct) |
| C13 | structured classical front ends (rotations, WHT, halves, coset splits) | 72 to 192 | CHECKED | #717 D4 |

**Case (a) of #2724** means "rank 8 and every factor in col M". It is #2338's Theorem P (row C3), not a separate PR.
#2724 closes case (b) and nothing else.

## 3. The open window

All of it is under the dense B at h = 9 unless a row says otherwise. Nobody has exhibited a decoder costing 63 or
less in any of these shapes [read].

| # | shape | needed for ≤ 63 | junk / helpers | known constraints | priority |
|---|---|---|---|---|---|
| **O1** | **rank-8 realisation, m ≥ 10, factors outside col M** (#2338 §5.2 (d)) | one-sided cost ≥ **60** for every rank-8 decoder (section 4) [computed]; factor excess `s + t ≤ m − 8` (#1633); average term < 6.3 at m = 10; a term costs ≥ `9 − δ(R) − δ(A)` (#2338 Lemma 3.1), so cheap terms need supports that avoid ≥ 3 label classes between them | **unchanged** (Gram up to gauge) | normal form `D = RA + N1`, `T = RB + N2`, `ABᵀ = 2I`, `N1Bᵀ = 0`, `AN2ᵀ = 0`, `N1N2ᵀ = 0` (#1224 Thm A, merged); a two-sided ten-term witness exists at **82** (#1633); must leave C6 to C9 | first to price; window only 60-63 one-sided |
| **O2** | **m = 9, rank-9 realisation, square** (case (c), #2338 §5.2 (c)) | `Σδ ≥ 18` over the 18 factor supports (#2338) | **changed** (rank-one Schur defect): helpers must be re-priced | contains every co-star decoder; co-stars with ≤ 7 normals do not exist (C12), with 8 or 9 normals open; the frame-module theorem constrains every nine-term decoder; #2338's SAT runs timed out at 25-28 min | second |
| **O3** | m ≥ 10, realisation rank ≥ 9 | as O1 | changed | no tool | third |
| **O4** | label table of rank h ≥ 10 (non-dense B), outside the Gram case | **≤ 59** to tie at h = 10, **≤ 55** at h = 11 (the hub's pricing, not a theorem) | changed | open: non-Gram rank-8 realisations with col-M factors at h ≥ 10; joint tables; h ≥ 12 not attempted | low |
| **O5** | more than one layer (centres feeding centres, second descents, recombined copies) | — | — | outside every floor here, cx-medium's included; the shipped verifier allows descents only for the 8 barrier centres (#717 S4) | outside a single-layer search |

**The cost band, by accounting.**

| accounting | where it applies | lower edge for every single-layer decoder | window |
|---|---|---|---|
| **one-sided** (copy rank = source-frame dimension, copy descends to frame 0) | the hub's shipped checker and frontier words [read: #2338 §1, §3.7]; cx-medium and its skeptic read our template the same way | **52.36** for every realisation: proved by cx-medium over every subspace dimension 1 to 9 (Colab runs `cxm-clq-a`, `cxm-clq-b`), and **recomputed independently here** (`calc/onesided.c`, own C code and own frame list, every Tmax identical) | **53 to 63** in O2 and O3 (O4 uses other labels; not computed) |
| one-sided, **rank-8 realisation** (O1) | as above | **60** (59.07) [computed: `calc/lpmwu.c` weights, `calc/verify_mu3.c` exact integer re-check over all 8,283,457 subspaces] | **60 to 63** in O1 |
| **pairing rank** (projected or asymmetric copy) | SPEC-legal; "a checker convention forbids it" (#2338 row 22); readers must already contain the copy's frame, or pay climbs (one bit-side case priced at −3.57 % net) | clique counting gives **at most 16.9**, so no useful floor [computed] | anything up to 63, with reader costs unpriced |

## 4. Why 52.36 and 64 do not conflict, and why 52.36 fails under pairing rank

- **cx-medium's argument** (proved in its own scope, and upheld by its skeptic):
  - In each exact 8-frame the decoder must deliver rank 8, and each total adds at most 1 there.
  - Under the one-sided charge a total pays its whole gather span W.
  - The best W, of dimension 6, meets 1,155 of the 1,260 frames.
  - So the cost is at least `8·1260·6/1155 = 52.36` for every single-layer decoder, whatever its realisation or
    number of terms.
- **The hub's 64 is not a counting bound.** It uses the rank of the realisation and the factor spaces (Theorems R and
  P, case (b)). That is why it covers only the shapes in section 2 and does not contradict 52.36.
- **Under pairing rank a term pays only the part of its gather span that pairs with its scatter span**, so cheap terms
  can meet many frames [computed: `calc/a2clique.py`, `calc/a2clique_long.py`, `calc/check_rank1.py`]:
  - **Rank 1.** The term gathers 8 lines (P112 indices 48, 51, 60, 63, 96, 99, 108, 111) and scatters to 16. Its
    block `B[R,A]` has 8 equal rows, so its F2-rank is 1. Gather and scatter both meet **598** of the 1,260 frames.
    This was re-checked from root inner products alone, with no labels.
  - **Higher ranks.** Rank 2 meets 897 frames, rank 3 1,155, and rank 5 all 1,260.
  - So the counting bound is at most `8·1260/598 = 16.9`. #717 D1(e) had sampled 43.8; this search goes further.
  - This is a statement about the *proof method*, not a decoder. The terms are not part of any exact decoder.
- **Counting over every spanning 8-set does not rescue the pairing-rank count either.** This is valid in O1. But the
  same cheap terms meet a random basis 1.1 to 2.7 times more often than 1/8 of their rank [computed:
  `calc/a2basis.py`].

**The O1 floor of 60 in the one-sided accounting** [computed; the argument is ours and not yet checked by a skeptic].
- **Why bases count in O1.** A real realisation of rank 8 equals `D(2G)D⁻¹` up to row scaling (hub #1212 Theorem 1,
  merged: one torus orbit). So for any 8 lines S whose roots form a basis of R^8, the block `N[S,S]` has rank 8.
- **The demand.** Each term adds rank at most 1 to that block, and only if its gather set meets S. So every basis
  needs at least 8 terms whose gather sets meet it. The frames are the special bases with orthogonal roots, and they
  are the only ones forced for every realisation.
- **The bound.** Fix any weights μ on a pool of bases. A term whose gather labels span W costs at least dim W and meets
  the μ-share `μ{S : S ∩ L(W) ≠ ∅}`. So `cost ≥ 8 / max_W (μ-share / dim W)`, with the maximum over every subspace
  W.
- **The weights.**
  - Uniform weights on the frames give cx-medium's 52.36, which is the calibration.
  - Uniform weights on 1,280 random bases give 56.94.
  - Weights tuned by multiplicative updates on frames plus 6,800 random bases give 58.99.
  - Adding 2,500 bases built to avoid the worst gather sets (`calc/gen_bases.py`) gives **59.07**, with 12.8 % of the
    mass on frames.
  - The binding gather sets are 20-line sets with a 6-dimensional label span and 32- to 34-line sets with a
    7-dimensional span. Under frames alone, the binding set is a D5 (20 lines, root rank 5).
- **Exactness.**
  - The final weights were rounded to integers, and the maximum was taken exactly, by integer cross-multiplication,
    over all 8,283,457 subspaces (`calc/verify_mu3.c`; it gives `1073741831·8·6/872468135 = 59.0733`).
  - Every pool set was re-checked to be 8 distinct lines with a nonzero integer determinant, and the 1,260 frames to
    be orthogonal (`calc/check_pool.py`).
  - The weights are in `calc/runBig/lpmwu_mu.txt` and the pool in `calc/pool_big.txt`.
  - The 6,800-basis result, 58.99, is in `calc/runFB2/`.
- **Where the method stops.** It cannot exceed 64, because every basis meets each `N(f)`. Its true optimum over all
  bases lies between 59.07 and 64 and is not known. The gain is slowing: +0.08 from 2,500 targeted bases.
- **What it does not cover.** Realisations of rank 9 or more (O2, O3), where a basis block need not have rank 8. Also
  the pairing-rank accounting, multi-layer decoders, and complex coefficients beyond #1212's real statement.

## 5. Notes for ctr-search

- **Calibrations available from the hub** [read]. These are exact objects the search must reproduce or refuse:
  - the coordinate decoder: 64 in both accountings [computed];
  - #1180's nine-term identity: pairing cost exactly 65, a split of a proportional pair;
  - #1633's two-sided ten-term P112 witness: copy bill 82, rising to 90 under one gauge, so its normal form is not
    cost-preserving;
  - frame decoders: E8 72, E7 49, D6 30, D4 12;
  - #968's P104 frame-pair decoder: 56;
  - #717 D3: cost 63 on six single-line P111 punctures. This is a positive control that a sub-64 decoder can be found
    when one exists.
- **Already closed, skip these:** fewer than 10 terms with a rank-8 realisation (C3 to C5); co-stars with ≤ 7 normals
  (C12); fixed-reader extensions (C7 to C9); root-code factors (C2); the Gram case at h = 10 and 11 (C11).
- **Price.** One unit of centre per call is worth about +2.8 % (an-ours; the hub's #717 P1 gives +2.81 %). Use it for
  O1 only. O2 to O4 change the junk or the window, so their helper circuit must be priced first.
- **Adjacent lever, not 112 ports.** #968's P104 unit costs 56 on 104 ports (0.538 per port, against 0.571). The hub
  priced that lever at a ceiling of ×1.054 (#1619) and found it negative once built (#717 D3). It is not in this
  window.

## Change log

- v1 16:24 EDT: first map from #2338, #2724, #2228, #717, the scope-map doc and Discussion 1267.
- v4 16:50 EDT: the O1 floor rises to 60 (59.07, larger pool with targeted bases, exact re-check); the pool premises
  were checked; the O1 row and the O4 scope were corrected.
- v3 16:46 EDT:
  - New O1 floor of 59 in the one-sided accounting, from a weighted count over bases (exact integer re-check).
  - cx-medium's 52.36 recomputed independently.
  - ctr-search advice updated.
- v2 16:31 EDT:
  - The accounting is settled from #2338 §3.7: the shipped format is one-sided, so 52.36 holds there.
  - The pairing-rank counting bound drops from 43.8 to 16.9, with an explicit rank-1 term.
  - Basis counting fails too.
  - Added: #1224's normal form, #1633's witness at 82, #2639's rank-11 Gram closure, and the P104 note.
