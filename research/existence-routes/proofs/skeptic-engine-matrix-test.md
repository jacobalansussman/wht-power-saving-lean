<!-- header: existence-routes -->
# The skeptic's record: engine matrix test

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** This is an AI skeptic's record of its check, as it wrote it. Its verdict is the status this revision uses for the track it checked, unless a later verdict is named in the section README.

---

# v-e3-engineE: skeptic's check of e3-engineE

15 October 2026, 06:49 to 07:06 EDT. I am the skeptic for track e3-engineE of exist3 (`../BRIEF.md`). Nothing
was posted. Predictions are in `STATUS.txt` (written before I read anything), and every run is in `RUNLOG.txt`.
I imported none of e3-engineE's code. My programs rebuild the letters as actual N x N matrices and build the frames
from them. From e3-engineE I used only the saved numbers in its `.npz` files. The one exception is the block-4 case,
which I regenerated with its own search program and then checked with mine.

## Verdict

**Not refuted. The calibration holds as reported.** That calibration includes e3-engineE's own admission that the
refusal half fails at block 3. Every proof I checked line by line is correct. Every computed claim I rebuilt came
out the same. I found three points of wording, none of which changes the result:

1. The outcome label "calibrated-no-hit" is generous. No version of the test passes the refusal half of the
   calibration at block 3: engine 2 in 5 is not refused. e3-engineE says so plainly in its calibration field.
   "Calibration failed (refusal half), no hit" would be the more accurate label.
2. The "yes, it stops at a constant" answer is proved only for the bound e3-engineE derived, 2 moves per array.
   Upper limits on what the pinned component could ever prove are computed only at blocks 3 and 4. e3-engineE states
   this caveat itself (STATUS addendum, NOTE section 6).
3. "absorb4 and px-t514: these tests read only (m, W, M)" is loose. The pinned condition eZe = eKe is a property of a
   word, so it could be evaluated on a given program. These two items are untested rather than inapplicable. This
   does not matter for the verdict, because the refusal half already fails at engine 2 in 5.

## What I rebuilt, and what came out

| Claim | My check | Result |
|---|---|---|
| Conventions: dir e_1 acts on Pauli labels as the transvection T_(e_1,0), K as (x,b) -> (x+b,b), L0 = Z-type labels | conjugated every Pauli by the actual matrices (`v_frames.py`) | confirmed, m = 2, 3, 4 |
| Frames 15 / 135 / 2,295; the free group has exactly m+1 orbits, one per distance (Prop E4) | Lagrangians by brute force over isotropic subspaces, not by an orbit; orbits under the images of the actual free letters (shifts, phases, address maps); checked on every frame, not a sample | confirmed m = 2, 3, 4 (m = 4 in the Colab job) |
| Spherical values and the two ratios c_1/mu_1 = 2 - 2^(1-m), c_m/mu_m = (2/3)(1 - (-1/2)^m) (Theorem E2(ii)) | joint eigenvalues of the averaged operators E M(t) E and E M(K) E on the free-group-invariant vectors, built from the actual t and K | confirmed m = 2, 3, 4: 1.5 / 0.5, 1.75 / 0.75, 1.875 / 0.625; the vectors (mu_j) and (c_j) are not proportional |
| No extra double cosets from a scalar e^(i pi/4) outside B | every letter has entries in Q(i), so the scalars in G are the powers of i, which lie in B | confirmed (an argument; it means the frame module is the whole trivially pinned component) |
| Prop E3: K = h P_K h^-1 with h a product of letters | exact Gaussian-integer arithmetic: Had P_K Had = 2^m K and 2 S_j dir(e_j) S_j = (1+i) Had_j (`v_e3.py`) | exact, m = 1..4; rank equalities at f = 1..3 (m = 2) and f = 1, 2 (m = 3) |
| Lemma E0: dir z = g dir(e_1) g^-1 for an address map g | exact, every z != 0 | exact, m = 1..4 |
| Theorem E10 (transfer form and bordered rank identity) | my own random engine word (4 moves, 3 registers, free mixing by random combinations of free letters), m = 2, f = 1, 2 | transfer form 3e-13 / 1.3e-12; rank identity at 5 values of lambda: holds |
| Realised fakes 2 in 4 (block 2), 2 in 4 and 3 in 6 (block 3), two-sided | (a) every saved free element is in P: I enumerated all of P (48 and 10,752 elements) from the actual free letters; (b) the condition checked as an operator identity E M(Z) E = E M(K) E on the whole frame module, and on 5 fresh random complex vectors (`v_fakecheck.py`) | (a) yes; (b) max error 1.5e-15 / 7.8e-16 / 2.2e-15 on the whole module, at most 2.3e-14 on the fresh vectors; largest single term 57, so no cancellation of huge terms; a control target gives 0.07 |
| Exact solutions near the fakes | Newton-Kantorovich at the saved points: full-rank Jacobian (complex step), smallest singular value 0.075 / 0.099, sampled Lipschitz constant (`v_kant.py`) | h = 2.7e-11 and 9.6e-10, far below 1/2: an exact real zero lies within 1e-12 of each block-3 fake. The Lipschitz constant is sampled, not interval-bounded, so this is certified numerics in spirit. A crude a-priori bound would still pass by several orders of magnitude. |
| One-sided form admits the same block-3 fakes (Prop E6) | span of ALL of P applied to the move outputs (43,010 and 64,515 vectors), full SVD | rank 268 of 270 and 402 of 405, with a gap from about 0.09 down to about 1e-12; targets at distance 6.4e-13 and 6.9e-12, the same at every tolerance. The 2 (3) missing directions are the averaged ones, as Prop E6 says, so the pass carries over to the exact nearby point. |
| Block 4: 2 in 5 realised (the saved runs had no file) | regenerated with e3-engineE's `hecke_vp.py` (same arguments and seed as its jobF3, plus a save file) on Colab ssh4, then checked with my `v_vpcheck.py` | regenerated: distance 3.42e-14, identical to its log. My check: all 8 free elements fix L0 and are automorphisms of the dual polar graph (so they lie in P); outer couplings solved by least squares (residual 3.9e-14, largest coefficient 10.8); E M(Z) E = E M(K) E on all 5 invariant vectors and 3 fresh random ones, max error 2.2e-14: HOLDS |

## The proofs, line by line

- **Lemma E1, Theorem E10.** These are correct. The chain expansion keeps the order of factors. The Schur complement
  of the unit lower-triangular block I - G gives rank N(lambda) = Mk + rank(rho(Z) - lambda).
- **Proposition E3.** Correct. K is conjugate in G to the free letter P_K by h = prod_j S_j dir(e_j) S_j. So the
  zero-move "word" T_0 = P_K I_W has the similarity type of K I_W in every representation, through the same h each
  time. Any test that reads only similarity data therefore accepts it. A pinned test such as the trivial component
  does not accept it: E P_K E = E, not E K E.
- **Proposition E4.** Correct. For the trivial component: the kernel of G -> Sp(2m,2) is the Paulis times the powers
  of i, and that kernel lies in B. So B\G/B is P\Sp/P, which has m+1 classes (frame distance). Each class is closed
  under inversion, so the Hecke algebra is commutative. For the Weil component, the Mackey count is correct, and
  "at most m+1" is the right claim.
- **Theorem E2.**
  - (i) is correct. The distance-m coefficient needs chains of at least m moves, so its matrix has rank at most
    M - m + 1.
  - (ii) is correct.
    - The trivial character gives tau = I, because the trivial character of C[G] kills t - 1.
    - c_j is not 0 for j >= 1, so A and B have rank W.
    - The weights lambda exist because (mu_j) and (c_j) are not proportional (I recomputed both).
    - N = sum lambda_j R_j is the identity plus a strictly lower-triangular part, and A N B = 0. So W <= M - W.
  - Two facts make the compression well defined, and I checked both:
    - the diagonal of R_j is the same mu_j for every move, by Lemma E0;
    - e(...)e does not depend on how a frame permutation is lifted to C[B], because the kernel is normal in G and
      lies in B.
  - (iii) The abstract construction is correct algebra: 2 mu_j + (c_j - 2 mu_j) = c_j. It is labelled correctly as
    being in the abstract model only.
- **Proposition E6.** Correct. V_M is a C[P]-submodule, so the targets lie in it if and only if each isotypic part
  does. The averaged part is the two-sided condition.
- **"Theorem 10 is the slice's trivially pinned component"** (in the remark after E2(ii)). This is plausible. Theorem
  10's proof (`exist2/ex2-families/PROOF.md`) uses the same bordered Schur complement per character, with the shift
  average. But it is a structural analogy and is argued rather than written out. It belongs under [argued], not
  inside the [proved] theorem. This is minor.

## Coverage, checked against the files

- `out_t_e10.txt` ends with "ALL E10 CHECKS PASS". It also has the zero-move check at m = 2 (f = 1..4) and m = 3
  (f = 1..3).
- `out_t_sph_exact.txt` runs m = 2..14 and ends with "ALL EXACT CHECKS PASS".
- `out_sp.txt`: 15 / 135 / 2,295 frames, m+1 orbits, and the intersection arrays. Its orbit check samples about 50
  columns; mine checks every frame.
- `out_modgame_generic.txt`: at m = 3, W = 2, generic plays reach the targets at M = 6, but not at 5 or 4. One run (W
  = 3) crashed with "SVD did not converge", and the file shows it. The note quotes only the cases that finished.
- Colab outputs: jobF3 (block 4, 2 in 5, 3.42e-14), jobG (3 4 8, best 4.6e-4, not found) and jobA3 (block 6, abstract)
  agree with NOTE.md.
- No msolve or Singular was used, so the x*x check does not apply.
- e3-engineE says it ran on the Mac past the 2-minute rule twice. RUNLOG confirms both runs and says it stopped them.

## What stays open (as e3-engineE says)

- The realised reach of the pinned component from block 5 on.
- The exact identity in all m+1 frame representations (the whole frame module, not only its averaged vector or its
  first column).
- The one-sided form at blocks 4 and 5.

None of these is an existence or exclusion claim of e3-engineE.

## Sources

All paths are under `the archive's checks/wht43/`.

- **Read:** `exist3/BRIEF.md`; `exist3/e3-engineE/` (NOTE.md, PROOF.md, STATUS.txt, RUNLOG.txt, sp.py,
  check_fake.py, framefake.py, fake_to_column.py, hecke_vp.py, output files); `moonshot/theory-guide/RULES.md`
  section 0 and the hb.py docstring; `exist2/ex2-families/PROOF.md` Theorem 10; `exist2/NOTE.md` (the calibration
  list).
- **My programs and outputs** (this folder):
  - `v_frames.py` -> `out_v_frames.txt`;
  - `v_fakecheck.py` -> `out_v_fakecheck.txt`;
  - `v_e3.py` -> `out_v_e3.txt`;
  - `v_kant.py` -> `out_v_kant.txt`;
  - `v_vpcheck.py`, and `regen/` (job script `jobV.sh`; Colab output in `regen/remote-out/ssh4/out_jobV.txt`).
