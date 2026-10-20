<!-- header: existence-routes -->
# The skeptic's record: finite fields

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** This is an AI skeptic's record of its check, as it wrote it. Its verdict is the status this revision uses for the track it checked, unless a later verdict is named in the section README.

---

# Skeptic's check of ex2-modp (solutions modulo a prime to complex units)

exist2 round, 14 October 2026, about 20:31 to 20:56 EDT by the Mac clock. Nothing was posted. My predictions were
written first in `STATUS.txt`; the log is in `RUNLOG.txt`. Labels: [computed] (a program and its output in this
folder), [read], [argued].

## Verdict

**Not refuted. The calibration holds**, with the limits ex2-modp declared: the E8 unit and the 8-move control word
were not run, and the block-2 stand-ins and px-t514 are blind.

Every existence claim and every exclusion claim I could rebuild came out the same with my own code. I also closed a
gap the track had disclosed: its polynomial systems existed only modulo two primes. I rebuilt them exactly and they
are the same systems, with one exception, which I certified separately.

## What I checked, and how

1. **The existence points** [computed]. I wrote my own parser and exact arithmetic (`vlib.py`, `vpoints.py`). I
   evaluated every claimed point exactly, and separately modulo both primes on the raw export files:
   - the 15 rational relaxation points;
   - the two points over Q(i) and Q(sqrt(-3)), using my own arithmetic in Q(t);
   - the 7 joint plain-unit chart points, with both arrays sharing their coefficients (the 8th chart has an empty
     core);
   - the 3 special-position controls.

   All 27 are exact zeros.
2. **The impossibility certificates** [computed]. `vcert.py` and `vall.py` re-check all 240 certificate outputs
   against the raw exports. My code does the parsing, the rational reconstruction and the exact identity
   1 = sum g_i f_i over Q, plus a second check modulo p2 on the raw p2 file.
   - 240 of 240 pass.
   - The other files in the certificate folders are the documented failed attempts (timeouts, "not in module"),
     which were replaced by later runs.
   - In all 64 cases where rows were dropped as "extras", the dropped rows involve only the shared coefficients,
     as claimed.
3. **Coverage** [computed]. I recomputed decide.py's list of exclusion steps from its own verdict files (`vcover.py`,
   `vcover_hard.py`).
   - The 8 complete skeletons: 122 of 122 steps are certified. One of them is the constant core -1, which is trivial.
   - The snapshot of g27_Fany4_s17 and s19: 75 of 75. Each skeleton has one exported chart that had not been decided
     yet, as the note says.
   - ex2-modp's own derived systems, for the substitutions and the extras, are built exactly as described
     (`vderived.py`).
4. **Exact rebuild of the systems** [computed]. This is the link ex2-modp disclosed as an assumption. Its systems
   over Q were rebuilt from two primes, and nothing proved they are the true systems. I removed the primes from
   nag-rerun's builders and rebuilt everything in exact rational arithmetic (`fullsysQ.py`, `vrelaxQ.py`).
   - All 23 relaxations are identical to the exports. Every coefficient is plus or minus 1.
   - All 956 screen and chart-relaxation systems equal the relaxation plus decide.py's added rows.
   - 160 of 162 full-chart cores are identical (`out_vcompare*.txt`).
   - The 2 exceptions are f27_Fp2_m3_s39, chart 0_1_2_2, both arrays. That export was written at 19:50 by an older
     version of the builder, which picked pivots in a different order. Today's builder, run modulo p, reproduces my
     exact core and not the export.
   - I certified that chart's exact core myself. Singular on a Colab proposed a certificate, and Python checked it
     exactly (`out_s39_cert_check.txt`).
   - So the 8 exclusions now rest on exactly built systems. They still assume that the closure model, a2-nag's
     relaxation lemma and decide.py's split into charts are right.
5. **The theory** [read]. I checked PROOF.md's Lemma 0 and Propositions 1 to 6 line by line, including the Hensel
   steps in Propositions 2 and 3, Rabinowitsch in 6(f) and properness in 6(b). They are correct.
   - I checked the Krick-Pardo-Sombra bound on the paper itself (arXiv math/9911094, Theorem 1 on page 2). It has no
     extra hypotheses, and the constant is as quoted.
   - The many-primes count (2,636 primes) recomputes correctly.

## Minor findings (none changes a conclusion)

- **The special-position control is weaker than stated.** The note says random coefficients never find the control,
  which is solvable only on a4 = 3, and that the modular hints do. But the random search draws integers from
  [-2, 2], which cannot contain 3. With [-3, 3] it finds a4 = 3 at try 40 (`out_spu_box.txt`). The 17 real hits
  support the claim that hints matter: without hints the search finds 6 of 17 with box 2, and 2 of 15 with box 3.
- **"These exclusions hold over every field" is broader than what was shown.** It is true for the exported integer
  systems. Every elimination pivot is 1, so it also holds for the chart's closure system. But the unit problem
  itself involves i and 1/2, and the chart and relaxation arguments are stated over C.
- **The two irrational hits are not proved to need i or a cube root of unity.** That rests on a Groebner basis modulo
  p and on which primes carry solutions, which is evidence. The exact points over Q(i) and Q(sqrt(-3)) are proved.
  The absence of rational points is not.
- **Refusal is automatic here.** For a method whose only proof step is an exact check, calibration (b) cannot fail.
  Proposition 3's obstruction test never met a real stray solution, because none occurred. I confirmed by
  exhaustive search:
  - no point over F_2 for the border relaxation, but points over F_3;
  - nothing over F_2, F_3 or F_5 for engine 2 in 5.

  The plain relaxations do have points over F_2, which fits Proposition 7.3. That proposition still rests on the
  unchecked MAP Proposition G.

## Files

All in `the archive's checks/wht43/exist2/v-ex2-modp/`:
- **programs:** `vlib.py`, `vpoints.py`, `vcert.py`, `vall.py`, `vcover.py`, `vcover_hard.py`, `vderived.py`,
  `fullsysQ.py` (nag-rerun's fullsys.py in exact arithmetic), `vcompare.py`, `vcompare_hard.py`, `vrelaxQ.py`,
  `vscreens.py`, `vfp.py`, `mkmod.py`, `vcheckQ.py`;
- **outputs:** `out_*.txt`; the Colab run is in `sjob/remote-out/ssh4/`;
- **exact cores:** `fq/`.
