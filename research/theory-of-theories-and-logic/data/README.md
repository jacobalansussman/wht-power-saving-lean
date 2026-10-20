# Data of the abstract rounds

Each folder holds the outputs and certificates behind one program folder of the same name in
[../programs/](../programs/), as they came out of the original runs on 14 and 15 October 2026. Files whose names end in
`_rerun.txt` are new outputs of the published programs from 19 October 2026, for programs whose original output was not
saved. Job-runner queue lines and numpy warning lines were removed from eight output files; nothing else was changed.

| folder | what is in it |
|---|---|
| [frame-graph-eigenvalues](frame-graph-eigenvalues/) | the eigenvalue table and cut caps |
| [one-exponent](one-exponent/) | lift profiles at blocks 3 to 6, the exact identity, the engine-priced caps; the block-3 and block-4 lifts in `lifts/` |
| [sp62-mixed-unit-certificate](sp62-mixed-unit-certificate/) | the integer matrices of Theorem M's certificate (`exact_mats_lam3.json`), the exact and skeptic outputs, the cut lists of the 30 irreducibles, GAP outputs |
| [venn-certificate-spare-row](venn-certificate-spare-row/) | the 89/47 certificate's cap and demands, the 180-certificate sweep and the mixed room |
| [chi-game-ceiling](chi-game-ceiling/) | the 13-row certificate for 169/57 (`certificate_q4_m3_K5_u_1D.json`), the lift of excess 18, both checkers' outputs |
| [exactification](exactification/) | the toy unit composed exactly, and the counterexample to the factor 2 |
| [clifford-cohomology](clifford-cohomology/) | the randomized exact checks of Theorem E+ and the control |
| [lift-gluing](lift-gluing/) | the counts behind E ≥ 68, the exact lift with excess 3 (`e3_Y.json`), the LP dual certificates, κ_m |
| [angle-law](angle-law/) | Theorem F up to block 40, Theorem H on hooks, frame-family checks |
| [near-miss-certificate-search](near-miss-certificate-search/) | all orbit cuts of the block-4 hook module, the hook-relative families, the skeptic's logs (cut 1434) |
| [zx-sum-free](zx-sum-free/) | the PyZX runs of the ZX track and the bug reproduction |
| [symmetric-block4-skeletons](symmetric-block4-skeletons/) | the relaxation verdicts of all 678 systems with every coefficient free, the 17 hits, the skeptic's independent verdicts, the open charts of g27_Fany4_s7, and the recipes of every system (`recipes_all678.tar.gz`) |
| [every-f-and-logical-form](every-f-and-logical-form/) | the engine groups, r′, the transfer counterexample, block-3 traces |
| [two-certificate-theorem](two-certificate-theorem/) | the six commuting-matrix certificates (`certs/`) and the checkers' outputs |
| [border-curves](border-curves/) | the three relaxation systems with exact border curves (msolve format, as inputs), their Singular versions, the curve start and the 80-digit refinement |

msolve 0.6.5 reads `x*x` as `x` (section 13 of the section README). The one msolve output from before the bug that is published is the relaxation table (`data/symmetric-block4-skeletons/relaxation_table.txt`, two primes): the bug cannot touch it, because relaxation systems never repeat a variable, and a skeptic re-decided 677 of its 678 systems with Singular over the rationals, with no disagreement (the 678th it decided modulo two primes, with the same answer). No other msolve output from before the bug was found is included.
