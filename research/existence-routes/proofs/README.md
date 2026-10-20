# Working proofs of this section

These are the working proof files of the three existence rounds of 14 and 15 October 2026, as the AI agents wrote them,
with private file paths removed, and the records of the AI skeptics who checked them. Each file starts with a header
that gives its status now; where the header and the file disagree, the header (the later verdict) wins. No person has
read these proofs. Only the two Lean files of `../programs/lean-existence-counts/` are checked by a proof checker.
"I" inside a file is the agent that wrote it, not me.

| file | what it is about | status now |
|---|---|---|
| [what-the-challenge-accepts.md](what-the-challenge-accepts.md) | P1 scalar closure, P2 one program for all lengths, P3 the engine accepts an existential unit, P4 linear descent, P5 open-set solvability | P1 and P3 proved in Lean; P2, P4, P5 argued by hand (two skeptics); the "field" catch refuted |
| [dimension-counts.md](dimension-counts.md) | Lemmas A to F: why counting cannot prove existence here | argued by hand; two skeptics; whether a count can exclude is open |
| [tau-theorem-and-laser.md](tau-theorem-and-laser.md) | Theorem 1 (τ-transfer), Lemma 2 (lockstep), Theorem 3 (Φ barrier), Proposition 4 (Kronecker powers), Lemma 6 (truncation), Propositions 5 and 7 | Theorem 1 argued by hand, two skeptics on its extension; Theorem 3 argued, no dedicated skeptic, its broad corollary refuted; Lemma 6 argued, no skeptic |
| [finite-fields.md](finite-fields.md) | Propositions 1 to 6: one prime never decides; three routes that work; refusing stray solutions | argued by hand, one skeptic; certificates checked by computer |
| [certified-numerics.md](certified-numerics.md) | Krawczyk certificates, Theorem 9 in the block-2 stand-in, Lemmas 7, 8, 10, Computations 12 to 14 | argued by hand and checked by computer, one skeptic, two summary sentences corrected |
| [families-commutative-slice.md](families-commutative-slice.md) | Theorem 10 (the matrix E), Theorems 12, 17, 18, Lemma 11, Corollaries 13, 15, 19 | argued by hand with exact certificates, one skeptic; partly superseded by the star lemma |
| [limits.md](limits.md) | Theorems 1, 7, 8, 8′, 9: limits in f and in the number of address bits | argued by hand, one skeptic, two reasons corrected |
| [price-reformulation.md](price-reformulation.md) | OpenAI's price language, Propositions 1 and 3 | Proposition 3 refuted; Proposition 1 open |
| [duality.md](duality.md) | the pipelining lemma, Strassen's theorem applied, the boosting corollary, certificate A | argued by hand, one skeptic |
| [star-lemma-and-slice.md](star-lemma-and-slice.md) | Theorems A, C, D, E and the star lemma (its section 7) | argued by hand with exact checks, one skeptic |
| [engine-matrix-test.md](engine-matrix-test.md) | Lemmas E0, E1, Theorem E10, Propositions E3, E4, E6, Theorem E2 | argued by hand with exact computations, one skeptic |
| [no-phase-group-algebra.md](no-phase-group-algebra.md) | Lemmas S, R, C, V, N, Theorem A, Proposition P | argued by hand, one skeptic; one generalisation refuted |
| [boosters.md](boosters.md) | Theorem B1, Corollary B2, the walk condition, the survivors | argued by hand, one skeptic; the headline refuted as stated |
| [seven-in-23-exact-level.md](seven-in-23-exact-level.md) | the lift of the 7-in-23 survivor, the f = 2 blindness, the exact criterion | the exclusion refuted; the blindness checked by computer |
| `skeptic-*.md` | the skeptics' records for the tracks above | the verdicts this section uses |

## Working names used inside the files

| name inside the files | here |
|---|---|
| `ex-legit`, `ex-v-legit-compute` | [what-the-challenge-accepts.md](what-the-challenge-accepts.md), `../programs/lean-existence-counts/`, `../programs/challenge-root-field/`, `../programs/restriction-of-scalars-toy/` |
| `ex-dimension`, `ex-v-dimension-compute`, `ex-v-ex-dimension-logic` | [dimension-counts.md](dimension-counts.md), `../programs/dimension-counts/` (`skeptic-compute/`, `skeptic-logic/`) |
| `ex-laser` | [tau-theorem-and-laser.md](tau-theorem-and-laser.md), `../programs/tau-theorem-and-laser/` |
| `ex2-modp`, `v-ex2-modp` | [finite-fields.md](finite-fields.md), `../programs/finite-fields/` (`skeptic/`); its `data/` and `cert*/` folders are in `../data/finite-fields-systems-and-certificates.tar.gz` |
| `ex2-certnum`, `v-ex2-certnum` | [certified-numerics.md](certified-numerics.md), `../programs/certified-numerics/` |
| `ex2-families`, `v-ex2-families` | [families-commutative-slice.md](families-commutative-slice.md), `../programs/families-commutative-slice/` |
| `ex2-limits`, `v-ex2-limits` | [limits.md](limits.md), `../programs/limits/` |
| `ex2-reform`, `v-ex2-reform` | [price-reformulation.md](price-reformulation.md), `../programs/price-reformulation/` |
| `ex2-duality`, `v-ex2-duality` | [duality.md](duality.md), `../programs/duality/` |
| `e3-slice6`, `v-e3-slice6` | [star-lemma-and-slice.md](star-lemma-and-slice.md), `../programs/star-lemma-and-slice/`; `starlp2/remote-out/` is `theorem-E-certificate/` |
| `e3-engineE`, `v-e3-engineE` | [engine-matrix-test.md](engine-matrix-test.md), `../programs/engine-matrix-test/` |
| `e3-rmbuild`, `v-e3-rmbuild` | [no-phase-group-algebra.md](no-phase-group-algebra.md), `../programs/no-phase-group-algebra/`; the theory guide's `hb.py` is copied beside it |
| `e3-boost`, `v-e3-boost` | [boosters.md](boosters.md), `../programs/boosters/` |
| `e3-bpexact`, `v-e3-bpexact` | [seven-in-23-exact-level.md](seven-in-23-exact-level.md), `../programs/seven-in-23-exact-level/` |
| `nag-rerun`, `a2-nag`, `mw-sym4`, `meta2`, `meta3`, `c2-finite-group`, `px-t514` | rounds told in `../../theory-of-theories-and-logic/` and `../../block-bounds-and-new-rules/` |
| `rl-prove`, `rl-v-logic`, `rl-vprove-compute`, `bound-check`, the hub's T74 and #1042 | the record track, `../../record-track/` |
| `w5-blueprint`, `w8-bpgl14`, `moonshot` | `../../moonshot-block4-block5/` |
| `remote-out/sshN/`, "Colab #N" | job folders on Linux machines of Google's Colab service; not published |
| `BRIEF.md`, `STATUS.txt`, `RUNLOG.txt` | the rounds' instructions and logs; not published |
