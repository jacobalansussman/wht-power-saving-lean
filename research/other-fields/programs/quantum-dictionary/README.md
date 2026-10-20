# quantum-dictionary: frames as stabiliser states, the interval lemma, twisted cancellation (13 October 2026)

These programs back sections 6 and 7 of [../../README.md](../../README.md): the dictionary between the transform
game and the stabiliser formalism, the shortest-path ("interval") lemma, the twisted cancellation found in the
engine's wider game, the rank bound that makes a pair save nothing, and the twisted label families tried at E8
scale. Nearly everything is exact arithmetic over F_2, over the Gaussian rationals, or with exact matrix entries.

Run every command from the folder that holds the program, with Python 3.9 or later. Programs marked "numpy" need
numpy; `e8_joint_sat.py` needs the `python-sat` package (`pip install python-sat`); the others use only the standard
library. Several programs write their output into an `out/` folder next to themselves; the outputs of the original
runs are in [../../data/quantum-dictionary/](../../data/quantum-dictionary/), in a folder of the same name. Programs in a
folder whose name contains `skeptic` were written by an independent AI agent, with its own code, to break the claims.

`lib/` holds two helper modules that the programs import: `frames.py` (the frame graphs: all Lagrangian subspaces of
F_2^(2m), their distances and the subspace frames) and `copyfloor.py` (an exact small model of the copy floor, the
cost of the scratch arrays that collect and deliver terms in the community's decoders). Both come from earlier rounds
of this project and are copied here unchanged apart from comments. `python3 frames.py` in `lib/` runs its own checks.

"Re-run" says what happened when the program was run again for this revision (Apple Mac, Python 3.9.6, numpy 2.0.2).

## dictionary-check/ and dictionary-skeptic/ (the dictionary table)

| program | what it checks | command | time | re-run |
|---|---|---|---|---|
| `dictionary-check/dictcheck.py` | every row of the dictionary at small m: shifts and sign patterns are Pauli X and Z; the free steps generate the Hadamard-free Clifford group (32 and 3,072 elements at m = 1, 2; label images 2, 48, 10,752 at m ≤ 3); a paid direction is √X after a CNOT circuit; frames are Lagrangians (3, 15, 135, 2,295) and stabiliser states (6, 60, 1,080, 36,720); distances, overlaps, the E8 Pauli map | `python3 dictcheck.py` (sections A to G; `python3 dictcheck.py C` runs one) | about 3 minutes, 260 MB | same output |
| `dictionary-check/dictcheck_extra.py` | OpenAI's frames Φ_U have label L_U; every Lagrangian is a pair (U, β) with β a symmetric form on U; subspace frames are β = dot product, CSS states β = 0 | `python3 dictcheck_extra.py` | 7 s | same output |
| `dictionary-skeptic/recheck.py` | the skeptic's independent re-check of the table, written from the definitions in the Lean engine rather than from the first program (numpy) | `python3 recheck.py` | about 1 to 4 minutes | same output (109 checks, 0 failures) |

`run_first_had_4_test_bugs.out` in the data folder is the skeptic's first run, in which four of its own tests had
bugs; `run_final.out` is the run after it fixed them. Both are kept because the first run is part of the record.

## literature-checks/ (the published facts the table relies on)

| program | what it checks | command | time | re-run |
|---|---|---|---|---|
| `dictcheck.py` | the counts from the quantum literature at small m: the number of Lagrangians, their split by distance from the start (the Bruhat cells), the 2(2^m − 1) neighbours of a frame | `python3 dictcheck.py` | seconds | same output |
| `convexcheck.py` | two dual-polar facts at m = 3: every frame on a shortest path between L and L' contains L ∩ L'; and the count of frames through a given isotropic subspace | `python3 convexcheck.py` | seconds | same output |
| `lempelcheck.py` | the least number of commuting paid directions whose product is a given symmetric phase pattern, against its Clifford rank (rank, or rank + 1 for alternating forms) | `python3 lempelcheck.py` | seconds | same output |

## sandwich-and-twists/ and the two sandwich skeptics (the interval lemma, twisted cancellation)

| program | what it checks | command | time | re-run |
|---|---|---|---|---|
| `sandwich-and-twists/twist.py` | the interval lemma on every ordered pair of frames, the excess of every frame, and the search for twisted pairs (one scratch total serving two ports whose labels are not orthogonal) over every symplectic twist: all 720 at m = 2, all 1,451,520 at m = 3 | `python3 twist.py 2`, `python3 twist.py 3` | m = 2: seconds; m = 3: about 4.5 minutes | same output for both |
| `sandwich-and-twists/witness.py` | prints the explicit m = 3 witness: labels 001 and 111 share one total at copy cost 2, every frame given by a basis of Pauli vectors | `python3 witness.py` | about 4.5 minutes | same output (it first prints the four lines of twist.py's preamble, which the saved file omits) |
| `sandwich-skeptic-math/counts.py` | independent counts of Lagrangians and intervals, without the shared frame library | `python3 counts.py` | seconds | run for this revision (no saved output; new output in the data folder) |
| `sandwich-skeptic-math/interval_check.py` | the interval lemma on all ordered pairs at m = 2, 3 and from 41 sources at m = 4 (numpy) | `python3 interval_check.py` | about a minute | as above |
| `sandwich-skeptic-math/witness_check.py` | the m = 3 witness, typed in from the printed frames and checked from scratch | `python3 witness_check.py` | seconds | as above |
| `sandwich-skeptic-math/timed_oldmodel_m2.py` | the skeptic's counter-example to "with timing respected even the old rule cannot share one total": one total serves two orthogonal ports at copy cost 3 at m = 2 | `python3 timed_oldmodel_m2.py` | seconds | as above |
| `sandwich-skeptic-model/indep.py` | the second skeptic's independent checks: all Lagrangians by brute force, the interval lemma and the excess histogram; the witness with a symplectic twist found by its own backtracking; the size of the twisted relation for m = 2 to 9 | `python3 indep.py 1`, `python3 indep.py 2 3` | seconds to minutes | same output |
| `sandwich-skeptic-model/indep_count.py` | interval sizes summed over all pairs equal the expected sum of subspace counts (179,955) | `python3 indep_count.py` | seconds | same output |

## twisted-families-e8/ (whole label families, the E8 decoder)

| program | what it checks | command | time | re-run |
|---|---|---|---|---|
| `e8_twist.py` | can one twist for the whole family (any Lagrangian label family) make the community's complex-side frame decoder cheaper? Per total: 56 lines, a clique of 7 pairwise non-orthogonal lines, the hub labels' F_2-rank 8 | `python3 e8_twist.py` | seconds | same output |
| `e8_minrank7.py`, `e8_minrank_check.py` | is there a rank-7 completion on one total (one copy move fewer)? exact backtracking; run on every total, two cliques each, uniform and per-total twists, and a sanity instance where a completion is known to exist | `python3 e8_minrank_check.py` | seconds to minutes | same output |
| `sub_twist.py` | the same question on sub-systems: E7, E8 without the axes (the hub's 112 ports), D8 and D4 | `python3 sub_twist.py` | seconds | same output |
| `e8_joint_sat.py` | a different twist on every total: a SAT encoding of a generous relaxation. Copy cost 63 is unsatisfiable, 64 is reached with the hub's labels. With `--loose` (one modelling condition dropped) 62 becomes satisfiable, so the result rests on that condition | `python3 e8_joint_sat.py min 7` (63: UNSAT), `python3 e8_joint_sat.py min 8` (64: SAT); see the data file for every run | seconds to about a minute each | not re-run for this revision (needs `python-sat`) |
| `lagfloor.py`, `lagmitm.py`, `lagv4.py`, `lagv4b.py` | small families: the exact copy floor of every 4-port pattern (218 up to relabelling) with payload F_3 and F_5, by enumeration and meet in the middle; no pattern goes below one copy move per port | `python3 lagv4.py 3`, `python3 lagv4b.py 3`, `python3 lagv4b.py 5` | seconds; F_5 about 1.5 to 6 minutes | same output (all three) |
| `lagframes.py` | random and structured twisted families at h = 3 and 4 with the exact frame condition, against the vector condition the other programs use | `python3 lagframes.py 3 60 4 rand` (also `stab0`, `stabF`, and `3 60 5 rand`, `4 30 4 rand`) | seconds to a minute; h = 4 about 6 minutes | same output (all five) |

`pricing.txt` in the data folder holds the what-if prices quoted in section 7 of the README (one copy move fewer per
call; one whole total fewer). They were computed with this project's copy of the community's pricing program, whose
control run reproduces the community's record of the time to the last digit. That program is not included in this
section; the prices are what-ifs for constructions nobody has built.

## transfer/ (what the quantum literature gives back)

| program | what it checks | command | time | re-run |
|---|---|---|---|---|
| `frames_geom.py` | frame geometry: the neighbours of the start frame by kind (half subspace, half not), the excess of every frame by its (U, β) description, the subspace distance of Koetter and Kschischang on subspace frames | `python3 frames_geom.py` | seconds | same output |
| `segnorm.py` | the Bravyi-Maslov (parabolic Bruhat) normal form in the engine's own matrices at m = 1, 2: every element is free · one block of rank r · free | `python3 segnorm.py` | seconds | same output |
| `pathbound.py` | for every group element g at m = 1, 2, the moves needed before and after it add up to m plus the excess of its frame | `python3 pathbound.py` | seconds | same output |
| `pauli_and_f.py` | the 120 E8 root lines as the odd-Y four-qubit Paulis: the largest pairwise non-orthogonal set (8), the 2,025 octads; and at m = 1 why the "every f" rule makes single-copy quantum identities useless here | `python3 pauli_and_f.py` | seconds | same output |
| `graphs3.py` | the 8 graph states on 3 qubits: this project's classes against local-complementation classes (they differ) | `python3 graphs3.py` | seconds | same output |
| `avgfam.py` | two quantum-natural label families for the averaging lemma of the fourth revision's network-coding note: both are far weaker than the best known table (1.75 and 1.25 against 2.561) | `python3 avgfam.py` | seconds | same output |
