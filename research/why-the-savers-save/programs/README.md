# Programs of the section "why the record savers save"

Every folder holds the programs of one track and, beside them, the outputs as recorded during the rounds of 13 to 15
October 2026 (`*.out`, `*.json`, `*.txt`, and `out/` folders). Running a program overwrites its outputs. Run each
program from its own folder: some read sibling files by relative path. Each program's header says what it does and
how to call it.

**Requirements.** Python 3.9 or later. Most programs use only the standard library. Programs that import numpy or
scipy need them installed (`pip install numpy scipy`); the SAT programs of the copy-wall round need python-sat
(`pip install python-sat`); the `.c` programs compile with any C compiler (`cc -O2 -o onesided onesided.c -lm`).
Programs that price a construction import the exact pricer `pricer.py` (my agents' re-implementation of the
community's moment engine), which ships once, with the tools of the lever package that the record-track section
describes: set `PRICER_DIR` to the folder that holds it. Three copy-wall programs also import `levers.py` (exchange rates
built on the pricer), which is in `pricing-levers/` and is found automatically. Programs that read a file of the
community say so in their header;
[../data/hub-inputs/README.md](../data/hub-inputs/README.md) lists those files with where to get them and their sha256.
Shared inputs of this section are in [../data/inputs/](../data/inputs/).

Statuses of the findings are in [../README.md](../README.md) and [../results.json](../results.json); this file only
says what each folder computes.

## Anatomy of savers (15 October)

| folder | what it shows | main commands | time |
|---|---|---|---|
| `anatomy-openai-unit/` | OpenAI's network: its family of units (`family.py`), the exact paid-move rule checked by brute force at small sizes (`frames_check.py`, `labels_check.py 7`), prices per direction against per child (`moment_price.py`), caps by ground-set size (`caps.py`: 25 points give 4.18e-10) | `python3 family.py`; `python3 frames_check.py` | seconds |
| `anatomy-hub-record-unit/` | the community's T74 and its checkpoints: free directions, centre charge, grouping (`anatomy.py`), ceilings (`caps.py`), the net free fraction 4.8924e-3 at every checkpoint since atm3 (`lineage.py`), window shares (`window.py`) | `python3 anatomy.py out/anatomy.json > out/anatomy.txt`; `python3 lineage.py ../../data/inputs/hub-checkpoint-profiles.jsonl` | under 1 s |
| `anatomy-our-lean-words/` | my 11 Lean-checked words: D = 320 per replica on every word (`anatomy_ours.py`), the per-call identity (`identity.py`), knobs and what-ifs (`knobs.py`), the stage-count knob (`stages.py`), the U8 word's classes (`u8word.py ../../data/inputs/U8X-parsed.json.gz`) | `python3 anatomy_ours.py anatomy_ours.json > run1.out` | under 1 s |
| `anatomy-synthesis-centres/` | the centre as an index code: copy ratios of the ADE root systems and E8 (`ade.py ade.json`), junk-clique windows (`clique.py`), the first-order cap and the template table (`cap.py`), the two-sided line on my U8 lowering (`twosided.py`) | `python3 ade.py ade.json > ade.out` | seconds |
| `anatomy-synthesis-skeptic/` | the skeptic's rebuild of the design table (`roots.py`), the cap check that found the 5.6e-3 ceiling (`capcheck.py`; makes `calc/` first), direct x-y gates in my lowerings (`xyprobe.py`) | `python3 roots.py roots.json`; `mkdir -p calc; python3 capcheck.py` | seconds |
| `second-line-two-sided-build/` | two-sided slack built exactly on my U8 word: calibration (`calib.py`), the crossing count (`twosided_local.py`), dirt gates (`dirtjoin.py`), minimum transitions (`mintrans*.py`, numpy and scipy), prices (`price_twosided*.py`), the Barnes-Wall screen (`bw16.py`); the hosting tool `hosty.py` is not shipped, its outputs are | `python3 calib.py`; `python3 twosided_local.py` | seconds |
| `second-line-two-sided-build-skeptic/` | the skeptic's own replay of the word (`vreplay.py`), calibration (`vcalib.py`), identity, direct gates, ports, its Barnes-Wall screen (`vbw16.py`) | `python3 vcalib.py` | seconds |
| `second-line-single-crossing-word/` | the crossing lemma on small cases: exhaustive search (`toybfs.py`), the quotient count and prices (`quotient.py`), Coxeter-Todd (`ct12.py`), reachability (`reach.py WORD.json.gz OUT.json`); `rerun/` holds the build's programs re-run with identical output | `python3 toybfs.py toybfs.json` | seconds |
| `second-line-single-crossing-word-skeptic/` | the deletion step on 120,000 random steps (`vproj.py`), the tie prices (`vprice.py`), multi-call search (`vmulticall.py`), banks (`vbank.py`), Coxeter-Todd from another construction (`vct12.py`, `vct12_phase.py`) | `python3 vproj.py`; `python3 vprice.py vprice.json` | seconds |
| `second-line-mirrored/` | the mirrored line per family: lines (`lines.py`), crossing counts (`counts.py`), twists (`twist*.py`; `all-twists/twistall.py` ran on a Colab machine), design and price (`design.py`, `price.py`), cleaning chains (`chaincost*.py`), screen (`screen.py`) | `python3 lines.py`; `python3 price.py` | seconds to a minute |
| `second-line-mirrored-skeptic/` | port facts and run bounds (`v1.py`), calibration and room (`v2price.py`), lazy-drop chains (`v3order.py`), readers and random twists (`v4.py`), the 120-per-call cleaning scheme, checked exactly (`v5perx.py`) | `python3 v5perx.py` | 1 s to 70 s |
| `small-block-window-wall/` | designs and their windows (`designs.py`), the centre floor for one-layer centres (`classfloor.py`), early arrivals (`early.py`), ceilings (`ceiling.py`), raw-array frames and networks (`frames.py`, `net.py`, numpy), bank exchanges by shear calls (`shears.py`, killed at 100 s for 5 banks in the original run) | `python3 designs.py`; `python3 classfloor.py` | seconds |
| `small-block-window-wall-skeptic/` | the skeptic's designs, class floor, raw arrays (numpy), ceilings and the piece-size sensitivity (`vceil.py`) | `python3 vdesign.py vdesign.json` | seconds |
| `medium-block-networks/` | the template at blocks 13 to 25: clique floors (`cliquefloor.py`, `cliquefloor9.py`, the 9-bit runs in `colab-runs/`), prices (`price_p112.py`, `price_t74_p1.py`), the cost breakdown (`items.py`), shrinking the template (`shrink.py`), ceilings (`ceilings.py`), natural decoders (`natural.py`), E7 calls and rigidity (`e7*.py`, numpy) | `python3 price_p112.py`; `python3 shrink.py` | seconds |
| `medium-block-networks-skeptic/` | independent prices (`vprice.py`), clique floors (`vclique.py`), hosting against banks (`vhost.py`), the even-window points at blocks 16 and 22 (`vsub*.py`, `vfloor8.py`, `vceil22.py n sigma c t jp`), raw E7 replay (`vraw.py`, numpy) | `python3 vprice.py` | 1 s to 25 s |

## The centre (15 October)

| folder | what it shows | main commands | time |
|---|---|---|---|
| `centre-floor-search/` | the theorem "no centre that keeps today's junk costs less than 64": the model (`p112.py`, numpy), calibration (`calib.py`: 64 and 72 found), the refused 48 (`negctl.py`), searches (`dsearch.py`, `bresp.py`, `sas.py`), ports per label dimension (`nd.py`), slices (`coflat.py`, `coflats_all.py`), counting (`profiles.py`, cross-checked by `xcheck.py`), the last cases (`m9.py`, `cover3.py`, `cover3b.py`, `h2k4.py`, `h2k4sym.py 2`), soundness on all exact decoders (`soundness.py`) | `python3 calib.py`; `python3 coflats_all.py && python3 m9.py`; `python3 profiles.py 63` | seconds to a few minutes each |
| `centre-floor-window/` | the counting floor 52.36 over the 1,260 orthogonal frames (`onesided.c`), the weighted floor 59.0733 (`gen_bases.py`, `lpmwu*.c`, `verify_mu*.c`; result in `runBig/verify.out`), the pairing-rank bound 16.9 (`a2clique*.py`, `check_rank1.py`, `a2basis.py`) | `cc -O2 -o onesided onesided.c -lm && ./onesided` | seconds; the weighted floor ran for minutes |
| `centre-floor-skeptic/` | the skeptic's own enumeration (`enum.c`), the 9-direction check (`m9check.py`), triples and K4 (`triples.c`, `k4.py`), profiles, calibration and the weighted floor (`mu.c`); its large intermediate lists are not shipped (they are regenerated by `triples.c` and `build.py`) | see each header | seconds to minutes |

## The copy wall and stage sharing (13 October)

| folder | what it shows | main commands | time |
|---|---|---|---|
| `copy-wall-floor-checker/` | the floor checker: calibration on E8 (bit P120 72, bit P111 63, complex P112 64), exact star and pairing floors for small units, evaluation of any decoder | `python3 copyfloor.py e8`; `python3 copyfloor.py eval <instance> decoder.json` | 4 s |
| `copy-wall-our-unit/` | the octad count (any decoder on P112 at least 53) and natural decoders (64 on P112, 72 on P120), what-if prices; uses `../e8-labels/e8.py` and numpy | `python3 toy.py` | minutes |
| `copy-wall-relays/` | the relay lemma's SAT encodings (`relaysat.py`, python-sat), witnesses, controls, the driver `run.py`, the no-gain random runs and the relay price | `python3 run.py 600 'P111 8 3 none'` (see `out/`) | minutes per case |
| `copy-wall-relays-skeptic-rerun/`, `copy-wall-relays-skeptic-model/` | the two skeptics: a clean re-run (`rerun-out/`), numbering cross-check, mixed absorption; the true-puncture re-run (`relay_true.py`, `numcheck2.py`) and row-free variants | see headers | minutes |
| `copy-wall-bit-side/` | Theorem A: case-A covers (`caseA_cover.py`), SAT runs (`sat_*.py`, `match_sat.py`, python-sat; batch scripts `batch*.sh` with `PYTHON` set to a Python that has it), flats and cliques, the price | `python3 caseA_cover.py word` | minutes |
| `copy-wall-bit-side-skeptic-rerun/`, `copy-wall-bit-side-skeptic-model/` | the clean re-run outputs and independent covers (`indep/`); the model skeptic's exact covers and SAT (`cover_exact.py`, `sat9.py`) | see headers | minutes |
| `copy-wall-isotropic-cuts/` and its two skeptics | the isotropic-cut classes and the counting lemma with its certificate on all 255 cuts; the skeptics' partition check and anisotropic cuts (`mine/aniso.py`) | `python3 s1_classes.py` | seconds to minutes |
| `copy-wall-small-span-totals/` | the octad Hall test (`hallsat.py`, `hallshape.py`), the 90 least-squares runs (`als*.py`, numpy), sub-units and controls | `python3 hallsat.py 7,7,7,7,7,7,7,7,7 60` | minutes |
| `copy-wall-stage-hand-over/` and its skeptics | the stage map replayed from the community's schedule (`stagemap.py`), the hand-count of routes | `python3 stagemap.py` | seconds |
| `copy-wall-idea-probes/` | the generators' probes: field searches at rank 8 over F_3 to F_7 (`codes/`), the real minrank 15 (`small/`), punctures and Petersen probes (`break/`), twins, Lagrangian and relay toys (`recombine/`), isotropic cliques (`pick/`) | see headers | seconds to minutes |
| `stage-sharing-frame-count/` and skeptic | the corrected frame-sharing count (6,470 of 12,210 pairs) on the community's bit geometry; the skeptic's independent count and the escape-(b) algebra | `python3 recount.py` (with the geometry file beside it) | 40 s |
| `stage-sharing-hand-over-toy/` and skeptic | the toy coded total in 3 moves with and without hand-over (`handover_toy.py 3`); the skeptic's rider check | `python3 handover_toy.py 3` | 40 s |
| `stage-sharing-riders/` and skeptic | the exact two-stage replay (`replay.py`), the integer bound of 18 on every word (`lb.py A`, `lb.py B`, numpy and scipy), its certificate check and controls; the skeptic's independent program `indep.py` (12 minutes) | `python3 replay.py`; `python3 lb.py A` | 1 to 2 minutes |
| `stage-sharing-permutation-door/` and skeptic | the permutation door on a tiny unit (`swapdoor.py`) and on the community's bit geometry (`hubcheck.py 6470`); the skeptic's carry repair (`indep.py`) | `python3 swapdoor.py` | 8 s; 35 s |
| `e8-labels/` | the E8 label tables of 10 October, used by `copy-wall-our-unit/` and `copy-wall-small-span-totals/` | imported | |
| `pricing-levers/` | `levers.py`: exchange rates of the community's complex side b (Kothe's U8 word) and bit side c, built on `pricer.py`; imported by `copy-wall-our-unit/whatif_shapes.py`, `copy-wall-relays/price_relay.py` and `copy-wall-bit-side/price.py` | `python3 levers.py` (with `PRICER_DIR` set) | seconds |

## The helpers (13 October)

These programs read the community's U8 or k8 word (see [../data/hub-inputs/README.md](../data/hub-inputs/README.md)).
The table gives no commands: each folder's `COMMANDS.txt` lists the commands as they were run during the round (with
the paths of the round's working folders), and each program's header says how to call it. The M-copy skeptics' programs
run from `helpers-m-copies/`, where their inputs are. Two light pricing programs were re-run for this revision (see
the table at the end); the rest were not.

| folder | what it shows |
|---|---|
| `helpers-census/`, `helpers-census-ours/` | the census of the 504 helpers (392 plane visitors, 112 M copies), the exact toy of the visit rules, the prices of removals on Capalbo's packing of U8 |
| `helpers-narrower-ends/` and its two skeptics | the frame cones of full-width helpers (exact no for all 378), the rank-7 pairs, the word builder and its prices; the skeptics' independent cones and the k8 prices |
| `helpers-half-reads/` | half reads on both sides; `bit/` holds the retirement of carrier 470 and its check |
| `helpers-bit-exits/` and its two skeptics | exits (impossible) and entrances of bit helpers |
| `helpers-m-copies/` and its two skeptics | the M-copy groups, chain covers, value-free covers; the skeptics' reach with the cut and the probes that refuted the floor as stated |
| `helpers-one-helper-two-planes/` | the plane visitors and why no helper serves two planes |
| `helpers-shared-pair-rows/`, `-build/` and two skeptics | the pair-row visit model (448), its block integer programmes, the literal build of one block, and the skeptics' helper-role traces that refuted the gain |
| `helpers-idea-probes/` | the generators' probes (fan-out, timelines, idle merges, block geometry, the search programmes) |
| `helpers-k8-census/` and skeptic | the same census and model on the k8 word (480 helpers), the what-if prices, and the skeptic's census that refuted the model gain |
| `helpers-k8-rank7-pair/` and skeptic | the price of the both-pairs assembly (`price_twin.py`, needs `PRICER_DIR`) and the skeptic's admission re-run; the full assembly used the community's own admission pipeline, which is not copied |

## Symmetry and other fields (13 October)

As for the helpers, the commands as run are in each folder's `COMMANDS.txt` where the track left one.

| folder | what it shows |
|---|---|
| `symmetry-lattice-units/` | Barnes-Wall, Leech and Reed-Muller line systems as units: clique bounds, the BW16 bit decoder, ceilings and what-ifs (numpy) |
| `symmetry-invariant-decoders/` | Kramer-Mesner searches for decoders invariant under symmetries of order 7 and 8 (`b/km7*.py`, `b/kmg.py`), controls at 64, packet decompositions (`a/`); its 71 MB flats cache is not shipped (regenerated by `b/frames.py`) |
| `symmetry-word-rewrites/` | superoptimisation of the k8 word: local rewrites, annealing, exhaustive windows, integer programmes (`lmilp.py`, `ilp.py`); reads the k8 word |
| `symmetry-stage-equivariant/` | the written argument that no total crosses a stage boundary, checked against the community's boundary files (`stage_equivariant.py`, reads the community's clone and the k8 word), and its exhaustive small cases (`python3 small_case.py`, about 5 minutes and 13 MB; `python3 small_case.py A Cs` repeats the quicker recorded run) |
| `symmetry-idle-regrain/`, `symmetry-walsh-packets/`, `symmetry-open-tests/` | the three confirmations |

## Re-runs for this revision (19 October 2026)

I re-ran the light programs below in a scratch copy of this folder on an Apple M4 Mac mini (Python 3.9.6; numpy
2.0.2 and scipy 1.13.1 where needed), each from its own folder. "Same" means byte-identical standard output, or an
identical JSON file, compared with the output recorded during the rounds.

| folder | programs re-run | result |
|---|---|---|
| `anatomy-our-lean-words/` | `anatomy_ours.py`, `knobs.py`, `identity.py`, `stages.py` | same |
| `anatomy-openai-unit/` | `family.py`, `frames_check.py`, `caps.py`, `moment_price.py` | ran; no output had been recorded; the printed figures equal the note's (best ground set 25 points, 4.18e-10; every frame identity holds for m = 2 to 5) |
| `anatomy-hub-record-unit/` | `anatomy.py`, `caps.py`, `lineage.py`, `window.py` | same |
| `anatomy-synthesis-centres/` | `ade.py`, `clique.py`, `cap.py`, `twosided.py` | same |
| `anatomy-synthesis-skeptic/` | `roots.py`, `capcheck.py` | same |
| `second-line-two-sided-build/` and skeptic | `calib.py`, `twosided_local.py`, `bw16.py`; `vcalib.py`, `vbw16.py` | same |
| `second-line-single-crossing-word/` and skeptic | `toybfs.py`, `quotient.py`, `ct12.py`; `vproj.py`, `vprice.py`, `vct12.py`, `vbank.py`, `vmulticall.py` | same |
| `second-line-mirrored/` and skeptic | `lines.py`, `price.py`; `v2price.py`, `v5perx.py` | same; `v5perx.py` had no recorded output and re-checks its scheme exactly (every check true) |
| `small-block-window-wall/` and skeptic | `designs.py`, `classfloor.py`, `early.py`, `ceiling.py`; `vdesign.py`, `vclass.py`, `vraw.py` | same |
| `small-block-window-wall/` | `shears.py` | the cases for 2 to 4 banks print as recorded; stopped after about 5 minutes in the 5-bank case, as the original run was |
| `medium-block-networks/` | `price_p112.py`, `shrink.py`, `price_t74_p1.py`, `items.py`, `natural.py`, `ceilings.py`, `cliquefloor.py` | same (`cliquefloor.py` by default computes a subset of the recorded designs; those agree) |
| `medium-block-networks-skeptic/` | `vprice.py`, `vclique.py`, `vhost.py` | same |
| `centre-floor-search/` | `calib.py`, `nd.py`, `coflat.py`, `coflats_all.py` then `m9.py`, `cover3.py`, `cover3b.py`, `profiles.py`, `xcheck.py`, `negctl.py`, `soundness.py` | same (`profiles.py` by default stops at 20 directions; those agree); `negctl.py` refuses the 48 decoder with 2,333 violations; `soundness.py` accepts all 13 exact decoders |
| `centre-floor-search/` | `h2k4sym.py 2` (the cases with 11 and 12 directions) | not re-run to completion for this revision: it ran past the 5-minute limit for light jobs on the Mac, and its run in the queue for longer jobs had not started when this folder was staged; the recorded output (no four sets with all four triples, so both cases closed) is shipped as it was |
| `centre-floor-window/` | `onesided.c` (compiled with `cc -O2 -o onesided onesided.c -lm`) | same: the counting floor 52.3636 per call over the 1,260 orthogonal frames |
| `copy-wall-floor-checker/` | `copyfloor.py e8` | same |
| `stage-sharing-permutation-door/` | `swapdoor.py` | same |
| `stage-sharing-riders/` | `replay.py` | same |
| `symmetry-stage-equivariant/` | `small_case.py` | part A same; parts C (the whole cover, 1,451,520 invocations) and B (h = 4, m = 7) ran to completion with no counter-case, where the recorded output of C had covered a sample of 2,001 invocations |
| `pricing-levers/` with `copy-wall-our-unit/`, `copy-wall-relays/`, `copy-wall-bit-side/` | `levers.py`, `whatif_shapes.py`, `price_relay.py`, `price.py`, with `PRICER_DIR` set to the lever package's tools | same (re-run after the fact check that shipped `levers.py`) |
| `helpers-census/`, `helpers-shared-pair-rows/` | `price_toy.py`, `price_market.py`, with `PRICER_DIR` set | same |
| `second-line-mirrored/` | `chaincost2.py` | same, apart from the input path it records |
| `anatomy-hub-record-unit/`, `symmetry-stage-equivariant/` | `lineage.py`; `small_case.py` (whole run, about 5 minutes, 13 MB) | same; same result as above, re-run a second time |

Everything else in this folder was not re-run for this revision: the SAT and integer-programming runs, the
programs that need the community's files, and the runs that took more than a few minutes. Their recorded outputs
are shipped as they were.
