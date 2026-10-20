# entropy-bounds: network-coding bounds for common-frame units (15 October 2026)

These programs back sections 8 and 9 of [../../README.md](../../README.md): the entropy linear programs (Shannon,
Ingleton and the Dougherty-Freiling-Zeger inequality) for codes on the frame graphs Γ_2, Γ_3 and Γ_4, the search for
lumpings of the frame graph that an LP could handle, the cut LP of a fixed small layout of sessions, Lemma W and
Lemma W+, and the exact screen for three-array layouts. The proofs are in
[../../proofs/entropy-bounds-proofs.md](../../proofs/entropy-bounds-proofs.md); the skeptic's report is
[../../proofs/entropy-bounds-skeptic.md](../../proofs/entropy-bounds-skeptic.md).

Everything needs Python 3.9 or later with numpy and scipy (the LPs use the HiGHS solver that comes with scipy). The
first-order runs on all of Γ_4 (`layout_cut_pdlp.py` and the drivers that call it) also need `highspy`. Run every
command from the folder that holds the program.

**Where things ran.** The small calibrations and checks run on a laptop in seconds. The large LPs (millions of
variables on all of Γ_4, 10 to 60 minutes each, several GB) ran on Linux machines of Google's Colab service, eight
processors and about 50 GB; the drivers in `lp-and-layouts/drivers/` and `skeptic/vrun.py` are the scripts those runs
used, and their outputs are in [../../data/entropy-bounds/colab-runs/](../../data/entropy-bounds/colab-runs/) and
[../../data/entropy-bounds/skeptic-colab-runs/](../../data/entropy-bounds/skeptic-colab-runs/) (file names start with the
name of the run). The large runs were not repeated for this revision. Do not start a large one on a machine with
less than about 16 GB of memory: one small calibration (`calib_small.py`, the undirected butterfly) already reached
3.4 GB.

## lp-and-layouts/ (the author's programs)

| program | what it does |
|---|---|
| `dpg.py` | the frame graphs Γ_1 to Γ_4 exactly: the Lagrangian subspaces of F_2^(2m) (3, 15, 135, 2,295), adjacency and distances |
| `nclp.py` | the entropy LP for the cost of a network code: one variable per symmetry class of sets of sessions and arcs, the closure reduction, Shannon rows, Ingleton and DFZ rows as cutting planes, lazy rows, and a first-order solve with a rigorous lower bound from the approximate dual |
| `calib_small.py`, `test_cp.py`, `test_closed.py`, `test_pdlp.py`, `time1.py`, `time2.py` | calibration on networks whose answers are known (a path; the directed butterfly, where the LP finds the coding gain 8/7) and timing tests |
| `k23gap.py` | the smallest demand sets on K_{2,3} where cuts and routing differ (four sessions) |
| `star.py` | the five-group lumping of Γ_2 (one frame, the three lines through it, the eight frames opposite it): routing 128, cut LP 112, 26 LP variables, ρ = 1.6 |
| `lay5_g2.py` | the block-2 layout of five sessions where the cut LP gives 9 against routing 10 (and Theorem 2W says 10), and how far it can be lumped |
| `cells_lp.py` | the entropy LP on actual lumpings of Γ_2 and Γ_3 |
| `quot.py`, `quot2.py`, `quot3.py`, `quot4.py`, `quot5.py`, `small5.py` | the searches for lumpings with a high ceiling ρ at blocks 2, 3 and 4 (the table in section 8 of the README) |
| `layout_cut.py`, `layout_cut2.py`, `layout_cut_q.py`, `layout_cut_pdlp.py`, `layout_search.py` | the directed cut LP of one fixed layout of sessions: one big LP (too slow), constraint generation with max-flow separation (exact), a fast certified bound by contraction, the first-order version with a rigorous dual bound, and a search for layouts with a low cut value |
| `dump_layouts.py` | writes out the frames of the two block-3 five-pair layouts with cut values 13 and 14 (`layouts_block3.json`) |
| `w3cut.py`, `w3probe.py`, `w3ident.py`, `w3compat.py`, `stretch2.py`, `lemma_calib.py` | session 2: the undirected cut LP with its dual, which cuts the duals use, the compatibility test, Lemma W for any number of arrays, and its calibration |
| `edef.py` | the cross deficiency e of a three-pair layout (Lemma W+: cut value at least 3m − e) |
| `w3image.py`, `w3screen.py`, `wscreen.py`, `kingcfg.py` | the exact screen: map every frame to its distances from one end of each pair, solve the cut LP on that small image, round the dual to fractions and check it exactly; `kingcfg.py` is a coarser relaxation that turned out too weak |
| `w3rand.py`, `w3rand2.py` | cut gaps for three pairs in random small graphs (they exist, always on non-compatible layouts) |
| `drivers/` | the scripts of the Colab runs (hill-climbing to hard layouts, the full Γ_4 solves, the Okamura-Seymour calibrations, structured layouts). They import the modules of the folder above. They were run against the versions of those modules at the time; the modules here are the final versions, and the drivers were not re-run for this revision |

Commands for the light programs, each a few seconds and under 200 MB, all re-run for this revision with the same
results as recorded in the original run log:

```
python3 test_cp.py        # butterfly multicast: full LP 7.0, cutting-plane 7.0
python3 test_closed.py    # butterfly 7.0 (full and closed LP), path 3.0
python3 k23gap.py         # 16 four-session demand sets on K_{2,3} with cut LP below routing
python3 star.py           # routing 128, cut LP 112, rho = 1.6, 26 ground elements
python3 lay5_g2.py        # full-graph directed cut LP 9.0 against routing 10
python3 dump_layouts.py   # rewrites layouts_block3.json (identical to the saved file)
python3 lemma_calib.py 2 3 100 1      # Lemma W calibration, block 2, three arrays (m W samples seed): 97 stretchable
                                      # layouts, every cut LP = 6 = W*m
python3 edef.py 3 3000 1              # e at block 3 (m samples seed): {0: 2352, 1: 639, 2: 9}, i.e. 78 %, 21 %, 0.3 %
```

## skeptic/ (the skeptic's own programs)

Written by an independent AI agent with its own frame graph and its own LPs; it imports the author's modules only
where it tests them (`nclp_comb.py`, `vineq.py`, and the index map in `vlemmas.py`).

| program | what it does | command | re-run |
|---|---|---|---|
| `vgamma.py` | its own Γ_m (15, 135, 2,295 frames; degrees 6, 14, 30) | `python3 vgamma.py 2 3 4` | not needed (used by the others) |
| `vsym.py` | the symmetry facts behind Lemma 1: transvections are automorphisms, the group is transitive on ordered opposite pairs (120 of 120; 8,640 of 8,640) and on arcs; the uniform cut value 3/2 at block 2 from all cuts | `python3 vsym.py 2 3` | same output |
| `vlemmas.py` | the index map onto `dpg.py`; the near-polygon claim behind Lemma W+ by brute force; explicit Lemma W packings checked exactly on random layouts; the Lemma W+ construction | `python3 vlemmas.py 2`, `python3 vlemmas.py 3` (`4` takes much longer) | m = 2, 3 same output |
| `vineq.py` | the Ingleton and DFZ rows of `nclp.py`: never positive on random linear rank functions, +1 on the Vámos matroid and on the Fano matroid over GF(2), as they must be | `python3 vineq.py` | same output |
| `vshannon.py` | its own Shannon LP: butterfly 7 against 8, combination network 9 against 10, no false gain on two unicasts | `python3 vshannon.py` | same output |
| `nclp_comb.py` | the author's LP on the combination network: 9 | `python3 nclp_comb.py` | same output |
| `vmisc.py` | the known gap layouts are not two-class stretchable; stretchable shares with more samples | `python3 vmisc.py` | same output |
| `vcut.py`, `vcut_tests.py` | its own cut LPs (compact flow form with its own rigorous Lagrangian bound) and the recomputation of the author's quoted values | `python3 vcut_tests.py` | not re-run |
| `vprep.py`, `vprep4.py`, `vgen.py`, `vrun.py`, `vfull4.py`, `vseed18.py`, `vtri.py`, `vshannon2.py` | preparation and Colab runs: `vprep4.py` writes Γ_4 in the author's numbering as two `.npy` files (5.5 MB, not shipped; regenerate them), `vrun.py` runs the cut LPs on them, `vseed18.py` replays one of the author's random screens, `vtri.py` builds the four-array triangle-plus-pendant layouts, `vshannon2.py` the Shannon LP on the three-cell lumping of Γ_2 (80) | see each file's first lines | not re-run |
| `jobs/` | the copy of `vshannon2.py` that the Colab run used, where it differs from the one above | | |
