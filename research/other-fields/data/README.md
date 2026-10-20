# Data of this section

The outputs of the original runs of 10 to 15 October 2026, in folders named like the program folders in
[../programs/](../programs/). They are plain text or small JSON files; the largest is under 100 KB. Private file
paths, machine names and timing lines of the machines' wrappers were removed; nothing else was changed. `SHA256SUMS`
lists a checksum for every file here (`shasum -a 256 -c SHA256SUMS` from this folder checks them).

| folder | what it holds |
|---|---|
| [famous-problems/](famous-problems/) | the outputs of the 20 checks of sections 1 to 5, one file per program (`out_*` for the surveys' programs, `skeptic_out_*` for the skeptics'); `finite-geometry/out_c4_hirzebruch.txt` is the arithmetic of Hirzebruch's inequality, whose program was not kept |
| [quantum-dictionary/](quantum-dictionary/) | the dictionary checks and the skeptic's re-check (its first run, with four buggy tests of its own, is kept as `run_first_had_4_test_bugs.out`); the interval lemma and twist search at m = 2 and 3 and the explicit witness; the skeptics' counts (the files in `sandwich-skeptic-math/` were produced for this revision, because the original run kept no output); the E8 family searches, with `e8_joint_sat.txt` holding every SAT run (including the looser model in which 62 becomes satisfiable) and `pricing.txt` the what-if prices; the transfer checks |
| [entropy-bounds/](entropy-bounds/) | `lp-and-layouts/`: the JSON outputs of the author's lumping searches, screens and layouts (`layouts_block3.json` holds the frames of the two block-3 five-pair layouts with cut values 13 and 14); `colab-runs/`: the outputs of the author's runs on Colab machines, named `<run>__<file>` (for example `job_w3s4a__w3screen_m4_s11.json` is one seed of the random three-pair screen at block 4, `job_hardpd3__hardpd.out` one full Γ_4 solve); `skeptic/` and `skeptic-colab-runs/`: the skeptic's outputs, with each run's task list (`job_<x>__tasks.json`) and output (`job_<x>__out.txt`) |
| [preprint-reading/](preprint-reading/) | the five outputs of the preprint reading; re-runs for this revision reproduced all five byte for byte |

Not included: the frame graph Γ_4 as `.npy` files (5.5 MB; `programs/entropy-bounds/skeptic/vprep4.py` regenerates
them), copies of programs that the Colab runs carried with them (the final versions are in `programs/`), solver
caches, run logs and status files, and any text of the two preprints.
