# orbit-designs: closed cancellation patterns, and why orbits multiply junk

`orb.py` builds orbit designs in the frame game. Each repair adds a whole orbit of a free symmetry σ (σ applied 0, 1,
…, n − 1 times to one stored array, each copy with its own scalar) instead of one landing. Each design then goes
through the exact screens: the linearised position test, the position test at the drawn scalars, the top-state count,
the layered screens and `heightfilter.py`, and then the exact solve. These programs import `../frame-game/`
(`fpa.py`, `design.py`, `design_w.py`) and `../heights/heightfilter.py`.

Run every command from this folder with Python 3.9 or later and numpy. The outputs of the original runs are in `../../data/orbit-designs/` with the same names. Some programs write a JSON file next to themselves or into `out/`; run them from a copy if you want to keep the shipped files.

| program | what it does | command | time |
|---|---|---|---|
| `cal.py` | calibration: the plain schedule and the order-3 control pass and their cut versions fail; the parity law (the order-3 pattern works exactly when the order of σ is odd) | `python3 cal.py` | 1 s |
| `t_lemma6.py` | 1,000 exact landing tests at blocks 3 to 5 | `python3 t_lemma6.py` | 1 s |
| `orb.py` | the orbit design search | `python3 orb.py MODE m trials depth seed [kinds] [char]`, MODE in `two`, `two2`, `three`; the last wave used, for example, `python3 orb.py two 4 1000 4 151` | minutes to hours |
| `replay.py` | rebuilds one design of an `orb.py` run deterministically and decides it with a larger budget | `python3 replay.py mode m depth seed trial cap [kinds] [save]` | minutes |
| `growth.py`, `growth_summ.py` | the same skeletons repaired by single landings and by orbits, side by side: junk growth | `python3 growth.py 4 60 6 501` | minutes |
| `cycle_census.py` | how often a cycle orbit (wave 7's kind) fits after one landing | `python3 cycle_census.py m trials seed` | minutes |
| `bilin.py` | numerical search for special orbit scalars on designs that pass the linearised test (evidence only) | see its header | minutes |
| `summ.py`, `tally.py` | the verdict table over all outputs | `python3 tally.py ../../data/orbit-designs/wave*/*/out_*.txt` | 1 s |

`../../data/orbit-designs/waveA/`, `waveB/` and `waveC/` hold the design runs, one folder per machine job.
`growth-runs/` holds the growth comparison and `out_tally.txt` the final table.
