# frame-game: the frame game, its exact solver, the position test and the packing-first designs

`fpa.py` is the frame game: the group algebra of Sp(2m, 2) with the stabiliser of the raw frame free. It has exact
Sp(2m, 2) arithmetic, transvections (the paid moves), free charts, the exact end solve over the rationals (closure,
peel, solve), the position test ("psi screen") and its linearised form, and a constructive landing chart. It uses
the standard library only.

`design.py` builds the two-array trunk designs (both inputs share one paid move, each climbs alone before and
finishes alone after) and repairs their junk by landing stored arrays on it one at a time. `design_w.py` is the
version in which one trunk move packs W inputs (3 in 10 at block 4).

Run every command from this folder with Python 3.9 or later and numpy. The outputs of the original runs are in `../../data/frame-game/` with the same names. Some programs write a JSON file next to themselves or into `out/`; run them from a copy if you want to keep the shipped files.

| program | what it does | command | time |
|---|---|---|---|
| `t_cal.py` | calibration: symplectic checks, the plain units pass, their truncations fail | `python3 t_cal.py` | 1 s |
| `t_o3ctl.py` | the order-3 control passes as 2 in 2m, and fails when cut by one move | `python3 t_o3ctl.py` | 1 s |
| `design.py` | the packing-first search (modes `dup`, `t2`, `psilin`) | `python3 design.py m trials depth seed [mode]`, for example `python3 design.py 4 1500 12 61 psilin` | minutes to hours |
| `design_w.py` | the W-input trunk designs | `python3 design_w.py 3 4 500 10 91` | about 30 min |
| `summ.py`, `summ2.py` | tallies over design outputs | `python3 summ2.py ../../data/frame-game/design-runs/out_design_*.txt` | 1 s |
| `t_psi.py`, `t_psilin.py`, `t_psi_repair.py`, `t_split.py`, `t_oneyes.py` | the position test on calibration words and designs; its linearised form with the scalars as symbols; the split relaxation | `python3 t_psi.py` and so on | seconds to minutes |
| `t_psibil.py`, `t_bilin.py`, `t_hit.py`, `t_hit2.py` | numerical searches for special coefficients (evidence only, needs scipy); `t_hit2.py` continues the design that was reported as closest (`../../data/frame-game/design-runs/hit_m4_s102_4.pkl`) | see each file's header | minutes |
| `t_lift.py`, `t_lift2.py` | the lift of a frame-game program to the engine at f = 1 mod 4, on raw complex arrays (imports `../ledger/ledger.py`) | `python3 t_lift.py`; `python3 t_lift2.py` | seconds |
| `t_blocks.py` | the block structure behind the frame game at block 2 (imports `../ledger/ledger.py`) | `python3 t_blocks.py` | seconds |
| `st.py` | single-trunk designs with designed packing (an early variant) | `python3 st.py 4 50 1` | minutes |

The "closest design" that `t_hit.py` and `t_hit2.py` examine was refuted by the skeptic: it fails the position test
exactly for every choice of coefficients (`../skeptic-checks/frame-game/` and `../../data/skeptic-checks/frame-game/out_vhit.txt`).
`hit_m4_s102_4.pkl` is that design, saved by `t_hit.py` during the round (the matching numerical run is
`out_psibil_4_150_12_102_dup_starts8_maxit300.txt`). It is a Python pickle, so load it only if you trust its source. The design-run outputs in
`../../data/frame-game/design-runs/` are named `out_design_<block>_<trials>_<depth>_<seed>[_mode].txt`.
