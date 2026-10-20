# The engine version of the slice's matrix test

**What it shows.** `sp.py` builds the frames and the engine's letters at blocks 2 to 4; `engE.py` and `t_e10.py` the transfer form
(Theorem E10); `t_sph_exact.py` the spherical values exactly; `modgame.py` the pinned component in an abstract model;
`hecke_*.py` and `framefake*.py` search for realisable "fakes" that pass the pinned test (2 in 4 and 3 in 6 at block 3,
2 in 5 at block 4); `fake_h_*.npz` are the saved solutions, checked by `check_fake.py`. `skeptic/` rebuilds every letter
as an actual matrix (`v_frames.py`, `v_e3.py`, `v_fakecheck.py`, `v_kant.py` for the exact solutions near the fakes,
`v_vpcheck.py`; `regen/` the regenerated block-4 case).

**Commands.** `python3 t_e10.py`; `python3 t_sph_exact.py`; `python3 check_fake.py fake_h_3_2_4.npz`; in `skeptic/`: `python3 v_frames.py`.

**Needs.** numpy and scipy; torch for `framefake_t.py` and `hecke_abs_t.py` (exact Jacobians).

**Time and memory.** Seconds to minutes; the block-4 runs were done on Colab machines.

**Expected output.** Their `out_*.txt`.

**Re-run for this revision.** No.

Nothing in this folder was re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
