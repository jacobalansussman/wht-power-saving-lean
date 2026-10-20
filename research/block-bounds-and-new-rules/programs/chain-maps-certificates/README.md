# Chain maps: certificates, readings and the parity identity

**What it shows.**
- `block2-27-lines/lines27.py` (numpy): at block 2 the frames of the reading at alpha are the 27 lines on a cubic
  surface; the Fermat cubic gives lifts whose counts equal the frame distance on all 729 ordered pairs (over C, and
  modulo 7 and 13): a perfect certificate, 2 moves per array.
- `block2-readings/s1s2_check.py`: exact checks of the two lemmas about readings; `residue.py`: the next readings (the
  16-element field, 325 frames) and the structure that Hirzebruch's inequality rules out; `f8_pairs.py`: transitivity
  on the 4,745 planes of the 8-element field. `skeptic/`: the skeptic's re-runs.
- `block3-frames/hookframes.py`: the block-3 hook as a family on the 135 frames (counts 80, 144, 224 at distances 1,
  2, 3); `geom891.py`: the 891 frames of the reading at alpha at block 3; `fermat.py`, `fermat_local.py`: the 405
  Fermat planes (exploration modulo a prime); `cert_check.py`: the exact check over Q(ω) of the lift family that gives
  231/131 = 1.763 (it needs the input files written by the track's `prep_lift.py`, which are not published; see
  below). `skeptic/`: exact re-checks of the hook frame counts and of the 891-frame lift.
- `parity-identity/relcheck.py` (pure Python, exact): with `row ^= row << 1` free, i^(parity of row 1) equals (1+i)/2
  + (1-i)/2 × (the column-wide sign flip of row 1) at f = 1 to 6, and the f = 1 shortcut does not represent the
  relation.
- `refuted-banded-lemma/e2_counter.py` (pure Python): the skeptic's counterexample to the lemma on families with
  banded label action (two global-parity maps whose product leaves the class).
- `period-2-maps/hook4odd_b2.py K`, `demcalc_b2.py`: chain maps that repeat every 2 columns are the engine at block 4
  with paired directions; the best hook gives 1.7826 pairs per array.

**Commands.** Each as `python3 -I NAME.py` from its folder (numpy where imported); `lines27.py`, `relcheck.py` and
`e2_counter.py` take seconds; the 891-frame scripts minutes. Original outputs in `../../data/chain-maps/` (skeptics'
outputs in `skeptic/`).

**Not published.** The lift input files of `cert_check.py` (generated) and the track's annealing program.

**Re-run for this revision:** yes for `lines27.py` (ALL CHECKS PASS), `relcheck.py` and `e2_counter.py`
(COUNTEREXAMPLE TO E2 AS WRITTEN CONFIRMED), on 19 October 2026, with output identical to the original
(`../../data/chain-maps/rerun_2026-10-19.txt`). The reading checks, the 891-frame scripts, the period-2 scripts and
the skeptics' programs were not re-run.
