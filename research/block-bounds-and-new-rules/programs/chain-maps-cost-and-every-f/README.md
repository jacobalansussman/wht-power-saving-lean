# Chain maps: their cost in OpenAI's model, and correctness at every f

**What it shows.**
- `ram_letters.py`: chain maps (`row ^= row << 1`, `row0 ^= row1 << 1`, their inverses) written as programs in a
  simulator of OpenAI's RAM model, with every operation counted: the work per entry stays constant as the number of
  digit columns grows (block sizes 2 to 4, up to 7 columns).
- `every_f.py`: "correct at every f" decided exactly as the equivalence of two weighted automata over Q(i) (pure
  Python, exact fractions). Part A checks for every f, by that test, the identity that a paid move seen through `row
  ^= row << 1` is a free shift, a one-column move and a chain of CZ gates; Part B shows why testing small f proves
  nothing (the map `row0 ^= row1 << e` composed with the transform agrees with it up to f = e and first differs at e +
  1).
- `skeptic/ram_sim.py`, `skeptic/partA_raw.py`: the skeptic's own simulator and raw-array checks.

**Commands.** `python3 -I every_f.py` and `python3 -I ram_letters.py` (pure Python, under 50 MB, seconds to a minute);
the same for the skeptic's two programs. Expected outputs: `../../data/chain-maps/out_every_f.txt` and
`out_ram_letters.txt`.

**Re-run for this revision:** yes, on 19 October 2026: `every_f.py` (6 minutes under load) and `ram_letters.py` (3
seconds) give output identical to the original (`../../data/chain-maps/rerun_2026-10-19.txt`). The skeptic's two
programs were not re-run.
