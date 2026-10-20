# The Coxeter carousel and the busy-register wall

Programs of the "wildcard" track of 16 October. Python with numpy. Each program prints its result and also writes it
to `out_<name>.txt` in the folder it runs in; the saved copies are in `../../data/coxeter-carousel/`.

| Program | What it shows | Command |
|---|---|---|
| `coxeter.py` | the step c = (move along direction 1, then relabel the directions cyclically) has c^m = the whole transform and c^(2m) = 1 in Sp(2m,2) for m = 2 to 8, with the Petrie-polygon distances; the same on real arrays for m = 2, 3 and f = 1 to 3 | `python3 coxeter.py` |
| `carousel.py` | the carousel with r = m carriers computes the transform exactly on streams of 10 inputs (m = 2, 3, 4); the age residual with r = m − 1 | `python3 carousel.py` |
| `carousel2.py` | the border check: with r = m − 1 the read error stays of order 1 however large the scalars; the block-3 menu of one-move steps | `python3 carousel2.py` |
| `carousel3.py` | the 192 block-3 steps walk 3 distinct geodesics, 64 each | `python3 carousel3.py` |
| `busy_wall.py` | the busy-register wall on random scalar shadows: the Hankel rank never exceeds the number of registers; residuals at r = p(δ+1) and one below | `python3 busy_wall.py` (a few minutes) |
| `carry.py` | the carry identities behind the verdict on Hamming-weight phasing | `python3 carry.py` |

Re-run for this revision: all six, same output (`busy_wall.py` took about 10 minutes at background priority).
