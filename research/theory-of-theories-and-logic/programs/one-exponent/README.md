# The exact identity behind "one direction costs two block calls", and the lift profiles behind the growth exponents

Section 1. `liftprof.c` reads a lift file (frames with their rational subspaces) and sums the dimension jumps over all pairs at each distance; the lifts `g3c.txt` and `g4c.txt` are in `data/one-exponent/lifts/`. The block-5 lift `g5c.txt` (36.6 MB) is not shipped; `gen.py` writes it again (`python3 gen.py 5 g5c.txt`, the same file byte for byte; its SHA-256 is in the section README). `gen.py` also writes `g3c.txt` and `g4c.txt` exactly as shipped.

Outputs of the original runs: [data/one-exponent/](../../data/one-exponent/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `archimedes.py` | one direction = two block calls, exact in Q(i) | `python3 archimedes.py` | 1 s | under 100 MB | Python 3 only | identical |
| `liftprof.c` | exact lift profiles (pair sums by distance) of a lift file | `cc -O2 -o liftprof liftprof.c -lm -lpthread && ./liftprof ../../data/one-exponent/lifts/g4c.txt 4` | under 1 s at blocks 3, 4 | under 100 MB | C compiler | identical apart from the input path |
| `liftsample.py` | sampled profiles of start lifts | `python3 liftsample.py 4 3000 1` | 2 s | under 100 MB | Python 3 only | identical |
| `gen.py` | the start lifts: frames of the dual polar graph and a lift of each (writes `g3c.txt`, `g4c.txt`, `g5c.txt` for m = 3, 4, 5) | `python3 gen.py 5 g5c.txt` | 33 s at m = 5 | 85 MB | Python 3 only | m = 3, 4, 5: byte-for-byte equal to the files used |

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
