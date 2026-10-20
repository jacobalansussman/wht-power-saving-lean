# One shared move then two private layers: the f = 2 test and the exact decisions

**What it shows.** The shape is W private starts, one shared move, then two private moves per array (W in 3W + 1: 2 in
7, 3 in 10, 4 in 13 at block 4).
- `rep2.py` is the f = 2 block test: at two digit columns the widened free letters split the data space into three
  pieces (16, 120 and 120 dimensions), and "the program reaches the target" becomes a row-space test in each. Unlike the
  frame-level test it sees the zero-divisor failure. It is a necessary condition, and it runs in milliseconds.
- `exact_small.py`, `exact_sigma2.py`, `lemma5_check.py` and `exact_sweep.py` build small skeletons around the special
  positions that the track's lemmas force, write each one as a polynomial system over the rationals (program
  coefficients and readout coefficients as unknowns, both outputs), and hand it to Singular. A reduced Groebner basis
  {1} means no complex coefficients make the skeleton a unit. `exact_ctl.py` writes the controls.
- `sing*/` hold the Singular inputs; `out_singular*.txt`, `out_sweep22.txt` and `out_sweepsig.txt` hold the verdicts:
  25 + 26 + 7 + 418 + 279 = 755 systems INCONSISTENT, 0 consistent, the rest too large (CAP).

**Commands.** Build the small C library that `bm.py` loads before the first run: `cc -O2 -shared -fPIC -o libmp_mac.so mp.c` on a Mac, or `cc -O2 -shared -fPIC -o libmp_lin.so mp.c` on Linux, in the folder that holds `mp.c`. (in `lib/`). Then (numpy and scipy):
- `python3 exact_small.py 3 711 3000 sing` writes the 25 systems of `out_singular.txt` into `sing/` (8 seconds on an Apple M4 Mac mini (16 GB));
- `Singular -q sing/<file>.sing` decides one system (Singular 4; the scripts print INCONSISTENT or CONSISTENT);
- `python3 exact_sweep.py ...` runs the sweep and calls Singular itself (see its docstring for the arguments).

The Singular logs begin with an error about `LIB "poly.lib"`: that library was missing on the machine that ran them;
the systems use only `std`, which is built in, so the verdicts are unaffected.

**Re-run for this revision:** `exact_small.py 3 711 3000` on 19 October 2026: its log and all 25 Singular inputs are
byte-for-byte identical to the saved ones. Singular itself was not re-run (it is not installed on the Mac); the
verdicts are the ones computed on 14 October and re-run by a skeptic, also modulo a prime.
