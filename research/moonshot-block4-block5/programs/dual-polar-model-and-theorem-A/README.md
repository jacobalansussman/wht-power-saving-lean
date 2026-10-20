# The exact model of the block-4 engine game, and the selective-advance census

**What it shows.**
- `sp.py` is the exact model: the frames of the block-4 engine game are the 2,295 Lagrangian planes of an 8-dimensional
  binary symplectic space (the dual polar graph of Sp(8,2)); a paid move is one step to a neighbouring plane, and the
  depth of a term is its distance from the input plane.
- `check_model.py` compares the model with the real 16 x 16 matrices of the engine on 400 random words (0 mismatches).
- `selective_sp.py` counts, over all 100,352 single-letter shared moves of the relevant kind, those for which one
  letter advances one term towards K while holding the other back (a "selective advance"): 43,008. This refuted an
  earlier observation that no such letter exists.
- `skeptic/vcensus.py` is the skeptic's independent recount from random real letters (same 43,008 and 57,344).

Theorem A itself (no 4th and 5th moves can finish a 2 in 7 after the start u1 = D b1 y1, u2 = D b2 y2,
u3 = D(r u1 + c r' u2)) is a written proof, not a program.

**Commands.** Build the small C library that `bm.py` loads before the first run: `cc -O2 -shared -fPIC -o libmp_mac.so mp.c` on a Mac, or `cc -O2 -shared -fPIC -o libmp_lin.so mp.c` on Linux, in the folder that holds `mp.c`. Then `python3 selective_sp.py`, `python3 check_model.py`, and in `skeptic/`, `python3 vcensus.py`
(numpy). 1 to 4 seconds each on an Apple M4 Mac mini (16 GB). On macOS, numpy 2.0 built against Apple's Accelerate library prints spurious `RuntimeWarning: divide by zero encountered in matmul` lines; they do not change any result. Add `-W ignore` to silence them.

**Expected output.** `out_selective_sp.txt`, `out_check_model.txt`, `skeptic/out_vcensus.txt`.

**Re-run for this revision:** yes, on 19 October 2026; identical outputs (apart from the timing line of `vcensus.py`).
