# The 8-move program with an order-3 halving end

**What it shows.** A correct 8-move program for 2 arrays at block 4 whose end is an order-3 "halving" (the arrays are
finished from x1 + x2 and x1 - sigma x2). It is the plain cost, not a saving: in the input basis y1 = x1 + x2,
y2 = x1 - sigma x2 it is two private 4-move chains. The script solves the end exactly, re-checks the formal identity, and
tests the program on random complex arrays at f = 1 to 5. It also confirms that the same prefix does not finish with
private ends in the basis x.

**Command.** Build the small C library that `bm.py` loads before the first run: `cc -O2 -shared -fPIC -o libmp_mac.so mp.c` on a Mac, or `cc -O2 -shared -fPIC -o libmp_lin.so mp.c` on Linux, in the folder that holds `mp.c`. Then `python3 verify_ctl.py` (needs numpy). About a second on an Apple M4 Mac mini (16 GB), 0.4 GB.

**Expected output.** `out_verify_ctl.json` (the program, the raw-array errors, the private-end test); on screen: exact
solve yes, formal re-check True, denominators [2], raw relative errors 4e-16 to 7e-16 for f = 1 to 5, private end in the
basis x: no, no.

**Re-run for this revision:** yes, on 19 October 2026; the JSON output is identical to the saved one.
