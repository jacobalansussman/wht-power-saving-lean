# Two pieces suffice with an order-3 letter

**What it shows.** With sigma the cyclic shift of three address bits (a free letter of order 3), the two pieces
x1 + x2 and x1 - sigma x2 give back K x1 exactly: Z1 (x1 + x2) + Z2 (x1 - sigma x2) = K x1 with
Z2 = K (1 - sigma + sigma^2)/2 and Z1 = Z2 sigma. It checks this with explicit complex matrices at blocks 3 and 4, and it
checks that the same trick fails for an order-2 letter, where 1 + u is a zero divisor (part B).

**Command.** `python3 cex.py` (needs numpy). Under a second on an Apple M4 Mac mini (16 GB). On macOS, numpy 2.0 built against Apple's Accelerate library prints spurious `RuntimeWarning: divide by zero encountered in matmul` lines; they do not change any result. Add `-W ignore` to silence them.

**Expected output.** `out_cex.txt`: `x1-coeff == K: True; x2-coeff == 0: True` at m = 3 and m = 4, and rank 4 of 8
(m = 3) and 8 of 16 (m = 4) for the order-2 filter in part B.

**Re-run for this revision:** yes, on 19 October 2026; the output is identical to `out_cex.txt`.
