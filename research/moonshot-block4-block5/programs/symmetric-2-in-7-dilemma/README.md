# The unit/zero-divisor dilemma: the two group-algebra facts

**What it shows.** Two exact checks in the regular representation of a cyclic group.
1. The cycle lemma: if a filter phi = sum_j a_j w^j must reproduce its own cyclic shift up to a scalar (the junk
   equation of the symmetric near miss), then a_j = c^-j with c^n = 1, so phi is n times an idempotent and its image at
   the kernel frame is a zero divisor, for every order n, odd included.
2. The sigma separation: 1 + s is invertible in the group algebra of Z_n exactly when n is odd, with inverse
   (1 - s + s^2 - ... + s^(n-1))/2.

The proof that the symmetric 2 in 7 is not a unit (each array's single route into its target must be invertible and a
zero divisor at once) is a written argument; this program checks only these two algebraic facts.

**Command.** `python3 cyclelemma.py` (needs sympy and numpy). Under a second on an Apple M4 Mac mini (16 GB).

**Expected output.** `rerun_cyclelemma_2026-10-19.txt`: for n = 2 to 6 every solution of the junk equation has the
form a_j = c^-j with c^n = 1 and a singular block at the kernel frame; 1 + s is singular for n = 2, 4, 6 and invertible
with determinant 2 for n = 3, 5, 7, with (1 - s + ... + s^(n-1))/2 as its inverse.

**Re-run for this revision:** yes, on 19 October 2026 (sympy 1.14, numpy 2.0); the track kept no saved output, so the
re-run's output is the one published here.
