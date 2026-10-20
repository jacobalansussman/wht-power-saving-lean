# The parity lattice at block 2 (division-free units need 5/3 moves per array)

Every free step of the everything-free game keeps the lattice of Gaussian-integer vectors whose entries all agree
modulo 1 − i and whose entries sum to an even number. At f = 2, one paid move pushes out 6 units of it and each
array needs 10, so a unit that never divides by 2 needs at least 5/3 moves per array, and a 2 in 3 would have to
divide by 2. This is argued by hand and not checked by a skeptic; Theorem 2W later superseded it at block 2 for every
unit, divisions included. The method (2-adic lattices built from Reed-Muller codes) is what remains of interest.

| File | What it shows |
|---|---|
| `padic.py` | Exact 2-adic lattice potentials over Z_(2)[i] (uniformiser 1 − i), Smith valuations, the letters of the game. |
| `verify_lattice.py` | An independent re-check by a different code path (Hermite bases and membership by back substitution): the lattice is kept by every monomial generator, one move pushes out 6, the target K 10. |
| `montecarlo.py` | A third, unrelated check by sampling: membership tested from the lattice's definition, giving the same 10 and 6. |

Commands (run from this folder; Python 3 only, no numpy):

```
python3 verify_lattice.py     # seconds
python3 montecarlo.py         # about 2 minutes
```

Re-runs of 19 October 2026: `../../data/parity-lattice-block2/rerun-2026-10-19/` (there were no saved outputs of
these two checks; the agent's notes record the numbers 6 and 10). Comments inside the programs use the agent's track
label ("mw-xw2tight"); the section's Sources maps it.
