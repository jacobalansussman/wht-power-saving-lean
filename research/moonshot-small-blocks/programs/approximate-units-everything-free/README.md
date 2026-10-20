# Approximate (border) units in the everything-free game at block 2

These programs check which identities survive the recursion at block 2 (the *kept* identities), the claims and
refutations about approximate units made in waves 3, 4 and 7, and the search for a "triangular" 2 in 3.

| File | What it shows |
|---|---|
| `check_chain.py` | The wave-3 skeptic's own code: the chain map on bit row 1 commutes exactly with the paid move D1 on row 0 (f = 1 to 8); its order 1, 2, 4, 4, 8, 8, 8, 8 grows without bound; it mixes columns; and the wave-3 relation R7 is exactly 0 at f = 1 to 6. Exact integer arithmetic. |
| `check_relX.py` | The same skeptic's toy approximate family whose cancellation uses that commutation: under the relative recursion its error is exactly 0 at f = 2 to 5, while generic sub-calls give 40 to 2,881. This is what refuted "the window holds no useful approximate unit". |
| `xw2lib.py`, `keptlib.py` | The block-2 game, and the same game with paid moves replaced by *fibrewise sub-calls* (each paid move acts as an arbitrary operator on its own bit row and the identity on the other), which models "whatever the smaller calls return". |
| `keptrel.py` | Which identities a fibrewise sub-call keeps: D1 D2 = D2 D1, and D1 commuting with every row-1-only step, are kept; D1 commuting with a row-0 flip, D1² = X1 and others are not. |
| `r8check.py` | The kept identity R8 (moves in directions 1, 2, 1, packing through row-0 projections (1 ± X1)/2) and its mirror. |
| `lamplighter.py` | The row-1 steps that commute with D1 contain the lamplighter group (shift = running XOR, lamps = single flips moved along by it), checked at f = 4 to 8. |
| `braid.py` | The braid relation D S D = S D S holds for the exact transform but is not kept by a general sub-call. |
| `propS_check.py` | Proposition S: a family built on R8 with the X-type halving converges to an exact 3-move program in disguise. |
| `vcheck.py`, `vcheck7.py` | The wave-4 skeptic's own rebuild (dense matrices, no shared code): R8, Proposition S, the lamplighter facts, the mirror running XOR (commutes with D1; with the running XOR it generates all of GL(3, 2) at f = 3, 2,520 elements at f = 4, at least 3.08 million at f = 5), and in `vcheck7.py` the Z-type halving that escapes Proposition S. `vcheck.py` writes its report to `vcheck_out.txt` next to itself. |
| `dual_toy.py` | Proposition Λ0: over the dual numbers, the two layer equations behind Theorem 2W admit an approximate 2 in 3 (explicit solution, exact integers). |
| `tri_als.py`, `tnorm.py` | The triangular search: exact block least squares for 3-move programs with outputs K x + T y and K y, fitted at f = 2, 3, 4 and re-evaluated on fresh probes at f = 2 to 5; `MU` (an environment variable) penalises the cross-term T; `tnorm.py` reports the size of T for a saved solution. |
| `rescore.py` | The wave-7 skeptic's re-evaluation of saved solutions on new probe inputs at f = 2, 4, 6. |

Commands (run from this folder; numpy needed):

```
python3 check_chain.py ; python3 check_relX.py             # seconds
python3 keptrel.py ; python3 r8check.py ; python3 lamplighter.py ; python3 braid.py ; python3 propS_check.py
python3 dual_toy.py                                          # seconds
python3 vcheck.py ; python3 vcheck7.py                       # a minute or two; vcheck.py about 1.6 GB
python3 rescore.py ../../data/approximate-units-everything-free/saved-programs/cont121chain.npz 1,2,1 chain K 9191 2,4,6
MU=0.01 python3 tri_als.py 1,2,1 chain 2,3,4 100 1 K 6     # one penalised search run: minutes; ran on Colab machines
```

The saved solutions in `../../data/approximate-units-everything-free/saved-programs/` are the five the skeptic
re-evaluated: `cont121chain.npz` (residual 0.017 with a cross-term 37 times the signal, the false descent),
`cont1212.npz` (the solvable 4-move control) and three penalised runs. The table of every search run is
`triangular_search_summary.txt` there (it was made from the raw run logs, which are not shipped). The batch whose
menu was lost to a machine outage used `MU=0.001 python3 tri_als.py <pattern> wide 2,3,4 80 3 K 6` for the patterns
1,2,1, 2,1,2, 1,2,3, 1,1,2 and 1,2,2; it is the first experiment of open question 1 in the section README.

Re-runs of 19 October 2026: `../../data/approximate-units-everything-free/rerun-2026-10-19/`. Comments inside the
programs use the agents' track labels ("v-w3-border", "xw-border", "v-xw-border", "w7-xwborder") and refer to their
working notes; the section's Sources maps them.
