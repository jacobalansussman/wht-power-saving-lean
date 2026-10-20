# The two-layer bound in the everything-free game (W + m moves), and the floor of approximate units

The two-layer bound (Theorem T′: every unit with W ≥ 2 arrays at block m needs W + m moves, in every game with bounded
free steps) is a written proof (argued by hand, not checked; one AI skeptic upheld it). These programs check its
numerical ingredients at block 2 in the everything-free game.

| File | What it shows |
|---|---|
| `xw2lib.py` | The everything-free game at block m (default 2): monomial free letters (relabellings, phases, controlled-S, chain maps and running XORs, parity phases, boundary letters) and the paid moves, as exact operators on arrays of 4^f entries. |
| `sanity2.py` | Checks of the library at f = 2 to 6: the parity-phase identity P = (1+i)/2 + (1−i)/2·Z, D·D = X, K unitary and flat, and the conjugation rules of CS. |
| `percol.py` | The layer-mixing relation: sandwiching the parity-phase identity between column-uniform words never mixes layers with three paid letters (0 of 1,296 cases) and does with four (72 cases). |
| `orth.py` | Lemma O numerically: the target between free frames against every word with one paid letter has normalised overlap at most 2^(−f/2) (exact traces to f = 5, estimates at f = 6 and 7). Heavier than the rest (it was run under the memory guard). |
| `relfind.py` | The relation finder: every identity valid at all f = 2 to 6 among menus of words with at most two paid letters (`main`: 457 words, 437 distinct; `orbit`: plus truncated orbit words of high-order free words; `ctl`: positive control). Peak 0.9 GB. |
| `mini.py`, `lemO.py` | The skeptic's own block-2 mini-game, written from scratch, and its independent check of Lemma O at f = 2 to 5. |
| `floor.py` | The skeptic's test of Theorem B (the integer floor of approximate units): fit a 2-in-3 shape at small f, then evaluate the same coefficients at f = 2 to 8; shapes with one deep move must stay at or above √(1/2) = 0.707. |

Commands (run from this folder; numpy needed):

```
python3 sanity2.py                       # seconds
python3 percol.py                        # seconds
python3 lemO.py                          # about 40 s, but about 2.6 GB of memory
python3 orth.py                          # several minutes, about 1 GB; not re-run for this revision
python3 relfind.py main                  # and: relfind.py orbit, relfind.py ctl; minutes, under 1 GB; not re-run
python3 floor.py ST 1,2,1 3,4 15 1 8     # also SR and CHW shapes; ran on a Colab machine; not re-run
```

Expected output: the files of the same names in `../../data/two-layer-bound-everything-free/` (the skeptic's outputs
carry the prefix `skeptic_`). Re-runs of 19 October 2026: `rerun-2026-10-19/` there.

Comments inside the programs use the agents' track labels and refer to their working notes; the section's Sources
maps the labels.
