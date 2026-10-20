# Block 2 needs two moves per array (Theorem 2W): the numerical facts behind it

Theorem 2W is a written proof (argued by hand, not checked; one AI skeptic re-derived it). These programs check the
facts it rests on, and the smallest case of the rank budget for approximate units (Lemma Λ1).

| File | What it shows |
|---|---|
| `check_orth.py` | The skeptic's own code. Exact normalised traces: a paid move seen through any free frame is nearly orthogonal to every free word, and the target is nearly orthogonal to every word with one paid move, with overlaps falling like 2^(−(f−1)/2) at block 2; at block 3 the same for the two-row target R against one-move words and for K against two-move words. Its free menu includes chain maps and their inverses, Toffoli across columns, CS, CZ, S, parity phases and non-monomial one-column steps. Controls (the target against a genuine two-move word) stay at 1. |
| `lemD.py`, `xw2lib.py` | The author's check of the same fact (a paid move against free words) with the block-2 game library, at f = 2 to 7 (exact traces to f = 4, sampled upper bounds above). |
| `check_lambda1.py` | The skeptic's exhaustive check of Lemma Λ1 over the two-element field in the smallest case (3 moves, 2 arrays, first order): every solvable configuration has exploding-readout rank b = 1 and base-input rank g = 1, and the exact case has no solution. Also a sample of the 3-in-4 case. |

Commands (run from this folder; numpy needed):

```
python3 check_orth.py 2 6 60 7      # block 2, f = 1..6     about 4 minutes, under 300 MB
python3 check_orth.py 3 4 12 8      # block 3, f = 1..4     about 1 minute (the saved output; 40 samples: 3.5 minutes)
python3 check_lambda1.py            #                       about 1.5 minutes
python3 lemD.py                     #                       seconds
```

Expected output: `../../data/block2-two-moves-per-array/out_orth_m2.txt` (the column `O'` falls 1.0, 0.71, 0.5,
0.35, 0.25, 0.18 for f = 1..6; `ctl` stays 1.0), `out_orth_m3.txt`, `out_lambda1.txt` ("solutions = 0" for the exact
case; "(b,g) counts {(1, 1): 4032}"), `out_lemD.txt`. The re-runs of 19 October 2026 are in
`../../data/block2-two-moves-per-array/rerun-2026-10-19/`.

Comments inside the programs use the agents' track labels (for example "v-w7-xwborder", the skeptic of the wave-7
approximate-units track) and refer to their working notes; the section's Sources maps the labels.
