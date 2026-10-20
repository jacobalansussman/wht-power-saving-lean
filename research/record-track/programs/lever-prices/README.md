# Programs: the prices of the recursion lever

Python 3.9 or newer, standard library only; run each with `python3 -I -B` from this folder. Times are on an Apple M4
Mac mini; every run is under 100 MB. Every program here was re-run on 19 October 2026 for this revision, and its
output matched the saved one in `../../data/lever-prices/`.

All of them solve the same balance: b = 1 - τ, where τ is the root of
Σ (paid blocks) n_r r^τ = Σ (arrays) k_w^τ. With every job k_w = 73 this is the community's moment formula.

| program | what it does | command | time | output |
|---|---|---|---|---|
| `lever_table.py` | the table of the README, sections 4.6 and 4.8: uniform, half cut, full cut and the non-separated cut on T74, #2979 (published fill and every head bin empty), #3945 (both) and #3815, each root bracketed with exact rational interval arithmetic to width 1e-16; checks that every figure stated by the earlier rounds lies in its bracket, and finds each full-cut break-even by certified comparisons | `python3 -I -B lever_table.py out.json` | 25 s | `lever_table.txt`, `lever_table.json` |
| `iprice.py` | the rigorous root used by `lever_table.py` (logarithms by atanh series and exp by Taylor series, both with tail bounds; no floating point in any decision). Written by the gate skeptic of 14 October | library | - | - |
| `rescue.py` | the plan of 14 October on T74: controls (the community's figure, the all-cut and non-separated what-ifs), then every separated variant, including cutting the heads by d axes for every d | `python3 -I -B rescue.py` | 9 s | `rescue.txt` |
| `breakeven.py` | break-even head counts on T74 for the full cut (9,443) and the half cut (1) | `python3 -I -B breakeven.py` | 11 s | `breakeven.txt` |
| `rescue1042.py` | the separated variants on #1042 (it imports `rescue.py`, which prints its own table first; the #1042 rows are the last three) | `python3 -I -B rescue1042.py \| tail -3` | 9 s | `rescue1042.txt` |
| `gateprice.py` | the gate's own root, written from the formula: T74 prices against the number of heads, for the head-maximising fill proposed on 14 October | `python3 -I -B gateprice.py out.json` | 1 s | `gateprice.json` |

Other outputs in `../../data/lever-prices/`: `g1-prices.txt` and `g1-prices.json` (the gate skeptic's rigorous
brackets on T74, made from the community's T74 profile file, which is not copied here: `lever_table.py` repeats its
main T74 rows from the profile numbers), `half-cut-gate-g1.json` (that gate's verdict, including the check of
the 24 + 24 split on all 112 ports), `lever-price-3945.txt` and `hub-root-controls.txt` (the #3945 pricing of
19 October and its controls: the community's roots of its newest circuits reproduced to all 24 digits, made with the
lever package's `lever_price.py` and `controls.py`), and `price-2979.json` (the #2979 gate's prices).

The profiles inside `lever_table.py` are the community's published global profiles, quoted as numbers: T74 from
#1839's construction, #2979 from `research/rank7-joint-endframes/allocation/PROFILE.json` on the community's main
branch, #3945 from `research/staged544-terminal-gap-fill/owner/PROFILE.json` at #3945's head, #3815 from #3815's
staged package. The same profiles are reproduced from the fills by the head test, digit for digit
(`../../data/heads/`).
