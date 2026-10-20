# A-fc: a higher table on #2979's words with the full cut (checked by computer, not proved in Lean)

**Status: checked by computer, not proved in Lean.** The Lean result of this revision is the full cut on #3945 with
every head kept, 1.7806935e-3 (`package/leverF`). This folder holds a table for #2979's words on which the same full
cut prices higher, at **1.7887046e-3** (8.0e-6 above). It is not built in Lean.

## What it is

The words, order and dispatch are those of #2979 (Marius-Paul Dumitrean's joint rank-7 end frames), unchanged; only
the table of where helper storage goes (the fill) is new. I re-solved it for the full cut: first a class solution of a
mixed-integer program (`class-solution/SOL-MIP0.json`), then a literal placement with the construction of #3815 by
Erwan Secousses (ChaloupeMusic; its per-tick Hall program), run in a tolerant mode. The tolerant mode needed 9 spare roles (8 of width 6 at
call 553 and 1 of width 4 at call 0) and dropped 8 planned width-6 slots (`spare-roles.json`). However it was made,
the final table is checked as a table: it passes #2979's unmodified owner engine, and the engine rejects all seven of
its corruption controls.

- 16,576 literal heads (8,288 on each of children 8 and 10), counted by three separate programs.
- The full cut: the heads are removed (their arrays take a job of 25 = 73 - 48); the table has no other rank-48
  block, so nothing is split, and every paid rank is at most 23 < 25, so the recursion with two job sizes
  (Theorem 1G) applies. Theorem 1G itself and the general theorem of the full
  cut are proved in Lean (they are the ones `package/leverF` uses on #3945 with every head kept). What is not in Lean
  is this table's own instance: its heads-first certificate and the rate fact 17,887,045/10^10 are checked by
  computer; the cut fact on #2979's words is argued by hand from the template it shares with #3945 (where it is
  proved in Lean), not checked.
- The price: b = 1.7887045643e-3 by three routes (exact rational enclosures certified on a 1e-24 grid; 60 digits;
  an independent exact interval program at 70 digits), which agree to 3e-60. The grid figure 17,887,045/10^10 holds,
  and one step higher fails.
- The community's own figure for this table, with no lever, is 1.6897030e-3, below #2979's own 1.7290340e-3: this
  table pays off only with the full cut. That is why the community's admission rejects it, by design: its
  strict-gain steps compare that lower figure with higher ones.
- Kappa is not claimed: it is capped by the selected bit side (about 1.7334e-3).

## The community's checks, run unmodified on a fresh clone of their repository at 12fa8c5c

| check | result for this table |
|---|---|
| `verify.py` (three-stage local replay of the words, 86 controls) | pass |
| `check_allocation.py`: pins, words, providers, the PR2839 owner engine, the physical count 66,813, 5 rejection controls | pass |
| `check_allocation.py`: the exact moment auditor (rank mass, root bracket, exact bisection) | pass, root 1.6897030e-3 |
| `check_allocation.py`: strict gain over the earlier headline | rejected (strict gain) |
| `check_cap.py`: cap certificate, rank mass, local profile paid, root at most the cap | pass |
| `check_cap.py`: recorded gap (root within 1e-6 of the cap) | rejected: the same lower figure, measured against the cap |
| complete replay: exports, fresh retained bootstrap, the native stages, owners, parent lift, rebind, discovery (7 steps), outer ledger (63 positive checks) | pass |
| complete admission | rejected only at "candidate must strictly improve full exact pinned frontier headline" |
| release stage (run past that step by my own glue): its `verify.py` (139 strict controls rejected) and outer checks | pass; then rejected at "ledger kappa does not exceed PR2269" (strict gain) |

The same harness passed #2979's published table and design A completely (positive controls).

## Re-create and re-check

    # the fill, from the class solution (needs numpy and highspy; about 20 seconds):
    git clone https://github.com/CrocSwap/integer-mult-bounds.git hub && git -C hub checkout 12fa8c5c844be21bfefb38a19503fbd867a11244
    bash research/record-track/programs/a-fc/r2/run2979.sh hub hub out-afc

    # the price by three routes, from the head recount (standard library only; under a minute):
    python3 -I -B research/record-track/programs/a-fc/price/hx_price.py research/record-track/data/a-fc/price/hx_heads-Afc.json --label A-fc

The decompressed fill is byte-identical to `design/FILL.json.gz`; the compressed file differs only in the 4-byte
time stamp of its gzip header. `price/` holds the head recount (`hx_heads-Afc.json`, with its control), the three-route
price (`price-Afc.json`), the controls (design B full cut 1.7806935e-3, design A 1.7797095e-3, the community's figures
of #3945 and #2979 digit for digit) and the constants of Theorem 1G (`thm1g-constants.json`). `SHA256SUMS` lists every
file of this folder.

## Sources

Found by a re-optimisation round on 20 October 2026 (03:17-03:21 EDT); the community's checks re-run on a fresh clone
at 12fa8c5c844be21bfefb38a19503fbd867a11244 the same morning; the price re-derived independently the same morning.
Fill sha256 7fc766c6c3e1561cee60b258138a918020b96d85ea782a2b01ad0ee035fd082e (decompressed a79f48af...).
