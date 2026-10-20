# Programs shared with the community package

These three programs are byte for byte the ones in the research package that I prepared for the community's
repository (CrocSwap/integer-mult-bounds), `research/lever-mixed-targets-2979/tools/`. That package's registered
check re-derives the lever's price from the community's own files of #2979, and replays #3945 and #3945 with every
head kept the same way. The other programs it uses are in [`tools/lever/`](../lever/) (`g0_staged.py`, `pricer.py`,
`mixed.py`, `lever_price.py`, `price_ours.py`), also byte for byte.

| program | what it does |
|---|---|
| `transport.py` | my transport of a community word with borrowed helpers into the community's native canonical form, replayed exactly on all 632 columns in both orientations |
| `lower_w.py` | my lowering of that canonical word into my certificate format (frames, warps carried by their maps, residual widths equal to the word's widths helper by helper), with seven wrong versions of the word that it must reject |
| `rate_cert.py` | an exact rational certificate of the mixed rate inequality of Theorem 1G at a figure N/10^10, and its refusal at (N+1)/10^10, with the same enclosures as the Lean rate modules |

Python 3.9 or newer, standard library only. Each program's first lines give its usage. Status: *checked by
computer*; the Lean proofs of the same facts are in the lever packages under `package/`.
