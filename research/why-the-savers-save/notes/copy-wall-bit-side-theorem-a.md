> Working note written by an AI agent during the round of 13 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# d-bitspan (agent d-g-9hyp): the bit side's 63 against off-pattern decoders (13 October 2026)

Folder alias: `../d-g-9hyp -> d-bitspan`. Python: the g-p112-decoder venv (`$VENV -I`), numpy + python-sat.
Lines are numbered as in the hub's s2-decoder `hamming_e8` ("geo" numbering). Two punctures are checked:
- `s2`  : 39,52,57,71,76,78,103,109,114 (the hub's s2 REPORT "shipped bit puncture"; the one named in the brief);
- `word`: 8,17,39,51,53,58,63,72,114 (#589's word puncture 10,15,51,53,56,66,72,89,111 mapped by map_check.py).
Their orthogonality graphs are not isomorphic (bound-check/g-bit-puncture/iso_check.out).

## Files
- `e8geo.py`      E8 lines, exact flat enumeration (integer normals), bit cliques, F_2 labels (U U^T = B0).
- `t_flats.py`    all flats of rank <= 7 with the cliques each meets (ranks 1-5 finished; stopped in rank 6). `t_cliques.py` counts cliques.
- `codim2.py`     E7 hyperplanes and the 4,900 root-orthogonal 6-flats (E6 = A2^perp, D6 = A1A1^perp): cliques missed.
- `caseA_cover.py` the pairwise cover lemma for case A (B = B0, reads outside the label space).
- `sat_notW0.py`  flat-free SAT: a decoder with K totals whose read span misses the label space W0?
- `sample_caseB.py` samples such decoders (K = 10) and prices each read span exactly (minimum-weight basis).
- `sat_light.py`  K = 10 read spans containing independent light vectors (single lines, rank-2/3 flats, root-orthogonal flats).
- `sat_fixed.py`  XOR decoder search with fixed flats (not run).
- `sat_caseA.py`, `sat_rflats.py`, `match_sat.py`  SAT models with reads confined to root-orthogonal flats.
- `price.py`      bit what-ifs with m-anatomy levers (control c0 = 1.7334270e-3 reproduced).
- `out/`          all outputs.
