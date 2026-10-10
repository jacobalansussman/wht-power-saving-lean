These are the files for [../for-network-coding-readers.md](../for-network-coding-readers.md). They use the Python standard library only. Nothing here is part of the proof-assistant development, and no person has reviewed it.

- `two-arrays.md` is the written argument for statement C(ii) of that note and for the operator form of its statement A.
- `python3 -B check_two_arrays.py cert_E162_bitmodel.json` recomputes the table of `two-arrays.md` (135 vertices, three integer rows each) in about a second. It ends with a line that starts `PASS`.
- `python3 -B subspace_metrics.py` checks the K_2,3, K_3,3 and Γ_3,3 labels of statement A over several fields in about a second. It ends with `ALL CHECKS PASS`.
