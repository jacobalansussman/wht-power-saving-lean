# OpenAI's second mechanism, read as a price world

Programs of two tracks of 16 October that examined the second mechanism inside OpenAI's existence proof ("Finite
tensor savings and exact Fourier circuits", 25 September 2026). Pure Python with numpy; each runs in seconds and well
under 1 GB. Saved outputs: `../../data/openai-mechanism-price-world/`. Re-run for this revision: all six, same output.

| Program | What it shows | Command |
|---|---|---|
| `openai_objects.py` | Part A: every identity of the reflection argument at d = 1 to 4 (the involution, its nonzeros, the parity form, the pivot relation of the Jordan-Wigner matrices); Part B: priced at 2 per pair the chain is a tie, so the whole margin is the dense pair; Part C: the cell behind "price ≤ number of nonzeros"; Part D: the zero-overhead ceilings against the dense-pair price | `python3 openai_objects.py > out_openai_objects.txt` |
| `feedback_state.py` | the McMillan degree of (I + zJ)^(⊗f) is f·4^f (exact ranks mod a prime, f = 1 to 4); closing one loop around the whole power gives the wrong matrix | `python3 feedback_state.py > out_feedback_state.txt` |
| `feedback_exact.py` | the same closures in exact fractions at f = 2 (largest difference 27/8) | `python3 feedback_exact.py > out_feedback_exact.txt` |
| `clifford_rank.py` | in the engine game a Clifford operator on one array costs log2 of its nonzeros per row: free letters keep the count, a Hadamard changes its log by at most 1 (6,400 random letters, 160 circuits) | `python3 clifford_rank.py > out_clifford_rank.txt` |
| `ceilings.py` | the zero-overhead ceilings of a fixed-matrix version (1.58e-3 with OpenAI's overheads), the overhead and dense-pair variants, the Lean packing bound | `python3 ceilings.py > out_ceilings.txt` |
| `dense_pair.py` | the 4 × 4 dense-pair identity re-derived exactly; at 3 bits no column order makes B(u)^(⊗3) diagonal plus rank one | `python3 dense_pair.py > out_dense_pair.txt` |

The ceilings are prices computed under an idealisation (every overhead set to zero), not theorems.
