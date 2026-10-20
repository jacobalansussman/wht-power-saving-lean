# Other kernels and other free steps

Programs of the "playgrounds" track of 16 October: the hook certificate and its relatives for kernels other than the
Hadamard kernel. Python with numpy; ranks are exact over Q (`exact_hook.py`) or modulo one or two large primes (a rank
modulo a prime can only undercount). Saved outputs: `../../data/kernel-playgrounds/`.

| Program | What it shows | Command used for the saved output |
|---|---|---|
| `hookk.py` | the block-3 hook (f = 5) for seven 2 × 2 kernels at two primes: 14/5 for every dense kernel, 173/80 for the triangular class | `python3 hookk.py 3 4 > out_hook_m3_k4.txt` |
| `exact_hook.py` | the same over Q exactly: H 80/224; the shear and the binomial involution 80/173 | `python3 exact_hook.py > out_exact_hook.txt` |
| `scan_k.py` | other members of the hook family at blocks 2 and 3 | `python3 scan_k.py M K1 K2 ...` (saved: M = 2 with k = 1 to 3; M = 3 with k = 1 to 3 and 5 to 7) |
| `chi_scan.py` | dense kernels [[1,1],[1,c]] for several cross-ratios, three primes | `python3 chi_scan.py > out_chi_scan.txt` |
| `venn2.py`, `venn2big.py` | two-factor "Venn" certificates at blocks 2 and 3 | `python3 venn2.py 2`, `python3 venn2.py 3`; `python3 venn2big.py 3,3 3,4 3,5 4,4` |
| `dualcert.py` | certificates in Λ^a V ⊗ Λ^b V* at block 2 | `python3 dualcert.py 2 > out_dualcert_m2.txt` |
| `patterns.py` | every pattern subspace at block 2 with f = 3 and f = 4 (triangular best 1.673; Hadamard 1.827) | `python3 patterns.py 2 3`, `python3 patterns.py 2 4` |
| `allinv.py` | all eight monomial-invariant subspaces of the hook module | `python3 allinv.py > out_allinv.txt` |
| `mixcert.py` | symmetric-power certificates at block 3 (the saved run was stopped at 10 minutes; partial; nothing beats the hook) | `python3 mixcert.py` |
| `hook_q3.py`, `hook_block.py` | the qutrit kernel at block 2, and a memory-light hook for the qutrit kernel at block 3 (k = 2) and block 4 (k = 2, 3) | `python3 hook_q3.py 2 1 2 3 4 5 6 7`; `python3 hook_block.py WHICH M K` |

`out_chi1415.txt` (the cross-ratio 14/15 kernel [[2,3],[5,7]] at blocks 2 and 3) was made with `scan_k.py`'s
functions from the command line. Re-run for this revision: `exact_hook.py`, `hookk.py 3 4`, `allinv.py`, `patterns.py 2 3` and `patterns.py 2 4` re-run: same output (`patterns.py 2 4` took 16 minutes at background priority); the other scans not re-run.
