# The 89/47 certificate for rank-2 targets with one spare row, and the mixed-unit room

Section 5. `skeptic_exact_cap.py` computes the cap 752 exactly over the rationals; `indep_venn_check.py` computes the cap and both demands modulo two primes (a modular rank can only understate the rational one, so the demands are lower bounds and the exact cap is an upper bound, which is the direction the certificate needs).

Outputs of the original runs: [data/venn-certificate-spare-row/](../../data/venn-certificate-spare-row/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `skeptic_exact_cap.py` | exact cap 752 of the 89/47 certificate | `python3 skeptic_exact_cap.py 3 4 0,2` | 1 s | under 100 MB | Python 3 only | identical |
| `indep_venn_check.py` | cap and demands of the 89/47 certificate modulo two primes | `python3 indep_venn_check.py 3 4 0,2` | 37 s | 1.5 GB | numpy | identical |
| `mixed_room.py` | the mixed-unit room from the 180-certificate sweep | `python3 mixed_room.py ../../data/venn-certificate-spare-row/out_venn_sweep_m3.jsonl` | under 1 s | under 100 MB | Python 3 only | identical |

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
