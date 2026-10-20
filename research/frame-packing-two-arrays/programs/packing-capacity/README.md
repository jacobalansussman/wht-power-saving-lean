# packing-capacity: what one paid move can carry

These programs import `../ledger/ledger.py`.

Run every command from this folder with Python 3.9 or later and numpy. The outputs of the original runs are in `../../data/packing-capacity/` with the same names. Some programs write a JSON file next to themselves or into `out/`; run them from a copy if you want to keep the shipped files.

| program | what it does | command | time |
|---|---|---|---|
| `vflow.py` | (needs scipy) Theorem A on known units: frame classes per register, the live-class flow, Hu's two-commodity LP for two arrays, and the exact irreducible packing P* over all relative coordinates | `python3 vflow.py 1`; `python3 vflow.py 2` | minutes |
| `vflowB.py` | Theorem A′ with coarser classes (image frames on the level path), raw coordinates | `python3 vflowB.py 3`; `python3 vflowB.py 4` | minutes |
| `geom.py` | exact geometry of one paid move on the frame graph: how many moves raise, keep or lower a term at each level (Lemma S), blocking sets, and the rigid-pair tables (exact at blocks 2 and 3, 300 samples at block 4) | `python3 geom.py 3 300` | 10 s, 0.45 GB |
| `tables.py` | what the proved floors allow in packing terms, and what 0.1 and 0.01 need | `python3 tables.py` | 1 s |
