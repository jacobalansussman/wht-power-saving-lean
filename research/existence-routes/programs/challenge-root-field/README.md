# The smallest saving that meets the challenge, and the root of unity it supplies

**What it shows.** `theta.py` computes 1 − θ from OpenAI's constants (paperM, paperW, paperDelta) at 60 digits: any proved saving
above it meets the challenge. `root_order.py` is a port of OpenAI's `PrimeSelection.lean` and `LengthSelection.lean`:
for each length n it prints the order d of the supplied root, so which roots of unity are powers of it; from 2^20 on,
d is divisible by 2^20·3·5·7·11·13·17, so ζ_8, √2, √3, √5, √13 and √17 are available. `skeptic/ord.py` and
`skeptic/ports.py` are an AI skeptic's independent port and its check of the switch at 2^20.

**Commands.** `python3 theta.py`; `python3 root_order.py`; in `skeptic/`: `python3 ord.py`, `python3 ports.py`.

**Needs.** Python 3 only.

**Time and memory.** Each under 2 seconds.

**Expected output.** `theta.py`: 1 − θ = 2.10643843004018e-13. The others: their `out_*.txt`.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 root_order.py` (1.3 s): matches the shipped output. Output: `rerun_2026-10-19_root_order.txt`.
- `python3 theta.py` (0.1 s): 1 − θ = 2.10643843004018e-13. Output: `rerun_2026-10-19_theta.txt`.
- `python3 ports.py` (in `skeptic/`) (1.2 s): same result; the shipped file has one more line, the runner's "[exited with code 0]". Output: `skeptic/rerun_2026-10-19_ports.txt`.
- `python3 ord.py` (in `skeptic/`) (0.1 s): matches the shipped output; the two key lines agree (308 lengths checked, none bad; 2^20 − 1 keeps d = n); the 300 random lengths differ from run to run. Output: `skeptic/rerun_2026-10-19_ord.txt`.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
