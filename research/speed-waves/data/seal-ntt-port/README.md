# SEAL NTT benchmark tables (section 2)

The measured results of the NEON NTT port into Microsoft SEAL 4.6.2. The patch itself is not in this revision.

| file | what it holds |
|---|---|
| `apple-m4/m4-table.csv` | SEAL's own benchmark program on the Apple M4: n, benchmark, pairs, median µs with the setting off and on, change in %, the lowest and highest change within a pair, verdict (7 pairs; generic 64-bit Arm build) |
| `apple-m4/control-neon-alone.csv` | the same program timed with the NEON code switched off at run time (identical code layout and memory): the NEON code's own effect |
| `apple-m4/control-not-neon.csv`, `control-on-vs-off.csv` | everything except the NEON code (tables, checks, layout), and the two programs against each other |
| `apple-m4/layout-table.csv` | the code-layout check behind the CKKS decryption finding |
| `neoverse-n2/run1-gcc.csv`, `run1-clang.csv`, `run2-gcc.csv`, `run2-clang.csv` | GitHub's Arm machines (Arm Neoverse N2): job, operation, n, pairs, median µs off and on, change in %, the lowest and highest change within a pair; two runs on fresh machines, GCC 13.3 and Clang 18.1 |
| `neoverse-n2/sysinfo-arm64-gcc-run1.txt` | the machine description of one Arm job |
| `scripts/analyse.py`, `combine.py`, `derived.py` | the M4 analysis (medians over pairs, verdicts, the derived multiply-plus-relinearise rows) |
| `scripts/summarize.py`, `combine_runs.py` | the Arm-run summaries and the combined tables |

Change is (ON / OFF) - 1, so negative means faster. A row whose within-pair spread includes zero is marked "no
change". The raw benchmark files (Google Benchmark JSON, every repetition) are not in this revision.
