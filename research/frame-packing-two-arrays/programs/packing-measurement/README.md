# packing-measurement: how much packing the known units carry

`fm_ledger.py` replays one call of a local word in the checker's gcert/1 format exactly (with fractions) in the
frame model, and computes the raw-coordinate ledger: packing P, branching B, descents D, sideways steps Z and idle
moves. `bw_ledger.py` does the same for the backward calls. `oa_ledger.py` measures OpenAI's network exactly at small
ground sets h, and `oa_closed.py` turns the counts into closed forms in h, evaluated at the size in OpenAI's paper.

**Inputs.** The E8 word and the block-110 circuit are in this repository: `tools/certificate/gcert1-e8-r783.json.gz`
and `tools/certificate/gcert1-p11-pr193.json.gz`. The three community-derived words measured in the README (U8X, X4X
and T74g) are not included; their SHA-256 hashes are in the section's Sources. `labels_check.py` is the anatomy
round's brute-force check of OpenAI's labels, copied unchanged; `oa_ledger.py` imports it.

Run from this folder, with an `out/` subfolder. The programs write their JSON there (`mkdir -p out` first). The
shipped outputs are in `../../data/packing-measurement/out/`. `summary.txt` there is the full table that the README
quotes.

| program | what it does | command | time |
|---|---|---|---|
| `calib.py` | calibration with `../ledger/ledger.py`: the ledger of the small exact units, the frame-model level formula checked exactly in the engine, and the gauge-exchange unit replayed in the frame model | `python3 calib.py` | 1 s |
| `fm_ledger.py` | forward calls (fresh, and with the two kinds of old sources) under four helper models | `python3 fm_ledger.py ../../../../tools/certificate/gcert1-e8-r783.json.gz E8r783 init=F1,F2,F3 none ext-each gates gates-one` | 3 s (block-110 circuit: minutes, more than 3 GB) |
| `bw_ledger.py` | backward calls B1 and B2 | `python3 bw_ledger.py <word> <NAME> B1 B2` | seconds |
| `pos_run.py` | distinct positions among the labels on each move | `python3 pos_run.py E8r783=<word>` | seconds |
| `pos_check.py` | the position formula checked on every frame of a word | `python3 pos_check.py <U8X word>` | seconds |
| `oa_ledger.py` | OpenAI's network, exact, at h = 6 to 9 | `python3 oa_ledger.py 6 7 8 9` | seconds |
| `oa_rows_check.py` | one invocation of OpenAI's network per stage, replayed with symbolic values | `python3 oa_rows_check.py` | seconds |
| `oa_closed.py` | closed forms in h, checked against the exact counts and evaluated at h = 22, 25 and 100 | `python3 oa_closed.py` (after `oa_ledger.py 6 7 8 9`) | 1 s |
| `summary.py` | the full table from the JSON files in `out/` | `python3 summary.py > out/summary.txt` | 1 s |

The path in the `fm_ledger.py` command assumes this section sits at `research/frame-packing-two-arrays/` in the
repository.
