# Streaming cells

Programs of the streaming track of 16 October. Pure Python, exact arithmetic. Saved outputs:
`../../data/streaming-cells/`. Re-run for this revision: both, same output (the saved `out_calib.txt` lacks its first
three lines, the block-1 calibration, which the re-run prints).

| Program | What it shows | Command |
|---|---|---|
| `stream_calib.py` | the transfer-function evaluator of a streaming cell in the commutative slice: plain cells are exact at every quadratic form of blocks 1 to 5, a cell with one move deleted is wrong at exactly half the forms; the slice consistency test (no direction set below m moves passes at blocks 3 to 5); free scalar feedback and a hidden paid loop leave the cost unchanged | `python3 stream_calib.py > out_calib.txt` (about 15 s to 75 s) |
| `shapes.py` | the hook room of stage shapes one move below plain at blocks 3 and 4, their value streamed and as mixed units, and the count of unknowns against equations | `python3 shapes.py > out_shapes.txt` |
