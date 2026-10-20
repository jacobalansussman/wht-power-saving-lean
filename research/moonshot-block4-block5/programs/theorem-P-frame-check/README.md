# Theorem P: the skeptic's exhaustive frame check

**What it shows.** Theorem P says that private starts, private layer-2 moves, one shared layer-3 move and private tops
never make a unit, for any number of arrays, in the engine game with single letters. Its proof uses four facts about the
2,295 frames of the block-4 engine game and a lemma (D4). `myframes.py` is the skeptic's own code (it does not import
the track's): it checks the four facts and both cases of Lemma D4 over every frame, under weaker hypotheses than the
proof, and finds 0 survivors; dropping one layer condition gives 1,344 survivors, so the check has teeth.

**Command.** `python3 myframes.py` (standard library). 14 seconds on an Apple M4 Mac mini (16 GB).

**Expected output.** `out_myframes.txt`.

**Re-run for this revision:** yes, on 19 October 2026; identical output.
