# Which ratio groups can host a shared end

**What it shows.** 3 in 10 and 4 in 13 at block 4 must end on genuinely shared moves. Over an abelian group R of
ratio letters, the rational group algebra splits into fields, and a rationed shared end for W arrays needs at least W of
them. `classify.py` counts the classes of minimal and rationed (W + 1)-piece ends for the small groups (signs, order 3
with signs, order 4, ...) and checks one representative of each class exactly. `skeptic/my_classify.py` is the
skeptic's independent recount; `skeptic/coord_classes.py` counts the classes of coordinate end moves and refutes the
claim that they form one class (signs give 4 classes, order 3 gives 2).

**Commands.** `python3 classify.py` (standard library, 8 seconds); in `skeptic/`: `python3 my_classify.py` (2 seconds)
and `python3 coord_classes.py` (numpy, under a second).

**Expected output.** `out_classify.txt` (and `out_classify.json`), `skeptic/log_my_classify.txt`,
`skeptic/log_coord_classes.txt`.

**Re-run for this revision:** all three, on 19 October 2026; identical outputs.
