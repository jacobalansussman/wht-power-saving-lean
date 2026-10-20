# preprint-reading: the two network codes of 7 October 2026, read against the transform game (10 October 2026)

These programs back section 10 of [../../README.md](../../README.md). They rebuild the graph and the local rule of
Braverman and He's code from the data in their preprint, measure how many free symbols the rule keeps over the field
with 9 elements and over the Gaussian integers, and test whether the code's graph can be laid on a frame graph. They
say nothing about whether either preprint is right; they compute properties of one reading of its data.

**The preprint's data are not copied here.** The programs need the two lists in the appendix of M. Braverman and
Z. He, "The Multiple Unicast Conjecture is False", arXiv:2610.10108, version 1: the 175 deleted edge indices and the
processing order of the 182 vertices. Make them from your own copy of the paper:

```
pdftotext 2610.10108v1.pdf bh.txt             # pdftotext comes with poppler; -layout works too
python3 extract_lists.py bh.txt bh_lists.json   # checks 175 distinct indices and a permutation of 0..181
```

Two short lists from the same appendix (the 9 edges fixed to zero and the 9 unused edges, 18 numbers) are written
into the programs, as the preprint gives them.

Then run, from this folder, with Python 3.9 or later (numpy only for the first program):

| program | what it does | command | time, memory | re-run |
|---|---|---|---|---|
| `extract_lists.py` | reads the two lists out of a text copy and checks them | as above | 1 s | tested on a text copy and on fresh `pdftotext` output (plain and `-layout`): the same lists |
| `rebuild_and_field.py` | rebuilds the 182-vertex graph (the points and lines of the projective plane over the field with 9 elements), the 735 remaining edges and 157 sessions, checks the counts the preprint states (5 in and 5 out at every vertex, every session 5 apart), and computes the dimension of the space of edge assignments obeying all 182 local rules: over the field with 9 elements, without the cubing at lines, and with the same integer rule read in Z[i]/p for p = 10007 and 1000003 (two digit conventions) | `python3 rebuild_and_field.py bh_lists.json` | 16 to 51 s, 75 MB | same output |
| `frame_graphs.py` | builds the frame graphs at blocks 2, 3 and 4 from scratch (15, 135, 2,295 frames), prints the numbers a placement must respect, and the saving a translated unit would have | `python3 frame_graphs.py` | 36 s | same output |
| `placement_test.py` | the placement test at block 5: pairs of vertices forced to keep their distance, and the 1,524 forced cases that need 8 or 9 different first steps where every frame graph has 7 | `python3 placement_test.py bh_lists.json` | 1 s | same output |
| `first_steps_and_drops.py` | the count of 7 first steps between frames 3 apart at blocks 3 and 4, and a greedy drop of the 9 sessions involved in most violations (1,342 violations remain) | `python3 first_steps_and_drops.py bh_lists.json` | 42 s | same output |
| `skeptic_collapsed_edges.py` | the skeptic's program: the forced pairs and the 7-against-8 count with separate loops; the case the first reading missed (placements that put both ends of some edges on one frame, which costs no move): block 4 closed (the collapsed edges form a matching, so at least 644 moves against 628), block 3 not decided; the exponent arithmetic | `python3 skeptic_collapsed_edges.py bh_lists.json` | 3 s | same output |

The outputs of the original runs are in [../../data/preprint-reading/](../../data/preprint-reading/), with the same
names. Re-runs for this revision (Apple Mac, Python 3.9.6, numpy 2.0.2) reproduced all five byte for byte.
