# The floor theorem of 10 October: exact tests

Programs of the floor round of 10 October. `f-engine`, `f-model`, `f-break` and `f-note` are the three authors' and
the note writer's tests; `g-engine`, `g-model` and `g-break` are the three skeptics', each written with its own code.
Pure Python (`f-engine/floor.py` also uses numpy). Ranks are exact over Q(i) or modulo primes p = 1 mod 4, with i a
square root of −1 (a rank modulo a prime can only undercount). Saved outputs: `../../data/floor-theorem/<folder>/`.

| Folder, program | What it shows | Command used for the saved output |
|---|---|---|
| `f-engine/floor.py` | matrices typed from the Lean definitions: every free letter keeps the diagonal; a block of rank k raises it by 2^m − 2^(m−k) (m = 2, 3, 4); parts B and C test the counting on more columns and the engine corollary | `python3 floor.py A > out_A.txt`, `... B > out_B.txt`, `... C > out_C.txt` |
| `f-model/checks.py`, `hatch2.py` | the model's checks; one address map that mixes two columns serves the kernel at m = 1, f = 2 with no block (an escape hatch) | `python3 checks.py > checks.out`; `python3 hatch2.py > hatch2.out` |
| `f-break/diag.py`, `walk.py` | the transfer ranks behind the floor (two primes), and random exact walks testing the transfer lemma at f = 2 | `python3 diag.py > diag.out`; `python3 walk.py > walk.out` |
| `f-break/span.py`, `span_clo.py`, `span_big.py` | the finite-f test: is the transform in the span of what zero-block programs can do? | `python3 span.py 1 2 65537` (and f = 3, 4, 5; `2 2 65537 bim`) → `span.out`; `python3 span_clo.py 2 2 65537 > span_clo.out`; `python3 span_clo.py 1 6 65537 > m1f6.out`, `... 2 3 65537 > m2f3.out`, `... 3 2 65537 > m3f2.out`; `python3 span_big.py 2 2 65537 > span_big.out` |
| `f-note/eig.py`, `diagrem.py` | the eigenvector support bound at lengths 4 and 16; the first family of the half-touch experiment (all 2^16 sign patterns at length 16) | `python3 eig.py > eig.out`; `python3 diagrem.py > diagrem.out` |
| `g-engine/check.py`, `even.py` | the skeptic's own field arithmetic: the potential lemma on random words, the block weights, the even-weight certificate at block 2 | `python3 check.py 1 3 4 > out_134.txt`; `python3 check.py 5 2 > out_25.txt`; `python3 check.py 2b > out_2b.txt`; `python3 even.py > out_even.txt` |
| `g-model/gcheck.py`, `gcheck2.py`, `gcheck3.py` | the model re-derived: weights of arbitrary column-uniform paid operators; column permutations free at m = 3, f = 3 | `python3 gcheck.py > gcheck.out`; `python3 gcheck2.py > gcheck2.out`; `python3 gcheck3.py > gcheck3.out`; `python3 gcheck3.py 469762049 > gcheck3b.out` |
| `g-break/g1.py`, `g2.py` | exact checks over Q(i) of the counterexamples to the refuted statements; spans of zero-block steps modulo a prime | `python3 g1.py > g1.out`; `python3 g2.py 65537 none 2500` (and modes `both`, `fan`, `fold`) → `g2.out` |

`g-break/rediag.out` and `rewalk.out` are the skeptic's re-runs of `f-break/diag.py` and `walk.py`.
Re-run for this revision: 22 runs re-run with the same output; `span_clo.py 3 2 65537` was stopped after 5 minutes before its first result, as the original run was; not re-run: `span.py`, `g2.py`.
