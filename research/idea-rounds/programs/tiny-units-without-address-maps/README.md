# Tiny units without address maps (10 to 11 October): selected checks

A selection from the night of 10 to 11 October: the programs behind the statements that section 1.3 of the README
reports, not the searches (their outputs were about 560 MB of search states and logs). Python with numpy unless said;
the two C programs need a C compiler. Saved outputs: `../../data/tiny-units-without-address-maps/<folder>/`.

`regmachine.py` is a hand transcription of the game the Lean engine certifies, written to test the agents' reading of
the Lean definitions (each definition names the Lean file and line of the third revision it copies). It is **not** Lean
and not checked by Lean. It has two simulators: on raw arrays for a given number of digit columns, and a formal one in
which every array holds a rational combination of "group element applied to input l". Several programs below import it.

**joint-test/** (block 3, two digit columns, the label-routing picture)

| Program | What it shows | Command |
|---|---|---|
| `model_cert.py` | builds the capacities and demands of the 15 label groups from 8 × 8 matrices, exactly; writes `model_cert.txt` (a copy is in the data folder) | `python3 model_cert.py model_cert.txt` |
| `hand.py` | the six fixed inequalities of the hand argument reject every multiset of 8 directions for 3 arrays | `python3 hand.py > out_hand.txt` (reads `model_cert.txt` from the current folder: run `model_cert.py` first, or copy it from the data folder; so do `static.py` and `staticlp.py`) |
| `static.py`, `staticlp.py` | static cuts for every multiset of K directions; the best combination of rows, 13/5 moves per array | `python3 static.py 3 8`; `python3 staticlp.py 9` |
| `check_cert.py` | checks the static certificate `static_cert_W3_k8.json` (in the data folder) from the label formulas alone | `python3 check_cert.py > out_check_cert.txt` (reads the copy in the data folder unless the file is in the current folder or `CERT_JSON` names it) |
| `sk_six.py` | the skeptic's own check that the six inequalities alone exclude 3 in 8, by enumeration | `python3 sk_six.py > out_six.txt` |
| `certfind.c`, `certcheck.c` | the cut certificates for all 5,764,801 ordered patterns of 8 moves, by dynamic programming over nested sets, and a checker that only adds numbers | `cc -O2 -o certfind certfind.c`; `cc -O2 -o certcheck certcheck.c`; then for D = 1 to 7: `./certfind model_cert.txt 3 8 all D D cert_W3_k8_pD.bin cert_W3_k8_pD.stats`, and `./certcheck model_cert.txt 3 8 all cert_W3_k8_p1.bin ... cert_W3_k8_p7.bin` (the seven files are 23 MB each and are not published) |

**face-certificate/** (block 3 at three digit columns; block 4 at four)

| Program | What it shows | Command |
|---|---|---|
| `facet.py` | the antisymmetric part of f digit columns and the span of the facet tuples: cap and demand per direction | `python3 facet.py M F` (saved: (m, f) = (2,2), (2,3), (3,2), (3,3), (3,4), (4,3), (4,5) to (4,8)) |
| `sk_general.py` | the checking agent's own recomputation, from `regmachine.py`'s matrices: the facet numbers, Theorem 2 of the general note, and the 28 bases | `python3 sk_general.py facet`, `... thm2`, `... bases` |
| `sk_cert.py` | the four-column certificate at block 4 (`cert_4col.txt`) recomputed with exact Gaussian integers: 3728/1219 = 3.058 moves per array | `python3 sk_cert.py > out_sk_cert.txt` (writes `my_cert_rows.json`) |

**symmetry-tests/**

| Program | What it shows | Command |
|---|---|---|
| `symtest.py` | exact tests of the symmetry maps on real schedules, replayed by `regmachine.py` (formal sum and f = 1, 2, 3): only the 6 coordinate permutations of the 168 address maps are symmetries of the sub-game; time reversal is one | `python3 symtest.py > symtest.out` |
| `gaugetest.py` | the gauge freedom of the straight-line form, tested the same way | `python3 gaugetest.py > gaugetest.out` |

Re-run for this revision: `model_cert.py`, `hand.py`, `static.py 3 8`, `staticlp.py 9`, `check_cert.py`, `sk_six.py`, `facet.py` for six (m, f) pairs, `sk_general.py` in its three modes, `sk_cert.py`, `symtest.py` and `gaugetest.py` re-run: same output (`facet.py 2 3` stops with a division by zero, because its subspace is empty, exactly as the saved run did); `certfind.c` and `certcheck.c` not re-run.
