# Block 3 with everything free: the exact collapse (cross-talk) test, and its two re-checks

When two arrays share a paid move, each output receives a copy of the other array's work, which must cancel. These
programs test, by exact linear algebra on stacked digit-column sizes, whether any mix of free words around paid moves
equals something with fewer paid moves (a *collapse*). Result: only trivial collapses (checked by computer; one
skeptic re-ran it against the complete lower span and classified every flag; a second skeptic confirmed the parts it
finished).

| File | What it shows |
|---|---|
| `xw3lib.py` | The everything-free game at block 3 (Toffoli, controlled-S, CCZ, chain maps and running XORs on every row, mirror chains, parity phases, boundary letters) with paid moves D_z, as exact operators on arrays of 8^f entries. |
| `sanity3.py` | Letter identities at f = 1 to 4 (the parity-phase identity, Toffoli and CS conjugation rules). |
| `orders.py` | Orders, as permutations, of free words such as A·B (running XOR then its mirror) at f = 2 to 6: 3, 7, 7, 31, 63. |
| `annih4.py`, `annih6.py` | The author's collapse tests: 3-move words with mixes at both internal positions (shared first, middle and last moves), and two shared moves in a row (the 3-in-8 shape). |
| `ctl_annih.py` | The positive control: a planted pair of columns whose sum lies in the lower span must be reported. |
| `krigid.py` | The zero-divisor test for readouts: no mix ρ of free words puts ρ·K inside the span of words with at most two moves. |
| `als3.py`, `fam.py` | The least-squares searches (block-coordinate exact least squares over program shapes, with planted controls). |
| `xtalk_ext.py`, `path_ext.py` | The skeptic's versions of the test against the complete span of every term with fewer moves. |
| `fu.py` | The follow-up: extract every flagged collapse at f = 3, 4, re-test the whole collapse subspace at f = 5 and 6 on fresh inputs, and separate words that are reducible (a direction applied twice). Ran on Colab machines, up to about 50 minutes per configuration. |
| `frame_check.py`, `identity_check.py` | That the frame of one flagged configuration commutes with the move (so it duplicates another), and the explicit trivial identities behind the surviving flags, at f = 3 to 6. |
| `classify.py`, `terms.py` | Read the saved follow-up outputs and print the classification table (coincidence of small f, or trivial in disguise) and the words that carry each flag. |
| `second-skeptic/` | The second skeptic's own letters (`vlet.py`), its independent re-test of single flags (`vchk.py`) and of the identities at f = 3 to 7 (`vident.py`). |

Commands (run from this folder; numpy needed). The menus are lists of free words; the ones used in the runs were

```
R9=1,X1,X2,X4,Z1,Z2,Z4,P01,P12
K5="1,T2,A,A.T2.Ai,Pr0"
K9="1,T2,CS01,A,Ai,A.T2.Ai,Pr0,P01,X1"
H15="1,P01,P12,C01,S1,X4,Z1,A,Ai,B,lX1,fS1,Pr0,T2,CS01"
H30="1,P01,P12,C01,C12,C20,S1,S2,S4,X4,Z1,CZ01,A,Ai,a0,a0i,E01,B,Bi,lS1,lX1,fS1,Pr0,Pr2,T2,T0,CS01,CS12,A.T2.Ai,Ai.CS01.A"

python3 sanity3.py ; python3 orders.py ; python3 ctl_annih.py                 # seconds each
python3 annih4.py "1;2;4" "$R9" "$H30" "$K9" "$H30" 3,4 1                    # about 2 minutes, about 1.2 GB; one 3-move configuration
python3 annih6.py "1;2" 4 7 "$R9" "$H30" "$K5" "$H30" 3,4 1                  # about 35 s, about 1 GB; two shared moves (also: 4 4, 7 3, ...)
python3 path_ext.py "1;2;4;7" "$R9" "1" "$K5" "$H15" "$H15" 3,4 1 5         # one batch-B configuration, complete span
python3 fu.py B1 "1;2;4;7" "$R9" "1" "$K5" "$H15" "$H15" 3,4 1              # follow-up of flag B1; heavy (Colab)
python3 frame_check.py ; python3 identity_check.py                           # seconds to a minute
python3 classify.py ../../data/block3-collapse-test/followup                 # reads the saved outputs
FU_DIR=../../data/block3-collapse-test/followup python3 terms.py             # reads the saved fu_*.npz
cd second-skeptic && python3 vident.py                                       # seconds, about 1.7 GB
```

Expected output: `../../data/block3-collapse-test/` (the table of every original run is `out_crosstab.md`; the
follow-up's classification is `followup/out_classify.txt`, one row per flag, ending in class (a) or (b), never (c)).
The heavy runs (the full campaigns of `annih4` and `annih6`, and `path_ext`, `fu`, `vchk` and the least-squares
campaigns) were not re-run for this revision; their saved outputs are in `data/`. One configuration of `annih4`
("1;2;4") and two of `annih6` (top moves 4,4 and 4,7) were re-run; the current code also prints a count of
non-trivial collapses (1 and 4 there), which are the repeated-direction collapses the section README describes. Re-runs of the light programs: `rerun-2026-10-19/` there.

Comments inside the programs use the agents' track labels (for example "ms-xw3" for the author, "v-ms-xw3" and
"xw3-extract" for the skeptic and its follow-up) and refer to their working notes; the section's Sources maps the
labels.
