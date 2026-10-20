# Proofs and working reports of this section

These are the working files of the entropy-bounds round of 15 October 2026, as the AI agents wrote them, with private
file paths removed. Each file starts with a header that gives its status now, after the skeptic's check. Where the
header and the file disagree, the header (the later verdict) wins. No person has read these proofs, and nothing in
them is checked in Lean. "I" inside a file is the agent that wrote it, not me.

The other rounds of this section (the famous-problems survey, the quantum dictionary, the reading of the two
preprints) have no separate proof files: their arguments are short, and the section README gives each one with its
status.

| file | what it holds | status now |
|---|---|---|
| [entropy-bounds-proofs.md](entropy-bounds-proofs.md) | Lemma 1 (the uniform layout is the worst), the validity of every LP row, Lemma 2 (contraction), Lemma 3 (every LP that writes rows for only some node sets is capped by routing on the lumped network), Lemma 4 (four cells never beat cuts), the calibration, the ceiling curves, the cut LP of a fixed layout, Lemma W (two-class stretch), Lemma W+ (moving an end), the exact screen | Lemma W and W+: argued by hand, not checked, one skeptic re-derived them and checked their certificates exactly; Lemmas 1 to 4: argued by hand, not checked (the skeptic checked the symmetry facts and the values, not the lumping proofs); numbers: checked by computer |
| [entropy-bounds-report.md](entropy-bounds-report.md) | the round's report: the answer, the tables, the screens, what would change the answer | upheld by one skeptic with wording corrections (listed in its header) |
| [entropy-bounds-skeptic.md](entropy-bounds-skeptic.md) | the skeptic's record: its own frame graphs, LPs and Shannon LP; the 16 hard layouts solved; the four-array cut gaps at block 4 | a skeptic's check |

## Working names used inside the files

| name inside the files | here |
|---|---|
| `nc-gain`, `calc/` | the author of these files; its programs are in [../programs/entropy-bounds/lp-and-layouts/](../programs/entropy-bounds/lp-and-layouts/) |
| `v-nc-gain` | the skeptic; [entropy-bounds-skeptic.md](entropy-bounds-skeptic.md) and [../programs/entropy-bounds/skeptic/](../programs/entropy-bounds/skeptic/) |
| `job_<name>/` | a run on a Colab machine; its programs are in [../programs/entropy-bounds/lp-and-layouts/drivers/](../programs/entropy-bounds/lp-and-layouts/drivers/) and its outputs are the files `job_<name>__*` in [../data/entropy-bounds/colab-runs/](../data/entropy-bounds/colab-runs/) |
| `sl-combinatorics` (Lemma 4.1, Lemma 4.2, Remark 4.7, Theorem 4.3, Corollary 4.4) | the slack round's proof file, published as [../../frame-packing-two-arrays/proofs/ledger-for-every-unit.md](../../frame-packing-two-arrays/proofs/ledger-for-every-unit.md) |
| `sl-algebra` | [../../frame-packing-two-arrays/proofs/common-frame-accounting.md](../../frame-packing-two-arrays/proofs/common-frame-accounting.md) |
| `v-slack` (the cut constants 3/2, 7/4, 15/8) | [../../frame-packing-two-arrays/proofs/slack-skeptic.md](../../frame-packing-two-arrays/proofs/slack-skeptic.md) |
| `fp-capacity` (Theorem A(d), the coding-gain term) | [../../frame-packing-two-arrays/proofs/packing-capacity-and-rigidity.md](../../frame-packing-two-arrays/proofs/packing-capacity-and-rigidity.md) |
| `fp-build`, `prop6`, Proposition 6 | [../../frame-packing-two-arrays/](../../frame-packing-two-arrays/) |
| Theorem R, the certificates 2, 14/5, 41/13 | [../../block-bounds-and-new-rules/](../../block-bounds-and-new-rules/) |
| Theorem 2W | [../../moonshot-small-blocks/](../../moonshot-small-blocks/) |
| `../../proofs/duality/dfz_check.py` (the DFZ encoding) | the duality round, in [../../existence-routes/](../../existence-routes/) |
