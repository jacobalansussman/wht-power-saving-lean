# Programs: assembling the head test for #2979, and the warp scan

The community's circuit #2979 ships its words and its allocation but no native word, and the head test (the lever
package's `heads_fill.py`, [`tools/lever/`](../../../../tools/lever/)) reads one field of a native word: the list of helper origins. These
programs build that field from the words, fetch the community's files read-only, and compare the borrowed helpers
("warps") of the community's circuits. Python 3.9 or newer, standard library only, plus the GitHub command line tool
`gh` for the fetch. Every community file named below is Apache-2.0; none is copied into this repository, and the
programs fetch them by commit and sha256.

| program | what it does |
|---|---|
| `fetch.sh <owner/repo> <commit> <path> <folder>` | read-only download of one file at a pinned commit; prints its sha256 and size |
| `mk_native_stub.py <WORD.json.gz> <out.json.gz>` | writes a stub native word holding only `helper_origin` = 224..703 minus the word's 72 removed helpers. On #3945, where the community ships the real native word, the stub equals the real field exactly; on #2979 my own transport of the words later produced the real field, and it equals the stub |
| `warpscan.py <word files>` | for each word file: its kind, the number of warps, initial corrections, zero buffers and warp edges |
| `warpcmp.py <word files>` | for each borrowed word: a digest of its (lender, map) pairs, of its whole warp set without positions, of its initial corrections and of its removed helpers |

## The head test on #2979

1. Fetch from the community's repository at main commit bb4bccbf42844d2801487450727c2f63e10d1d0f, folder
   `research/rank7-joint-endframes/`: `allocation/FILL.json.gz` (sha256 c2136be9434cda511a1a8a44d98f9f72ea5376830a8426674253b4ee351c4c39), `allocation/ORDER.json`
   (519d1d5dbb997a219fcfcad34c0b03e7838bef24199784d7e964fa10090f66e1), `allocation/DISPATCH.json` (65e435bd68e97243c003677271b88485a29d4d0fae2dff1354100404a145e74c), `data/WORD-base.json.gz` (abaf67bfcb257528658c906b612f57025617592bd7889828f92b97a5f4f4300d) and
   `data/WORD-extra.json.gz` (e0eaf64c4d5d324693d43889bce1b9e53e8da00402105c8c978fd102f208bed7). Fetch `PLAN.json.gz` (d566cfb87b71639d4138e1f89b543e2dd6a7772070226a2a14c799a6835bf21e) and `TEMPLATE.json.gz`
   (21297b20f1f2dc76eb9dd8124416068c3cb096bab9df17fa4da19bc68dfe082e) from #3945's `research/staged544-terminal-gap-fill/owner/inputs/` at its head
   91f81d131252db704f945891da1ea54e29d46b8b (repository Dugongue/integer-mult-bounds); #2979's lineage pins the same
   two files.
2. Lay them out as the head test expects: `pkg2979/FILL.json.gz` and `pkg2979/inputs/{ORDER.json, DISPATCH.json,
   PLAN.json.gz, TEMPLATE.json.gz, WORD-base.json.gz, WORD-extra.json.gz}`, and write the stubs:

       python3 -I -B mk_native_stub.py pkg2979/inputs/WORD-base.json.gz  pkg2979/inputs/old-NATIVE-base.json.gz
       python3 -I -B mk_native_stub.py pkg2979/inputs/WORD-extra.json.gz pkg2979/inputs/old-NATIVE-extra.json.gz

3. Run the head test, its control and the recount without chronology (the programs of the lever package):

       python3 -I -B <lever tools>/heads_fill.py pkg2979                          > heads-2979.json
       python3 -I -B <lever tools>/heads_fill.py pkg2979 --control-bridge-first   > heads-2979-control.json
       python3 -I -B <lever tools>/recount.py pkg2979/FILL.json.gz 74             > recount-2979.txt

Expected (re-run on 19 October 2026 for this revision, about 1 second and 0.27 GB; identical to the first run apart
from the folder names): 16,086 literal heads of rank 48 (child 8: 8,064; child 10: 8,022), 12,104 neither (children
12 and 13), the idle profile {1: 17,310, 3: 4,956, 5: 1,728, 48: 28,190} and the full paid profile equal to the
community's `allocation/PROFILE.json`; the control gives 0 heads; the recount gives 8,064 / 8,022 / 6,512 / 5,592.
The outputs are `../../data/heads/2979-head-test.json`, `2979-head-test-control.json` and `2979-head-recount.txt`.
The independent count in my own lowering and Lean network (the same 16,086) is
`../../data/heads/2979-heads-in-our-lowering.json`; its programs are the lever package's `g0_staged.py` and
`lower_b.py`.

## The warp scan

    python3 -I -B warpcmp.py <the base words of #2979, #3945, #3815, #3682, #3642, #3148, #3605> <their extra words>

Re-run on 19 October 2026 for this revision (output `../../data/heads/warp-comparison.txt`): every one of these
words has 72 warps, all non-identity, from 72 lenders, 36 initial corrections and the same 72 removed helpers; the
(lender, map) digest is the same on every base word and the same on every extra word. Only the warps' positions and
frame labels move from word to word (#3653 uses #3642's words). So each of these circuits needs exactly the warp work
#3945 needs, and no more of it. The best circuit without warps is T74, whose half cut (1.6968621e-3) is
below #3945's own figure (*checked by computer*, one program; the scan for circuits without warps in the band
1.6950e-3 to 1.7310e-3 was by keywords over the pull requests' titles and bodies, so it is a search, not a proof).
