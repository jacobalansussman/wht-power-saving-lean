# Data

## `native/`: the community's two words in canonical form

Each `CANONICAL-WORD.json.gz` is the output of #3945's own program `native/run_native.py`, run unchanged by the
agents on a Linux machine from the pull request's word files; all five of its stages passed for each word. The
content hash (sha256 of the uncompressed JSON) is the one the community's receipts record.

| folder | word | file sha256 | content sha256 |
|---|---|---|---|
| `base` | #3945, provider gauge44 (`inputs/NEW-base.json.gz`) | `2868441d...` | `938f0924e78233bd...` |
| `extra` | #3945, provider extra (`inputs/NEW-extra.json.gz`) | `37f6dfa5...` | `6df1ec686b7ea392...` |
| `old-base` | #3815 (the words of #3734), gauge44 | `c6dafae9...` | `0db11248b311e5ae...` |
| `old-extra` | #3815, extra | `81fd6cf1...` | `a82077d5c0245d67...` |

For the two old words the canonical content equals, as JSON, the community's stored `owner/inputs/old-NATIVE-*`
(in `third-party/`).

## `gcert/`: my conversion of those words

`tools/lever/lower_b.py` converts each canonical word into this project's certificate format, extended for the
community's borrowed helpers (the format is called *gcert/1b*). The files without `-copies` are the ones the
Lean data is generated from: no source copies (the sources are written to and restored). The `-copies` files are
the same words before the copies are folded back, kept for the checks.

What the format holds, in short:

- 112 sources, 112 targets and 408 helper slots; frames are subspaces of GF(2)^9 given by reduced bases.
- Each slot's start and end frame; 134 slots start above the empty frame and 91 end at a hyperplane (the
  helpers the community calls zero buffers).
- The gates of the two phases of a run: additions at a frame, frame climbs, and **warps**, which move a register
  from one frame to another for free by an invertible linear map G of the 9 label bits (the community's
  borrowed-helper device; 72 maps per word, listed with their columns).
- The 12 targets that are lifted to a one-dimensional frame in the first phase, and the community's
  precondition gates, kept in the community's order.
- Each slot's *residual*, the width of the climb the undo step will need; it equals, helper by helper, the
  community's width.

What the conversion drops, as for the earlier words: the community's compensation, its 36 twisted initial
corrections and its full cleanup (with the 91 zero buffers). The engine derives the dirt and the undo step from
the clean product, and the checks in `results/g0/lower-*.json` show that this is exact: prepending the
community's own dirt phase gives exactly target + source on every target.

| file | sha256 |
|---|---|
| `S3945-base.gcert1b.json.gz` | `c453b298f77a1a8c...` |
| `S3945-extra.gcert1b.json.gz` | `1fae6030e4873b5f...` |
| `S3815-base.gcert1b.json.gz` | `09c8433ce591138...` |
| `S3815-extra.gcert1b.json.gz` | `d37cda4e3426ad52...` |

`residuals-3945.json` and `residuals-3815.json` give, per provider, the 408 residual widths, the 91 slots that end
at a hyperplane, and the original helper number of each slot.

## `t74ctl/`: the regression data

The community's earlier T74 words and their owner fill, repackaged by `tools/lever/t74_control.py` in the format
of #3945's package. The head test of `tools/lever/g0_staged.py` must find on it the 11,010 heads (6,150 + 4,860)
that the first head test found on 14 October 2026; it does (`results/g0/g0-T74-control.json`).
