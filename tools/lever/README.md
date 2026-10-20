# The lever's search and pricing programs

Everything here is Python 3.9 or newer with the standard library only, and runs on a laptop. Run each one with
`python3 -I -B` from the root of the repository. `bash tools/reproduce-search.sh` runs them all and compares
their outputs with `results/` ([VERIFY.md](../../VERIFY.md), section 12.6). [PATCHES.txt](PATCHES.txt) lists the only changes made to
them when they were copied here: absolute import paths replaced by the program's own folder, and an option
`--xhost` added to `g0_staged.py`.

Paths below: `P = third-party/CrocSwap-integer-mult-bounds/staged544-terminal-gap-fill/owner` (the community's
owner data for #3945) and `W = $P/inputs`.

## 1. The head test on the community's own data

    python3 -I -B tools/lever/heads_fill.py $P                          > heads-3945.json
    python3 -I -B tools/lever/heads_fill.py $P --control-bridge-first   > heads-3945-control.json
    python3 -I -B tools/lever/recount.py $P/FILL.json.gz 360            > recount-3945.txt

`heads_fill.py` follows every array of #3945's unit through the community's schedule: the local sweeps, the
bridges between arrays, the idle sweeps, the hosted helper roles of the owner fill, and the terminal events. A
block is a *literal head* if nothing touches its array before it. It re-implements the window rules of the
community's `owner/verify.py`; it does not import it. Result: 57,312 literal heads of rank 48 (28,672 in child 8,
28,640 in child 10), 80,640 rank-48 blocks that are neither heads nor tails, and the community's idle and full
profiles reproduced exactly. The control moves the bridge that reads the head arrays before their idle sweep: the
heads drop to 0. `recount.py` counts the fully idle bins with no chronology at all and gives the same numbers.

## 2. The prices

    python3 -I -B tools/lever/controls.py      > controls.txt
    python3 -I -B tools/lever/lever_price.py   > lever-price.txt

`controls.py` reproduces the community's complex roots of its newest words to all 24 digits, with two
independent formulas: `pricer.py`, an exact re-implementation of the community's moment engine (certified by
rational enclosures on the grid of 10^-24), and `mixed.py`, the balance of the README ("The lever, with the numbers") solved with
60-digit decimals. `mixed.py` takes targets of several sizes; with the single target 73 it is the community's
formula.

`lever_price.py` first reproduces the figures of the plan of 14 October 2026 on T74 (full cut 1.6948358e-3, half
cut 1.6968621e-3, break-even 9,443 heads), then prices every variant on #3945 and the half cut on each of the
community's newest words (the README's table of lever prices).

## 3. My conversion of the community's words

    python3 -I -B tools/lever/lower_b.py data/native/base/CANONICAL-WORD.json.gz S3945-base.gcert1b.json.gz lower-base.json \
        --word $W/WORD-base.json.gz --with-copies S3945-base-copies.gcert1b.json.gz
    python3 -I -B tools/lever/mutants.py data/native/base/CANONICAL-WORD.json.gz $W/WORD-base.json.gz mutants-base.json

(and the same with `extra`; for #3815, `data/native/old-*` and `$W/old-WORD-*`). `lower_b.py` converts a canonical
word into gcert/1b ([data/README.md](../../data/README.md)) and checks the result: frames nested at every climb and
carried through every warp; the block histogram equal to the community's native profile; the cost identity; the
residuals equal to the community's widths, helper by helper; the exit conditions and Shea's scalar facts in the free
algebra of warps; and the dirt check. `mutants.py` makes seven wrong versions of the word (a frame lowered, two
columns of a warp swapped, an initial correction dropped, a coefficient doubled, an exit helper widened, a cut
coefficient changed, a lifted addition moved) and checks that `lower_b.py` rejects each one. It does.

## 4. The head test in my lowering and my network

    python3 -I -B tools/lever/g0_staged.py . $P $P/FILL.json.gz data/gcert/S3945-base.gcert1b.json.gz \
        data/gcert/S3945-extra.gcert1b.json.gz $P/PROFILE.json g0-3945.json --label "#3945" \
        --gates <Shea's tree>/lean/complex588/Work/Bridge401/Gates.lean --xhost package/u8x0hq/overlay/Work/XHost
    python3 -I -B tools/lever/price_ours.py g0-3945.json 0.001730663236534549309146 0.001730663236534549309147

`g0_staged.py` repeats the head test with my own data instead of the community's: the widths are my residuals,
the order of steps is that of the Lean network (read from the Lean files of `Work/XHost` and from Shea's
`Gates.lean`), and every run is cut into segments at its call times. It checks that every one of the 1,321,920
hosted roles lies inside a window of my network and that the global profile equals the community's digit for
digit, and it finds the same 57,312 heads. `price_ours.py` prices that profile with the community's engine (it
equals the community's stated root to 24 digits) and prices the lever with my head count.

`t74_control.py` made `data/t74ctl/` from the earlier T74 files (which are not in this repository);
`g0_staged.py` on that data must find the 11,010 heads of the first head test, and does.

## 5. The generators of the Lean data

They ship inside each lever package, in `package/<name>/generators/`, with `REGENERATE.sh`, which makes every generated Lean file of the package again from the community's files and compares it with the package's `SHA256SUMS` (see each package's `generators/README.txt`).
