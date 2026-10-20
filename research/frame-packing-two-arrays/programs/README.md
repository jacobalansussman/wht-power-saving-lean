# Programs of this section

Each folder holds the programs behind one group of findings, a README with the commands, and its outputs in
`../data/<same folder name>/`. The programs import their neighbours by folder: `heights/` imports `../ledger/`,
`orbit-designs/` imports `../frame-game/` and `../heights/`, `skeptic-checks/single-term-game/vcross.py` and
`vsym.py` import `../../single-term-game/`, `skeptic-checks/frame-game/vselfinj.py` and `vpsilin.py` import
`../../frame-game/`, and `skeptic-checks/frame-game-completion/export_progs.py` imports
`../../frame-game-completion/`. Keep the folder layout as it is.

Requirements: Python 3.9 or later; numpy for most programs; scipy for `ledger/cutlp.py`,
`packing-capacity/vflow.py`, `single-term-game/t_phi.py`, `frame-game/t_bilin.py`, `t_psibil.py`, `t_hit.py`, `t_hit2.py` and
`skeptic-checks/packing-capacity/vhu.py`. Nothing else. Everything is exact (fractions, Gaussian rationals, exact group
elements) unless a README says a search is numerical.

| folder | what it shows |
|---|---|
| [ledger/](ledger/) | an exact engine at blocks 1 to 4 with the full ledger M = W·m − 2N + X + 2D + Z + B − P; the exhaustive identity checks; the cut constants |
| [common-frame-accounting/](common-frame-accounting/) | frame geometry from actual matrices; what a 0.1 or 0.01 unit needs |
| [heights/](heights/) | `heightfilter.py` (Theorem C as an exact filter for any word); the named units in the true model; the `ledger.py` scalar warning |
| [packing-measurement/](packing-measurement/) | the packing ledger of one call of each local word, and of OpenAI's network |
| [packing-capacity/](packing-capacity/) | Theorem A on actual units; what one move can serve; the rigidity tables |
| [frame-game/](frame-game/) | the frame game with its exact solver and the position test; the packing-first designs |
| [frame-game-completion/](frame-game-completion/) | the completion rule and `lift.py`, which turns a frame-game program into an engine word |
| [orbit-designs/](orbit-designs/) | orbit designs with the exact screens; the parity-law calibration |
| [single-term-game/](single-term-game/) | the exact single-term game solver; the invariant tests for Conjecture T |
| [skeptic-checks/](skeptic-checks/) | the skeptics' own programs, written without importing the tracks' code |

## Re-run for this revision

On 19 October 2026 I re-ran the programs below from a clean copy of these folders, on an Apple M4 Mac mini. The
comparison is line by line against the shipped output, ignoring blank lines. "Same" means equal apart from timing
lines, the paths printed in warnings, and where stderr lines fall.

| folder | command | time | result |
|---|---|---|---|
| ledger | `python3 t_known.py` | 29 s | same |
| ledger | `python3 t_pauliblind.py` | under 1 s | same |
| ledger | `python3 calc_needs.py` | under 1 s | same |
| ledger | `python3 t_strict.py` | 30 s | same, plus a final summary line (34 cases, no identity violations, no non-units) that the shipped output lacks |
| ledger | `python3 cutconst.py 2` | under 1 s | same (c_2 = 1.5 over all 16,383 frame sets) |
| ledger | `python3 t_exhaust.py 2 2 3` | 2 s | 28,850 runs, 0 identity violations, as in the run record |
| ledger | `python3 t_exhaust.py 2 2 4` | 345 s | same: 692,402 runs, 0 identity violations, 24 units, all the plain words |
| common-frame-accounting | `python3 frames.py` | under 1 s | same |
| common-frame-accounting | `python3 needs.py` | under 1 s | same |
| heights | `python3 checks.py` | under 1 s | same |
| heights | `python3 heightfilter.py` | under 1 s | same |
| heights | `python3 t_scalar.py` | under 1 s | same (the word `ledger.py` accepts fails in the true model) |
| heights | `python3 t_chain.py 1`, `2`, `3` | 2 s each | same (the shipped file holds the three seeds) |
| packing-measurement | `python3 calib.py` | under 1 s | same |
| packing-measurement | `python3 oa_ledger.py 7` | under 1 s | same as the h = 7 block of the shipped `oa_ledger_7_8_9.out` |
| packing-measurement | `python3 fm_ledger.py <tools/certificate/gcert1-e8-r783.json.gz> E8r783 init=F1,F2,F3 none ext-each gates gates-one` | 3 s | same: every call type and helper model of the E8 word |
| packing-capacity | `python3 tables.py` | under 1 s | same |
| packing-capacity | `python3 geom.py 3 300` | 10 s, 0.45 GB | same |
| frame-game | `python3 t_cal.py` | under 1 s | same |
| frame-game | `python3 t_o3ctl.py` | under 1 s | same |
| frame-game-completion | `python3 t_units.py 2 8 6` | under 1 s | same |
| frame-game-completion | `python3 t_formal.py` | 17 s | same |
| frame-game-completion | `python3 t_hom.py 2 8` | under 1 s | same |
| frame-game-completion | `python3 t_dsize.py` | 4 s | same |
| orbit-designs | `python3 cal.py` | under 1 s | same |
| orbit-designs | `python3 t_lemma6.py` | under 1 s | same |
| single-term-game | `python3 t_stg_toys.py D5 0 1500000` | under 1 s | same as the D5 lines of `../data/single-term-game/run10/outA_final.txt` |
| skeptic-checks/slack | `python3 vframes.py` | under 1 s | same |
| skeptic-checks/slack | `python3 vstandin.py` | under 1 s | same |
| skeptic-checks/heights | `python3 vthmA.py` | 2 s | same |
| skeptic-checks/packing-capacity | `python3 vsharp.py` | under 1 s | same |
| skeptic-checks/packing-measurement | `python3 vcal.py` | 4 s | same |
| skeptic-checks/packing-measurement | `python3 voa.py` | under 1 s | same |
| skeptic-checks/frame-game-completion | `python3 vexh2.py` | 39 s | same: all 10,800 products at block 2 exact |
| skeptic-checks/closed-patterns | `python3 v_lemma1.py` | 18 s | same: the parity law on 40 cases |
| skeptic-checks/frame-game | `python3 vtally.py` | under 1 s | same: every design tally recounted from the raw outputs |
| skeptic-checks/frame-game | `python3 vhit.py` | under 1 s | same, except that the shipped file repeats one line |
| skeptic-checks/frame-game | `python3 vselfinj.py 4 30 12 61` | 17 s | same |
| skeptic-checks/frame-game | `python3 vpsilin.py 4 30 12 61` | 16 s | same |
| skeptic-checks/single-term-game | `python3 vb4.py D7all:0 5` | under 1 s | same as its line in `prop45_local.txt` |
| skeptic-checks/single-term-game | `python3 vinv.py D5 4` | under 1 s | Conjecture S holds on all 487 states (no shipped file for this case) |

Not re-run for this revision, because each takes from many minutes to hours or needs a machine with more memory:
the block-3 exhaustive ledger check (`t_exhaust.py 3 2 4`, 23 minutes in the original run), the block-1 census
(`t_census.py 1 2 2 3`, 13 minutes), the design waves of `frame-game/` and `orbit-designs/`, the block-3 and block-4
cube and Hamming runs of `single-term-game/` (up to 111 million states), the skeptic's block-4 cube run
(`vb4.py B4 7 4 k 6` for k = 0 to 5, about 10 minutes each on a Colab machine), and the large raw-array checks
(`skeptic-checks/heights/vraw.py`, `skeptic-checks/packing-capacity/vraw_big.py`). Their outputs are shipped as
they came out.
