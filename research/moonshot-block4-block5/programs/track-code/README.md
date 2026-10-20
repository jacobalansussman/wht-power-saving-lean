# Track code, as it ran

These folders hold the programs of the tracks whose findings are told in the section README but which have no curated
folder above. They are copied as they ran, with private paths replaced by relative ones; their outputs, logs and
machine job folders are not published. None of them was re-run for this revision. Most need numpy. A script that
imports `bm.py` needs the small C library next to that `bm.py`: `cc -O2 -shared -fPIC -o libmp_mac.so mp.c` on a Mac,
or `cc -O2 -shared -fPIC -o libmp_lin.so mp.c` on Linux. Each script's first lines say what it does and how it was
called.

| folder | what it is |
|---|---|
| `common/` | `bm.py` and `mp.c` (the exact game engine with closures and mod-p ranks), `blk9.py` and `nullsp.c` (hook blocks), shared by the folders below |
| `fano-blueprint/` | the Fano-gadget blueprint at block 4: Lemma P and S checks, the hook depth profile (`delta3.py`, `hookk.py`), the private-ending count, the template searches; `skeptic/` holds the order-3 counterexample's companions `frust.py` and `hk.py` |
| `symmetric-block4-skeletons/` | the symmetric block-4 searches that produced the swap-symmetric 2-in-7 near miss (`symsolve.py`, `symsearch.py`, `symlib.py`) |
| `symmetric-2-in-7-repairs/` | the near miss pushed: normalisation, the full test, the epsilon family, the junk absorber, the twisted gadgets, order-3 twists |
| `singer-template/` | the Singer-cycle template for 2 in 7, 3 in 10 and 4 in 13 (`t_full.py`, `tops.py`, `s4lib.py`); `skeptic/` holds the joint solver `jointcheck.py` and `pooltest.py` |
| `two-shared-layer-2-moves/` | Theorem P's track: the frame calculus `sp.py`, the letter searches, the exact sweeps and the second-round checks; `skeptic/` holds `robust.py` and the two jobs `gapcheck.py` and `pools.py` |
| `block5-early-separation/` | block 5, wave 2: Corollary B++, the top geometry and the censuses (`diag*/`); `skeptic/chk_bf7d.py` refuted "leftovers find no partner" |
| `block5-two-branch-trunk/` | block 5, wave 4: the two-branch trunk (`split.py`, the old `peel.py`, the exact probes); `skeptic/` holds the twin and other-branch checks |
| `block5-split-trunk-orbits/` | block 5, wave 7: Lemma E1-sigma, Lemma T, the collision censuses (`e1census.py`, `cover.py`, `lcensus.py`, `orbstat.py`), the orbit designs (`jdesign.py`) and `famcheck.py`; `skeptic/` holds the cycle-orbit runs |

The block-5 folders before wave 8 used `peel.py`, which missed forced zeros (see `../block5-sparse-designs/`). Their
"alive" layer-3 closures are not evidence of anything; wave 8 decided every rebuildable one again, exactly, with
`peel2.py` and exact tests.
