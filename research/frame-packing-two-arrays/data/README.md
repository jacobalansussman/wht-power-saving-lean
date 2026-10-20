# Data of this section

The outputs of the original runs of 15 October 2026, under the names the programs gave them. Each folder matches a
folder of [../programs/](../programs/), whose README says which command wrote which file. Everything here is small
text or JSON (about 2.8 MB in all). The one exception is a 4 KB Python pickle, `frame-game/design-runs/hit_m4_s102_4.pkl`.

| folder | what is in it |
|---|---|
| `ledger/` | the ledger checks: named units, random words and units, the exhaustive word enumerations at blocks 2 and 3, the block-1 census, the cut constants and the cut LP |
| `common-frame-accounting/` | frame facts at blocks 1 to 3, what 0.1 and 0.01 need, the order-3 control check |
| `heights/` | the named units in the true model, the height filter's self-test, the `ledger.py` scalar warning, the token games and beam searches |
| `packing-measurement/out/` | the packing ledger of every call measured and of OpenAI's network; `summary.txt` is the full table |
| `packing-capacity/` | Theorem A on known units (`vflow_*`, `vflowB_*`), the geometry and rigidity tables (`geom.*`), the needs tables |
| `frame-game/` | calibration outputs; `design-runs/` holds the packing-first design runs and the numerical position-test searches |
| `frame-game-completion/` | the lifted test programs on raw arrays, the exact block-2 identities, the homomorphism checks, the correction sizes |
| `orbit-designs/` | `waveA/`, `waveB/`, `waveC/` (the orbit design runs, one folder per machine job), `growth-runs/`, and the tallies |
| `single-term-game/runN/` | the machine runs of the single-term game: run 10 toy games, runs 12 to 14 the cube and Hamming sub-games, run 15 and 17 invariants, run 16 the block-4 cube 2 in 7 (six parts), run 18 the Φ tests, run 19 the 2 in 5 decisions |
| `skeptic-checks/` | the skeptics' outputs, in folders matching `../programs/skeptic-checks/`; `single-term-game/block4-cube-2in7/` is the independent 107-million-state run |
