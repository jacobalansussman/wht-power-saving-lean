fill/made-by/: the programs and the record of the run that made ../FILL.json.gz

  r2/run2979.sh             the command (one run, about 1 minute on a cloud machine with 8 processors):
                            bash r2/run2979.sh u1-2979sep runs/u2979-sep/SOL-MIP0.json 0.0017884
                            with the community repository's #3815 construction (research/staged544-clock-ceiling)
                            and #2979's package (research/rank7-joint-endframes) as the hub inputs
  tools/hall_fc.py          #3815's per-tick Hall LP (construction/hall_lp.py, loaded from a copy and re-parameterised
                            in memory; the hub file is not changed) with the objective of the full cut: every head
                            removed, every other paid block above 24 split
  tools/build_fc.py, fc_cost.py, offsets3945.py   round 1's builder (LP, MIP on the generated columns, placement) and
                            its helpers; r2/tools/build_fc2.py is its copy with one more placer
  r2/tools/fallback_patch.py  the tolerant placement: #3815's reverse-time placement (construction/build_fill.py),
                            patched in memory; a role with no planned slot left at its call goes to a bin open at that
                            call that is not a head bin, not yet used at that call, with room in its final capacity
  r2/tools/fillmodel.py     the fill model the builder reads
  SOL-MIP0.json             the class solution (the MIP point) the placement started from
  u1-2979sep.fallback.json  the spare places used: 8 roles of width 6 at call 553 and 1 of width 4 at call 0; 8 planned
                            width-6 slots left idle
  BUILD.json                the run's record: the PR2839 owner engine's verdict (PASS) with its 7 rejection controls,
                            the idle remainders and the paid profile
  logs/u1.out, eval-u1.out  the run's log and the evaluation's log
