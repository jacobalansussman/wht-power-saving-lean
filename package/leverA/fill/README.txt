fill/: design A's fill (the community's #2979 circuit with every head kept) and the commands that made it

FILL.json.gz   sha256 764395fcc2d883e9dc1ecfedce2131c745dc7ad1b9c41d7b52ef54da16d7f71c, 271,728 rows.
               It replaces only allocation/FILL.json.gz of #2979. The words, ORDER, DISPATCH, PROVIDERS, TEMPLATE, PLAN
               and PROFILE are #2979's files, byte for byte (generators/tree/wht45/allheads/aA-data/in/alloc2979H/
               holds the same FILL beside #2979's other four allocation files, and generators/REGENERATE.sh checks it).

What changed and why
  In #2979's allocation model the empty rank-48 bins can sit in children 8/10 (where a bin that stays empty is a
  literal head: the first operation on its own array) or in children 12/13 (never heads), at the same price; #2979's
  solver split them arbitrarily (16,086 heads in its published fill). I added one constraint per owner: every bin
  of the child-8/10 group is empty. The solver then finds a fill with the same paid profile, the same idle mix and
  the same price for #2979 itself (its b stays 1.7290340e-3), and 16,576 literal heads (8,288 in child 8 and 8,288
  in child 10), the most possible: 224 per copy.

How it was made (on a Colab; #2979's own code from the community repository's main branch at 20ecf2891, unchanged)
  0. environment: CPython 3.12, numpy 2.4.5, scipy 1.16.1, ortools 9.15.6755 (#2979's recorded versions)
       uv venv -p /usr/bin/python3.12 v2979
       uv pip install numpy==2.4.5 scipy==1.16.1 ortools==9.15.6755 mpmath==1.3.0 sympy==1.14.0 gmpy2==2.3.2
  1. construction/setup2979.sh: copies #2979's construction/allocation/*.py (fill17.py, double.py, ...) and
     tools/relax.py, root.py beside the scratch tree they expect (only symbolic links are added).
  2. construction/mkpatch.py code/fill17.py code/fill17_heads.py: fill17.py plus the option --force-heads (the
     constraint above) and --max-heads; nothing else in fill17.py changes.
  3. construction/run_probe2979.sh --force-heads force --seconds 2400, which runs
       fill17_heads.py --work runs/force-t37 --chronology pr2833-t37 --providers allocation/PROVIDERS.json
         --mix-file construction/allocation/T37-FACE-MIX.json --buffers 91 --seed 1 --workers 1 --catalogue tight
         --force-heads --seconds 2400
       double.py runs/force-t37 runs/force-t74 --providers allocation/PROVIDERS.json --buffers 91
       groupstats.py runs/force-t74/FILL.json.gz ORDER.json DISPATCH.json PROVIDERS.json TEMPLATE.json.gz
     CP-SAT (one worker, seed 1): OPTIMAL in 298.8 s; the owner check inside fill17 passes on the T37 fill; double.py's
     owner engine passes on the T74 fill (271,728 roles, idle {1: 17,310, 3: 4,956, 5: 1,728, 48: 28,190}); total
     heads 16,576; FILL sha256 764395fc... (logs/force-heads-run.out). The T37 model, solver record, patterns, owner
     check and flow are in construction/T37-*.

Checks of the fill (checks/)
  - #2979's registered check_allocation.py on a package copy with only allocation/FILL.json.gz and its pin replaced:
    PASS_RANK7_LITERAL_OWNERS_AND_EXACT_MOMENT, five rejection controls rejected; the unchanged package passes too
    (checkalloc-o-model.out; again, independently, ah-gate/regA/).
  - #2979's complete admission replay (reproduce.py) on a fresh copy where only the FILL differs, with fresh
    natives: PASS_RANK7_JOINT_ENDFRAMES_COMPLETE_REPLAY, 139 of 139 strict controls rejected, kappa = the claim,
    PASS_RELEASE (ah-gate/rpA/, ah-gate/GATE.md).
  - heads counted five ways (16,576 each), the paid profile equal to #2979's PROFILE.json digit for digit, and the
    half-cut price 1.7564708e-3 by two routes (eval-o-eval.json, heads_fill.json, GROUPSTATS.json, ah-gate/GATE.md).
