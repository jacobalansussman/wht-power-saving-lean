# Candidate (A): U8 word with exit-frame helpers, `GXD.U8X`, saving 1.2704023e-3

DFT of every length (OpenAI family 130 statement), exponent `1 - 12704023/10^10`. The kernel interface is at
`12704024/10^10`.

| | saving | status |
|---|---|---|
| previous Lean-checked record (Ryan Shea, Kothe's #588 circuit with 112 source copies) | 1.2331685e-3 | Lean-checked |
| **this package (A), certificate `GXD.U8X`** | **1.2704023e-3** (+3.02 %) | **Lean-checked on the Colab by k-lower, 13 Oct 2026.** Standard axioms, Compare.lean PASS, the one-unit-stronger challenge FAILS. The official comparator (nanoda and the Lean kernel) accepted it; the negative config was rejected |

The construction is the hub's U8 word (Jeremy Kothe #869/#870, Chafik Boukhalfa #890, Capalbo #942; canonical word
sha256 `24db880b...`, rebuilt by k-word with the hub's own `reproduce_barch.py`). It has 126 exit-frame helpers, and the
hub's 171 tail gates are lowered as the last gates of phase B. It keeps Shea's 112 source copies. The bank table is
pad-minimal: 4,721 banks, one pad of 41, W = 11889.

## What the package contains

- `overlay/`: the NEW files only. Paths are relative to Ryan Shea's `lean/complex588`
  (shea256/fourier-transform-below-nlogn, commit 5239484; `lean/` is byte-identical at 0e61af6).
  - `XExit/{Pack,Halves,Compat,Hub,End}.lean`: x-exit-lean's exit-frame theory (5 modules).
  - `Work/KLow/{ExitEnd,ZCheck,LabelsX,NoY,GlueX}.lean`: k-lower's engine.
  - 68 generated data, check and instance modules: `Work/GCert/Data/{U8X*,RateU8Xx,InstU8X}.lean`,
    `Work/GCert/Data/Gen/U8X*`, `Work/FourierU8X/*` and the challenge files `Work/FourierU8X{,Neg}/UniformFourierChallenge.lean`.
  - `comparator/fourier-u8x.json` and `comparator/fourier-u8x-neg.json`.
  These are k-lower's package files, byte for byte. `package_a.py` checks them against k-lower's `SHA256SUMS.overlay`.
- `MODULES.txt`: the 73 modules in build order (topological). It leaves out the challenges and `Axioms`, which
  `checks.sh` builds, and the XExit modules, which belong to the base tree.
- `SHA256SUMS.overlay`, `SHA256SUMS` (every file of the package); the tarball's own hash is in `k-combine-u8x.tar.gz.sha256`.
- `scripts/`: `setup.sh` (apply the overlay), `build-direct.sh` (one module at a time), `checks.sh` (k-lower's:
  axioms, sorry grep, both challenges, Compare.lean positive and negative) and `cmp.sh` (the official comparator, positive then negative).

## Build (Colab: 8 cores, 50 GB; at most 2 heavy Lean processes at once on the machine)

```sh
. /content/work/env.sh
mkdir -p /content/work/A && tar xzf k-combine-u8x.tar.gz -C /content/work/A && cd /content/work/A/k-combine-u8x
sha256sum -c --quiet SHA256SUMS && echo PACKAGE-OK
# base tree: your own COPY of a built Shea tree with the XExit library, e.g. colab-build's c-xexit
cp -a /content/work/c-xexit /content/work/A-tree
T=/content/work/A-tree/lean/complex588
bash scripts/setup.sh $PWD $T     # refuses to overwrite a different file; identical files are accepted
#   on a plain Shea tree instead: PATCH_LAKEFILE=1 bash scripts/setup.sh $PWD $T
#   then: (cd $T && lake build XExit.Pack XExit.Halves XExit.Compat XExit.Hub XExit.End)
tmux -L A new -d -s A "bash scripts/build-direct.sh $PWD $T /content/work/logs/A > /content/work/logs/A.out 2>&1; \
  bash scripts/checks.sh U8X $T /content/work/logs/A > /content/work/logs/A-checks.out 2>&1"
```

k-lower's run took 23 minutes for the 73 modules (16:58–17:21 UTC). Every module returned rc 0, with no `sorryAx`.

## Checks at the exact exponent (in `checks.sh`)

```sh
cd $T
lake env lean Work/FourierU8X/Axioms.lean
#   transform_mainU8X, convolution_mainU8X, RAM.hillsU8X_program: [propext, Classical.choice, Quot.sound]
lake env lean --run tools/Compare.lean Work.FourierU8X.UniformFourierChallenge Work.FourierU8X.Main \
  OAI.PowerSaving.transform_mainU8X OAI.PowerSaving.convolution_mainU8X        # RESULT: PASS (1 - 12704023/10^10)
lake env lean --run tools/Compare.lean Work.FourierU8XNeg.UniformFourierChallenge Work.FourierU8X.Main \
  OAI.PowerSaving.transform_mainU8X OAI.PowerSaving.convolution_mainU8X        # NEGATIVE (12704024): must FAIL
```

## Official comparator

These are c-base's tools: comparator d03acab, lean4export 076e8e5 and nanoda 3a24072. `jail2.sh` stands in for
landrun and systemd-run, which this container lacks.

```sh
tmux -L Acmp new -d -s c "bash scripts/cmp.sh U8X $T /content/work/logs/A-cmp > /content/work/logs/A-cmp.out 2>&1"
#   positive: 'nanoda kernel accepts the solution', 'Lean default kernel accepts the solution', 'Your solution is okay!', rc 0
#   negative (comparator/fourier-u8x-neg.json, 12704024): rc != 0, "Const does not match ... decimalExponentU8X"
```

Shea's direct form (RUNLOG.md) is
`COMPARATOR_LANDRUN=... COMPARATOR_LEAN4EXPORT=... COMPARATOR_NANODA=... lake env comparator comparator/fourier-u8x.json`.

## Verdicts already obtained (k-lower, 13 Oct 2026, Colab)

- Positive comparator: "nanoda kernel accepts the solution", "Lean default kernel accepts the solution",
  "Your solution is okay!", rc 0. Wall time 4,609 s; Lake rebuilt the 73 modules inside the jail, one at a time.
  Maximum resident memory 9.4 GB.
- Negative comparator: rc 1, `uncaught exception: Const does not match between challenge and target
  'OAI.PowerSaving.decimalExponentU8X'`, at 19:35 UTC.
- Logs are under `~/misc/fft/verify/checks/wht43/colab-run/logs/k-lower-u8x/` (on the Colab: `/content/work/logs/k-lower-u8x/`).

## Credits

- Jeremy Kothe (#869/#870, #588, #799), Chafik Boukhalfa (#890, #741) and Capalbo (#942): the U8 word and its
  lineage, with the hub authors (exit frames and tail #591/#593/#599) and FanPUYen/Hao Yan (owners).
- Ryan Shea: the Lean instantiation, the generators and the bank/seam pipeline.
- Nil Mamano: GInvE, Gauge and the record files.
- Dugongue (#401) and jacklightChen (#377): the bank geometry.
- OpenAI: the model and the Fourier reduction.
- Jacob Sussman: the engine and the checkers.
- Agents: x-exit-lean (XExit), k-word (the reproduced word) and k-lower (the lowering, KLow and the build).
