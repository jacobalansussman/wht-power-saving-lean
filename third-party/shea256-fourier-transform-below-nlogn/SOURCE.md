# Files from Ryan Shea's repository

Unmodified copies of files of github.com/shea256/fourier-transform-below-nlogn at commit
52394849c7ecf3412d39b925d75fc7fdb98994d8, at the same paths as in that repository. Licence: Apache-2.0
(`LICENSE`); notices: `NOTICE`, `lean/complex588/NOTICE`, `lean/complex588/NOTICE.complex588`.

* `research/complex-transfer-588-596/experiments/` and `verification/`: Shea's generator and verification
  scripts, which the generators in `tools/community/c588m-gen/` drive or adapt (each adapted file says what
  it changes).
* `lean/complex588/artifacts/`: his published certificate of Kothe's #588 circuit with 112 source copies
  (`candidate.gcert.json.gz`), its bank schedule and its price; the inputs of `tools/community/c588m-gen/`.
* `lean/complex588/third-party/record-templates/`: the record files of the Fourier chain (Nil Mamano's
  submission to bigospeedrun, as Shea's tree carries them); the templates of the generated Fourier files.

Each file was compared by sha256 with the clone at that commit on the machine where the packages of
`community/` were built.
