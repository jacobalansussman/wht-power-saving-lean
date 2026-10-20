# Data and certificates

Each folder here belongs to the program folder of the same name in `../programs/`, and holds that folder's recorded
outputs, certificates and input lists exactly as the rounds produced them (private paths removed). `sh
../programs/workdir.sh` merges each pair into one working folder.

The compressed files hold the outputs of the long runs on remote machines, one folder per run, as the programs that
read them expect (`workdir.sh` unpacks each `job-outputs.tar.gz` into `jobs/` and `systems.tar.gz` into `gbin/`):

- `natural-form-star-tool/job-outputs.tar.gz`: every per-class verdict of the block-5 scans (10, 9 and 8 in 31), the
  deep passes, the repeated-target sample, the block-4 rank cover and the certificate checks (37 MB unpacked).
- `eight-in-31-rank-bounds/job-outputs.tar.gz` and its `independent-check/` twin: the stage-by-stage verdicts of
  round 3 and the skeptic's cascade.
- `three-in-14-groebner/systems.tar.gz`: the 140 local systems (Singular over Q and over a prime, msolve over two
  primes, and JSON) of the 3-in-14 sets and their controls; `job-outputs.tar.gz`: the Gröbner results and cofactor
  checks, and the four-array probe; `independent-check/job-outputs.tar.gz`: the skeptic's runs.

Key files: `natural-form-star-tool/all8.txt` (the first-pass verdict for each of the 50,895 classes of 8 in 31; 1,669 pass
there, and a deeper pass with every unneeded form as a test form excluded 403 of them),
`natural-form-star-tool/surv8_final.txt` (the 1,266 survivors of round 2), `eight-in-31-rank-bounds/verdicts.txt`,
`undecided.txt` (the 127), `hyperplane_classes.txt` (the 24), `eight-in-31-rank-bounds/independent-check/
extra_exclusions.txt` (the skeptic's 31), `three-in-14-groebner/certs.json` (zero-block certificates for the 44
sets), `natural-form-solver/independent-check/survivors_b5w4.txt` (the nine four-array sets at block 5 next to
plain), `one-move-block3-symmetry/independent-check/out_v4_gauge_table.txt` (the gauge-symmetry screens).

`SHA256SUMS` lists the sha256 of every file in this folder (`shasum -a 256 -c SHA256SUMS`).
