# Toeplitz-tile FIR filter on the M4's matrix unit, fp32 (section 4)

The filter as outer products: the input in 16 polyphase rows, a Toeplitz table of the taps, 4 accumulator tiles of
the matrix unit, each vertical slice 16 contiguous outputs. Surviving result: **-59 % / -53 % / -50 %** at 16 / 64 /
256 taps against Accelerate's `sgemm` on block-Toeplitz panels, one thread, 2^20 samples, in cache; it loses on short
calls. Status: checked by computer, two critics.

| file | what it does |
|---|---|
| `fir.c` | the matrix-unit kernel and its baselines (vDSP `conv`, vDSP and PFFFT overlap-save); `check`, `time`, `olsscan` |
| `olsparts.c` | times the parts of the overlap-save baseline |
| `critic/alt.c` | the critic's stronger baselines (includes `../fir.c` unchanged): `sgemm` block-Toeplitz, multi-threaded vDSP and overlap-save |
| `results/`, `critic/results/` | the original outputs (`time_2e20_same.txt` is the headline table; `critic/results/alt_2e20_same.txt` the critic's) |

```sh
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$PFFFT/include fir.c $PFFFT/lib/libpffft.a -framework Accelerate -o fir
for L in 16 64 256 1024 4096; do ./fir check $L 20; done        # correctness: every line PASS
./fir time 64 20 same                                             # one row of the headline table
(cd critic && clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$PFFFT/include alt.c $PFFFT/lib/libpffft.a -framework Accelerate -o alt)
```

Needs: an Apple M4 (SME2), PFFFT. Time: the checks take seconds; a timing row about a minute. Memory under 1 GB.

Re-run for this revision: `fir check` at 16, 64 and 256 taps gave output identical to `results/check_2e20.txt`.
