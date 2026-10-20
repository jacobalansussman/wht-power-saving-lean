# One-pass double-precision Welch spectrum on the M4's matrix unit (section 3c)

Window, FFT (the batched kernel of [../sme-fp64-batched-fft](../sme-fp64-batched-fft/), with the window folded into
its first pass) and the |X|² accumulation in one pass over the signal. Surviving result, one P-core, 2^22 samples,
50 % Hann, the whole call: **-33 % / -27 % / -11 % / -41 %** at segment lengths 512 / 1,024 / 2,048 / 4,096 against
FFTW double with NEON and the real-input routes; on four P-cores no size clears (the claimed four-core gains are refuted). Status: checked by computer, one
critic and the fairness critic.

| file | what it does |
|---|---|
| `as0.c` | the Welch pipeline and its baselines (vDSP, FFTW double, PFFFT double, KFR) in one C file |
| `smegi.h`, `smegp.h`, `ap0lib.c`, `ap0lib2.c`, `ap0p.c` | this folder's copy of the kernel (with the window additions in `smegi.h`) |
| `wbench.c`, `cluster.c` | one-core and four-core timing (`wbench N NOVERLAP LOG2L C`) |
| `check.py`, `check_big.py` | 239 checks against `scipy.signal.welch`, every bin; `AS0_BROKEN=1` must fail all |
| `as0py.py`, `welch_sme.py`, `e2e.py` | the Python entry points and the end-to-end comparison with scipy and MNE |
| `critic/au1bench.c`, `critic/au1cluster.c` | the critic's harness with a correctly built FFTW double and P-core plans |
| `critic-fairness/welch_alt.py` | the 6-line batched `rfft` Welch in Python that took most of the Python headline |
| `results/`, `critic/results/` | original outputs (`t12.txt` one core, `t13.txt` four cores; `critic/results/r_t12.txt` the critic's) |

```sh
# kfrd.o from ../baseline-wrappers (KFR), then:
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$PFFFT/include -I$FFTW_DOUBLE_NEON/include -shared -fPIC as0.c \
  kfrd.o <KFR objects> $FFTW_DOUBLE_NEON/lib/libfftw3.a $PFFFT/lib/libpffft.a -framework Accelerate -lc++ -o libas0.dylib
python3 check.py                       # 239 PASS
AS0_BROKEN=1 python3 check.py          # 239 FAIL, as it must
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$PFFFT/include -I$FFTW_DOUBLE_NEON/include wbench.c \
  kfrd.o <KFR objects> $FFTW_DOUBLE_NEON/lib/libfftw3.a $PFFFT/lib/libpffft.a -framework Accelerate -lc++ -o wbench
./wbench 4096 2048 22 1                # one row of the one-core table
```

`welch_sme.py` loads `libap0v2.dylib` from [../sme-fp64-batched-fft](../sme-fp64-batched-fft/); build it there from
`ap0lib2.c`. Needs: an Apple M4, PFFFT, FFTW double with NEON, KFR, numpy and scipy (MNE for the MNE rows). Time:
checks about a minute; timing tables a few minutes. Memory under 1 GB.

Not re-run for this revision.
