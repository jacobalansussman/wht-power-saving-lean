# Batched fp32 FFTs on the M4's matrix unit (section 4)

A 256-point complex FFT as 16 × 16 with each pass a DFT matrix applied with outer products, batched, with the
drain of one batch overlapped with the next. Surviving result: **-28 % per core** against vDSP's DFT at batch 64; on
the whole chip only -5 % to -12 % in cache and nothing from DRAM; it loses at batch 1. Status: checked by computer, one critic.

| file | what it does |
|---|---|
| `gen.py`, `gen9.py` | write the generated kernels `fft256_gen.h`, `fft256_gen2.h`, `fft256_gen3.h` and `fft256_gen9.h` |
| `fft256.c` | the kernels (versions v1 to v9) against vDSP and PFFFT; `check` and `time` |
| `fft400.c` | the 400-point (Whisper frame) variant, with FFTW as the baseline |
| `contend.c` | shows that vDSP itself uses the matrix unit (it slows down when the unit is busy) |
| `results/` | the original outputs (`final_same.txt`, `final_separate.txt`) |

```sh
python3 gen.py fft256_gen.h && python3 gen9.py fft256_gen9.h
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$PFFFT/include fft256.c $PFFFT/lib/libpffft.a -framework Accelerate -o fft256
./fft256 check                                                   # every line PASS
PROGS=0,4,6,7,9,10,11,12 ./fft256 time same 64 256
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$PFFFT/include -I$FFTW/include fft400.c \
  $PFFFT/lib/libpffft.a $FFTW/lib/libfftw3f.a -framework Accelerate -o fft400 && ./fft400 check
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$PFFFT/include contend.c $PFFFT/lib/libpffft.a -framework Accelerate -o contend && ./contend
```

Program indices in `PROGS`: 0 v1, 4 v5, 6 v7, 7 v8, 8 v9, 9 `vDSP_DFT_Execute`, 10 `vDSP_fft_zop`, 11 `vDSP_fftm_zop`,
12 PFFFT. Needs: an Apple M4, PFFFT, FFTW single for `fft400`. Time: seconds to a minute. Memory under 1 GB.

Re-run for this revision: `fft256 check` gave output identical to `results/check256.txt`.
