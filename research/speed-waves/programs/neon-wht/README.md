# The NEON Walsh-Hadamard baseline and the FFT library baselines (section 4)

`wht_neon.h` is a header-only NEON Walsh-Hadamard transform (fp32, in place, unnormalised, natural order), written as
the shared WHT baseline in wave 11. It takes 8 % to 24 % less time than FFHT's NEON code at each of the six sizes
tested from 2^8 to 2^22, one thread (`results/wht_same.txt`, `results/wht_separate.txt`). Status: checked by computer, one program, not
repeated by anyone.

| file | what it does |
|---|---|
| `wht_neon.h` | the transform: a 64-float register kernel, radix-16/8/4/2 register passes, cache blocking above 2^14 |
| `wht_bench.c` | checks FFHT and `wht_neon.h` against a double-precision WHT for 2^1 to 2^24, then times both |
| `fft_bench.c` | the FFT baselines of wave 11: vDSP, FFTW single with NEON, PFFFT, complex and real |
| `accel_example.c` | small checked examples of vDSP calls (FFTs, `conv`, `desamp`) |
| `results/` | the original outputs |

```sh
clang -O3 -mcpu=apple-m4 -I$FFHT/include wht_bench.c $FFHT/lib/libffht.a -o wht_bench
./wht_bench check                                  # every size: PASS; last line ALL PASS
./wht_bench time same 8 10 12 16 20 22             # -> compare with results/wht_same.txt
./wht_bench time separate 8 10 12 16 20 22
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$FFTW/include -I$PFFFT/include fft_bench.c \
  $PFFFT/lib/libpffft.a $FFTW/lib/libfftw3f.a -framework Accelerate -o fft_bench
./fft_bench c same 256 1024 4096 65536 1048576
./fft_bench r same 1024 65536
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations accel_example.c -framework Accelerate -o accel_example && ./accel_example
```

Needs: any arm64 Mac; FFHT, PFFFT and FFTW (single) for the baselines. Time: the check takes seconds; each timing
table about a minute. Memory under 1 GB.

Re-run for this revision: `wht_bench check` from this folder gave output identical to `results/wht_check.txt`.
