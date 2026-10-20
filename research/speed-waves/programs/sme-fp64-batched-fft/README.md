# Batched complex128 FFTs on the M4's matrix unit (section 3a)

A transform of 8^P points as P radix-8 passes, each pass a complex 8 × 8 DFT matrix applied with double-precision
outer products (32 per group of 64 points); two groups accumulate while the previous job's tiles drain; passes of
different transforms are interleaved. 1,024 and 2,048 points add one NEON pass. numpy's memory layout (interleaved
complex128 rows). Surviving result, one P-core, batch 64, in cache: **-41 % at 512 points** against KFR and **-44 % at
4,096** against FFTW double with NEON (the -52 % first measured at 4,096 was against an FFTW plan that was 17-20 %
slow; `baseline-recheck/` holds the program that found it). Status: checked by computer, one critic plus the baseline
recheck.

| file | what it does |
|---|---|
| `sme64.h` | the 512-point kernel (three passes), the first version |
| `smeg.h` | the same plan for N = 8^P (512, 4,096, 32,768) with prefetch |
| `smegi.h` | `smeg.h` with interleaved complex input and output (numpy layout, no conversion pass) |
| `smegp.h`, `ap0p.c` | 2 · 8^P points by decimation in time (1,024, 8,192), the last butterfly in NEON |
| `ap0lib.c`, `ap0lib2.c` | the numpy-layout front end for every size, and the library built for Python |
| `fftnp.c`, `fftnp2.c` | check and time against PFFFT double, pocketfft and vDSP (`fftnp N check B`, `fftnp N same B`) |
| `fft512.c`, `fftg.c` | the split-complex kernels alone |
| `cluster.c` | the whole P-cluster: one matrix-unit thread plus NEON threads |
| `probe.c` | the probe of the unit's double-precision outer-product rate |
| `e2e.py`, `ap0_backend.py` | from Python through ctypes; a sketch of a `scipy.fft` backend |
| `pocketd.cc` | a C entry point to pocketfft (the baseline numpy and scipy use) |
| `baseline-recheck/smevsfftw.c` | the kernel against a correctly built FFTW double with NEON and P-core wisdom |
| `results/` | original outputs: `tnp_512_4096_v2.txt` (one core), `tcluster.txt` (cluster), `te2e.txt` (Python), `baseline-recheck/results/t_sme.txt` |

```sh
clang++ -O3 -mcpu=apple-m4 -std=c++17 -I$POCKETFFT -c pocketd.cc -o pocketd.o
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$PFFFT/include fftnp.c pocketd.o $PFFFT/lib/libpffft.a \
  -framework Accelerate -lc++ -o fftnp
./fftnp 512 check 64       # CHECK SME fp64 interleaved n=512 B=64 max|err|/max|X| 4.769e-16
./fftnp 4096 check 64      # 5.338e-16
BROKEN=1 ./fftnp 2048 check 64    # a flipped twiddle: the SME line shows an error near 6.6e-4 (it must fail)
./fftnp 4096 same 64       # timing, against the libraries
# for Python: build libap0.dylib from ap0lib.c (and libap0v2.dylib from ap0lib2.c) as shared libraries:
clang -O3 -mcpu=apple-m4 -shared -fPIC ap0lib.c -o libap0.dylib && python3 e2e.py 512 64 4096 64
```

Needs: an Apple M4, PFFFT, pocketfft; KFR and a correctly built FFTW double with NEON for the strongest baselines
(see [../README.md](../README.md)). Time: the checks take seconds; each timing table a few minutes. Memory under 1 GB
(DRAM runs use 16-32 MiB arrays).

Re-run for this revision: the three check commands above gave the errors shown, identical to the original binary.
