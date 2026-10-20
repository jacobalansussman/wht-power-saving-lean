# A complex128 Fourier-optics step on the M4's matrix unit (section 3b)

One angular-spectrum propagation step, `ifft2(fft2(E) * H)` in complex128: row passes with the batched kernel of
[../sme-fp64-batched-fft](../sme-fp64-batched-fft/), column passes through vertical slices of the matrix unit's tiles
(`colfft.h`), the transfer function multiplied in between, output in index-reversed order. Surviving result, one
P-core, from DRAM: **-38 % / -47 % / -62 %** at 1,024² / 2,048² / 4,096² against FFTW double with NEON using its own
2-D plan; on the whole P-cluster a tie at 1,024² (the claimed cluster gain there is refuted). Status: checked by computer, one critic and the fairness critic.

| file | what it does |
|---|---|
| `colfft.h` | the column FFT through the matrix unit's vertical slices |
| `prop.c` | the step and its baselines (PFFFT or KFR rows plus transpose, FFTW 2-D); modes `check`, `step`, `fft2`, `parts`, and the cluster modes; settings by environment variables (`SW` strip width, `KFRW=1` KFR workers, `TAIL`, `NOFUSE`, `NOFFTW`, `BROKEN`) |
| `as1lib.c`, `prop_py.py` | the step as a library for Python, against numpy and scipy |
| `pocket2.cc` | a C entry point to pocketfft |
| `critic/prop_au.c` | the critic's version with FFTW's own 2-D plans (1 and 4 threads, PATIENT) |
| `critic-fairness/prop_t.c`, `mlx_asp.py` | the fairness critic's FFTW timings and the GPU (MLX, complex64) comparison |
| `results/`, `critic/results/`, `critic-fairness/results/` | original outputs (`tfinal.txt`, `td.txt`, `tj.txt` one core; `tmtk.txt`, `tn.txt` cluster; `critic/results/t_c*.txt` the critic's tables) |

```sh
clang++ -O3 -mcpu=apple-m4 -std=c++17 -I$POCKETFFT -c pocket2.cc -o pocket2.o
# libkfrb.dylib from ../baseline-wrappers (KFR) in this folder, then:
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -Wno-uninitialized -I$PFFFT/include -I$FFTW_DOUBLE_NEON/include \
  prop.c pocket2.o $PFFFT/lib/libpffft.a $FFTW_DOUBLE_NEON/lib/libfftw3.a -L. -lkfrb -Wl,-rpath,@executable_path \
  -framework Accelerate -lc++ -o prop
./prop 1024 check; ./prop 2048 check                  # every SME CHECK line PASS
BROKEN=1 ./prop 1024 check                            # one flipped column twiddle: the SME CHECK lines must FAIL
./prop 1024 step                                      # one P-core timing table
SW=16 NOFFTW=1 ./prop 4096 step
SW=16 KFRW=1 TAIL=3 ./prop 2048 mt                    # whole P-cluster
```

FFTW plans: the program loads FFTW wisdom from `fftw-d-wisdom.txt` (or `FFTW_WF`); plan on a performance core.
Needs: an Apple M4, PFFFT, pocketfft, FFTW double with NEON, KFR. Time: checks under a minute (4,096² takes longer);
timing tables a few minutes. Memory: 1,024² and 2,048² well under 1 GB; 4,096² about 1-2 GB.

Not re-run for this revision.
