# Programs of the speed-first waves

Each folder holds the programs behind one group of findings in [../README.md](../README.md), their result files
(`results/`), and a README with the commands. The file headers keep the names of the rounds that wrote them; the
section's Sources map those names to these folders.

## What you need

- **Processor programs**: an Apple silicon Mac. Everything named `sme-*` needs an **Apple M4** (or a later chip with
  SME); the NEON programs run on any arm64 Mac. Apple clang (the Xcode command-line tools) and the Accelerate
  framework. Build SME code with `-mcpu=apple-m4`; `-march=armv9*` builds crash with an illegal instruction.
- **Python programs**: Python 3.12 with numpy 2.x and scipy (the original runs used numpy 2.5.4 and scipy 1.18.1,
  some numpy 2.0.2 and scipy 1.13.1).
- **FPGA programs**: the open FPGA tools. The runs used Yosys 0.69, nextpnr 0.11, Icarus Verilog and Verilator 5 from
  YosysHQ's OSS CAD Suite. Set `OSS_CAD_SUITE` to its folder, or put the tools on your `PATH`.

## Third-party libraries (not included; build them yourself)

The programs compare against these libraries. None is copied here. The commands below set environment variables
that the folder READMEs use.

| variable | library | how it was built |
|---|---|---|
| `PFFFT` (folder with `include/` and `lib/libpffft.a`) | PFFFT, the marton78/pffft fork (BSD-like) | `clang -O3 -mcpu=apple-m4 -DPFFFT_ENABLE_NEON -c pffft.c pffft_common.c pffft_double.c && ar rcs libpffft.a pffft.o pffft_common.o pffft_double.o`. **Without `-DPFFFT_ENABLE_NEON` it is scalar.** |
| `POCKETFFT` (folder with `pocketfft_hdronly.h`) | pocketfft, C++ header (Martin Reinecke, BSD-3), the FFT inside numpy and scipy | header only |
| `FFHT` (folder with `include/fht.h` and `lib/libffht.a`) | FFHT (MIT), the copy inside PyTorch's ExecuTorch (`extension/llm/custom_ops/spinquant/third-party/FFHT`) | `clang -O3 -mcpu=apple-m4 -std=c99 -c fht.c fast_copy.c && ar rcs libffht.a fht.o fast_copy.o`; float only on Arm |
| `FFTW` (install prefix) | FFTW 3.3.10, single precision with NEON (GPL) | `CC=clang CFLAGS="-O3 -mcpu=apple-m4" ./configure --prefix=$FFTW --enable-static --disable-shared --disable-fortran --enable-armv8-cntvct-el0 --enable-float --enable-neon` |
| `FFTW_DOUBLE_NEON` (install prefix) | FFTW 3.3.10, double precision with NEON (GPL) | `CC=clang CFLAGS="-O3 -mcpu=apple-m4" ./configure --prefix=$FFTW_DOUBLE_NEON --build=aarch64-apple-darwin --host=aarch64-apple-darwin --enable-neon --enable-armv8-cntvct-el0 --disable-fortran --disable-shared --enable-static`. **Without the host triple, double precision has no NEON; without `--enable-armv8-cntvct-el0` FFTW has no timer on macOS and its MEASURE / PATIENT plans are untuned.** Plan wisdom on a performance core. |
| `KFR` | KFR (GPL-2 or commercial), commit 7c1e55c of kfrlib/kfr | its own CMake build; [baseline-wrappers](baseline-wrappers/) holds the two small C wrappers the programs call (`kfrd.cc`, `kfrb.cc`) |

FFTW and KFR are GPL. They were used only for measurement, and nothing here links them into anything distributed.

## Folders

| folder | section | contents |
|---|---|---|
| [timing-harness](timing-harness/) | method | `bench.h`, the shared timing harness |
| [baseline-wrappers](baseline-wrappers/) | method | C wrappers around KFR for the double-precision baselines |
| [fpga-measure](fpga-measure/) | method | synthesis and place-and-route measurement, equivalence checking, UP5K timing |
| [neon-wht](neon-wht/) | 4 | the NEON Walsh-Hadamard baseline and the FFT library baselines |
| [sme-fp32-fir](sme-fp32-fir/) | 4 | Toeplitz-tile FIR on the matrix unit |
| [sme-fp32-fft](sme-fp32-fft/) | 4 | batched 256- and 400-point fp32 FFTs on the matrix unit |
| [sme-fp64-batched-fft](sme-fp64-batched-fft/) | 3a | batched complex128 FFTs on the matrix unit |
| [sme-optics-step](sme-optics-step/) | 3b | the angular-spectrum optics step |
| [sme-welch](sme-welch/) | 3c | one-pass Welch spectrum |
| [sme-fp64-convolve](sme-fp64-convolve/) | 3d | numpy-compatible convolve and correlate |
| [sme-cluster-trylock](sme-cluster-trylock/) | 3e | per-cluster try-lock for multi-process code |
| [fftw-neon-probe](fftw-neon-probe/) | 5 | the shipped-FFTW and pocketfft probes |
| [packed-pfb](packed-pfb/) | 6a | the telescope filter bank with two signals per DSP |
| [carry-folded-sdf-fft](carry-folded-sdf-fft/) | 6b | the streaming FFT with multiplexers in carry cells |
| [up5k-folded-fft](up5k-folded-fft/) | 6c | folded FFTs on the iCE40-UP5K |
| [up5k-streaming-fft-and-wht](up5k-streaming-fft-and-wht/) | 6c | streaming 4,096-point FFT and WHT on the iCE40-UP5K |

## Timing

The timings in `results/` ran on one P-core of an Apple M4 Mac mini while other jobs used the efficiency cores. Each
result file starts with lines such as `[timing lock] got lock 23:08:40; load: { 22.28 20.39 14.99 }`, which record
the machine-wide timing lock and the load average at the time. To time on your own machine, run the same commands on
a quiet machine; on macOS the harness asks for the user-interactive quality-of-service class, which keeps the thread on
a performance core.
