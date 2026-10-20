# numpy's convolve and correlate on the M4's matrix unit (section 3d, wave 15)

A drop-in for `np.convolve` and `np.correlate` in float64 and complex128 that computes the direct sum with Toeplitz
tiles on the matrix unit (all 8 double-precision tiles; the input transposed through the unit in streaming mode),
NEON for short calls, and scipy's `oaconvolve` above a measured crossover. Same modes and argument rules as numpy.
Result: **-68 % to -87 %** against the best of numpy, scipy and vDSP from 9 to 1,025 taps (whole Python call, 10^6
samples); against a tuned C FFT overlap-save it wins only below about 450 taps. Status: checked by computer, one
program, not repeated by anyone (no critic read it).

| file | what it does |
|---|---|
| `fc64.c` | the kernels: matrix-unit Toeplitz tiles, blocked accumulation (the accurate default), NEON and scalar paths, multi-thread split, NaN check inside the transpose |
| `_fc64ext.c` | an optional CPython extension (lower call overhead) |
| `fc64.py` | the drop-in (`fc64.convolve`, `fc64.correlate`) with the dispatcher; loads `libfc64.dylib` |
| `check.py` | correctness against numpy and a correctly rounded reference (`quick`, `full`, `broken`) |
| `check_nan3.py`, `check_auto.py` | NaN and infinity patterns; the dispatcher's choices |
| `bench.py`, `summarize.py` | the timing grid against numpy, scipy and vDSP; the tables |
| `olsd.c` | the C-level comparison with a PFFFT double overlap-save (the strongest FFT route) |
| `err_table.py`, `errb5.py` | the error tables |
| `results/` | original outputs (`t3_same_1e6a.txt` and `t3_same_1e6b.txt` are the headline grid; `t17_olsd_*.txt` the C comparison) |

```sh
clang -O3 -mcpu=apple-m4 -shared fc64.c -framework Accelerate -o libfc64.dylib
PYINC=$(python3 -c "import sysconfig; print(sysconfig.get_paths()['include'])")
EXT=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
clang -O3 -mcpu=apple-m4 -shared -undefined dynamic_lookup -I$PYINC _fc64ext.c fc64.c -framework Accelerate -o _fc64ext$EXT
python3 check.py quick        # SUMMARY 534 of 534 PASS
python3 check.py broken       # BROKEN COPY: 4 of 4 checks FAIL (all must FAIL)
python3 check_nan3.py         # NaN / inf pattern = numpy's (compare with results/check_nan3.txt)
VECLIB_MAXIMUM_THREADS=1 python3 bench.py same f64 1e6 3,9,33,129 fc_sme,fc_smeb8,np,sp_auto,sp_oa,sp_fft,vdsp
clang -O3 -mcpu=apple-m4 -I$PFFFT/include olsd.c fc64.c $PFFFT/lib/libpffft.a -framework Accelerate -o olsd
./olsd scan 257 20; ./olsd time 257 20 2048 55
```

**The multi-thread path leaked memory before a fix made after the crash of 12 October** (each helper thread's
scratch buffers were never freed). The copy here is the fixed one: scratch now lives in a heap block freed when the
thread ends, and a helper that cannot be started has its share computed by the caller. Outputs were bit-identical to
the old builds in all 3,471 cases tried. The per-cluster try-lock of [../sme-cluster-trylock](../sme-cluster-trylock/)
is not wired in, so do not call this from many processes at once.

Needs: an Apple M4, Python 3.12 with numpy and scipy, PFFFT for `olsd`. Time: `check.py quick` about a minute; the
timing grid several minutes. Memory: under 1 GB up to 10^6 samples, about 1 GB at 10^7.

Re-run for this revision: `check.py quick` 534 of 534 PASS, `check.py broken` 4 of 4 FAIL, and `check_nan3.py` output
identical to `results/check_nan3.txt`.
