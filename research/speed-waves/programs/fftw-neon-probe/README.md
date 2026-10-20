# The FFTW that people install on Arm has no vector code (section 5)

| file | what it does |
|---|---|
| `shipfftw.c` | loads conda-forge's FFTW 3.3.10 for osx-arm64 with `dlopen` and times it against FFTW 3.3.10 built with NEON, same API, same plan flag, taking turns in one process |
| `fourstep.py` | numpy's single long complex128 FFT against a four-step built from numpy's own batched calls |
| `pfvec.cpp` | pocketfft's batched rows with and without its vector path (`-DPOCKETFFT_NO_VECTORS`) |
| `results/` | `t_shipfftw.txt` (the 1.27-2.6 times), `t_fourstep.txt`, `t_pfvec.txt` (vectors 30 % slower at 64 × 4,096) |

To count the NEON codelets in any FFTW library: `strings libfftw3.dylib | grep -c _neon` (0 in conda-forge's
osx-arm64 and linux-aarch64 builds of both precisions, 583 in a correct NEON build).

```sh
# download conda-forge's fftw package for osx-arm64 (3.3.10, build 112) and unpack it into ./sf; then:
clang -O3 -mcpu=apple-m4 -I$FFTW_DOUBLE_NEON/include shipfftw.c $FFTW_DOUBLE_NEON/lib/libfftw3.a $FFTW/lib/libfftw3f.a -o shipfftw
./shipfftw sf/lib/libfftw3.dylib c2c 512      # LIB KIND N; KIND c2c or r2c; libfftw3f.dylib for single precision
./shipfftw sf/lib/libfftw3.dylib c2c 4096     # compare with results/t_shipfftw.txt
python3 fourstep.py
clang++ -std=c++17 -O3 -mcpu=apple-m4 -I$POCKETFFT pfvec.cpp -o pf_vec
clang++ -std=c++17 -O3 -mcpu=apple-m4 -DPOCKETFFT_NO_VECTORS -I$POCKETFFT pfvec.cpp -o pf_novec
./pf_vec; ./pf_novec
```

The first argument is the downloaded library; it is opened with `dlopen`, so its symbols do not clash with the NEON
builds linked into the program.
Status: checked by computer, one program, not repeated by anyone. Needs: an arm64 Mac, both FFTW builds, pocketfft.
Time: a few minutes. Not re-run for this revision.
