# Wrappers around KFR

Two small C entry points around KFR's double-precision complex DFT, used as the strongest double-precision baseline
in several folders (sme-optics-step, sme-welch, sme-cluster-trylock).

- `kfrd.cc`: `kd_new(n)`, `kd_exec(h, in, out)`: out-of-place, interleaved complex128 (numpy layout), forward,
  unnormalised.
- `kfrb.cc`: `kb_batch(h, n, x, y, B)`: B rows through `kd_exec`, for calls from Python.

KFR is GPL-2 or commercial and is not included. It was used only for measurement. Build KFR (commit 7c1e55c of
kfrlib/kfr) with its own CMake files, then:

```sh
clang++ -std=c++20 -O3 -mcpu=apple-m4 -I$KFR/include -c kfrd.cc -o kfrd.o
# libkfrb.dylib (used by sme-optics-step with KFRW=1):
clang++ -std=c++20 -O3 -mcpu=apple-m4 -shared -fPIC -I$KFR/include kfrd.cc kfrb.cc <KFR's dft objects or library> -o libkfrb.dylib
```
