// au-h0 probe: pocketfft (the C++ header numpy 2.x and scipy vendor) with and without its vector path.
// Built twice: default (clang: VLEN double = 2 on NEON) and -DPOCKETFFT_NO_VECTORS (what an MSVC build gets,
// by the header's own compiler test). Batched complex128 / real fp64 rows, the numpy call shape fft(x, axis=-1).
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <time.h>
#include "pocketfft_hdronly.h"
using namespace pocketfft;
static double now() { return clock_gettime_nsec_np(CLOCK_UPTIME_RAW) * 1e-9; }
int main(int argc, char **argv) {
  size_t rows = atoi(argv[1]), n = atoi(argv[2]); int r2c_ = atoi(argv[3]);
  std::vector<std::complex<double>> a(rows * n), b(rows * n);
  std::vector<double> ra(rows * n);
  srand(1);
  for (size_t i = 0; i < rows * n; i++) { a[i] = {rand() / (double)RAND_MAX - .5, rand() / (double)RAND_MAX - .5}; ra[i] = a[i].real(); }
  shape_t shape{rows, n}; shape_t axes{1};
  stride_t sc{(ptrdiff_t)(n * sizeof(std::complex<double>)), (ptrdiff_t)sizeof(std::complex<double>)};
  stride_t sr{(ptrdiff_t)(n * sizeof(double)), (ptrdiff_t)sizeof(double)};
  stride_t so{(ptrdiff_t)((n / 2 + 1) * sizeof(std::complex<double>)), (ptrdiff_t)sizeof(std::complex<double>)};
  auto run = [&]() {
    if (r2c_) r2c(shape, sr, so, axes, FORWARD, ra.data(), b.data(), 1.0, 1);
    else c2c(shape, sc, sc, axes, FORWARD, a.data(), b.data(), 1.0, 1);
  };
  run();
  // checksum so both builds can be compared
  double cs = 0; for (size_t i = 0; i < rows * (r2c_ ? n / 2 + 1 : n); i++) cs += std::abs(b[i]);
  std::vector<double> t;
  for (int r = 0; r < 9; r++) {
    long it = 1; double dt;
    for (;;) { double t0 = now(); for (long i = 0; i < it; i++) run(); dt = now() - t0; if (dt > 0.2) break; it *= 2; }
    t.push_back(dt / it / (rows * n) * 1e9);
  }
  std::sort(t.begin(), t.end());
#ifdef POCKETFFT_NO_VECTORS
  const char *v = "NO_VECTORS(msvc-like)";
#else
  const char *v = "vectors(clang/gcc)";
#endif
  printf("RESULT %-22s %s rows=%zu n=%zu best %.4f med %.4f ns/pt checksum %.10e\n", v, r2c_ ? "r2c" : "c2c", rows, n, t[0], t[4], cs);
}
