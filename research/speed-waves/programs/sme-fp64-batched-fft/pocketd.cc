// pocketd.cc -- ap-0: C wrapper around pocketfft (C++ header-only, BSD-3; the FFT inside numpy >= 1.17 and scipy.fft),
// complex128. Same plan class scipy.fft uses for c2c (pocketfft_c<double>). In place, interleaved.
#include "pocketfft_hdronly.h"
using namespace pocketfft::detail;
extern "C" {
void *pk_cd_new(int n) { return new pocketfft_c<double>((size_t)n); }
void pk_cd_exec(void *p, double *c) { ((pocketfft_c<double> *)p)->exec((cmplx<double> *)c, 1.0, true); }
}
