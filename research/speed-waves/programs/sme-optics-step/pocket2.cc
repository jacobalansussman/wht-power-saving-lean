// pocket2.cc -- as-1: C wrapper around pocketfft (header-only, BSD-3; the FFT inside numpy and scipy.fft), complex128,
// in place, interleaved, forward or backward, with a scale factor.
#include "pocketfft_hdronly.h"
using namespace pocketfft::detail;
extern "C" {
void *pk2_new(int n) { return new pocketfft_c<double>((size_t)n); }
void pk2_exec(void *p, double *c, int fwd, double fct) { ((pocketfft_c<double> *)p)->exec((cmplx<double> *)c, fct, fwd != 0); }
}
