// kfrd.cc -- ar-0: C wrapper around KFR (kfrlib/kfr, GPL-2 / commercial; private measurement only) complex<double> DFT,
// out of place, interleaved (complex128 = numpy layout), forward, unnormalised.
#include <kfr/dft.hpp>
#include <cstdio>
using namespace kfr;
struct kd { dft_plan<double> *p; u8 *tmp; };
extern "C" {
void *kd_new(long n) { kd *h = new kd; h->p = new dft_plan<double>((size_t)n); h->tmp = (u8 *)aligned_alloc(64, ((h->p->temp_size + 63) / 64) * 64 + 64); return h; }
void kd_exec(void *v, const double *in, double *out) { kd *h = (kd *)v; h->p->execute((complex<double> *)out, (const complex<double> *)in, h->tmp, false); }
void kd_dump(void *v) { ((kd *)v)->p->dump(); }
}
