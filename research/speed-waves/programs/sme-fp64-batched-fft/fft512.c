// fft512.c -- ap-0: check and time the SME fp64 512-point batched FFT against vDSP double, PFFFT double, pocketfft.
// usage: fft512 check B | fft512 time same|separate B
#include <Accelerate/Accelerate.h>
#include "../timing-harness/bench.h"
#include "pffft_double.h"
#include "sme64.h"
#define NN 512
#define LOGN 9
void *pk_cd_new(int n); void pk_cd_exec(void *p, double *c);
typedef struct { bench_pool pool; long B; vDSP_DFT_SetupD dft; FFTSetupD fs; PFFFTD_Setup *pf; double *work; void *pk; } ctx_t;
static ctx_t C;
// pool buffer layout: split view [xr BN][xi BN][yr BN][yi BN]; interleaved view [x 2BN][y 2BN]
static void r_sme(void *c, long i) { double *b = bench_pool_get(&C.pool, i); long M = C.B * NN; fft512_sme(b, b + M, b + 2 * M, b + 3 * M, C.B); }
static void r_vdft(void *c, long i) { double *b = bench_pool_get(&C.pool, i); long M = C.B * NN;
  for (long t = 0; t < C.B; t++) vDSP_DFT_ExecuteD(C.dft, b + t * NN, b + M + t * NN, b + 2 * M + t * NN, b + 3 * M + t * NN); }
static void r_vfft(void *c, long i) { double *b = bench_pool_get(&C.pool, i); long M = C.B * NN;
  for (long t = 0; t < C.B; t++) { DSPDoubleSplitComplex a = {b + t * NN, b + M + t * NN}, o = {b + 2 * M + t * NN, b + 3 * M + t * NN};
    vDSP_fft_zopD(C.fs, &a, 1, &o, 1, LOGN, FFT_FORWARD); } }
static void r_pffft(void *c, long i) { double *b = bench_pool_get(&C.pool, i); long M = C.B * NN;
  for (long t = 0; t < C.B; t++) pffftd_transform_ordered(C.pf, b + 2 * t * NN, b + 2 * M + 2 * t * NN, C.work, PFFFT_FORWARD); }
static void r_pocket(void *c, long i) { double *b = bench_pool_get(&C.pool, i); long M = C.B * NN;  // scipy's c2c: copy in->out, exec in place
  memcpy(b + 2 * M, b, sizeof(double) * 2 * M); for (long t = 0; t < C.B; t++) pk_cd_exec(C.pk, b + 2 * M + 2 * t * NN); }

// reference: DFT with compensated (double-double) accumulation, exact index reduction
static double cw[NN], sw[NN];
static void two_sum(double a, double b, double *s, double *e) { *s = a + b; double z = *s - a; *e = (a - (*s - z)) + (b - z); }
static void ref_bin(const double *xr, const double *xi, int k, double *Xr, double *Xi) {
  double sr = 0, er = 0, si = 0, ei = 0;
  for (int j = 0; j < NN; j++) { int m = (int)(((long)j * k) % NN); double c = cw[m], s = sw[m];  // w = c - i s
    double terms[4] = {xr[j] * c, xi[j] * s, xi[j] * c, -xr[j] * s};
    double errs[4] = {fma(xr[j], c, -terms[0]), fma(xi[j], s, -terms[1]), fma(xi[j], c, -terms[2]), fma(-xr[j], s, -terms[3])};
    double t, e; two_sum(sr, terms[0], &t, &e); sr = t; er += e + errs[0]; two_sum(sr, terms[1], &t, &e); sr = t; er += e + errs[1];
    two_sum(si, terms[2], &t, &e); si = t; ei += e + errs[2]; two_sum(si, terms[3], &t, &e); si = t; ei += e + errs[3]; }
  *Xr = sr + er; *Xi = si + ei;
}
// error of one program's output (split or interleaved) on transforms tsel; returns max err / max|X| and rel RMS
static void check_prog(const char *name, double *b, int inter, double *maxrel, double *rms) {
  long M = C.B * NN; double emax = 0, xmax = 0, se = 0, sx = 0;
  long ts[3] = {0, C.B / 2, C.B - 1};
  for (int q = 0; q < 3; q++) { long t = ts[q];
    const double *xr = b + t * NN, *xi = b + M + t * NN;  // split input copy is always kept in the split slots (see fill)
    for (int k = 0; k < NN; k++) {
      double Rr, Ri; ref_bin(xr, xi, k, &Rr, &Ri);
      double gr = inter ? b[2 * M + 2 * (t * NN + k)] : b[2 * M + t * NN + k], gi = inter ? b[2 * M + 2 * (t * NN + k) + 1] : b[3 * M + t * NN + k];
      double e = hypot(gr - Rr, gi - Ri), x = hypot(Rr, Ri); if (!(e <= emax)) emax = e; if (x > xmax) xmax = x; se += e * e; sx += x * x; } }
  *maxrel = emax / xmax; *rms = sqrt(se / sx);
  printf("CHECK %-22s B=%ld  max|err|/max|X| %.3e  rel RMS %.3e\n", name, C.B, *maxrel, *rms);
}
int main(int argc, char **argv) {
  bench_init();
  const char *mode = argc > 1 ? argv[1] : "check"; long B = argc > 2 ? atol(argv[2]) : 64; int sep = argc > 1 && !strcmp(mode, "separate");
  C.B = B; for (int m = 0; m < NN; m++) { cw[m] = cos(2 * M_PI * m / NN); sw[m] = sin(2 * M_PI * m / NN); }
  if (getenv("BROKEN")) sme64_broken = 1;
  sme64_init();
  long M = B * NN;
  bench_pool_init(&C.pool, sizeof(double) * 4 * M, sep ? BENCH_SEPARATE : BENCH_SAME);
  C.dft = vDSP_DFT_zop_CreateSetupD(NULL, NN, vDSP_DFT_FORWARD); C.fs = vDSP_create_fftsetupD(LOGN, kFFTRadix2);
  C.pf = pffftd_new_setup(NN, PFFFT_COMPLEX); C.work = pffftd_aligned_malloc(2 * NN * 8); C.pk = pk_cd_new(NN);
  for (int p = 0; p < C.pool.count; p++) { double *b = C.pool.buf[p]; bench_fill_d(b, 2 * M); }
  if (!strcmp(mode, "check")) {
    double *b = C.pool.buf[0]; double mr, rr;
    // split input in [0, 2M); make an interleaved copy for PFFFT/pocket in a separate buffer
    double *bi = bench_alloc(sizeof(double) * 4 * M); memcpy(bi, b, sizeof(double) * 2 * M);
    for (long j = 0; j < M; j++) { bi[2 * j] = b[j]; bi[2 * j + 1] = b[M + j]; }
    r_sme(NULL, 0); check_prog("SME fp64", b, 0, &mr, &rr);
    r_vdft(NULL, 0); check_prog("vDSP_DFT_ExecuteD", b, 0, &mr, &rr);
    r_vfft(NULL, 0); check_prog("vDSP_fft_zopD", b, 0, &mr, &rr);
    // interleaved programs: point the pool at bi temporarily
    void *save = C.pool.buf[0]; C.pool.buf[0] = bi;
    double *keep = bench_alloc(sizeof(double) * 2 * M); memcpy(keep, bi, sizeof(double) * 2 * M);
    r_pffft(NULL, 0); memcpy(bi, b, sizeof(double) * 2 * M);  // check_prog reads split input from b's layout
    { double *tmp = bench_alloc(sizeof(double) * 4 * M); memcpy(tmp, b, sizeof(double) * 2 * M); memcpy(tmp + 2 * M, bi + 2 * M, sizeof(double) * 2 * M); check_prog("PFFFT-d ordered", tmp, 1, &mr, &rr); free(tmp); }
    memcpy(bi, keep, sizeof(double) * 2 * M); r_pocket(NULL, 0);
    { double *tmp = bench_alloc(sizeof(double) * 4 * M); memcpy(tmp, b, sizeof(double) * 2 * M); memcpy(tmp + 2 * M, bi + 2 * M, sizeof(double) * 2 * M); check_prog("pocketfft c2c", tmp, 1, &mr, &rr); free(tmp); }
    C.pool.buf[0] = save;
    return 0;
  }
  bench_prog progs[5] = {{"SME fp64 (split)", r_sme, NULL, NULL}, {"vDSP_DFT_ExecuteD (split)", r_vdft, NULL, NULL},
                         {"vDSP_fft_zopD (split)", r_vfft, NULL, NULL}, {"PFFFT-d ordered (inter)", r_pffft, NULL, NULL},
                         {"pocketfft c2c copy+exec", r_pocket, NULL, NULL}};
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)M; bench_stat st[5]; char title[128];
  snprintf(title, sizeof title, "c2c fp64 n=512 B=%ld %s", B, sep ? "separate" : "same");
  int np = getenv("ONLY_SME") ? 2 : 5;
  bench_compare(title, progs, np, &o, st);
  return 0;
}
