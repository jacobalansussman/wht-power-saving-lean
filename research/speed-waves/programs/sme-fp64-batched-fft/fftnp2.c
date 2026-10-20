// fftnp.c -- ap-0: numpy layout (complex128 interleaved rows). SME front end (ap0lib.c) vs PFFFT-d ordered,
// pocketfft (scipy's c2c: copy in->out then exec in place), vDSP interleaved DFT (double).
// usage: fftnp N check B | fftnp N same|separate B
#include <Accelerate/Accelerate.h>
#include "../timing-harness/bench.h"
#include "pffft_double.h"
#include "ap0lib.c"
#include "ap0p.c"
static ap0p_t *H2;
void *pk_cd_new(int n); void pk_cd_exec(void *p, double *c);
static long NN; static long BB; static bench_pool pool; static ap0_t *H; static PFFFTD_Setup *pf; static double *work; static void *pk; static vDSP_DFT_Interleaved_SetupD vi;
static void r_sme(void *c, long i) { double *b = bench_pool_get(&pool, i); ap0_fft(H, b, b + 2 * BB * NN, BB); }
static void r_sme2(void *c, long i) { double *b = bench_pool_get(&pool, i); ap0p_fft(H2, b, b + 2 * BB * NN, BB); }
static void r_pffft(void *c, long i) { double *b = bench_pool_get(&pool, i); for (long t = 0; t < BB; t++) pffftd_transform_ordered(pf, b + 2 * t * NN, b + 2 * BB * NN + 2 * t * NN, work, PFFFT_FORWARD); }
static void r_pocket(void *c, long i) { double *b = bench_pool_get(&pool, i); memcpy(b + 2 * BB * NN, b, 16 * BB * NN); for (long t = 0; t < BB; t++) pk_cd_exec(pk, b + 2 * BB * NN + 2 * t * NN); }
static void r_vint(void *c, long i) { double *b = bench_pool_get(&pool, i); for (long t = 0; t < BB; t++) vDSP_DFT_Interleaved_ExecuteD(vi, (const DSPDoubleComplex *)(b + 2 * t * NN), (DSPDoubleComplex *)(b + 2 * BB * NN + 2 * t * NN)); }
static double *cw, *sw;
static void two_sum(double a, double b, double *s, double *e) { *s = a + b; double z = *s - a; *e = (a - (*s - z)) + (b - z); }
static void ref_bin(const double *x, long k, double *Xr, double *Xi) {
  double sr = 0, er = 0, si = 0, ei = 0;
  for (long j = 0; j < NN; j++) { long m = (j * k) % NN; double c = cw[m], s = sw[m], xr = x[2 * j], xi = x[2 * j + 1];
    double t0 = xr * c, t1 = xi * s, t2 = xi * c, t3 = -xr * s, t, e;
    two_sum(sr, t0, &t, &e); sr = t; er += e + fma(xr, c, -t0); two_sum(sr, t1, &t, &e); sr = t; er += e + fma(xi, s, -t1);
    two_sum(si, t2, &t, &e); si = t; ei += e + fma(xi, c, -t2); two_sum(si, t3, &t, &e); si = t; ei += e + fma(-xr, s, -t3); }
  *Xr = sr + er; *Xi = si + ei;
}
static void check(const char *name, const double *b) {
  double emax = 0, xmax = 0, se = 0, sx = 0; long ts[3] = {0, BB / 2, BB - 1}; long ks = NN > 2048 ? 37 : 1;
  for (int q = 0; q < 3; q++) for (long k = 0; k < NN; k += ks) { double Rr, Ri; ref_bin(b + 2 * ts[q] * NN, k, &Rr, &Ri);
    const double *g = b + 2 * BB * NN + 2 * (ts[q] * NN + k); double e = hypot(g[0] - Rr, g[1] - Ri), x = hypot(Rr, Ri);
    if (!(e <= emax)) emax = e; if (x > xmax) xmax = x; se += e * e; sx += x * x; }
  printf("CHECK %-26s n=%ld B=%ld  max|err|/max|X| %.3e  rel RMS %.3e\n", name, NN, BB, emax / xmax, sqrt(se / sx));
}
int main(int argc, char **argv) {
  bench_init(); if (getenv("SKIP")) ap0_skip = atoi(getenv("SKIP")); NN = atol(argv[1]); const char *mode = argv[2]; BB = atol(argv[3]); int sep = !strcmp(mode, "separate");
  cw = malloc(8 * NN); sw = malloc(8 * NN); for (long m = 0; m < NN; m++) { cw[m] = cos(2 * M_PI * m / NN); sw[m] = sin(2 * M_PI * m / NN); }
  H = ap0_new(NN); H2 = ap0p_new(NN); pf = pffftd_new_setup((int)NN, PFFFT_COMPLEX); work = pffftd_aligned_malloc(16 * NN); pk = pk_cd_new((int)NN);
  vi = vDSP_DFT_Interleaved_CreateSetupD(NULL, NN, vDSP_DFT_FORWARD, 0);
  bench_pool_init(&pool, 32 * BB * NN, sep ? BENCH_SEPARATE : BENCH_SAME);
  for (int p = 0; p < pool.count; p++) bench_fill_d(pool.buf[p], 2 * BB * NN);
  if (!strcmp(mode, "check")) {
    double *b = pool.buf[0];
    if (getenv("BROKEN")) { H->wi[3] = -H->wi[3]; H2->wi[37] = -H2->wi[37]; }
    r_sme2(NULL, 0); check("SME DIT-2 (LD4D, no pre-pass)", b);
    r_sme(NULL, 0); check(H->R > 1 && getenv("BROKEN") ? "SME (BROKEN radix-R twiddle)" : "SME fp64 interleaved", b);
    r_pffft(NULL, 0); check("PFFFT-d ordered", b); r_pocket(NULL, 0); check("pocketfft c2c", b);
    if (vi) { r_vint(NULL, 0); check("vDSP_DFT_Interleaved_ExecD", b); } else printf("vDSP interleaved double: no setup at n=%ld\n", NN);
    return 0;
  }
  bench_prog progs[5] = {{"SME DIT-2 (LD4D)", r_sme2, NULL, NULL}, {"SME fp64 interleaved", r_sme, NULL, NULL}, {"PFFFT-d ordered", r_pffft, NULL, NULL},
                         {"pocketfft c2c copy+exec", r_pocket, NULL, NULL}, {"vDSP_DFT_Interleaved_D", r_vint, NULL, NULL}};
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)(BB * NN); bench_stat st[5]; char title[128];
  snprintf(title, sizeof title, "c2c fp64 interleaved n=%ld B=%ld %s", NN, BB, sep ? "separate" : "same");
  bench_compare(title, progs, getenv("ONLY_SME") ? 2 : (vi ? 5 : 4), &o, st);
  return 0;
}
