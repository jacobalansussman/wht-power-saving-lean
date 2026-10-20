// smevsfftw.c -- inf-strongest: ap-0's SME complex128 kernel (wave 13, read-only include) against the timed FFTW-d NEON
// (inf-strongest/fftwd), same process, programs taking turns, ar-0's batches. usage: smevsfftw N same|separate B
#include <fftw3.h>
#include "../../timing-harness/bench.h"
#include "../ap0lib.c"
static long NN, BB; static bench_pool pool; static ap0_t *H; static fftw_plan pr, pm, pp;
static void r_sme(void *c, long i) { double *b = bench_pool_get(&pool, i); ap0_fft(H, b, b + 2 * BB * NN, BB); }
static void r_row(void *c, long i) { double *b = bench_pool_get(&pool, i); for (long t = 0; t < BB; t++) fftw_execute_dft(pr, (fftw_complex *)(b + 2 * t * NN), (fftw_complex *)(b + 2 * BB * NN + 2 * t * NN)); }
static void r_many(void *c, long i) { double *b = bench_pool_get(&pool, i); fftw_execute_dft((fftw_plan)c, (fftw_complex *)b, (fftw_complex *)(b + 2 * BB * NN)); }
int main(int argc, char **argv) {
  bench_init(); NN = atol(argv[1]); int sep = !strcmp(argv[2], "separate"); BB = atol(argv[3]);
  H = ap0_new(NN);
  double *a = bench_alloc(16 * BB * NN), *o = bench_alloc(16 * BB * NN); int n = (int)NN;
  uint64_t t0 = bench_now_ns();
  pr = fftw_plan_dft_1d(n, (fftw_complex *)a, (fftw_complex *)o, FFTW_FORWARD, FFTW_MEASURE);
  pm = fftw_plan_many_dft(1, &n, (int)BB, (fftw_complex *)a, NULL, 1, n, (fftw_complex *)o, NULL, 1, n, FFTW_FORWARD, FFTW_MEASURE);
  pp = sep ? pm : fftw_plan_many_dft(1, &n, (int)BB, (fftw_complex *)a, NULL, 1, n, (fftw_complex *)o, NULL, 1, n, FFTW_FORWARD, FFTW_PATIENT);  // DRAM batch: PATIENT skipped (planning time)
  printf("# FFTW-d (timer build) plans MEASURE row, MEASURE many, PATIENT many: %.2f s on this P-core\n", (bench_now_ns() - t0) * 1e-9);
  bench_pool_init(&pool, 32 * BB * NN, sep ? BENCH_SEPARATE : BENCH_SAME);
  for (int p = 0; p < pool.count; p++) bench_fill_d(pool.buf[p], 2 * BB * NN);
  // check: SME and the three FFTW plans against each other (all rows), and FFTW against a direct DFT on row 0
  { double *b = pool.buf[0], *y[4]; void (*f[4])(void *, long) = {r_sme, r_row, r_many, r_many}; void *cx[4] = {0, 0, pm, pp};
    for (int k = 0; k < 4; k++) { y[k] = malloc(16 * BB * NN); memset(b + 2 * BB * NN, 0xff, 16 * BB * NN); f[k](cx[k], 0); memcpy(y[k], b + 2 * BB * NN, 16 * BB * NN); }
    static double ref[2 * 16384], wc[16384], ws[16384]; for (long j = 0; j < NN; j++) { wc[j] = __cospi(2.0 * j / NN); ws[j] = -__sinpi(2.0 * j / NN); } for (long kk = 0; kk < NN; kk++) { double sr = 0, si = 0; for (long j = 0; j < NN; j++) { long q = (j * kk) % NN; sr += b[2 * j] * wc[q] - b[2 * j + 1] * ws[q]; si += b[2 * j] * ws[q] + b[2 * j + 1] * wc[q]; } ref[2 * kk] = sr; ref[2 * kk + 1] = si; }
    const char *nm[4] = {"SME ap-0", "FFTW-d MEASURE row", "FFTW-d MEASURE many", "FFTW-d PATIENT many"};
    for (int k = 0; k < 4; k++) { double se = 0, sr = 0, e2 = 0, s2 = 0; for (long i = 0; i < 2 * NN; i++) { se += (y[k][i] - ref[i]) * (y[k][i] - ref[i]); sr += ref[i] * ref[i]; }
      for (long i = 0; i < 2 * BB * NN; i++) { e2 += (y[k][i] - y[1][i]) * (y[k][i] - y[1][i]); s2 += y[1][i] * y[1][i]; }
      printf("CHECK %-22s n=%ld B=%ld row0 relRMS vs direct DFT %.2e; all rows vs FFTW row plan %.2e %s\n", nm[k], NN, BB, sqrt(se / sr), sqrt(e2 / s2), sqrt(e2 / s2) < 1e-14 ? "PASS" : "FAIL"); } }
  bench_prog progs[4] = {{"SME fp64 (ap-0)", r_sme, NULL, NULL}, {"FFTW-d MEASURE row", r_row, NULL, NULL}, {"FFTW-d MEASURE many", r_many, NULL, pm}, {"FFTW-d PATIENT many", r_many, NULL, pp}};
  bench_opts op = BENCH_OPTS_DEFAULT; op.units = (double)(BB * NN); bench_stat st[4]; char title[96];
  snprintf(title, sizeof title, "smevsfftw n=%ld B=%ld %s", NN, BB, sep ? "separate" : "same");
  bench_compare(title, progs, 4, &op, st);
  return 0;
}
