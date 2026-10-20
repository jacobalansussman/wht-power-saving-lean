// au1bench.c -- wave 14 critic au-1 of claim T2 (as-0's SME Welch). One whole fp64 Welch call, one P-core.
// Uses as-0's as0.c unchanged (read only, included), but links the TIMER-ENABLED FFTW 3.3.10 double NEON build of
// inf-strongest ($FFTW_DOUBLE_NEON) and imports its P-core PATIENT wisdom, then plans
// PATIENT in this process (run under the timing lock, so on a P-core). Adds the real-input routes a practitioner would write:
//   FFTW-d r2c per frame (PATIENT, timer build) and PFFFT-d real (ordered) per frame; same NEON prep and |X|^2 code style.
// usage: au1bench N NOVERLAP LOG2L [hann|hamming]
#include "../as0.c"
#include "../../timing-harness/bench.h"
static as0_t *H; static double *X, *PSD, *REF; static long L;
static fftw_plan FR; static double *RIN; static fftw_complex *ROUT; static double *ACCR;
static PFFFTD_Setup *PR; static double *PIN, *POUT, *PWORK;
static int SME_ENGINE;

static inline void prep_real(const double *a, const double *w, double *r, long N) {
  float64x2_t ma = vdupq_n_f64(fmean(a, N));
  for (long n = 0; n < N; n += 4) {
    vst1q_f64(r + n, vmulq_f64(vsubq_f64(vld1q_f64(a + n), ma), vld1q_f64(w + n)));
    vst1q_f64(r + n + 2, vmulq_f64(vsubq_f64(vld1q_f64(a + n + 2), ma), vld1q_f64(w + n + 2)));
  }
}
static void finish_r(long N, double scale, const double *acc, long nseg, double *psd) {
  double sc = scale / (double)nseg;
  for (long k = 0; k <= N / 2; k++) psd[k] = acc[k] * sc * ((k == 0 || k == N / 2) ? 1.0 : 2.0);
}
static void welch_fftw_r2c(const double *x, double *psd) {
  long N = H->N, st = H->step, nseg = as0_nseg(H, L); double *acc = ACCR; memset(acc, 0, (N / 2 + 2) * sizeof(double));
  double ny = 0;
  for (long s = 0; s < nseg; s++) {
    prep_real(x + s * st, H->win, RIN, N); fftw_execute(FR);
    acc_int(acc, (const double *)ROUT, N / 2);  // bins 0..N/2-1
    double a = ROUT[N / 2][0], b = ROUT[N / 2][1]; ny += a * a + b * b;
  }
  acc[N / 2] = ny; finish_r(N, H->scale, acc, nseg, psd);
}
static void welch_pffft_real(const double *x, double *psd) {
  long N = H->N, st = H->step, nseg = as0_nseg(H, L); double *acc = ACCR; memset(acc, 0, (N / 2 + 2) * sizeof(double));
  double dc = 0, ny = 0;
  for (long s = 0; s < nseg; s++) {
    prep_real(x + s * st, H->win, PIN, N);
    pffftd_transform_ordered(PR, PIN, POUT, PWORK, PFFFT_FORWARD);  // [X0, X_{N/2}, re1, im1, ...]
    dc += POUT[0] * POUT[0]; ny += POUT[1] * POUT[1];
    acc_int(acc, POUT, N / 2);  // acc[0] wrong (DC^2 + Ny^2), fixed below
  }
  acc[0] = dc; acc[N / 2] = ny; finish_r(N, H->scale, acc, nseg, psd);
}
static void r_sme(void *c, long i) { as0_welch_rows(H, X, 1, L, PSD, SME_ENGINE); }
static void r_fftw_c2c(void *c, long i) { as0_welch_rows(H, X, 1, L, PSD, 7); }
static void r_kfr(void *c, long i) { as0_welch_rows(H, X, 1, L, PSD, 5); }
static void r_pfd(void *c, long i) { as0_welch_rows(H, X, 1, L, PSD, 3); }
static void r_fftw_r2c(void *c, long i) { welch_fftw_r2c(X, PSD); }
static void r_pffft_real(void *c, long i) { welch_pffft_real(X, PSD); }

int main(int argc, char **argv) {
  bench_init();
  long N = atol(argv[1]), nov = atol(argv[2]); L = 1L << atol(argv[3]);
  const char *wn = argc > 4 ? argv[4] : "hann";
  double *w = malloc(N * sizeof(double));
  for (long n = 0; n < N; n++) { double c = cos(2 * M_PI * n / N); w[n] = strcmp(wn, "hamming") ? 0.5 - 0.5 * c : 0.54 - 0.46 * c; }
  int wis = fftw_import_wisdom_from_filename("infra/inf-strongest/wisdom_d_P.txt");
  printf("# P-core wisdom import: %s\n", wis ? "ok" : "FAILED");
  setenv("AS0_FFTW_PATIENT", "1", 1);
  uint64_t t0 = bench_now_ns();
  H = as0_new(N, nov, w, 1.0);
  RIN = fftw_malloc(sizeof(double) * (N + 4)); ROUT = fftw_malloc(sizeof(fftw_complex) * (N / 2 + 2));
  FR = fftw_plan_dft_r2c_1d((int)N, RIN, ROUT, FFTW_PATIENT);
  printf("# planning (c2c + r2c PATIENT, in process) %.2f s\n", (bench_now_ns() - t0) * 1e-9);
  printf("# FFTW c2c plan: "); fftw_print_plan(H->fw); printf("\n# FFTW r2c plan: "); fftw_print_plan(FR); printf("\n");
  ACCR = bench_alloc(sizeof(double) * (N / 2 + 8));
  PR = pffftd_new_setup((int)N, PFFFT_REAL); PIN = pffftd_aligned_malloc(8 * N); POUT = pffftd_aligned_malloc(8 * N); PWORK = pffftd_aligned_malloc(8 * N);
  X = bench_alloc(sizeof(double) * L); for (long j = 0; j < L; j++) X[j] = bench_rand_pm1() + 0.3 * sin(0.05 * (double)j) + 0.1;
  long F = N / 2 + 1; PSD = bench_alloc(sizeof(double) * F); REF = bench_alloc(sizeof(double) * F);
  int fused = (H->sme && ((ap0v2_t *)H->sme)->a && ((ap0v2_t *)H->sme)->a->R == 1);
  SME_ENGINE = fused ? 1 : 4;
  // exact reference: plain long-double DFT Welch on the first frames only would be slow; use vDSP pipeline (scipy-checked by check.py)
  as0_welch_rows(H, X, 1, L, REF, 2);
  bench_prog all[6] = {{"FFTW-d c2c pair (timer,PAT)", r_fftw_c2c, NULL, NULL}, {"FFTW-d r2c/frame (timer,PAT)", r_fftw_r2c, NULL, NULL},
                       {"PFFFT-d real/frame ordered", r_pffft_real, NULL, NULL}, {"KFR-d c2c pair", r_kfr, NULL, NULL},
                       {"PFFFT-d c2c pair z-order", r_pfd, NULL, NULL}, {fused ? "SME-F" : "SME-D", r_sme, NULL, NULL}};
  int np = 6;
  for (int k = 0; k < np; k++) { memset(PSD, 0, F * sizeof(double)); all[k].run(NULL, 0); bench_check_d(all[k].name, PSD, REF, F, 1e-13); }
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)L; o.unit_name = "input sample"; bench_stat st[6]; char title[160];
  snprintf(title, sizeof title, "au1 welch fp64 N=%ld nov=%ld %s L=2^%s nseg=%ld", N, nov, wn, argv[3], as0_nseg(H, L));
  bench_compare(title, all, np, &o, st);
  int bi = 0; for (int k = 1; k < np - 1; k++) if (st[k].median < st[bi].median) bi = k;
  for (int k = 0; k < np; k++) printf("MS %-30s best %8.3f median %8.3f ms/call\n", st[k].name, st[k].best * L * 1e-6, st[k].median * L * 1e-6);
  printf("VS SME vs best library (%s): best %+6.1f %%  median %+6.1f %%\n", st[bi].name,
         100 * (st[np - 1].best / st[bi].best - 1), 100 * (st[np - 1].median / st[bi].median - 1));
  return 0;
}
