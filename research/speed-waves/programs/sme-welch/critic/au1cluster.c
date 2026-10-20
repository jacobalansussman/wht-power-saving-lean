// au1cluster.c -- critic au-1: as-0's 4-P-core Welch (cluster.c logic, as0.c unchanged) re-run with the TIMER FFTW-d
// (inf-strongest build, P-core wisdom, PATIENT planned in process under the timing lock) and the real-input routes per frame.
// usage: au1cluster N NOVERLAP LOG2L C [hann|hamming]
#include "../as0.c"
#include "../../timing-harness/bench.h"
#include <pthread.h>
#include <stdatomic.h>
static double *X, *PSD, *REF; static long C, L, F;
static as0_t *HS[4]; static int ENG[4];
static fftw_plan FR[4]; static double *RIN[4]; static fftw_complex *ROUT[4]; static double *ACCR[4];
static PFFFTD_Setup *PR; static double *PIN[4], *POUT[4], *PWORK[4];
static inline void prep_real(const double *a, const double *w, double *r, long N) {
  float64x2_t ma = vdupq_n_f64(fmean(a, N));
  for (long n = 0; n < N; n += 2) vst1q_f64(r + n, vmulq_f64(vsubq_f64(vld1q_f64(a + n), ma), vld1q_f64(w + n)));
}
static void fin_r(const as0_t *h, const double *acc, long nseg, double *psd) {
  long N = h->N; double sc = h->scale / (double)nseg;
  for (long k = 0; k <= N / 2; k++) psd[k] = acc[k] * sc * ((k == 0 || k == N / 2) ? 1.0 : 2.0);
}
static void welch_r2c(int w, const double *x, double *psd) {
  as0_t *h = HS[w]; long N = h->N, st = h->step, nseg = as0_nseg(h, L); double *acc = ACCR[w], ny = 0;
  memset(acc, 0, (N / 2 + 2) * sizeof(double));
  for (long s = 0; s < nseg; s++) {
    prep_real(x + s * st, h->win, RIN[w], N); fftw_execute(FR[w]);
    acc_int(acc, (const double *)ROUT[w], N / 2); double a = ROUT[w][N / 2][0], b = ROUT[w][N / 2][1]; ny += a * a + b * b;
  }
  acc[N / 2] = ny; fin_r(h, acc, nseg, psd);
}
static void welch_pr(int w, const double *x, double *psd) {
  as0_t *h = HS[w]; long N = h->N, st = h->step, nseg = as0_nseg(h, L); double *acc = ACCR[w], dc = 0, ny = 0;
  memset(acc, 0, (N / 2 + 2) * sizeof(double));
  for (long s = 0; s < nseg; s++) {
    prep_real(x + s * st, h->win, PIN[w], N); pffftd_transform_ordered(PR, PIN[w], POUT[w], PWORK[w], PFFFT_FORWARD);
    dc += POUT[w][0] * POUT[w][0]; ny += POUT[w][1] * POUT[w][1]; acc_int(acc, POUT[w], N / 2);
  }
  acc[0] = dc; acc[N / 2] = ny; fin_r(h, acc, nseg, psd);
}
static atomic_long nextc; static atomic_int gen, done;
static void work(int w) {
  long c;
  while ((c = atomic_fetch_add_explicit(&nextc, 1, memory_order_relaxed)) < C) {
    if (ENG[w] == 8) welch_r2c(w, X + c * L, PSD + c * F);
    else if (ENG[w] == 9) welch_pr(w, X + c * L, PSD + c * F);
    else as0_welch_rows(HS[w], X + c * L, 1, L, PSD + c * F, ENG[w]);
  }
}
static void *worker(void *arg) {
  int w = (int)(long)arg; pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
  int seen = 0;
  for (;;) {
    int g; while ((g = atomic_load_explicit(&gen, memory_order_acquire)) == seen) __asm__ volatile("yield");
    if (g < 0) return NULL; seen = g; work(w); atomic_fetch_add_explicit(&done, 1, memory_order_acq_rel);
  }
}
typedef struct { int eng[4]; } cfg_t;
static void run_cfg(void *c, long i) {
  cfg_t *cf = c; for (int w = 0; w < 4; w++) ENG[w] = cf->eng[w];
  atomic_store(&nextc, 0); atomic_store(&done, 0); atomic_fetch_add_explicit(&gen, 1, memory_order_acq_rel);
  work(0); while (atomic_load_explicit(&done, memory_order_acquire) < 3) __asm__ volatile("yield");
}
int main(int argc, char **argv) {
  bench_init();
  long N = atol(argv[1]), nov = atol(argv[2]); L = 1L << atol(argv[3]); C = atol(argv[4]);
  const char *wn = argc > 5 ? argv[5] : "hann"; F = N / 2 + 1;
  double *w = malloc(N * sizeof(double));
  for (long n = 0; n < N; n++) { double cc = cos(2 * M_PI * n / N); w[n] = strcmp(wn, "hamming") ? 0.5 - 0.5 * cc : 0.54 - 0.46 * cc; }
  int wis = fftw_import_wisdom_from_filename("infra/inf-strongest/wisdom_d_P.txt");
  printf("# P-core wisdom import: %s\n", wis ? "ok" : "FAILED");
  setenv("AS0_FFTW_PATIENT", "1", 1);
  for (int k = 0; k < 4; k++) {
    HS[k] = as0_new(N, nov, w, 1.0);
    RIN[k] = fftw_malloc(sizeof(double) * (N + 4)); ROUT[k] = fftw_malloc(sizeof(fftw_complex) * (N / 2 + 2));
    FR[k] = fftw_plan_dft_r2c_1d((int)N, RIN[k], ROUT[k], FFTW_PATIENT); ACCR[k] = bench_alloc(sizeof(double) * (N / 2 + 8));
    PIN[k] = pffftd_aligned_malloc(8 * N); POUT[k] = pffftd_aligned_malloc(8 * N); PWORK[k] = pffftd_aligned_malloc(8 * N);
  }
  PR = pffftd_new_setup((int)N, PFFFT_REAL);
  printf("# FFTW c2c plan: "); fftw_print_plan(HS[0]->fw); printf("\n# FFTW r2c plan: "); fftw_print_plan(FR[0]); printf("\n");
  int S = (HS[0]->sme && ((ap0v2_t *)HS[0]->sme)->a && ((ap0v2_t *)HS[0]->sme)->a->R == 1) ? 1 : (HS[0]->plM ? 4 : 0);
  X = bench_alloc(sizeof(double) * C * L); for (long j = 0; j < C * L; j++) X[j] = bench_rand_pm1() + 0.3 * sin(0.05 * (double)(j % L)) + 0.1;
  PSD = bench_alloc(sizeof(double) * C * F); REF = bench_alloc(sizeof(double) * C * F);
  as0_welch_rows(HS[0], X, C, L, REF, 2);
  pthread_t th[3]; for (long k = 1; k < 4; k++) pthread_create(&th[k - 1], NULL, worker, (void *)k);
  cfg_t cfgs[7] = {{{7, 7, 7, 7}}, {{8, 8, 8, 8}}, {{9, 9, 9, 9}}, {{5, 5, 5, 5}}, {{S, S, S, S}}, {{S, S, 7, 7}}, {{S, 7, 7, 7}}};
  const char *names[7] = {"4 x FFTW-d c2c pair (timer)", "4 x FFTW-d r2c/frame", "4 x PFFFT-d real/frame", "4 x KFR-d c2c pair",
                          "4 SME threads", "2 SME + 2 FFTW-d c2c", "1 SME + 3 FFTW-d c2c"};
  int NC = 7; bench_prog progs[7];
  for (int k = 0; k < NC; k++) {
    progs[k] = (bench_prog){names[k], run_cfg, NULL, &cfgs[k]};
    memset(PSD, 0, sizeof(double) * C * F); run_cfg(&cfgs[k], 0);
    char lab[96]; snprintf(lab, sizeof lab, "%s (all %ld channels)", names[k], C); bench_check_d(lab, PSD, REF, C * F, 1e-13);
  }
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)(C * L); o.unit_name = "input sample"; bench_stat st[7]; char title[160];
  snprintf(title, sizeof title, "au1 cluster welch fp64 N=%ld nov=%ld %s C=%ld L=2^%s", N, nov, wn, C, argv[3]);
  bench_compare(title, progs, NC, &o, st);
  int bl = 0; for (int k = 1; k < 4; k++) if (st[k].median < st[bl].median) bl = k;
  int bs = 4; for (int k = 5; k < NC; k++) if (st[k].median < st[bs].median) bs = k;
  for (int k = 0; k < NC; k++) printf("MS %-30s best %9.3f median %9.3f ms/call\n", names[k], st[k].best * C * L * 1e-6, st[k].median * C * L * 1e-6);
  printf("VS best SME mix (%s) vs best 4-thread library (%s): median %+6.1f %%  best %+6.1f %%\n", names[bs], names[bl],
         100 * (st[bs].median / st[bl].median - 1), 100 * (st[bs].best / st[bl].best - 1));
  atomic_store(&gen, -1); return 0;
}
