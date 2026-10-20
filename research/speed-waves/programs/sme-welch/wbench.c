// wbench.c -- wave 14 as-0: time one whole float64 Welch call (C channels x L samples) with each engine.
// usage: wbench N NOVERLAP LOG2L C [window: hann|hamming]   (all programs take turns, bench.h, ns per input sample)
#include "as0.c"
#include "../timing-harness/bench.h"
static as0_t *H; static double *X, *PSD, *REF; static long C, L;
static void r_smeA(void *c, long i) { as0_welch_rows(H, X, C, L, PSD, 0); }
static void r_smeF(void *c, long i) { as0_welch_rows(H, X, C, L, PSD, 1); }
static void r_vdsp(void *c, long i) { as0_welch_rows(H, X, C, L, PSD, 2); }
static void r_pfd(void *c, long i) { as0_welch_rows(H, X, C, L, PSD, 3); }
static void r_smeD(void *c, long i) { as0_welch_rows(H, X, C, L, PSD, 4); }
static void r_kfr(void *c, long i) { as0_welch_rows(H, X, C, L, PSD, 5); }
static void r_smeW(void *c, long i) { as0_welch_rows(H, X, C, L, PSD, 6); }
static void r_fftw(void *c, long i) { as0_welch_rows(H, X, C, L, PSD, 7); }
int main(int argc, char **argv) {
  bench_init();
  long N = atol(argv[1]), nov = atol(argv[2]); L = 1L << atol(argv[3]); C = atol(argv[4]);
  const char *wn = argc > 5 ? argv[5] : "hann";
  double *w = malloc(N * sizeof(double));
  for (long n = 0; n < N; n++) { double c = cos(2 * M_PI * n / N); w[n] = strcmp(wn, "hamming") ? 0.5 - 0.5 * c : 0.54 - 0.46 * c; }  // periodic
  H = as0_new(N, nov, w, 1.0);
  X = bench_alloc(sizeof(double) * C * L); for (long j = 0; j < C * L; j++) X[j] = bench_rand_pm1() + 0.3 * sin(0.05 * (double)(j % L)) + 0.1;
  long F = N / 2 + 1; PSD = bench_alloc(sizeof(double) * C * F); REF = bench_alloc(sizeof(double) * C * F);
  int fused = (H->sme && ((ap0v2_t *)H->sme)->a && ((ap0v2_t *)H->sme)->a->R == 1);
  as0_welch_rows(H, X, C, L, REF, 2);  // vDSP as the in-C reference (check.py checks every engine against scipy)
  bench_prog all[5] = {{"vDSP-D split pipeline", r_vdsp, NULL, NULL}, {"PFFFT-d z-order pipeline", r_pfd, NULL, NULL},
                       {"KFR-d pipeline", r_kfr, NULL, NULL}, {"SME-A (NEON acc)", r_smeA, NULL, NULL}, {"SME-F (acc in last pass)", r_smeF, NULL, NULL}};
  if (!fused && H->plM) all[4] = (bench_prog){"SME-D (DIF prep, acc in pass)", r_smeD, NULL, NULL};
  if (getenv("WITH_FFTW")) all[1] = (bench_prog){"FFTW-d NEON pipeline", r_fftw, NULL, NULL};
  if (fused && getenv("WITH_W")) all[3] = (bench_prog){"SME-W (window in pass 1)", r_smeW, NULL, NULL};
  int np = (fused || H->plM) ? 5 : 4;
  if (getenv("ONLY_SME")) { all[0] = all[np - 1]; np = 1; }
  for (int k = (getenv("ONLY_SME") ? 0 : 1); k < np; k++) { all[k].run(NULL, 0); bench_check_d(all[k].name, PSD, REF, C * F, 1e-13); }
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)(C * L); o.unit_name = "input sample"; bench_stat st[5]; char title[160];
  long nseg = as0_nseg(H, L);
  snprintf(title, sizeof title, "welch fp64 N=%ld nov=%ld %s C=%ld L=2^%s (%.0f MiB) nseg=%ld", N, nov, wn, C, argv[3], 8.0 * C * L / 1048576, nseg);
  bench_compare(title, all, np, &o, st);
  for (int k = 0; k < np; k++) printf("MS %-28s best %8.3f median %8.3f ms/call\n", st[k].name, st[k].best * C * L * 1e-6, st[k].median * C * L * 1e-6);
  for (int k = 3; k < np; k++) printf("VS %-28s vs KFR-d: best %+6.1f %%  median %+6.1f %%\n", st[k].name, 100 * (st[k].best / st[2].best - 1), 100 * (st[k].median / st[2].median - 1));
  for (int k = 1; k < np; k++) printf("VS %-28s vs vDSP: best %+6.1f %%  median %+6.1f %%\n", st[k].name, 100 * (st[k].best / st[0].best - 1), 100 * (st[k].median / st[0].median - 1));
  return 0;
}
