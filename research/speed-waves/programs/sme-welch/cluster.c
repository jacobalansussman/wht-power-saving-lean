// cluster.c -- wave 14 as-0: multi-channel float64 Welch (C channels x L samples) on the 4 performance cores.
// 4 threads (main + 3 spinning workers); channels handed out dynamically (atomic counter), so faster threads take more
// channels (no share tuning). Each thread has its own as0 handle (own scratch). Engines per thread slot from a config.
// usage: cluster N NOVERLAP LOG2L C [hann|hamming]
#include "as0.c"
#include "../timing-harness/bench.h"
#include <pthread.h>
#include <stdatomic.h>
static double *X, *PSD, *REF; static long C, L, F;
static as0_t *HS[4]; static int ENG[4];
static atomic_long nextc; static atomic_int gen, done;
static void work(int w) { long c; while ((c = atomic_fetch_add_explicit(&nextc, 1, memory_order_relaxed)) < C) as0_welch_rows(HS[w], X + c * L, 1, L, PSD + c * F, ENG[w]); }
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
  for (int k = 0; k < 4; k++) HS[k] = as0_new(N, nov, w, 1.0);
  int S = (HS[0]->sme && ((ap0v2_t *)HS[0]->sme)->a && ((ap0v2_t *)HS[0]->sme)->a->R == 1) ? 1 : (HS[0]->plM ? 4 : 0);
  X = bench_alloc(sizeof(double) * C * L); for (long j = 0; j < C * L; j++) X[j] = bench_rand_pm1() + 0.3 * sin(0.05 * (double)(j % L)) + 0.1;
  PSD = bench_alloc(sizeof(double) * C * F); REF = bench_alloc(sizeof(double) * C * F);
  as0_welch_rows(HS[0], X, C, L, REF, 2);  // single-thread vDSP reference (every engine is checked against scipy in check.py)
  pthread_t th[3]; for (long k = 1; k < 4; k++) pthread_create(&th[k - 1], NULL, worker, (void *)k);
  cfg_t cfgs[7] = {{{2, 2, 2, 2}}, {{5, 5, 5, 5}}, {{S, 5, 5, 5}}, {{S, S, 5, 5}}, {{S, S, S, S}}, {{S, 2, 2, 2}}, {{3, 3, 3, 3}}};
  const char *names[7] = {"4 x vDSP-D pipeline", "4 x KFR-d pipeline", "1 SME + 3 KFR-d", "2 SME + 2 KFR-d", "4 SME threads", "1 SME + 3 vDSP-D", "4 x PFFFT-d pipeline"};
  if (getenv("WITH_FFTW")) {  // second set: FFTW-d NEON (ar-f's build) as the NEON library
    cfg_t c2[7] = {{{2, 2, 2, 2}}, {{7, 7, 7, 7}}, {{S, 7, 7, 7}}, {{S, S, 7, 7}}, {{S, S, S, S}}, {{5, 5, 5, 5}}, {{S, 5, 5, 5}}};
    const char *n2[7] = {"4 x vDSP-D pipeline", "4 x FFTW-d NEON pipeline", "1 SME + 3 FFTW-d", "2 SME + 2 FFTW-d", "4 SME threads", "4 x KFR-d pipeline", "1 SME + 3 KFR-d"};
    memcpy(cfgs, c2, sizeof cfgs); memcpy(names, n2, sizeof names);
  }
  int NC = getenv("NCFG") ? atoi(getenv("NCFG")) : 7;
  bench_prog progs[7];
  for (int k = 0; k < NC; k++) {
    progs[k] = (bench_prog){names[k], run_cfg, NULL, &cfgs[k]};
    memset(PSD, 0, sizeof(double) * C * F); run_cfg(&cfgs[k], 0);
    char lab[96]; snprintf(lab, sizeof lab, "%s (all %ld channels)", names[k], C); bench_check_d(lab, PSD, REF, C * F, 1e-13);
  }
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)(C * L); o.unit_name = "input sample"; bench_stat st[7]; char title[160];
  snprintf(title, sizeof title, "cluster welch fp64 N=%ld nov=%ld %s C=%ld L=2^%s (%.0f MiB)", N, nov, wn, C, argv[3], 8.0 * C * L / 1048576);
  bench_compare(title, progs, NC, &o, st);
  for (int k = 0; k < NC; k++) printf("MS %-24s best %9.3f median %9.3f ms/call   vs 4 x vDSP median %+6.1f %%  vs slot-1 library x 4 median %+6.1f %%\n", names[k], st[k].best * C * L * 1e-6, st[k].median * C * L * 1e-6, 100 * (st[k].median / st[0].median - 1), 100 * (st[k].median / st[1].median - 1));
  atomic_store(&gen, -1); return 0;
}
