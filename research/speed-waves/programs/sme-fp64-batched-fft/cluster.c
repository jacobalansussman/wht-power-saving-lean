// cluster.c -- ap-0: whole performance cluster (4 P-cores, one shared SME unit), numpy layout (complex128 rows).
// Configurations: 4 threads of a baseline vs 1 SME thread + 3 baseline threads (SME gets fraction f of the rows).
// Spinning pool (main thread = worker 0). One call = the whole array (B rows of N). Outputs checked against
// single-thread pocketfft. usage: cluster N B
#include <Accelerate/Accelerate.h>
#include <stdatomic.h>
#include "../timing-harness/bench.h"
#include "pffft_double.h"
#include "ap0lib.c"
void *pk_cd_new(int n); void pk_cd_exec(void *p, double *c);
static long NN, BB; static double *X, *Y; static ap0_t *H; static PFFFTD_Setup *pf; static void *pk; static vDSP_DFT_Interleaved_SetupD vi;
enum { L_PFFFT, L_VDSP, L_POCKET, L_SME };
typedef struct { int lib[4]; double f; int nthr; } cfg_t;  // lib per worker; f = SME fraction when worker 0 is SME
static cfg_t CUR; static _Atomic long gen = 0; static _Atomic int done = 0;
static void run_rows(int lib, long r0, long r1, double *work) {
  if (r1 <= r0) return;
  switch (lib) {
  case L_SME: ap0_fft(H, X + 2 * r0 * NN, Y + 2 * r0 * NN, r1 - r0); break;
  case L_PFFFT: for (long t = r0; t < r1; t++) pffftd_transform_ordered(pf, X + 2 * t * NN, Y + 2 * t * NN, work, PFFFT_FORWARD); break;
  case L_VDSP: for (long t = r0; t < r1; t++) vDSP_DFT_Interleaved_ExecuteD(vi, (const DSPDoubleComplex *)(X + 2 * t * NN), (DSPDoubleComplex *)(Y + 2 * t * NN)); break;
  case L_POCKET: memcpy(Y + 2 * r0 * NN, X + 2 * r0 * NN, 16 * (r1 - r0) * NN); for (long t = r0; t < r1; t++) pk_cd_exec(pk, Y + 2 * t * NN); break;
  }
}
static void range(int w, long *r0, long *r1) {
  int n = CUR.nthr;
  if (CUR.lib[0] == L_SME && n > 1) { long s = (long)(CUR.f * BB + 0.5); if (w == 0) { *r0 = 0; *r1 = s; return; }
    long rest = BB - s; *r0 = s + rest * (w - 1) / (n - 1); *r1 = s + rest * w / (n - 1); return; }
  *r0 = BB * w / n; *r1 = BB * (w + 1) / n;
}
static void *worker(void *arg) {
  int w = (int)(long)arg; pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
  double *work = pffftd_aligned_malloc(16 * NN); long seen = 0;
  for (;;) { long g; while ((g = atomic_load_explicit(&gen, memory_order_acquire)) == seen) ; seen = g; if (g < 0) break;
    if (w < CUR.nthr) { long r0, r1; range(w, &r0, &r1); run_rows(CUR.lib[w], r0, r1, work); }
    atomic_fetch_add_explicit(&done, 1, memory_order_acq_rel); }
  return NULL;
}
static double *work0;
static void run_cfg(void *c, long i) {
  CUR = *(cfg_t *)c; atomic_store(&done, 0); atomic_fetch_add_explicit(&gen, 1, memory_order_acq_rel);
  long r0, r1; range(0, &r0, &r1); run_rows(CUR.lib[0], r0, r1, work0);
  while (atomic_load_explicit(&done, memory_order_acquire) < 3) ;
}
int main(int argc, char **argv) {
  bench_init(); NN = atol(argv[1]); BB = atol(argv[2]);
  H = ap0_new(NN); pf = pffftd_new_setup((int)NN, PFFFT_COMPLEX); pk = pk_cd_new((int)NN); vi = vDSP_DFT_Interleaved_CreateSetupD(NULL, NN, vDSP_DFT_FORWARD, 0);
  X = bench_alloc(16 * BB * NN); Y = bench_alloc(16 * BB * NN); bench_fill_d(X, 2 * BB * NN); work0 = pffftd_aligned_malloc(16 * NN);
  pthread_t th[3]; for (long w = 1; w < 4; w++) pthread_create(&th[w - 1], NULL, worker, (void *)w);
  // reference: single-thread pocketfft
  double *ref = bench_alloc(16 * BB * NN); memcpy(ref, X, 16 * BB * NN); for (long t = 0; t < BB; t++) pk_cd_exec(pk, ref + 2 * t * NN);
  cfg_t cfgs[9] = {{{L_PFFFT, L_PFFFT, L_PFFFT, L_PFFFT}, 0, 4}, {{L_SME, L_PFFFT, L_PFFFT, L_PFFFT}, 0.40, 4},
                   {{L_SME, L_PFFFT, L_PFFFT, L_PFFFT}, 0.50, 4}, {{L_VDSP, L_VDSP, L_VDSP, L_VDSP}, 0, 4},
                   {{L_SME, L_VDSP, L_VDSP, L_VDSP}, 0.45, 4}, {{L_POCKET, L_POCKET, L_POCKET, L_POCKET}, 0, 4},
                   {{L_SME, L_PFFFT, L_PFFFT, L_PFFFT}, 1.0, 1}, {{L_PFFFT, L_PFFFT, L_PFFFT, L_PFFFT}, 0, 1},
                   {{L_SME, L_PFFFT, L_PFFFT, L_PFFFT}, 0.30, 4}};
  const char *names[9] = {"4 x PFFFT-d", "1 SME (40%) + 3 PFFFT-d", "1 SME (50%) + 3 PFFFT-d", "4 x vDSP interleaved D",
                          "1 SME (45%) + 3 vDSP-D", "4 x pocketfft", "1 SME thread alone", "1 PFFFT-d thread alone", "1 SME (30%) + 3 PFFFT-d"};
  bench_prog progs[9];
  for (int k = 0; k < 9; k++) {
    progs[k] = (bench_prog){names[k], run_cfg, NULL, &cfgs[k]};
    memset(Y, 0, 16 * BB * NN); run_cfg(&cfgs[k], 0);
    double emax = 0, xmax = 0; for (long j = 0; j < 2 * BB * NN; j++) { double e = fabs(Y[j] - ref[j]); if (!(e <= emax)) emax = e; if (fabs(ref[j]) > xmax) xmax = fabs(ref[j]); }
    printf("CHECK %-28s all %ld rows vs 1-thread pocketfft: max|diff|/max|X| %.3e %s\n", names[k], BB, emax / xmax, emax / xmax < 1e-14 ? "PASS" : "FAIL");
  }
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)(BB * NN); bench_stat st[9]; char title[128];
  snprintf(title, sizeof title, "cluster c2c fp64 interleaved n=%ld rows=%ld (%.0f MiB in)", NN, BB, 16.0 * BB * NN / 1048576);
  bench_compare(title, progs, 9, &o, st);
  atomic_store(&gen, -1); return 0;
}
