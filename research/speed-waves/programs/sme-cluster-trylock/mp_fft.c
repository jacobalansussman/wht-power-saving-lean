/* mp_fft.c -- as-3: multi-process SME-owner study, kernel = ap-0's fp64 complex row FFTs (wave 13, smegi.h via
 * ap0lib.c, included unchanged). numpy layout (complex128 interleaved), forward, out of place.
 * NEON path = PFFFT-double ordered (ap-0's best baseline at 512; vDSP interleaved double ties).
 * One call = B rows of N points (env N=512, B=64). Check: every output value of the last call of each path against
 * a plain fp64 O(N^2) DFT (exact index reduction), max|err|/max|X| < 1e-13. usage: mp_fft W:POL[:E=n] ... */
#include <Accelerate/Accelerate.h>
#include "../timing-harness/bench.h"
#include "pffft_double.h"
#include "../sme-fp64-batched-fft/ap0lib.c"
#include <pthread.h>
#include "fftw3.h"
static int g_fftwd; static fftw_plan g_fp; /* NEONLIB=fftwd: FFTW 3.3.10 double built WITH NEON by ar-f (wave 13), PATIENT plan */
static long NN = 512, BB = 64;
void *kd_new(long n); void kd_exec(void *v, const double *in, double *out);
static int g_kfr; /* NEONLIB=kfr: KFR complex<double> (ar-0 found it the strongest fp64 NEON library at 512 in cache) */
typedef struct { double *x, *yN, *yS, *px, *py; PFFFTD_Setup *pf; double *work; ap0_t *H; void *kd; } fw;
static void *k_ctx[16];
static double K_UNITS;
static const char *K_UNIT = "point", *K_TITLE = "fft64";
static size_t K_ARENA_PER_WORKER;
static double *X0, *REF; static int g_copy_mode; static PFFFTD_Setup *g_pf;
static void *arena_get(size_t);
static void k_global(void) {
  if (getenv("N")) NN = atol(getenv("N"));
  if (getenv("B")) BB = atol(getenv("B"));
  K_UNITS = (double)NN * BB; K_ARENA_PER_WORKER = 2 * 16 * NN * BB + 4096;
  g_copy_mode = getenv("COPY") ? atoi(getenv("COPY")) : 0;
  X0 = bench_alloc(16 * NN * BB); bench_rng_state = 0x5151ull; bench_fill_d(X0, 2 * NN * BB);
  REF = bench_alloc(16 * NN * BB);
  double *cw = malloc(8 * NN), *sw = malloc(8 * NN);
  for (long m = 0; m < NN; m++) { cw[m] = cos(2 * M_PI * m / NN); sw[m] = sin(2 * M_PI * m / NN); }
  for (long t = 0; t < BB; t++) {
    const double *x = X0 + 2 * t * NN;
    for (long k = 0; k < NN; k++) {
      double re = 0, im = 0;
      for (long j = 0; j < NN; j++) { long m = (j * k) % NN; double c = cw[m], s = sw[m]; re += x[2 * j] * c + x[2 * j + 1] * s; im += x[2 * j + 1] * c - x[2 * j] * s; }
      REF[2 * (t * NN + k)] = re; REF[2 * (t * NN + k) + 1] = im;
    }
  }
  g_pf = pffftd_new_setup((int)NN, PFFFT_COMPLEX);
  g_kfr = getenv("NEONLIB") && !strcmp(getenv("NEONLIB"), "kfr");
  if (g_kfr) K_TITLE = "fft64kfr";
  g_fftwd = getenv("NEONLIB") && !strcmp(getenv("NEONLIB"), "fftwd");
  if (g_fftwd) {
    K_TITLE = "fft64fftwd";
    fftw_complex *a = fftw_malloc(16 * NN), *b = fftw_malloc(16 * NN);
    g_fp = fftw_plan_dft_1d((int)NN, a, b, FFTW_FORWARD, FFTW_PATIENT);
    if (!g_fp) { fprintf(stderr, "fftw plan failed\n"); exit(1); }
  }
}
static void *k_make(int w) {
  fw *f = calloc(1, sizeof *f);
  f->x = arena_get(16 * NN * BB); memcpy(f->x, X0, 16 * NN * BB);
  f->yS = arena_get(16 * NN * BB);
  f->yN = bench_alloc(16 * NN * BB); f->px = bench_alloc(16 * NN * BB); memcpy(f->px, X0, 16 * NN * BB); f->py = bench_alloc(16 * NN * BB);
  f->pf = g_pf; f->work = pffftd_aligned_malloc(16 * NN);
  f->H = ap0_new(NN);
  if (g_kfr) f->kd = kd_new(NN);
  if (getenv("BROKEN")) { if (f->H->R == 1) f->H->pl->Fr[2][3] = -f->H->pl->Fr[2][3]; else f->H->wi[3] = -f->H->wi[3]; }
  return f;
}
static void k_owner_init(void) {}
static void k_worker_init(void *c) { (void)c; }
static void k_neon(void *c) {
  fw *f = c; const double *x = g_copy_mode ? f->px : f->x;
  if (g_kfr) { for (long t = 0; t < BB; t++) kd_exec(f->kd, x + 2 * t * NN, f->yN + 2 * t * NN); return; }
  if (g_fftwd) { for (long t = 0; t < BB; t++) fftw_execute_dft(g_fp, (fftw_complex *)(x + 2 * t * NN), (fftw_complex *)(f->yN + 2 * t * NN)); return; }
  for (long t = 0; t < BB; t++) pffftd_transform_ordered(f->pf, x + 2 * t * NN, f->yN + 2 * t * NN, f->work, PFFFT_FORWARD);
}
static void k_sme_part(void *c) { fw *f = c; ap0_fft(f->H, f->x, f->yS, BB); }
static void k_tail(void *c) { (void)c; }
static void k_copy_in(void *c) { fw *f = c; memcpy(f->x, f->px, 16 * NN * BB); }
static void k_copy_out(void *c) { fw *f = c; memcpy(f->py, f->yS, 16 * NN * BB); }
static int k_check(void *c, int path, const char *lab) {
  fw *f = c; const double *y = path ? (g_copy_mode ? f->py : f->yS) : f->yN;
  double emax = 0, xmax = 0;
  for (long i = 0; i < 2 * NN * BB; i++) { double e = fabs(y[i] - REF[i]); if (!(e <= emax)) emax = e; if (fabs(REF[i]) > xmax) xmax = fabs(REF[i]); }
  int ok = emax / xmax < 1e-13;
  printf("CHECK %-10s n=%ld B=%ld all %ld values vs fp64 DFT: max|err|/max|X| %.3e %s\n", lab, NN, BB, 2 * NN * BB, emax / xmax, ok ? "PASS" : "FAIL");
  return ok;
}
#include "mpcore.h"
int main(int argc, char **argv) { setvbuf(stdout, NULL, _IOLBF, 0); bench_init(); return mp_main(argc, argv); }
