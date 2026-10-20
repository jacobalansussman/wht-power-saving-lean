// au-h0 probe: FFTW 3.3.10 as shipped by conda-forge for osx-arm64 (no NEON codelets, found by strings)
// against FFTW 3.3.10 built with NEON (our builds), same API, same plan flag (MEASURE), one process, taking turns.
// The conda library is dlopen'ed (RTLD_LOCAL); the NEON builds are linked statically.
#include <dlfcn.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

static double now(void) { return clock_gettime_nsec_np(CLOCK_UPTIME_RAW) * 1e-9; }

typedef void *(*pd_t)(int, void *, void *, int, unsigned);
typedef void *(*pr_t)(int, void *, void *, unsigned);
typedef void (*ex_t)(void *);

static double timeit(ex_t ex, void *p, int n) {
  // loop of >= 0.2 s
  long it = 1; double t;
  for (;;) { double t0 = now(); for (long i = 0; i < it; i++) ex(p); t = now() - t0; if (t > 0.2) break; it *= 2; }
  return t / it / n * 1e9; // ns per point
}
static int cmpd(const void *a, const void *b) { double x = *(double *)a, y = *(double *)b; return (x > y) - (x < y); }

int main(int argc, char **argv) {
  const char *lib = argv[1]; // conda libfftw3.dylib or libfftw3f.dylib
  int fl = strstr(lib, "fftw3f") != NULL;
  const char *kind = argv[2]; int n = atoi(argv[3]);
  void *h = dlopen(lib, RTLD_NOW | RTLD_LOCAL);
  if (!h) { printf("dlopen %s\n", dlerror()); return 1; }
  const char *pfx = fl ? "fftwf_" : "fftw_";
  char nm[64];
  snprintf(nm, 64, "%splan_dft_1d", pfx); pd_t cpd = (pd_t)dlsym(h, nm);
  snprintf(nm, 64, "%splan_dft_r2c_1d", pfx); pr_t cpr = (pr_t)dlsym(h, nm);
  snprintf(nm, 64, "%sexecute", pfx); ex_t cex = (ex_t)dlsym(h, nm);
  int r2c = !strcmp(kind, "r2c");
  size_t es = fl ? sizeof(float) : sizeof(double);
  void *in = calloc(2 * n + 16, es), *out = calloc(2 * n + 16, es), *in2 = calloc(2 * n + 16, es), *out2 = calloc(2 * n + 16, es);
  void *pc, *pn; ex_t nex;
  if (fl) {
    pn = r2c ? (void *)fftwf_plan_dft_r2c_1d(n, in2, out2, FFTW_MEASURE) : (void *)fftwf_plan_dft_1d(n, in2, out2, FFTW_FORWARD, FFTW_MEASURE);
    nex = (ex_t)fftwf_execute;
  } else {
    pn = r2c ? (void *)fftw_plan_dft_r2c_1d(n, in2, out2, FFTW_MEASURE) : (void *)fftw_plan_dft_1d(n, in2, out2, FFTW_FORWARD, FFTW_MEASURE);
    nex = (ex_t)fftw_execute;
  }
  pc = r2c ? cpr(n, in, out, FFTW_MEASURE) : cpd(n, in, out, FFTW_FORWARD, FFTW_MEASURE);
  // same random input in both, check outputs agree
  srand(1); int m = r2c ? n : 2 * n;
  for (int i = 0; i < m; i++) { double v = rand() / (double)RAND_MAX - 0.5; if (fl) { ((float *)in)[i] = v; ((float *)in2)[i] = v; } else { ((double *)in)[i] = v; ((double *)in2)[i] = v; } }
  cex(pc); nex(pn);
  int mo = r2c ? 2 * (n / 2 + 1) : 2 * n; double md = 0, mx = 0;
  for (int i = 0; i < mo; i++) { double a = fl ? ((float *)out)[i] : ((double *)out)[i], b = fl ? ((float *)out2)[i] : ((double *)out2)[i]; md = fmax(md, fabs(a - b)); mx = fmax(mx, fabs(b)); }
  printf("%s %s n=%d  max|conda-neon|/max = %.2e\n", fl ? "fp32" : "fp64", kind, n, md / mx);
  double tc[9], tn[9];
  for (int r = 0; r < 9; r++) { tc[r] = timeit(cex, pc, n); tn[r] = timeit(nex, pn, n); printf("  rep %d conda %.4f  neon %.4f ns/pt\n", r, tc[r], tn[r]); }
  qsort(tc, 9, 8, cmpd); qsort(tn, 9, 8, cmpd);
  printf("RESULT %s %s n=%d conda(no NEON) best %.4f med %.4f | NEON build best %.4f med %.4f | conda/neon best %.2fx\n", fl ? "fp32" : "fp64", kind, n, tc[0], tc[4], tn[0], tn[4], tc[0] / tn[0]);
  return 0;
}
