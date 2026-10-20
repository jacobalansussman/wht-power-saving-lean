/* olsd.c -- av-1 (wave 15): kernel-level (C, no Python) fp64 comparison of the SME Toeplitz-tile direct FIR against
 * the strongest fp64 FFT route I could build in the time: overlap-save with PFFFT double (NEON), real transforms in
 * z-order, filter spectrum pre-scaled, block size N scanned per tap count. "valid" correlation form:
 *   x has n + L - 1 samples, y[m] = sum_k g[k] x[m + k], m < n.
 * Build: clang -O3 -mcpu=apple-m4 -I$B/pffft/include olsd.c fc64.c $B/pffft/lib/libpffft.a -framework Accelerate -o olsd
 *   (B=$BASELINES)
 * Use: olsd check L LOGN | olsd scan L LOGN | olsd time L LOGN N_ols [progs bitmask: 1 sme, 2 sme_b8, 4 ols, 8 vdsp]
 */
#include "../timing-harness/bench.h"
#include "pffft_double.h"

int fc_corr_d(const double *x, long nx, const double *g, int L, long o0, long nout, double *y, int path);
int fc_vdsp_d(const double *x, long nx, const double *g, int L, long o0, long nout, double *y);
extern int fc_nancheck;

typedef struct {
  const double *x, *g;
  double *y;
  long n;
  int L, N, path, nc;
  PFFFTD_Setup *pf;
  double *X, *Y, *H, *w, *t, *seg;
} ctx_t;

static void ols_setup(ctx_t *c, int N) {
  c->N = N;
  c->pf = pffftd_new_setup(N, PFFFT_REAL);
  if (!c->pf) { fprintf(stderr, "pffftd setup failed N=%d\n", N); exit(1); }
  c->X = bench_alloc(8 * (size_t)N); c->Y = bench_alloc(8 * (size_t)N); c->H = bench_alloc(8 * (size_t)N);
  c->w = bench_alloc(16 * (size_t)N); c->t = bench_alloc(16 * (size_t)N); c->seg = bench_alloc(8 * (size_t)N);
  double *hp = bench_alloc(8 * (size_t)N);
  for (int k = 0; k < c->L; k++) hp[k] = c->g[c->L - 1 - k] / (double)N; /* h = g reversed */
  pffftd_transform(c->pf, hp, c->H, c->w, PFFFT_FORWARD);
  free(hp);
}
static void run_ols(void *vc, long i) {
  ctx_t *c = vc;
  int N = c->N, L = c->L;
  long n = c->n, xlen = n + L - 1, step = N - L + 1, m = 0;
  const double *x = c->x;
  double *y = c->y;
  for (; m + step <= n && m + N <= xlen; m += step) {
    pffftd_transform(c->pf, x + m, c->X, c->w, PFFFT_FORWARD);
    pffftd_zconvolve(c->pf, c->X, c->H, c->Y);
    pffftd_transform(c->pf, c->Y, c->t, c->w, PFFFT_BACKWARD);
    memcpy(y + m, c->t + L - 1, 8 * step);
  }
  while (m < n) {
    long avail = xlen - m, cp = avail < N ? avail : N;
    memcpy(c->seg, x + m, 8 * cp); memset(c->seg + cp, 0, 8 * (N - cp));
    long k = n - m < step ? n - m : step;
    pffftd_transform(c->pf, c->seg, c->X, c->w, PFFFT_FORWARD);
    pffftd_zconvolve(c->pf, c->X, c->H, c->Y);
    pffftd_transform(c->pf, c->Y, c->t, c->w, PFFFT_BACKWARD);
    memcpy(y + m, c->t + L - 1, 8 * k);
    m += k;
  }
}
static void run_sme(void *vc, long i) {
  ctx_t *c = vc;
  fc_nancheck = c->nc;
  fc_corr_d(c->x, c->n + c->L - 1, c->g, c->L, c->L - 1, c->n, c->y, c->path);
}
static void run_vdsp(void *vc, long i) {
  ctx_t *c = vc;
  fc_vdsp_d(c->x, c->n + c->L - 1, c->g, c->L, c->L - 1, c->n, c->y);
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: olsd check|scan|time L LOGN [N_ols] [progs]\n"); return 1; }
  bench_init();
  int L = atoi(argv[2]), logn = atoi(argv[3]);
  long n = 1L << logn;
  double *x = bench_alloc(8 * (size_t)(n + L - 1)), *g = bench_alloc(8 * (size_t)L), *y = bench_alloc(8 * (size_t)n);
  bench_fill_d(x, n + L - 1);
  bench_fill_d(g, L);
  ctx_t base = {x, g, y, n, L, 0, 1, 3, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  if (!strcmp(argv[1], "check")) {
    int Nc = 1024; while (Nc < 2 * L) Nc *= 2; Nc *= 2;
    ctx_t c = base; ols_setup(&c, Nc);
    long ns = 2000;
    double *ref = malloc(8 * ns), *got = malloc(8 * ns);
    long *pos = malloc(sizeof(long) * ns);
    for (long s = 0; s < ns; s++) pos[s] = s < 500 ? s : (s < 1000 ? n - 1 - (s - 500) : (long)(bench_rand_u64() % n));
    for (long s = 0; s < ns; s++) { double a = 0; for (int k = 0; k < L; k++) a = __builtin_fma(g[k], x[pos[s] + k], a); ref[s] = a; }
    const char *nm[4] = {"sme", "sme_b8", "ols_pffftd", "vdsp_convD"};
    for (int p = 0; p < 4; p++) {
      for (long m = 0; m < n; m++) y[m] = NAN;
      if (p == 0) { c.path = 1; run_sme(&c, 0); } else if (p == 1) { c.path = 5; run_sme(&c, 0); }
      else if (p == 2) run_ols(&c, 0); else run_vdsp(&c, 0);
      for (long s = 0; s < ns; s++) got[s] = y[pos[s]];
      char lab[64]; snprintf(lab, 64, "%s L=%d n=2^%d", nm[p], L, logn);
      bench_check_d(lab, got, ref, ns, 1e-12);
    }
    return 0;
  }
  if (!strcmp(argv[1], "scan")) { /* time OLS block sizes quickly (3 repeats of 0.2 s; scan only) */
    for (int N = 1024; N <= 65536; N *= 2) {
      if (N < 2 * L) continue;
      ctx_t c = base; ols_setup(&c, N);
      bench_prog pr[1] = {{"ols", run_ols, NULL, &c}};
      bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)n; o.unit_name = "output"; o.print_reps = 0;
      bench_stat st[1]; char t[64]; snprintf(t, 64, "olsscan L=%d n=2^%d N=%d", L, logn, N);
      bench_compare(t, pr, 1, &o, st);
    }
    return 0;
  }
  int Nols = argc > 4 ? atoi(argv[4]) : 4096, mask = argc > 5 ? atoi(argv[5]) : 15;
  ctx_t c1 = base, c2 = base, c3 = base, c4 = base, c5 = base, c6 = base;
  c2.path = 5;
  c5.nc = 0; c6.nc = 1; /* sme with no non-finite check / with the output scan on the core (c1, c2: check inside the streaming transpose, the default) */
  ols_setup(&c3, Nols);
  bench_prog all[6] = {{"sme", run_sme, NULL, &c1}, {"sme_b8", run_sme, NULL, &c2}, {"ols_pffftd", run_ols, NULL, &c3},
                       {"vdsp_convD", run_vdsp, NULL, &c4}, {"sme_noscan", run_sme, NULL, &c5}, {"sme_outscan", run_sme, NULL, &c6}};
  bench_prog pr[6]; int np = 0;
  for (int p = 0; p < 6; p++) if (mask & (1 << p)) pr[np++] = all[p];
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)n; o.unit_name = "output";
  bench_stat st[6]; char t[96]; snprintf(t, 96, "valid fp64 L=%d n=2^%d same-array (OLS N=%d)", L, logn, Nols);
  bench_compare(t, pr, np, &o, st);
  return 0;
}
