/* alt.c -- critic al-1 (wave 11): stronger / practitioner baselines for claim C2 (SME2 Toeplitz FIR).
 * Reuses aj-1's fir.c unchanged (included with main renamed) so the SME kernel and the baselines are byte-identical.
 * Adds:  vdsp_conv_pos  (vDSP_conv in correlation form, positive filter stride, pre-reversed taps)
 *        vdsp_conv_mt   (NT threads, GCD dispatch_apply, user-interactive QoS, one slice of outputs each)
 *        pffft_ols_mt   (NT threads of aj-1's PFFFT overlap-save, own work buffers per thread)
 *        ols2_vdsp_mt   (NT threads of aj-1's vDSP overlap-save)
 *        sme_mt         (NTS threads of aj-1's SME kernel: does the shared unit scale?)
 *        blas_toep      (Accelerate cblas_sgemm on block-Toeplitz panels: Y += X_j * T_j, no copies)
 * Build: clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$B/pffft/include alt.c $B/pffft/lib/libpffft.a -framework Accelerate -o alt
 * Use:   alt check L LOGN ; alt time L LOGN same|separate MASK
 * env: NT (default 4), NTS (default 2), BLAS_B (default 64), BLAS_CH (default 16384), OLS_N, OLS_PF_N as in fir.c
 */
#define main fir_main
#include "../fir.c"
#undef main
#include <dispatch/dispatch.h>

#define MAXT 8
static int NT = 4, NTS = 2, BB = 64; static long BCH = 16384;
static ctx_t CT[MAXT], CS[MAXT];
static bench_pool XPT[MAXT], YPT[MAXT], XPS[MAXT], YPS[MAXT];
static float *GREV; /* g[k] = h[L-1-k] */
static float *TJ;   /* BLAS Toeplitz panels: J panels of B x B, row-major [d][i] */
static int JJ;

static void run_vdsp_pos(void *vc, long i) {
  ctx_t *c = vc;
  const float *x = bench_pool_get(c->xp, i);
  float *y = bench_pool_get(c->yp, i);
  vDSP_conv(x, 1, GREV, 1, y, 1, c->n, c->L);
}
static int VPOS = 0;
static void run_vdsp_mt(void *vc, long i) {
  (void)vc;
  dispatch_apply(NT, dispatch_get_global_queue(QOS_CLASS_USER_INTERACTIVE, 0), ^(size_t t) {
    if (VPOS) run_vdsp_pos(&CT[t], i); else run_vdsp(&CT[t], i);
  });
}
static void run_pf_mt(void *vc, long i) {
  (void)vc;
  dispatch_apply(NT, dispatch_get_global_queue(QOS_CLASS_USER_INTERACTIVE, 0), ^(size_t t) { run_ols_pf(&CT[t], i); });
}
static void run_ols2_mt(void *vc, long i) {
  (void)vc;
  dispatch_apply(NT, dispatch_get_global_queue(QOS_CLASS_USER_INTERACTIVE, 0), ^(size_t t) { run_ols2(&CT[t], i); });
}
static void run_sme_mt(void *vc, long i) {
  (void)vc;
  dispatch_apply(NTS, dispatch_get_global_queue(QOS_CLASS_USER_INTERACTIVE, 0), ^(size_t t) { run_sme(&CS[t], i); });
}

/* block Toeplitz with BLAS: y[B r + i] = sum_j sum_d x[B (r + j) + d] * T_j[d][i],  T_j[d][i] = g[B j + d - i] */
static void blas_setup(int L) {
  int B = BB;
  JJ = (B + L - 2) / B + 1;
  TJ = bench_alloc(4 * (size_t)JJ * B * B);
  for (int j = 0; j < JJ; j++)
    for (int d = 0; d < B; d++)
      for (int ii = 0; ii < B; ii++) {
        int k = B * j + d - ii;
        TJ[(size_t)j * B * B + d * B + ii] = (k >= 0 && k < L) ? GREV[k] : 0.f;
      }
}
static void run_blas(void *vc, long i) {
  ctx_t *c = vc;
  const float *x = bench_pool_get(c->xp, i);
  float *y = bench_pool_get(c->yp, i);
  int B = BB;
  long n = c->n;
  for (long m = 0; m < n; m += BCH) {
    long cn = BCH < n - m ? BCH : n - m;
    int M = (int)(cn / B);
    for (int j = 0; j < JJ; j++)
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, B, B, 1.0f, x + m + (long)B * j, B,
                  TJ + (size_t)j * B * B, B, j ? 1.0f : 0.0f, y + m, B);
  }
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: alt check|time L LOGN [same|separate] [mask]\n"); return 1; }
  const char *mode = argv[1];
  int L = atoi(argv[2]), logn = atoi(argv[3]);
  long n = 1L << logn;
  if (getenv("NT")) NT = atoi(getenv("NT"));
  if (getenv("NTS")) NTS = atoi(getenv("NTS"));
  if (getenv("BLAS_B")) BB = atoi(getenv("BLAS_B"));
  if (getenv("BLAS_CH")) BCH = atol(getenv("BLAS_CH"));
  if (getenv("VPOS")) VPOS = atoi(getenv("VPOS"));
  if (getenv("OLS_PF_N")) OLS_PF_N = atoi(getenv("OLS_PF_N"));
  if (getenv("SME_CH")) SME_CH = atol(getenv("SME_CH"));
  bench_init();
  float *h = bench_alloc(4 * (size_t)L);
  bench_fill_f(h, L);
  GREV = bench_alloc(4 * (size_t)L);
  for (int k = 0; k < L; k++) GREV[k] = h[L - 1 - k];
  size_t pad = 4 * (size_t)(BB * 4 + 4096); /* BLAS panels read up to B*(J-1) - (L-1) zero floats past x */
  size_t xbytes = 4 * (size_t)(n + L - 1) + pad, ybytes = 4 * (size_t)n;
  int sep = argc > 4 && !strcmp(argv[4], "separate");
  bench_pool xp, yp;
  bench_pool_init(&xp, xbytes, (!strcmp(mode, "time") && sep) ? BENCH_SEPARATE : BENCH_SAME);
  bench_pool_init(&yp, ybytes, (!strcmp(mode, "time") && sep) ? BENCH_SEPARATE : BENCH_SAME);
  for (int b = 0; b < xp.count; b++) bench_fill_f(xp.buf[b], n + L - 1);
  int N = getenv("OLS_N") ? atoi(getenv("OLS_N")) : default_ols_N(L);
  ctx_t c;
  make_ctx(&c, L, n, &xp, &yp, h, N);
  blas_setup(L);
  /* per-thread slices */
  for (int pass = 0; pass < 2; pass++) {
    int T = pass ? NTS : NT;
    bench_pool *XP = pass ? XPS : XPT, *YP = pass ? YPS : YPT;
    ctx_t *C = pass ? CS : CT;
    long ns = n / T;
    for (int t = 0; t < T; t++) {
      XP[t].count = xp.count; XP[t].bytes = xp.bytes; XP[t].buf = malloc(sizeof(void *) * xp.count);
      YP[t].count = yp.count; YP[t].bytes = yp.bytes; YP[t].buf = malloc(sizeof(void *) * yp.count);
      for (int b = 0; b < xp.count; b++) XP[t].buf[b] = (float *)xp.buf[b] + t * ns;
      for (int b = 0; b < yp.count; b++) YP[t].buf[b] = (float *)yp.buf[b] + t * ns;
      make_ctx(&C[t], L, ns, &XP[t], &YP[t], h, N);
    }
  }
  printf("# alt L=%d n=2^%d NT=%d NTS=%d BLAS_B=%d J=%d (useful MAC fraction %.3f) BLAS_CH=%ld ols_N=%d pf_N=%d VPOS=%d\n", L,
         logn, NT, NTS, BB, JJ, (double)L / (JJ * BB), BCH, N, c.pN, VPOS);
  bench_prog all[10] = {{"vdsp_conv", run_vdsp, NULL, &c},       {"vdsp_conv_pos", run_vdsp_pos, NULL, &c},
                        {"vdsp_conv_mt", run_vdsp_mt, NULL, &c}, {"ols_pffft", run_ols_pf, NULL, &c},
                        {"pffft_ols_mt", run_pf_mt, NULL, &c},   {"ols2_vdsp", run_ols2, NULL, &c},
                        {"ols2_vdsp_mt", run_ols2_mt, NULL, &c}, {"sme_toeplitz", run_sme, NULL, &c},
                        {"sme_mt", run_sme_mt, NULL, &c},        {"blas_toep", run_blas, NULL, &c}};
  if (!strcmp(mode, "check")) {
    const float *x = xp.buf[0];
    long nchk = 0, *idx = malloc(sizeof(long) * 40000);
    for (long m = 0; m < 4096; m++) idx[nchk++] = m;
    for (long m = n - 4096; m < n; m++) idx[nchk++] = m;
    for (long m = n / 4 - 2048; m < n / 4 + 2048; m++) idx[nchk++] = m; /* thread-slice boundaries */
    for (long m = n / 2 - 2048; m < n / 2 + 2048; m++) idx[nchk++] = m;
    for (int k = 0; k < 16384; k++) idx[nchk++] = (long)(bench_rand_u64() % (uint64_t)n);
    double *ref = malloc(8 * nchk);
    for (long k = 0; k < nchk; k++) ref[k] = ref_one(x, h, L, idx[k]);
    float *got = malloc(4 * nchk);
    char lab[96];
    for (int f = 0; f < 10; f++) {
      memset(yp.buf[0], 0xff, ybytes);
      all[f].run(all[f].ctx, 0);
      const float *y = yp.buf[0];
      for (long k = 0; k < nchk; k++) got[k] = y[idx[k]];
      snprintf(lab, sizeof lab, "%s L=%d n=2^%d (%ld outs)", all[f].name, L, logn, nchk);
      bench_check_f(lab, got, ref, nchk, 1e-5);
    }
    return 0;
  }
  int mask = argc > 5 ? atoi(argv[5]) : 1023;
  bench_prog pr[10]; int np = 0;
  for (int f = 0; f < 10; f++) if (mask & (1 << f)) pr[np++] = all[f];
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)n; o.unit_name = "output";
  bench_stat st[10];
  char t[128]; snprintf(t, sizeof t, "alt L=%d n=2^%d %s NT=%d B=%d", L, logn, sep ? "separate" : "same", NT, BB);
  bench_compare(t, pr, np, &o, st);
  for (int f = 0; f < np; f++)
    printf("GMACS %-40s %-14s best %8.4f ns/out = %7.1f GMAC/s (median %8.4f)\n", t, st[f].name, st[f].best, L / st[f].best,
           st[f].median);
  return 0;
}
