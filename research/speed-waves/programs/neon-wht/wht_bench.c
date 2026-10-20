/* wht_bench.c -- check FFHT-NEON (executorch copy) and our NEON textbook WHT against a plain
 * reference for 2^1..2^24, then time both (fp32, in place, single thread).
 * Usage: wht_bench check            -> correctness only
 *        wht_bench time MODE LOGN... -> MODE = same | separate
 * Build: see baselines/README.md */
#include "../timing-harness/bench.h"
#include "wht_neon.h"
#include "fht.h"

typedef struct {
  int log_n;
  bench_pool pool;
  float *init; /* pristine input, copied in by reset */
  int which;   /* 0 = FFHT, 1 = ours */
} wctx;

static void w_reset(void *v) {
  wctx *c = v;
  size_t n = (size_t)1 << c->log_n;
  for (int i = 0; i < c->pool.count; i++) memcpy(c->pool.buf[i], c->init, n * sizeof(float));
}
static void w_run(void *v, long i) {
  wctx *c = v;
  float *b = bench_pool_get(&c->pool, i);
  if (c->which == 0) fht_float(b, c->log_n);
  else wht_neon_f32(b, c->log_n);
}

static int check_all(void) {
  int ok = 1;
  for (int L = 1; L <= 24; L++) {
    size_t n = (size_t)1 << L;
    float *x = bench_alloc(n * 4), *y = bench_alloc(n * 4);
    double *r = malloc(n * 8);
    bench_fill_f(x, n);
    for (size_t i = 0; i < n; i++) r[i] = x[i];
    memcpy(y, x, n * 4);
    bench_ref_wht_d(r, L);
    fht_float(x, L);
    wht_neon_f32(y, L);
    char lab[64];
    snprintf(lab, sizeof lab, "FFHT-NEON  2^%d", L);
    ok &= bench_check_f(lab, x, r, n, 1e-5);
    snprintf(lab, sizeof lab, "ours-NEON  2^%d", L);
    ok &= bench_check_f(lab, y, r, n, 1e-5);
    free(x); free(y); free(r);
  }
  printf("ALL %s\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}

int main(int argc, char **argv) {
  bench_init();
  if (getenv("WHT_STRIPS")) wht_neon_strips = atoi(getenv("WHT_STRIPS"));
  if (getenv("WHT_B")) wht_neon_block_log = atoi(getenv("WHT_B"));
  if (getenv("WHT_MAXR")) wht_neon_high_maxr = atoi(getenv("WHT_MAXR"));
  if (wht_neon_high_maxr < 1 || wht_neon_high_maxr > 4) wht_neon_high_maxr = 3; /* 0 would loop forever */
  printf("# ours: strips=%d B=%d high maxr=%d\n", wht_neon_strips, wht_neon_block_log, wht_neon_high_maxr);
  if (argc < 2 || !strcmp(argv[1], "check")) return check_all();
  int mode = !strcmp(argv[2], "separate") ? BENCH_SEPARATE : BENCH_SAME;
  for (int a = 3; a < argc; a++) {
    int L = atoi(argv[a]);
    size_t n = (size_t)1 << L;
    wctx c[2];
    for (int w = 0; w < 2; w++) {
      c[w].log_n = L;
      c[w].which = w;
      bench_pool_init(&c[w].pool, n * 4, mode);
      c[w].init = bench_alloc(n * 4);
      if (w == 0) { if (getenv("WHT_ZERO")) memset(c[w].init, 0, n * 4); else bench_fill_f(c[w].init, n); }
      else memcpy(c[w].init, c[0].init, n * 4); /* same input for both */
    }
    /* quick check on this exact path before timing */
    {
      double *r = malloc(n * 8);
      for (size_t i = 0; i < n; i++) r[i] = c[0].init[i];
      bench_ref_wht_d(r, L);
      for (int w = 0; w < 2; w++) {
        w_reset(&c[w]);
        w_run(&c[w], 0);
        bench_check_f(w ? "ours-NEON pre-timing" : "FFHT-NEON pre-timing", c[w].pool.buf[0], r, n, 1e-5);
      }
      free(r);
    }
    bench_prog p[2] = {{"FFHT-NEON", w_run, w_reset, &c[0]}, {"ours-NEON-textbook", w_run, w_reset, &c[1]}};
    bench_opts o = BENCH_OPTS_DEFAULT;
    o.units = (double)n;
    bench_stat st[2];
    char title[64];
    snprintf(title, sizeof title, "WHT fp32 2^%d %s", L, mode ? "separate" : "same");
    bench_compare(title, p, 2, &o, st);
    for (int w = 0; w < 2; w++) {
      bench_pool_free(&c[w].pool);
      free(c[w].init);
    }
  }
  return 0;
}
