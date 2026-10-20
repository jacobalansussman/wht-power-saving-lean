/* bench_example.c -- the smallest use of bench.h: two in-place fp32 WHT programs (a plain loop and
 * a 2-stage-unrolled loop), checked against the double reference, then timed in turns, in the
 * same-array mode and the separate-arrays mode.
 * Build: clang -O3 -mcpu=apple-m4 bench_example.c -o bench_example
 * Run:   ./bench_example 10            (log2 of the size) */
#include "bench.h"

typedef struct { int log_n; bench_pool pool; const float *init; } ctx_t;

static void wht_plain(float *x, int L) {
  size_t n = (size_t)1 << L;
  for (size_t h = 1; h < n; h <<= 1)
    for (size_t i = 0; i < n; i += 2 * h)
      for (size_t j = i; j < i + h; j++) { float a = x[j], b = x[j + h]; x[j] = a + b; x[j + h] = a - b; }
}
static void wht_radix4(float *x, int L) { /* two stages per sweep */
  size_t n = (size_t)1 << L, h = 1;
  for (; 4 * h <= n; h *= 4)
    for (size_t i = 0; i < n; i += 4 * h)
      for (size_t j = i; j < i + h; j++) {
        float a = x[j], b = x[j + h], c = x[j + 2 * h], d = x[j + 3 * h];
        float s0 = a + b, s1 = a - b, s2 = c + d, s3 = c - d;
        x[j] = s0 + s2; x[j + h] = s1 + s3; x[j + 2 * h] = s0 - s2; x[j + 3 * h] = s1 - s3;
      }
  if (h < n)
    for (size_t j = 0; j < h; j++) { float a = x[j], b = x[j + h]; x[j] = a + b; x[j + h] = a - b; }
}
static void reset(void *v) {
  ctx_t *c = v;
  for (int i = 0; i < c->pool.count; i++) memcpy(c->pool.buf[i], c->init, sizeof(float) << c->log_n);
}
static void run_plain(void *v, long i) { ctx_t *c = v; wht_plain(bench_pool_get(&c->pool, i), c->log_n); }
static void run_r4(void *v, long i) { ctx_t *c = v; wht_radix4(bench_pool_get(&c->pool, i), c->log_n); }

int main(int argc, char **argv) {
  int L = argc > 1 ? atoi(argv[1]) : 10;
  size_t n = (size_t)1 << L;
  bench_init();
  float *init = bench_alloc(n * 4);
  bench_fill_f(init, n);
  /* 1. reference check BEFORE timing */
  double *ref = malloc(n * 8);
  for (size_t i = 0; i < n; i++) ref[i] = init[i];
  bench_ref_wht_d(ref, L);
  float *t = bench_alloc(n * 4);
  memcpy(t, init, n * 4); wht_plain(t, L); bench_check_f("plain", t, ref, n, 1e-5);
  memcpy(t, init, n * 4); wht_radix4(t, L); bench_check_f("radix4", t, ref, n, 1e-5);
  /* 2. time, both modes */
  for (int mode = BENCH_SAME; mode <= BENCH_SEPARATE; mode++) {
    ctx_t a = {L, {0}, init}, b = {L, {0}, init};
    bench_pool_init(&a.pool, n * 4, mode);
    bench_pool_init(&b.pool, n * 4, mode);
    bench_prog p[2] = {{"plain", run_plain, reset, &a}, {"radix4", run_r4, reset, &b}};
    bench_opts o = BENCH_OPTS_DEFAULT;
    o.units = (double)n;
    bench_stat st[2];
    bench_compare(mode == BENCH_SAME ? "example same-array" : "example separate-arrays", p, 2, &o, st);
    bench_pool_free(&a.pool); bench_pool_free(&b.pool);
  }
  return 0;
}
