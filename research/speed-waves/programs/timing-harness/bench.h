/* bench.h -- header-only timing harness for the speed-first waves.
 *
 * What it does
 *   - monotonic clock: clock_gettime_nsec_np(CLOCK_UPTIME_RAW) on macOS (CLOCK_MONOTONIC_RAW elsewhere)
 *   - calibrates, per program, a call count so that ONE repeat lasts >= 0.2 s (opts.min_rep_s)
 *   - runs >= 9 repeats (opts.min_reps); in every round the programs TAKE TURNS inside this one
 *     process (the starting program rotates each round, so no program is always first)
 *   - prints every repeat, then best and median, in ns per unit (unit = point by default)
 *   - buffer pools for the two modes: BENCH_SAME (one array, again and again: warm cache) and
 *     BENCH_SEPARATE (a ring of arrays totalling >= 64 MiB, so each call meets cold-ish data)
 *   - reference-check helpers (max abs error, relative error) to call BEFORE timing
 *   - bench_init(): asks macOS for the user-interactive QoS class (keeps the thread on a P-core)
 *
 * Use (see bench_example.c next to this file):
 *     bench_prog progs[2] = { {"mine", run_mine, reset_mine, &ctx1}, {"lib", run_lib, reset_lib, &ctx2} };
 *     bench_opts o = BENCH_OPTS_DEFAULT; o.units = n;          // ns per point
 *     bench_stat st[2]; bench_compare("wht n=1024 same-array", progs, 2, &o, st);
 *   run(ctx, i) performs ONE transform (i = call number within the repeat; use bench_pool_get(pool, i)
 *   to walk a ring of separate arrays). reset(ctx) (optional, may be NULL) runs before each repeat,
 *   outside the timed region.
 *   Run timing binaries on a quiet machine, on a performance core (the
 *   measurements here held a timing lock on the M4's performance cores).
 */
#ifndef WHT42_BENCH_H
#define WHT42_BENCH_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef __APPLE__
#include <pthread.h>
#include <sys/sysctl.h>
#endif

static inline uint64_t bench_now_ns(void) {
#ifdef __APPLE__
  return clock_gettime_nsec_np(CLOCK_UPTIME_RAW);
#else
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC_RAW, &t);
  return (uint64_t)t.tv_sec * 1000000000ull + (uint64_t)t.tv_nsec;
#endif
}

typedef struct {
  const char *name;
  void (*run)(void *ctx, long i); /* ONE transform per call */
  void (*reset)(void *ctx);       /* optional; before each repeat, not timed */
  void *ctx;
} bench_prog;

typedef struct {
  double min_rep_s;      /* seconds per repeat, >= 0.2 */
  int min_reps;          /* repeats, >= 9 */
  double units;          /* units per call (e.g. n points) for the ns-per-unit figure */
  const char *unit_name; /* "point", "transform", ... */
  int print_reps;        /* 1 = print every repeat (house rule: keep 1) */
} bench_opts;
#define BENCH_OPTS_DEFAULT {0.2, 9, 1.0, "point", 1}

#define BENCH_MAX_REPS 64
typedef struct {
  const char *name;
  long calls_per_rep;
  int nrep;
  double rep[BENCH_MAX_REPS]; /* ns per unit, in the order run */
  double best, median;
} bench_stat;

static void bench_init(void) {
#ifdef __APPLE__
  pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
  char brand[128] = {0};
  size_t len = sizeof brand;
  sysctlbyname("machdep.cpu.brand_string", brand, &len, NULL, 0);
  printf("# bench.h: cpu=%s qos=user-interactive clock=CLOCK_UPTIME_RAW\n", brand);
#endif
}

static int bench_cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

static double bench_time_calls(bench_prog *p, long k) {
  if (p->reset) p->reset(p->ctx);
  uint64_t t0 = bench_now_ns();
  for (long i = 0; i < k; i++) p->run(p->ctx, i);
  return (double)(bench_now_ns() - t0);
}

static long bench_calibrate(bench_prog *p, double min_s) {
  long k = 1;
  bench_time_calls(p, 1); /* warm up: first touch, page faults, lazy setup */
  for (;;) {
    double ns = bench_time_calls(p, k);
    if (ns >= min_s * 1e9) return k;
    if (ns < min_s * 1e9 / 16) k *= 8;
    else k = (long)ceil((double)k * min_s * 1e9 * 1.1 / ns);
  }
}

/* Run all programs in turns; fills st[0..n-1]; returns 0. */
static int bench_compare(const char *title, bench_prog *progs, int n, const bench_opts *o, bench_stat *st) {
  int reps = o->min_reps < 9 ? 9 : o->min_reps;
  if (reps > BENCH_MAX_REPS) reps = BENCH_MAX_REPS;
  double min_s = o->min_rep_s < 0.2 ? 0.2 : o->min_rep_s;
  printf("## %s  (%d repeats of >= %.2f s each, programs take turns; ns per %s)\n", title, reps, min_s, o->unit_name);
  for (int j = 0; j < n; j++) {
    st[j].name = progs[j].name;
    st[j].nrep = reps;
    st[j].calls_per_rep = bench_calibrate(&progs[j], min_s);
  }
  for (int r = 0; r < reps; r++) {
    for (int jj = 0; jj < n; jj++) {
      int j = (jj + r) % n;
      double ns = bench_time_calls(&progs[j], st[j].calls_per_rep);
      st[j].rep[r] = ns / ((double)st[j].calls_per_rep * o->units);
      if (o->print_reps) printf("rep %2d  %-28s %12.4f\n", r, progs[j].name, st[j].rep[r]);
    }
  }
  for (int j = 0; j < n; j++) {
    double tmp[BENCH_MAX_REPS];
    memcpy(tmp, st[j].rep, sizeof(double) * reps);
    qsort(tmp, reps, sizeof(double), bench_cmp_d);
    st[j].best = tmp[0];
    st[j].median = (reps & 1) ? tmp[reps / 2] : 0.5 * (tmp[reps / 2 - 1] + tmp[reps / 2]);
  }
  for (int j = 0; j < n; j++)
    printf("RESULT %-40s %-28s best %10.4f  median %10.4f  ns/%s  (calls/rep %ld)\n", title, st[j].name,
           st[j].best, st[j].median, o->unit_name, st[j].calls_per_rep);
  fflush(stdout);
  return 0;
}

/* ---------- buffers ---------- */
enum { BENCH_SAME = 0, BENCH_SEPARATE = 1 };
typedef struct {
  int count;
  size_t bytes;
  void **buf;
} bench_pool;

static void *bench_alloc(size_t bytes) { /* 128-byte aligned (M4 cache line), zeroed */
  void *p = NULL;
  if (posix_memalign(&p, 128, bytes ? bytes : 128)) {
    fprintf(stderr, "bench_alloc failed\n");
    exit(1);
  }
  memset(p, 0, bytes);
  return p;
}

/* SAME: 1 array.  SEPARATE: enough arrays to total >= 64 MiB (at least 2), past the M4's 16 MB L2. */
static void bench_pool_init(bench_pool *p, size_t bytes, int mode) {
  size_t target = (size_t)64 << 20;
  p->bytes = bytes;
  p->count = (mode == BENCH_SAME) ? 1 : (int)((target + bytes - 1) / bytes);
  if (mode != BENCH_SAME && p->count < 2) p->count = 2;
  p->buf = (void **)malloc(sizeof(void *) * p->count);
  for (int i = 0; i < p->count; i++) p->buf[i] = bench_alloc(bytes);
}
static inline void *bench_pool_get(const bench_pool *p, long i) { return p->buf[i % p->count]; }
static void bench_pool_free(bench_pool *p) {
  for (int i = 0; i < p->count; i++) free(p->buf[i]);
  free(p->buf);
}

/* ---------- inputs and reference checks ---------- */
static uint64_t bench_rng_state = 0x9E3779B97F4A7C15ull;
static inline uint64_t bench_rand_u64(void) { /* xorshift64* */
  bench_rng_state ^= bench_rng_state >> 12;
  bench_rng_state ^= bench_rng_state << 25;
  bench_rng_state ^= bench_rng_state >> 27;
  return bench_rng_state * 0x2545F4914F6CDD1Dull;
}
static inline double bench_rand_pm1(void) { return (double)(bench_rand_u64() >> 11) * (2.0 / 9007199254740992.0) - 1.0; }
static void bench_fill_f(float *x, size_t n) { for (size_t i = 0; i < n; i++) x[i] = (float)bench_rand_pm1(); }
static void bench_fill_d(double *x, size_t n) { for (size_t i = 0; i < n; i++) x[i] = bench_rand_pm1(); }

/* Prints "CHECK label: max abs err E, rel err R (tol T) PASS|FAIL"; returns 1 on pass.
 * rel err = max|got-ref| / max|ref|. */
static int bench_check_f(const char *label, const float *got, const double *ref, size_t n, double tol) {
  double emax = 0, rmax = 0;
  for (size_t i = 0; i < n; i++) {
    double e = fabs((double)got[i] - ref[i]);
    if (!(e <= emax)) emax = e; /* catches NaN */
    if (fabs(ref[i]) > rmax) rmax = fabs(ref[i]);
  }
  double rel = rmax > 0 ? emax / rmax : emax;
  int ok = rel <= tol;
  printf("CHECK %-40s max abs err %.3e  rel err %.3e  (tol %.1e) %s\n", label, emax, rel, tol, ok ? "PASS" : "FAIL");
  return ok;
}
static int bench_check_d(const char *label, const double *got, const double *ref, size_t n, double tol) {
  double emax = 0, rmax = 0;
  for (size_t i = 0; i < n; i++) {
    double e = fabs(got[i] - ref[i]);
    if (!(e <= emax)) emax = e;
    if (fabs(ref[i]) > rmax) rmax = fabs(ref[i]);
  }
  double rel = rmax > 0 ? emax / rmax : emax;
  int ok = rel <= tol;
  printf("CHECK %-40s max abs err %.3e  rel err %.3e  (tol %.1e) %s\n", label, emax, rel, tol, ok ? "PASS" : "FAIL");
  return ok;
}

/* Plain reference WHT (natural/Hadamard order, unnormalised) in double. */
static void bench_ref_wht_d(double *x, int log_n) {
  size_t n = (size_t)1 << log_n;
  for (size_t h = 1; h < n; h <<= 1)
    for (size_t i = 0; i < n; i += 2 * h)
      for (size_t j = i; j < i + h; j++) {
        double a = x[j], b = x[j + h];
        x[j] = a + b;
        x[j + h] = a - b;
      }
}

#endif /* WHT42_BENCH_H */
