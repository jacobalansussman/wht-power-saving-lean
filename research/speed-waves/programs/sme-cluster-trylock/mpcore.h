/* mpcore.h -- as-3 (wave 14): multi-PROCESS driver for "who gets the SME unit" policies on the M4.
 *
 * W worker processes (fork) each run calls of one kernel in a loop; the parent measures aggregate throughput over
 * R windows and every call's latency (only calls that start inside the measured windows).
 * Policies (POL env or argv):
 *   neon      every worker runs the NEON path
 *   sme       every worker runs the SME path itself (each enters streaming mode; the unit is shared in hardware)
 *   static:K  workers 0..K-1 run SME, the rest NEON (the wave-13 split)
 *   token:K   K tokens in shared memory; a call try-acquires one: success -> SME part (token released before the
 *             caller's NEON tail), failure -> NEON path at once.
 *   owner     one resident owner PROCESS (extra, user-interactive QoS) runs the SME part of submitted jobs; a caller
 *             try-claims the owner's single slot (CAS); success -> submit, wait (spin SPIN_US then sleep on the
 *             address), then its NEON tail; failure -> NEON path at once. Buffers live in one MAP_SHARED arena created
 *             before fork, so a job is a worker index (zero copy). COPY=1: the caller also copies its input into the
 *             slot's staging buffer and the output back (the cost of a drop-in with private buffers).
 * The kernel file defines the hooks below before including this header.
 */
#include <os/os_sync_wait_on_address.h>
#include <signal.h>
#include <stdatomic.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>

#define MAXW 16
#define MAXLAT 65536

typedef struct {
  _Atomic int ready, go, stop, measuring;
  _Atomic long units[MAXW];
  _Atomic long nsme[MAXW], nneon[MAXW];
  _Atomic int tokens;
  _Atomic int tokens_e; /* ctoken: the E-cluster pool */
  _Atomic int holders_p; /* cauto: SME holders on the P cluster */
  _Atomic uint64_t last_ns[MAXW]; _Atomic int last_p[MAXW]; /* cauto: when each worker last started a call, and on which cluster */
  _Atomic long sme_on_e; /* SME calls made from an E-core (cpu number < 6) */
  /* owner */
  _Atomic uint32_t slot;   /* 0 free, w+1 claimed by worker w */
  _Atomic uint32_t state;  /* 0 idle, 1 submitted, 2 done */
  _Atomic int owner_sleeping, caller_sleeping, owner_exit, owner_ready;
  _Atomic long owner_jobs;
  int job_w;
  /* results */
  int nlat[MAXW];
  long ncall_meas[MAXW];
  int ok[MAXW];
  double neon_lat_sum[MAXW]; long neon_lat_n[MAXW];
  double sme_lat_sum[MAXW]; long sme_lat_n[MAXW];
  float lat[MAXW][MAXLAT]; /* microseconds */
} shm_t;

static shm_t *S;
static int g_spin_us = 3;
static int g_copy = 0;

static void *shm_alloc(size_t bytes) {
  void *p = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);
  if (p == MAP_FAILED) { perror("mmap"); exit(1); }
  memset(p, 0, bytes);
  return p;
}
/* arena allocator for buffers that the owner must see (pre-fork, same address in every process) */
static char *g_arena; static size_t g_arena_used, g_arena_size;
static void *arena_get(size_t bytes) {
  size_t a = (g_arena_used + 127) & ~(size_t)127;
  if (a + bytes > g_arena_size) { fprintf(stderr, "arena full\n"); exit(1); }
  g_arena_used = a + bytes; return g_arena + a;
}

static inline void wait_while_eq(_Atomic uint32_t *a, uint32_t v, _Atomic int *sleeping_flag) {
  uint64_t t0 = bench_now_ns(), lim = (uint64_t)g_spin_us * 1000;
  while (atomic_load_explicit(a, memory_order_acquire) == v) {
    if (bench_now_ns() - t0 > lim) {
      atomic_store(sleeping_flag, 1);
      while (atomic_load(a) == v) os_sync_wait_on_address((void *)a, v, 4, OS_SYNC_WAIT_ON_ADDRESS_SHARED);
      atomic_store(sleeping_flag, 0);
      return;
    }
  }
}
static inline void wake(_Atomic uint32_t *a, _Atomic int *sleeping_flag) {
  if (atomic_load(sleeping_flag)) os_sync_wake_by_address_all((void *)a, 4, OS_SYNC_WAKE_BY_ADDRESS_SHARED);
}

enum { P_NEON, P_SME, P_STATIC, P_TOKEN, P_OWNER, P_SPIN, P_CTOKEN, P_CAUTO };
static int g_W; static uint64_t g_win_ns = 5000000;
static int g_relaxed; /* RELAXED=1: token / cauto acquire and release with memory_order_relaxed */
static int g_KE;
/* the CPU number: on this macOS (26, M4) the low bits of TPIDR_EL0 hold it [ran: cpunum.c]: P-cores 6..9, E-cores 0..5 */
static inline int cur_cpu(void) { uint64_t p; __asm__ volatile("mrs %0, TPIDR_EL0" : "=r"(p)); return (int)(p & 0xff); }
static inline int try_take(_Atomic int *pool) {
  int t = atomic_load(pool);
  while (t > 0) if (atomic_compare_exchange_weak(pool, &t, t - 1)) return 1;
  return 0;
}
static int g_nbr_us = 400; /* spin:K neighbours: scalar busy loop of this many us per call (no units counted) */
static int g_pol, g_K;

/* ---- hooks the kernel file provides ----
 * static void *k_ctx[MAXW];                       worker contexts made pre-fork by k_make(w)
 * static void k_global(void);                     tables, reference (parent, pre-fork)
 * static void *k_make(int w);
 * static void k_neon(void *c);                    whole call on NEON
 * static void k_sme_part(void *c);                the part that needs the SME unit (locally streaming)
 * static void k_tail(void *c);                    the caller's NEON tail after the SME part (may be empty)
 * static void k_copy_in(void *c), k_copy_out(void *c);   COPY=1 staging copies
 * static int k_check(void *c, int path, const char *lab);  1 = PASS
 * static const double K_UNITS; static const char *K_UNIT, *K_TITLE;
 */

static void owner_loop(void) {
  pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
  k_owner_init();
  atomic_store(&S->owner_ready, 1);
  for (;;) {
    wait_while_eq(&S->state, 0, &S->owner_sleeping);
    if (atomic_load(&S->owner_exit)) break;
    if (atomic_load_explicit(&S->state, memory_order_acquire) != 1) continue;
    int w = S->job_w;
    k_sme_part(k_ctx[w]);
    atomic_fetch_add(&S->owner_jobs, 1);
    atomic_store_explicit(&S->state, 2, memory_order_release);
    wake(&S->state, &S->caller_sleeping);
  }
  _exit(0);
}

static int do_call(int w, void *c) { /* returns 1 if SME was used */
  switch (g_pol) {
    case P_NEON: k_neon(c); return 0;
    case P_SME: k_sme_part(c); k_tail(c); return 1;
    case P_STATIC: if (w < g_K) { k_sme_part(c); k_tail(c); return 1; } k_neon(c); return 0;
    case P_TOKEN: {
      if (g_relaxed) {
        int t = atomic_load_explicit(&S->tokens, memory_order_relaxed);
        while (t > 0) if (atomic_compare_exchange_weak_explicit(&S->tokens, &t, t - 1, memory_order_relaxed, memory_order_relaxed)) break;
        if (t > 0) { k_sme_part(c); atomic_fetch_add_explicit(&S->tokens, 1, memory_order_relaxed); k_tail(c); return 1; }
        k_neon(c); return 0;
      }
      int t = atomic_load(&S->tokens);
      while (t > 0) if (atomic_compare_exchange_weak(&S->tokens, &t, t - 1)) break;
      if (t > 0) { if (cur_cpu() < 6) atomic_fetch_add(&S->sme_on_e, 1); k_sme_part(c); atomic_fetch_add(&S->tokens, 1); k_tail(c); return 1; }
      k_neon(c); return 0;
    }
    case P_CTOKEN: {
      _Atomic int *pool = cur_cpu() >= 6 ? &S->tokens : &S->tokens_e;
      if (try_take(pool)) { if (pool == &S->tokens_e) atomic_fetch_add(&S->sme_on_e, 1); k_sme_part(c); atomic_fetch_add(pool, 1); k_tail(c); return 1; }
      k_neon(c); return 0;
    }
    case P_CAUTO: { /* per-cluster token with an automatic limit: E-cores never take SME; on the P cluster at most
                     1 holder while <= 2 workers are active there (started a call in the last 5 ms), else 2 */
      int onP = cur_cpu() >= 6; uint64_t now = bench_now_ns();
      atomic_store_explicit(&S->last_ns[w], now, memory_order_relaxed); atomic_store_explicit(&S->last_p[w], onP, memory_order_relaxed);
      if (onP) {
        int A = 0;
        for (int v = 0; v < g_W; v++) if (atomic_load_explicit(&S->last_p[v], memory_order_relaxed) && now - atomic_load_explicit(&S->last_ns[v], memory_order_relaxed) < g_win_ns) A++;
        int lim = A >= 3 ? 2 : 1, h = atomic_load(&S->holders_p);
        if (g_relaxed) {
          h = atomic_load_explicit(&S->holders_p, memory_order_relaxed);
          while (h < lim) if (atomic_compare_exchange_weak_explicit(&S->holders_p, &h, h + 1, memory_order_relaxed, memory_order_relaxed)) { k_sme_part(c); atomic_fetch_sub_explicit(&S->holders_p, 1, memory_order_relaxed); k_tail(c); return 1; }
        } else
        while (h < lim) if (atomic_compare_exchange_weak(&S->holders_p, &h, h + 1)) { k_sme_part(c); atomic_fetch_sub(&S->holders_p, 1); k_tail(c); return 1; }
      }
      k_neon(c); return 0;
    }
    case P_SPIN: if (w < g_K) { uint64_t t0 = bench_now_ns(); volatile uint64_t x = 1; while (bench_now_ns() - t0 < (uint64_t)g_nbr_us * 1000) for (int i = 0; i < 64; i++) x = x * 6364136223846793005ull + 1; return 2; }
      k_neon(c); return 0;
    case P_OWNER: {
      uint32_t z = 0;
      if (atomic_load_explicit(&S->slot, memory_order_relaxed) == 0 && atomic_compare_exchange_strong(&S->slot, &z, (uint32_t)w + 1)) {
        if (g_copy) k_copy_in(c);
        S->job_w = w;
        atomic_store(&S->state, 1);
        wake(&S->state, &S->owner_sleeping);
        wait_while_eq(&S->state, 1, &S->caller_sleeping);
        atomic_store(&S->state, 0);
        atomic_store_explicit(&S->slot, 0, memory_order_release);
        if (g_copy) k_copy_out(c);
        k_tail(c); return 1;
      }
      k_neon(c); return 0;
    }
  }
  return 0;
}

static void worker_loop(int w, int qos_e) {
  pthread_set_qos_class_self_np(qos_e ? QOS_CLASS_BACKGROUND : QOS_CLASS_USER_INTERACTIVE, 0);
  void *c = k_ctx[w];
  k_worker_init(c);
  /* warm both paths once (first touch, TLS) */
  k_neon(c);
  atomic_fetch_add(&S->ready, 1);
  while (!atomic_load(&S->go)) ;
  int used_sme = 0, used_neon = 0;
  while (!atomic_load(&S->stop)) {
    int meas = atomic_load_explicit(&S->measuring, memory_order_relaxed);
    uint64_t t0 = bench_now_ns();
    int s = do_call(w, c);
    uint64_t t1 = bench_now_ns();
    if (s == 2) continue;
    atomic_fetch_add_explicit(&S->units[w], (long)K_UNITS, memory_order_relaxed);
    if (s) { used_sme = 1; atomic_fetch_add_explicit(&S->nsme[w], 1, memory_order_relaxed); }
    else { used_neon = 1; atomic_fetch_add_explicit(&S->nneon[w], 1, memory_order_relaxed); }
    if (meas) {
      double us = (t1 - t0) / 1e3;
      S->lat[w][S->nlat[w] % MAXLAT] = s ? -(float)us : (float)us; S->nlat[w]++;
      S->ncall_meas[w]++;
      if (s) { S->sme_lat_sum[w] += us; S->sme_lat_n[w]++; } else { S->neon_lat_sum[w] += us; S->neon_lat_n[w]++; }
    }
  }
  int ok = 1; char lab[64];
  if (used_neon) { snprintf(lab, sizeof lab, "w%d neon", w); ok &= k_check(c, 0, lab); }
  if (used_sme) { snprintf(lab, sizeof lab, "w%d sme", w); ok &= k_check(c, 1, lab); }
  S->ok[w] = ok;
  _exit(0);
}

static int cmpf(const void *a, const void *b) { float x = *(const float *)a, y = *(const float *)b; return (x > y) - (x < y); }

/* run one configuration: W workers, policy string, nE of them at background QoS (E-cores) */
static void run_config(int W, const char *pol, int nE) {
  if (!strcmp(pol, "neon")) g_pol = P_NEON;
  else if (!strcmp(pol, "sme")) g_pol = P_SME;
  else if (!strncmp(pol, "static:", 7)) { g_pol = P_STATIC; g_K = atoi(pol + 7); }
  else if (!strncmp(pol, "token:", 6)) { g_pol = P_TOKEN; g_K = atoi(pol + 6); }
  else if (!strncmp(pol, "owner", 5)) { g_pol = P_OWNER; g_spin_us = pol[5] == ':' ? atoi(pol + 6) : (getenv("SPIN_US") ? atoi(getenv("SPIN_US")) : 3); }
  else if (!strncmp(pol, "ctoken:", 7)) { g_pol = P_CTOKEN; g_K = atoi(pol + 7); const char *q = strchr(pol + 7, ':'); g_KE = q ? atoi(q + 1) : 0; }
  else if (!strcmp(pol, "cauto")) g_pol = P_CAUTO;
  else if (!strncmp(pol, "spin:", 5)) { g_pol = P_SPIN; g_K = atoi(pol + 5); }
  else { printf("bad policy %s\n", pol); return; }
  memset(S, 0, sizeof *S); g_W = W;
  S->tokens = (g_pol == P_TOKEN || g_pol == P_CTOKEN) ? g_K : 0;
  S->tokens_e = g_pol == P_CTOKEN ? g_KE : 0;
  pid_t opid = 0;
  if (g_pol == P_OWNER) { opid = fork(); if (opid == 0) owner_loop(); while (!atomic_load(&S->owner_ready)) ; }
  pid_t pid[MAXW];
  for (int w = 0; w < W; w++) { pid[w] = fork(); if (pid[w] == 0) worker_loop(w, w >= W - nE); }
  while (atomic_load(&S->ready) < W) ;
  atomic_store(&S->go, 1);
  struct timespec warm = {0, 300000000}; nanosleep(&warm, NULL);
  int R = getenv("R") ? atoi(getenv("R")) : 7; double res[64];
  atomic_store(&S->measuring, 1);
  uint64_t tm0 = bench_now_ns();
  for (int r = 0; r < R; r++) {
    long u0 = 0, u1 = 0;
    for (int w = 0; w < W; w++) u0 += atomic_load(&S->units[w]);
    uint64_t t0 = bench_now_ns();
    struct timespec ws = {0, 400000000}; nanosleep(&ws, NULL);
    uint64_t t1 = bench_now_ns();
    for (int w = 0; w < W; w++) u1 += atomic_load(&S->units[w]);
    res[r] = (double)(t1 - t0) / (double)(u1 - u0);
  }
  atomic_store(&S->measuring, 0);
  double tmeas = (bench_now_ns() - tm0) / 1e9;
  atomic_store(&S->stop, 1);
  for (int w = 0; w < W; w++) waitpid(pid[w], NULL, 0);
  if (opid) { atomic_store(&S->owner_exit, 1); atomic_store(&S->state, 3); os_sync_wake_by_address_all((void *)&S->state, 4, OS_SYNC_WAKE_BY_ADDRESS_SHARED); waitpid(opid, NULL, 0); }
  /* latency */
  long nl = 0; for (int w = 0; w < W; w++) nl += S->nlat[w] < MAXLAT ? S->nlat[w] : MAXLAT;
  float *all = malloc(sizeof(float) * (nl + 1)); long k = 0;
  for (int w = 0; w < W; w++) { int n = S->nlat[w] < MAXLAT ? S->nlat[w] : MAXLAT; memcpy(all + k, S->lat[w], sizeof(float) * n); k += n; }
  float *an = malloc(sizeof(float) * (nl + 1)), *as = malloc(sizeof(float) * (nl + 1)); long nn_ = 0, ns_ = 0;
  for (long i = 0; i < nl; i++) { if (all[i] < 0) { as[ns_++] = -all[i]; all[i] = -all[i]; } else an[nn_++] = all[i]; }
  qsort(all, nl, sizeof(float), cmpf); qsort(an, nn_, sizeof(float), cmpf); qsort(as, ns_, sizeof(float), cmpf);
  double p50 = nl ? all[nl / 2] : 0, p99 = nl ? all[(long)(nl * 0.99)] : 0, pmax = nl ? all[nl - 1] : 0;
  double n50 = nn_ ? an[nn_ / 2] : 0, n99 = nn_ ? an[(long)(nn_ * 0.99)] : 0, s50 = ns_ ? as[ns_ / 2] : 0, s99 = ns_ ? as[(long)(ns_ * 0.99)] : 0;
  free(all); free(an); free(as);
  double sres[64]; memcpy(sres, res, sizeof(double) * R); qsort(sres, R, sizeof(double), bench_cmp_d);
  long ns = 0, nn = 0; double sl = 0, nlsum = 0; long sln = 0, nln = 0; int ok = 1;
  for (int w = 0; w < W; w++) { ns += S->nsme[w]; nn += S->nneon[w]; sl += S->sme_lat_sum[w]; sln += S->sme_lat_n[w]; nlsum += S->neon_lat_sum[w]; nln += S->neon_lat_n[w]; ok &= S->ok[w]; }
  printf("reps %s W=%d E=%d %s:", K_TITLE, W, nE, pol);
  for (int r = 0; r < R; r++) printf(" %.3f", res[r]);
  printf("\n");
  printf("PATHS %-10s W=%d E=%d %-9s neon calls p50 %.1f p99 %.1f us (n=%ld) | sme calls p50 %.1f p99 %.1f us (n=%ld) | owner jobs %ld | SME calls from E-cores %ld\n", K_TITLE, W, nE, pol, n50, n99, nn_, s50, s99, ns_, (long)atomic_load(&S->owner_jobs), (long)atomic_load(&S->sme_on_e));
  if (g_relaxed && (g_pol == P_TOKEN || g_pol == P_CAUTO)) printf("NOTE relaxed token ops\n");
  printf("RESULT %-10s W=%d E=%d %-9s%s median %8.3f best %8.3f ns/%s agg  = %8.3f M%s/s | call p50 %8.1f p99 %8.1f max %8.1f us | SME calls %5.1f %% | mean us: sme %.1f neon %.1f | calls %ld | CHECK %s\n",
         K_TITLE, W, nE, pol, g_copy && g_pol == P_OWNER ? "+copy" : "", sres[R / 2], sres[0], K_UNIT, 1e3 / sres[R / 2], K_UNIT, p50, p99, pmax,
         100.0 * ns / (ns + nn > 0 ? ns + nn : 1), sln ? sl / sln : 0, nln ? nlsum / nln : 0, (long)(tmeas > 0 ? sln + nln : 0), ok ? "PASS" : "FAIL");
  fflush(stdout);
}

static int mp_main(int argc, char **argv) {
  /* usage: prog W:POL[:E=n] ...  e.g. 4:neon 4:static:1 4:token:1 4:owner 8:owner:E=4 */
  if (getenv("SPIN_US")) g_spin_us = atoi(getenv("SPIN_US"));
  if (getenv("COPY")) g_copy = atoi(getenv("COPY"));
  if (getenv("NBR_US")) g_nbr_us = atoi(getenv("NBR_US"));
  if (getenv("RELAXED")) g_relaxed = atoi(getenv("RELAXED"));
  S = shm_alloc(sizeof(shm_t));
  int maxW = 0;
  for (int a = 1; a < argc; a++) { int W = atoi(argv[a]); if (W > maxW) maxW = W; }
  if (maxW > MAXW) maxW = MAXW;
  g_arena_size = (size_t)(maxW + 1) * K_ARENA_PER_WORKER + (64 << 20);
  g_arena = shm_alloc(g_arena_size);
  k_global();
  for (int w = 0; w < maxW; w++) k_ctx[w] = k_make(w);
  for (int a = 1; a < argc; a++) {
    char buf[64]; snprintf(buf, sizeof buf, "%s", argv[a]);
    int W = atoi(buf); char *p = strchr(buf, ':'); if (!p) continue; p++;
    int nE = 0; char *e = strstr(p, ":E="); if (e) { nE = atoi(e + 3); *e = 0; }
    run_config(W, p, nE);
  }
  return 0;
}
