/* mp_lm.c -- as-3: multi-process SME-owner study, kernel = ap-2's Whisper log-mel (wave 13, lm.c, included unchanged).
 * NEON path = B-gmel-flog (FFTW 3.3.10 two-frame STFT + generated NEON mel + fast log: ap-2's strongest baseline).
 * SME part = sme_lm400 (STFT + power + mel in ZA, one streaming call); tail = NEON fast log + clamp (flog_all).
 * One call = one 30 s chunk = 3000 frames. usage: mp_lm W:POL[:E=n] ...   (see mpcore.h) */
#define main lm_main
#include "logmel/lm.c"
#undef main

#define FR 3000
typedef struct { lctx n, s; float *priv_sig, *priv_out; } lmw;
static void *k_ctx[16];
static double K_UNITS = FR;
static const char *K_UNIT = "frame", *K_TITLE = "logmel";
static size_t K_ARENA_PER_WORKER = (size_t)4 * (160L * FR + 600) + (size_t)4 * NM * (FR + 40) + 4096;
static int g_copy_mode;
static void k_global(void) {
  make_tables_400();
  BLASSetThreading(BLAS_THREADING_SINGLE_THREADED);
  fftwf_import_wisdom_from_filename("$BASELINES/fftw/wisdom-m4-fftwf-measure.txt");
  make_mel_tables(0);
  lctx c0; setup(&c0, K_B_NMF, FR);
  make_ref(FR);
  if (getenv("BREAK") && atoi(getenv("BREAK")) == 1) make_mel_tables(1); /* broken filterbank: every check must FAIL */
  g_copy_mode = getenv("COPY") ? atoi(getenv("COPY")) : 0;
}
static void *arena_get(size_t);
static void *k_make(int w) {
  lmw *m = calloc(1, sizeof *m);
  if (!setup(&m->n, K_B_GMF, FR) || !setup(&m->s, K_S_ZAF, FR)) { fprintf(stderr, "setup failed\n"); exit(1); }
  long need = 160L * (FR - 1) + 400 + 64;
  float *sig = arena_get(sizeof(float) * need); memcpy(sig, g_sig, sizeof(float) * need);
  m->n.sig = m->s.sig = sig;
  m->s.out = arena_get(sizeof(float) * NM * (FR + 34));
  m->priv_sig = bench_alloc(sizeof(float) * need); memcpy(m->priv_sig, g_sig, sizeof(float) * need);
  m->priv_out = bench_alloc(sizeof(float) * NM * (FR + 34));
  if (g_copy_mode) m->n.sig = m->priv_sig; /* the NEON path reads the private copy */
  return m;
}
static void k_owner_init(void) { tls_touch(); }
static void k_worker_init(void *c) { (void)c; tls_touch(); }
static void k_neon(void *c) { run_core(&((lmw *)c)->n); }
static void k_sme_part(void *c) { lmw *m = c; tls_touch(); sme_lm400(m->s.sig, 160, m->s.win, m->s.out, FR / 2); }
static void k_tail(void *c) {
  lmw *m = c;
  if (g_copy_mode) { lctx t = m->s; t.out = m->priv_out; flog_all(&t); } else flog_all(&m->s);
}
static void k_copy_in(void *c) { lmw *m = c; memcpy((void *)m->s.sig, m->priv_sig, sizeof(float) * (160L * (FR - 1) + 400)); }
static void k_copy_out(void *c) { lmw *m = c; memcpy(m->priv_out, m->s.out, sizeof(float) * NM * FR); }
static int k_check(void *c, int path, const char *lab) {
  lmw *m = c;
  const float *o = path ? (g_copy_mode ? m->priv_out : m->s.out) : m->n.out;
  return check_out(lab, o, FR);
}
#include "mpcore.h"
int main(int argc, char **argv) { setvbuf(stdout, NULL, _IOLBF, 0); bench_init(); return mp_main(argc, argv); }
