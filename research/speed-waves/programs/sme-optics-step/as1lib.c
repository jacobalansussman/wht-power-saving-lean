// as1lib.c -- as-1: the SME 2-D FFT / propagation step as a shared library for Python (ctypes).
// as1_new(N, H) copies H / N^2 (H = transfer function, complex128 N x N, C order); as1_prop(h, E, out): out = propagated
// field at -n (reversed order, E unchanged); as1_prop_nat(h, E, out): natural order (adds one reversal pass);
// as1_step(h, U): U <- propagated field at -n, in place (chain use); as1_fft2(h, x, y): y = fft2(x).
#define main prop_main
#include "prop.c"
#undef main
typedef struct { prop_t *P; double *W, *W2, *G; long N; } as1_t;
void *as1_new(long N, const double *H) {
  as1_t *h = calloc(1, sizeof *h); h->N = N; h->G = col_alloc(16 * N * N); h->W = col_alloc(16 * N * N); h->W2 = col_alloc(16 * N * N);
  double s = 1.0 / ((double)N * (double)N); for (long i = 0; i < 2 * N * N; i++) h->G[i] = H[i] * s;
  col_pfd = 4; h->P = prop_new(N, h->G, 0); return h;
}
void as1_set_pfd(int d) { col_pfd = d; }
void as1_set_fused(int f) { use_fused = f; }
void as1_prop(void *hv, const double *E, double *out) { as1_t *h = hv; ap0v2_fft(h->P->row, E, h->W, h->N); if (use_fused && h->P->R > 1) prop_cols_fused(h->P, h->W); else prop_cols(h->P, h->W); ap0v2_fft(h->P->row, h->W, out, h->N); }
static void reverse2d(const double *x, double *y, long N) {  // y[a][b] = x[-a][-b]
  for (long a = 0; a < N; a++) {
    const double *s = x + 2 * ((N - a) % N) * N; double *d = y + 2 * a * N;
    vst1q_f64(d, vld1q_f64(s));
    for (long b = 1; b < N; b++) vst1q_f64(d + 2 * b, vld1q_f64(s + 2 * (N - b)));
  }
}
void as1_prop_nat(void *hv, const double *E, double *out) { as1_t *h = hv; ap0v2_fft(h->P->row, E, h->W, h->N); if (use_fused && h->P->R > 1) prop_cols_fused(h->P, h->W); else prop_cols(h->P, h->W); ap0v2_fft(h->P->row, h->W, h->W2, h->N); reverse2d(h->W2, out, h->N); }
void as1_step(void *hv, double *U) { as1_t *h = hv; prop_step(h->P, U, h->W); }
void as1_fft2(void *hv, const double *x, double *y) { as1_t *h = hv; sme_fft2(h->P, x, y); }

// whole cluster from Python: 1 SME thread (the caller) + 3 NEON threads (PFFFT-d); E unchanged, out = field at -n
static int as1_mt_ready = 0;
void as1_mt_init(void *hv) { as1_t *h = hv; if (as1_mt_ready) return; if (getenv("KFRW")) mt_kfr = 1; if (getenv("TAIL")) mt_tail = atoi(getenv("TAIL")); NN = h->N; PP = h->P; Gt = h->G; GT = col_alloc(16 * h->N * h->N); transpose(Gt, GT, h->N);
  PF.N = h->N; PF.s = pffftd_new_setup((int)h->N, PFFFT_COMPLEX); PF.work = pffftd_aligned_malloc(16 * h->N); mt_init(3); as1_mt_ready = 1; }
void as1_prop_mt(void *hv, const double *E, double *out, int sme_on) { as1_t *h = hv;
  MT.U = (double *)E; MT.O = out; MT.W = h->W; long N = h->N;
  mt_phase(0, N / 8, sme_on); mt_phase(1, N / PP->SW, sme_on); mt_phase(2, N / 8, sme_on); }
