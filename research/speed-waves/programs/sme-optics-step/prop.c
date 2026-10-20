// prop.c -- as-1 (wave 14): complex128 2-D FFT and angular-spectrum propagation step on the M4's SME unit.
// Rows: ap-0's kernel (wave-13/ap-0/ap0lib2.c, numpy layout). Columns: colfft.h (no transpose: 8-column row pieces are
// the data operand). Propagation step with FORWARD transforms only:  U -> rows -> [per strip of SW columns, in L2:
// cols, x G, cols] -> rows.  Since FFT2(FFT2(y))[n] = N^2 y[-n], FFT2(G * FFT2(U)) with G = H/N^2 equals the propagated
// field at index -n (both axes reversed). H of the angular spectrum is symmetric (H[-k] = H[k]), so a chain of steps
// alternates natural / reversed order and is exact after an even number of steps (a real-space screen V between steps
// must then be applied as V[-n] on odd steps: a second precomputed table). Optional REV=1: the last row pass of the step
// takes ... (not built; see report).
// Baselines (C, one thread): PFFFT-double rows + blocked transpose + rows (+ G, rows, transpose, rows); FFTW 3.3.10
// double (scalar) 2-D plans, FFTW_MEASURE (time limit 60 s).
// usage: prop N check | prop N step | prop N fft2 | prop N parts
#include "../sme-fp64-batched-fft/ap0lib2.c"
#include "colfft.h"
#include "../timing-harness/bench.h"
#include <arm_neon.h>
#include <pffft_double.h>
#include <fftw3.h>
void *pk2_new(int n);
void *kd_new(long n); void kd_exec(void *h, const double *in, double *out);  // KFR complex<double> DFT (ar-0's wrapper, libkfrb.dylib; GPL-2, private)
void pk2_exec(void *p, double *c, int fwd, double fct);

typedef struct { long N, SW; void *row; colplan *cp, *cpL; double *s0, *s1, *s2, *s3, *twr, *twi; int R; const double *G; double *Gs; } prop_t;
static prop_t *prop_new(long N, const double *G, int broken) {
  prop_t *P = calloc(1, sizeof *P); P->N = N; P->G = G;
  P->SW = getenv("SW") ? atol(getenv("SW")) : (N <= 1024 ? 64 : 32);
  P->row = ap0v2_new(N); P->cp = col_new(N, broken);
  P->s0 = col_alloc(16 * N * P->SW); P->s1 = col_alloc(16 * N * P->SW); P->s2 = col_alloc(16 * N * P->SW); P->s3 = col_alloc(16 * N * P->SW);
  P->Gs = col_alloc(16 * N * N);  // G re-laid strip by strip (strip c0: rows of SW values, contiguous), so its reads stream
  for (long c0 = 0; c0 < N; c0 += P->SW) for (long r = 0; r < N; r++) memcpy(P->Gs + 2 * (c0 * N + r * P->SW), G + 2 * (r * N + c0), 16 * P->SW);
  P->cpL = col_new_o(N, broken, 1); P->R = P->cp->R[0] == 8 ? 1 : P->cp->R[0];
  if (P->R > 1) { long M = N / P->R; P->twr = col_alloc(8 * N); P->twi = col_alloc(8 * N);
    for (int i = 1; i < P->R; i++) for (long k = 0; k < M; k++) { long e = ((long)i * k) % N; double th = 2.0 * M_PI * (double)e / (double)N;
      P->twr[(i - 1) * M + k] = cos(th); P->twi[(i - 1) * M + k] = -sin(th); } }
  return P;
}
static void hmul(const double *Gq, long N, double *x, long sw) {  // Gq: this strip's G block (row stride sw)
  for (long r = 0; r < N; r++) {
    const double *g = Gq + 2 * r * sw; double *a = x + 2 * r * sw;
    for (long c = 0; c < sw; c += 2) {
      float64x2x2_t v = vld2q_f64(a + 2 * c), w = vld2q_f64(g + 2 * c), o;
      o.val[0] = vfmsq_f64(vmulq_f64(v.val[0], w.val[0]), v.val[1], w.val[1]);
      o.val[1] = vfmaq_f64(vmulq_f64(v.val[0], w.val[1]), v.val[1], w.val[0]);
      vst2q_f64(a + 2 * c, o);
    }
  }
}
static void prop_cols(prop_t *P, double *W) {
  long N = P->N, SW = P->SW;
  for (long c0 = 0; c0 < N; c0 += SW) {
    col_fft_ls(P->cp, W + 2 * c0, N, P->s2, SW, SW, P->s0, P->s1);
    hmul(P->Gs + 2 * c0 * N, N, P->s2, SW);
    col_fft_ls(P->cp, P->s2, SW, W + 2 * c0, N, SW, P->s0, P->s1);
  }
}
// Fused middle pass (NEON, strip in L2) for N = R * 8^a, R in {2, 4}: last radix-R stage of column FFT #1 (plan cpL, small
// radix last: y[k + o M] = sum_i w_R^{io} w_N^{ik} x[k + i M]), times G, first radix-R stage of column FFT #2 (plan cp:
// z[R k + o'] = sum_o w_R^{o o'} y'[k + o M]).  X, Z: strips (row stride sw); G rows at column offset c0.
#define CMULR(ar, ai, br, bi, orr, oi) do { orr = vfmsq_f64(vmulq_f64(ar, br), ai, bi); oi = vfmaq_f64(vmulq_f64(ar, bi), ai, br); } while (0)
static void mid_pass(const prop_t *P, const double *X, double *Z, long sw, long c0) {
  long N = P->N, M = N / P->R; const double *Gq = P->Gs + 2 * c0 * N;
  if (P->R == 2) {
    for (long k = 0; k < M; k++) {
      float64x2_t wr = vdupq_n_f64(P->twr[k]), wi = vdupq_n_f64(P->twi[k]);
      const double *x0 = X + 2 * k * sw, *x1 = X + 2 * (k + M) * sw, *g0 = Gq + 2 * k * sw, *g1 = Gq + 2 * (k + M) * sw;
      double *z0 = Z + 2 * (2 * k) * sw, *z1 = Z + 2 * (2 * k + 1) * sw;
      for (long c = 0; c < sw; c += 2) {
        float64x2x2_t a = vld2q_f64(x0 + 2 * c), b = vld2q_f64(x1 + 2 * c), ga = vld2q_f64(g0 + 2 * c), gb = vld2q_f64(g1 + 2 * c), o0, o1;
        float64x2_t tr, ti; CMULR(b.val[0], b.val[1], wr, wi, tr, ti);
        float64x2_t y0r = vaddq_f64(a.val[0], tr), y0i = vaddq_f64(a.val[1], ti), y1r = vsubq_f64(a.val[0], tr), y1i = vsubq_f64(a.val[1], ti);
        float64x2_t u0r, u0i, u1r, u1i; CMULR(y0r, y0i, ga.val[0], ga.val[1], u0r, u0i); CMULR(y1r, y1i, gb.val[0], gb.val[1], u1r, u1i);
        o0.val[0] = vaddq_f64(u0r, u1r); o0.val[1] = vaddq_f64(u0i, u1i); o1.val[0] = vsubq_f64(u0r, u1r); o1.val[1] = vsubq_f64(u0i, u1i);
        vst2q_f64(z0 + 2 * c, o0); vst2q_f64(z1 + 2 * c, o1);
      }
    }
  } else {  // R == 4
    for (long k = 0; k < M; k++) {
      float64x2_t w1r = vdupq_n_f64(P->twr[k]), w1i = vdupq_n_f64(P->twi[k]), w2r = vdupq_n_f64(P->twr[M + k]), w2i = vdupq_n_f64(P->twi[M + k]);
      float64x2_t w3r = vdupq_n_f64(P->twr[2 * M + k]), w3i = vdupq_n_f64(P->twi[2 * M + k]);
      const double *x[4], *g[4]; double *z[4];
      for (int i = 0; i < 4; i++) { x[i] = X + 2 * (k + i * M) * sw; g[i] = Gq + 2 * (k + i * M) * sw; z[i] = Z + 2 * (4 * k + i) * sw; }
      for (long c = 0; c < sw; c += 2) {
        float64x2x2_t a0 = vld2q_f64(x[0] + 2 * c), a1 = vld2q_f64(x[1] + 2 * c), a2 = vld2q_f64(x[2] + 2 * c), a3 = vld2q_f64(x[3] + 2 * c);
        float64x2_t t0r = a0.val[0], t0i = a0.val[1], t1r, t1i, t2r, t2i, t3r, t3i;
        CMULR(a1.val[0], a1.val[1], w1r, w1i, t1r, t1i); CMULR(a2.val[0], a2.val[1], w2r, w2i, t2r, t2i); CMULR(a3.val[0], a3.val[1], w3r, w3i, t3r, t3i);
        // forward DFT4: y0 = s02 + s13, y2 = s02 - s13, y1 = d02 - i d13, y3 = d02 + i d13
        float64x2_t s02r = vaddq_f64(t0r, t2r), s02i = vaddq_f64(t0i, t2i), d02r = vsubq_f64(t0r, t2r), d02i = vsubq_f64(t0i, t2i);
        float64x2_t s13r = vaddq_f64(t1r, t3r), s13i = vaddq_f64(t1i, t3i), d13r = vsubq_f64(t1r, t3r), d13i = vsubq_f64(t1i, t3i);
        float64x2_t yr[4], yi[4];
        yr[0] = vaddq_f64(s02r, s13r); yi[0] = vaddq_f64(s02i, s13i); yr[2] = vsubq_f64(s02r, s13r); yi[2] = vsubq_f64(s02i, s13i);
        yr[1] = vaddq_f64(d02r, d13i); yi[1] = vsubq_f64(d02i, d13r); yr[3] = vsubq_f64(d02r, d13i); yi[3] = vaddq_f64(d02i, d13r);
        for (int o = 0; o < 4; o++) { float64x2x2_t gg = vld2q_f64(g[o] + 2 * c); float64x2_t ur, ui; CMULR(yr[o], yi[o], gg.val[0], gg.val[1], ur, ui); yr[o] = ur; yi[o] = ui; }
        s02r = vaddq_f64(yr[0], yr[2]); s02i = vaddq_f64(yi[0], yi[2]); d02r = vsubq_f64(yr[0], yr[2]); d02i = vsubq_f64(yi[0], yi[2]);
        s13r = vaddq_f64(yr[1], yr[3]); s13i = vaddq_f64(yi[1], yi[3]); d13r = vsubq_f64(yr[1], yr[3]); d13i = vsubq_f64(yi[1], yi[3]);
        float64x2x2_t o0, o1, o2, o3;
        o0.val[0] = vaddq_f64(s02r, s13r); o0.val[1] = vaddq_f64(s02i, s13i); o2.val[0] = vsubq_f64(s02r, s13r); o2.val[1] = vsubq_f64(s02i, s13i);
        o1.val[0] = vaddq_f64(d02r, d13i); o1.val[1] = vsubq_f64(d02i, d13r); o3.val[0] = vsubq_f64(d02r, d13i); o3.val[1] = vaddq_f64(d02i, d13r);
        vst2q_f64(z[0] + 2 * c, o0); vst2q_f64(z[1] + 2 * c, o1); vst2q_f64(z[2] + 2 * c, o2); vst2q_f64(z[3] + 2 * c, o3);
      }
    }
  }
}
static int use_fused = 1;
static void prop_cols_fused(prop_t *P, double *W) {
  long N = P->N, SW = P->SW; int S = P->cp->S;
  for (long c0 = 0; c0 < N; c0 += SW) {
    col_fft_rls(P->cpL, 0, S - 1, W + 2 * c0, N, P->s2, SW, SW, P->s0, P->s1);
    mid_pass(P, P->s2, P->s3, SW, c0);
    col_fft_rls(P->cp, 1, S, P->s3, SW, W + 2 * c0, N, SW, P->s0, P->s1);
  }
}
static void prop_step(prop_t *P, double *U, double *W) { ap0v2_fft(P->row, U, W, P->N); if (use_fused && P->R > 1) prop_cols_fused(P, W); else prop_cols(P, W); ap0v2_fft(P->row, W, U, P->N); }
static void sme_cols(prop_t *P, double *W) { long N = P->N, SW = P->SW; for (long c0 = 0; c0 < N; c0 += SW) col_fft_ls(P->cp, W + 2 * c0, N, W + 2 * c0, N, SW, P->s0, P->s1); }
static void sme_fft2(prop_t *P, const double *U, double *W) { ap0v2_fft(P->row, U, W, P->N); sme_cols(P, W); }

// ---- PFFFT-double baseline ----
typedef struct { long N; PFFFTD_Setup *s; double *work; } pf_t;
static void pf_rows(pf_t *F, const double *x, double *y, int fwd) {
  for (long r = 0; r < F->N; r++) pffftd_transform_ordered(F->s, x + 2 * r * F->N, y + 2 * r * F->N, F->work, fwd ? PFFFT_FORWARD : PFFFT_BACKWARD);
}
static int tr_naive = 0;
static void transpose(const double *x, double *y, long N) {  // complex128 N x N; 32 x 32 blocks through a local tile (no L1 set aliasing)
  if (tr_naive == 1) { for (long i0 = 0; i0 < N; i0 += 32) for (long j0 = 0; j0 < N; j0 += 32) for (long i = i0; i < i0 + 32; i++) {
      const double *s = x + 2 * (i * N + j0); for (long j = 0; j < 32; j++) vst1q_f64(y + 2 * ((j0 + j) * N + i), vld1q_f64(s + 2 * j)); } return; }
  if (tr_naive == 2) {  // previous version: 32 x 32 blocks through a local tile
    static double tile[32 * 32 * 2] __attribute__((aligned(128)));
    for (long i0 = 0; i0 < N; i0 += 32) for (long j0 = 0; j0 < N; j0 += 32) {
      for (long i = 0; i < 32; i++) { const double *s = x + 2 * ((i0 + i) * N + j0); for (long j = 0; j < 32; j++) vst1q_f64(tile + 2 * (j * 32 + i), vld1q_f64(s + 2 * j)); }
      for (long j = 0; j < 32; j++) { double *d = y + 2 * ((j0 + j) * N + i0); for (long i = 0; i < 32; i++) vst1q_f64(d + 2 * i, vld1q_f64(tile + 2 * (j * 32 + i))); }
    }
    return;
  }
  // 8 source rows x 64 columns: the 8 rows fill the 8 ways of the L1 sets they share; each destination line (8 complex =
  // 128 B) is written whole, once.
  for (long i0 = 0; i0 < N; i0 += 8) for (long j0 = 0; j0 < N; j0 += 64) {
    const double *s0 = x + 2 * (i0 * N + j0);
    for (long j = 0; j < 64; j++) {
      double *d = y + 2 * ((j0 + j) * N + i0); const double *s = s0 + 2 * j;
      float64x2_t a0 = vld1q_f64(s), a1 = vld1q_f64(s + 2 * N), a2 = vld1q_f64(s + 4 * N), a3 = vld1q_f64(s + 6 * N);
      float64x2_t a4 = vld1q_f64(s + 8 * N), a5 = vld1q_f64(s + 10 * N), a6 = vld1q_f64(s + 12 * N), a7 = vld1q_f64(s + 14 * N);
      vst1q_f64(d, a0); vst1q_f64(d + 2, a1); vst1q_f64(d + 4, a2); vst1q_f64(d + 6, a3); vst1q_f64(d + 8, a4); vst1q_f64(d + 10, a5); vst1q_f64(d + 12, a6); vst1q_f64(d + 14, a7);
    }
  }
}
static void pf_step(pf_t *F, const double *G, double *U, double *W, double *T) {
  long N = F->N; pf_rows(F, U, W, 1); transpose(W, T, N); pf_rows(F, T, W, 1);
  for (long i = 0; i < N * N; i += 2) {  // W *= G^T (= G: square grid, equal spacing)
    float64x2x2_t v = vld2q_f64(W + 2 * i), w = vld2q_f64(G + 2 * i), o;
    o.val[0] = vfmsq_f64(vmulq_f64(v.val[0], w.val[0]), v.val[1], w.val[1]); o.val[1] = vfmaq_f64(vmulq_f64(v.val[0], w.val[1]), v.val[1], w.val[0]);
    vst2q_f64(W + 2 * i, o);
  }
  pf_rows(F, W, T, 0); transpose(T, W, N); pf_rows(F, W, U, 0);
}
static void pf_fft2(pf_t *F, const double *U, double *W, double *T) { long N = F->N; pf_rows(F, U, W, 1); transpose(W, T, N); pf_rows(F, T, W, 1); transpose(W, T, N); }

// ---- KFR baseline: rows + transpose, forward transforms only (same reversed-order trick as the SME step: same work) ----
static void *KD; static void kf_rows(const double *x, double *y, long N) { for (long r = 0; r < N; r++) kd_exec(KD, x + 2 * r * N, y + 2 * r * N); }
static void kf_step(const double *G, double *U, double *W, double *T, long N) {
  kf_rows(U, W, N); transpose(W, T, N); kf_rows(T, W, N);
  for (long i = 0; i < N * N; i += 2) { float64x2x2_t v = vld2q_f64(W + 2 * i), w = vld2q_f64(G + 2 * i), o;
    o.val[0] = vfmsq_f64(vmulq_f64(v.val[0], w.val[0]), v.val[1], w.val[1]); o.val[1] = vfmaq_f64(vmulq_f64(v.val[0], w.val[1]), v.val[1], w.val[0]); vst2q_f64(W + 2 * i, o); }
  kf_rows(W, T, N); transpose(T, W, N); kf_rows(W, U, N);
}
static void kf_fft2(const double *U, double *W, double *T, long N) { kf_rows(U, W, N); transpose(W, T, N); kf_rows(T, W, N); transpose(W, T, N); }
// ---- globals for the bench callbacks ----
/* GT = G transposed: H is NOT exactly symmetric under transpose ((c - fx^2) - fy^2 rounds differently) */ static long NN; static double *U, *W, *T, *Gt, *GT; static prop_t *PP; static pf_t PF; static fftw_plan FWf, FWb;
static void gmul(double *x, const double *g, long n) {
  for (long i = 0; i < n; i += 2) {
    float64x2x2_t v = vld2q_f64(x + 2 * i), w = vld2q_f64(g + 2 * i), o;
    o.val[0] = vfmsq_f64(vmulq_f64(v.val[0], w.val[0]), v.val[1], w.val[1]); o.val[1] = vfmaq_f64(vmulq_f64(v.val[0], w.val[1]), v.val[1], w.val[0]);
    vst2q_f64(x + 2 * i, o);
  }
}
static void r_sme_step(void *c, long i) { prop_step(PP, U, W); }
static void r_pf_step(void *c, long i) { pf_step(&PF, GT, U, W, T); }
static void r_fw_step(void *c, long i) { fftw_execute(FWf); gmul(W, Gt, NN * NN); fftw_execute(FWb); }
static void r_sme_fft2(void *c, long i) { sme_fft2(PP, U, W); }
static void r_kf_step(void *c, long i) { kf_step(GT, U, W, T, NN); }
static void r_kf_fft2(void *c, long i) { kf_fft2(U, W, T, NN); }
static void r_kf_rows(void *c, long i) { kf_rows(U, W, NN); }
static void r_pf_fft2(void *c, long i) { pf_fft2(&PF, U, W, T); }
static void r_fw_fft2(void *c, long i) { fftw_execute(FWf); }
static void r_sme_rows(void *c, long i) { ap0v2_fft(PP->row, U, W, NN); }
static void r_sme_cols(void *c, long i) { sme_cols(PP, W); }
static void r_sme_colsH(void *c, long i) { prop_cols(PP, W); }
static void r_sme_colsF(void *c, long i) { if (PP->R > 1) prop_cols_fused(PP, W); else prop_cols(PP, W); }
static void r_pf_rows(void *c, long i) { pf_rows(&PF, U, W, 1); }
static void r_transpose(void *c, long i) { transpose(U, W, NN); }
static void r_copy(void *c, long i) { memcpy(W, U, 16 * NN * NN); }
static void r_strip_l2(void *c, long i) { for (long q = 0; q < NN / PP->SW; q++) col_fft_ls(PP->cp, PP->s2, PP->SW, PP->s2, PP->SW, PP->SW, PP->s0, PP->s1); }


// ---- whole performance cluster: 1 SME thread + 3 NEON threads (PFFFT-d), dynamic work sharing (step 4 of the brief) ----
// Phases: rows U->W (units of 8 rows; the SME thread takes 8 units = 64 rows at a time), column strips on W (unit = strip of
// SW columns; NEON threads: gather the strip transposed, PFFFT-d forward rows, x G row (G symmetric), PFFFT-d forward rows,
// scatter back -- the same arithmetic as the SME strip), rows W->U. sme_on = 0: the main thread does NEON work too
// (= a 4-thread NEON baseline with the same plan).
#include <pthread.h>
#include <stdatomic.h>
typedef struct { _Atomic long next; long total; int kind, sme; _Atomic int gen, done; int nw, quit; double *U, *W, *O; double *wk[4], *t1[4], *t2[4]; } mt_t;
static mt_t MT; static int mt_sme_share_units = 8; static int mt_kfr = 0; static void *KDT[4];  // mt_kfr: NEON threads use KFR-d (one plan per thread)
static inline void nrow(int t, const double *x, double *y, double *work) { if (mt_kfr) kd_exec(KDT[t], x, y); else pffftd_transform_ordered(PF.s, x, y, work, PFFFT_FORWARD); }
static void neon_rows_unit(int t, long u, const double *x, double *y, double *work) {
  long N = NN; for (long r = u * 8; r < u * 8 + 8 && r < N; r++) nrow(t, x + 2 * r * N, y + 2 * r * N, work);
}
static void neon_strip(int t, long q, double *T1, double *T2, double *work) {
  long N = NN, SW = PP->SW, c0 = q * SW; double *W_ = MT.W;
  for (long r0 = 0; r0 < N; r0 += 16) for (long i = 0; i < SW; i++) for (long r = r0; r < r0 + 16; r++) vst1q_f64(T1 + 2 * (i * N + r), vld1q_f64(W_ + 2 * (r * N + c0 + i)));
  for (long i = 0; i < SW; i++) {
    nrow(t, T1 + 2 * i * N, T2 + 2 * i * N, work);
    gmul(T2 + 2 * i * N, GT + 2 * (c0 + i) * N, N);  // column c0+i of G = row c0+i of G^T
    nrow(t, T2 + 2 * i * N, T1 + 2 * i * N, work);
  }
  for (long r0 = 0; r0 < N; r0 += 16) for (long i = 0; i < SW; i++) for (long r = r0; r < r0 + 16; r++) vst1q_f64(W_ + 2 * (r * N + c0 + i), vld1q_f64(T1 + 2 * (i * N + r)));
}
static int mt_tail = 0;  // column phase: NEON threads leave the last mt_tail strips to the (2x faster) SME thread
static void mt_neon_work(int t) {
  for (;;) {
    if (t > 0 && MT.sme && MT.kind == 1 && mt_tail && atomic_load(&MT.next) >= MT.total - mt_tail) break;  // only when an SME thread is there to take them
    long u = atomic_fetch_add(&MT.next, 1); if (u >= MT.total) break;
    if (MT.kind == 0) neon_rows_unit(t, u, MT.U, MT.W, MT.wk[t]); else if (MT.kind == 2) neon_rows_unit(t, u, MT.W, MT.O, MT.wk[t]); else neon_strip(t, u, MT.t1[t], MT.t2[t], MT.wk[t]);
  }
}
static void *mt_worker(void *arg) {
  int t = (int)(long)arg; pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0); int seen = 0;
  for (;;) { while (atomic_load(&MT.gen) == seen) { if (MT.quit) return NULL; __asm__ volatile("yield"); } seen = atomic_load(&MT.gen);
    mt_neon_work(t); atomic_fetch_add(&MT.done, 1); }
}
static void mt_init(int nw) {
  long N = NN; MT.nw = nw; for (int t = 0; t < 4; t++) KDT[t] = kd_new(N);
  for (int t = 0; t < 4; t++) { MT.wk[t] = col_alloc(16 * N); MT.t1[t] = col_alloc(16 * N * PP->SW); MT.t2[t] = col_alloc(16 * N * PP->SW); }
  for (int t = 1; t <= nw; t++) { pthread_t th; pthread_create(&th, NULL, mt_worker, (void *)(long)t); }
}
static void mt_phase(int kind, long total, int sme_on) {
  MT.kind = kind; MT.sme = sme_on; MT.total = total; atomic_store(&MT.done, 0); atomic_store(&MT.next, 0); atomic_fetch_add(&MT.gen, 1);
  long N = NN, SW = PP->SW;
  if (!sme_on) mt_neon_work(0);
  else for (;;) {
    long rem = total - atomic_load(&MT.next), take = kind == 1 ? 1 : mt_sme_share_units;
    if (mt_tail && kind != 1) { long t2 = rem / 8; if (t2 < take) take = t2 < 1 ? 1 : t2; }
    long u = atomic_fetch_add(&MT.next, take); if (u >= total) break; long ue = u + take < total ? u + take : total;
    if (kind == 0) ap0v2_fft(PP->row, MT.U + 2 * u * 8 * N, MT.W + 2 * u * 8 * N, (ue - u) * 8);
    else if (kind == 2) ap0v2_fft(PP->row, MT.W + 2 * u * 8 * N, MT.O + 2 * u * 8 * N, (ue - u) * 8);
    else { long c0 = u * SW;
      if (use_fused && PP->R > 1) { int S = PP->cp->S; col_fft_rls(PP->cpL, 0, S - 1, MT.W + 2 * c0, N, PP->s2, SW, SW, PP->s0, PP->s1); mid_pass(PP, PP->s2, PP->s3, SW, c0); col_fft_rls(PP->cp, 1, S, PP->s3, SW, MT.W + 2 * c0, N, SW, PP->s0, PP->s1); }
      else { col_fft_ls(PP->cp, MT.W + 2 * c0, N, PP->s2, SW, SW, PP->s0, PP->s1); hmul(PP->Gs + 2 * c0 * N, N, PP->s2, SW); col_fft_ls(PP->cp, PP->s2, SW, MT.W + 2 * c0, N, SW, PP->s0, PP->s1); } }
  }
  while (atomic_load(&MT.done) < MT.nw) __asm__ volatile("yield");
}
static void mt_step(double *U_, double *W_, int sme_on) {
  MT.U = U_; MT.O = U_; MT.W = W_; long N = NN;
  mt_phase(0, N / 8, sme_on); mt_phase(1, N / PP->SW, sme_on); mt_phase(2, N / 8, sme_on);
}
static void r_mt_sme(void *c, long i) { mt_step(U, W, 1); }
static void r_mt_neon(void *c, long i) { mt_step(U, W, 0); }
// ---- reference and checks ----
static void ref_fft2(double *x, long N, int fwd, double fct) {  // pocketfft rows then columns
  void *pk = pk2_new((int)N); double *col = col_alloc(16 * N);
  for (long r = 0; r < N; r++) pk2_exec(pk, x + 2 * r * N, fwd, 1.0);
  for (long c = 0; c < N; c++) {
    for (long r = 0; r < N; r++) { col[2 * r] = x[2 * (r * N + c)]; col[2 * r + 1] = x[2 * (r * N + c) + 1]; }
    pk2_exec(pk, col, fwd, fct);
    for (long r = 0; r < N; r++) { x[2 * (r * N + c)] = col[2 * r]; x[2 * (r * N + c) + 1] = col[2 * r + 1]; }
  }
  free(col);
}
static void two_sum(double a, double b, double *s, double *e) { *s = a + b; double z = *s - a; *e = (a - (*s - z)) + (b - z); }
typedef struct { double r, re, i, ie; } ddc;
static void dd_acc(ddc *a, double xr, double xi, double c, double s) {  // a += (xr + i xi)(c + i s), compensated
  double t0 = xr * c, t1 = -xi * s, t2 = xr * s, t3 = xi * c, t, e;
  two_sum(a->r, t0, &t, &e); a->r = t; a->re += e + fma(xr, c, -t0);
  two_sum(a->r, t1, &t, &e); a->r = t; a->re += e + fma(-xi, s, -t1);
  two_sum(a->i, t2, &t, &e); a->i = t; a->ie += e + fma(xr, s, -t2);
  two_sum(a->i, t3, &t, &e); a->i = t; a->ie += e + fma(xi, c, -t3);
}
static void dd_bin(const double *x, long N, long k1, long k2, double *Xr, double *Xi) {  // 2-D DFT bin, compensated sums
  double *rr = malloc(8 * N), *ri = malloc(8 * N);
  for (long n1 = 0; n1 < N; n1++) {
    ddc a = {0, 0, 0, 0};
    for (long n2 = 0; n2 < N; n2++) { long e = (n2 * k2) % N; double th = -2.0 * M_PI * (double)e / (double)N; dd_acc(&a, x[2 * (n1 * N + n2)], x[2 * (n1 * N + n2) + 1], cos(th), sin(th)); }
    rr[n1] = a.r + a.re; ri[n1] = a.i + a.ie;
  }
  ddc b = {0, 0, 0, 0};
  for (long n1 = 0; n1 < N; n1++) { long e = (n1 * k1) % N; double th = -2.0 * M_PI * (double)e / (double)N; dd_acc(&b, rr[n1], ri[n1], cos(th), sin(th)); }
  *Xr = b.r + b.re; *Xi = b.i + b.ie; free(rr); free(ri);
}
static double relerr(const double *a, const double *b, long n2, double *rms) {  // max|a-b| / max|b|, and rel RMS
  double em = 0, rm = 0, se = 0, sr = 0;
  for (long i = 0; i < n2; i++) { double e = fabs(a[i] - b[i]); if (!(e <= em)) em = e; if (fabs(b[i]) > rm) rm = fabs(b[i]); se += e * e; sr += b[i] * b[i]; }
  *rms = sqrt(se / sr); return em / rm;
}
static void make_G(double *G, long N) {  // angular spectrum, as scout-c0's pipes2.py: dx = 10 um, lambda = 633 nm, z = 0.1 m; G = H / N^2
  double lam = 633e-9, dx = 1e-5, z = 0.1;
  for (long a = 0; a < N; a++) for (long b = 0; b < N; b++) {
    double fx = (double)(a < (N + 1) / 2 ? a : a - N) / (N * dx), fy = (double)(b < (N + 1) / 2 ? b : b - N) / (N * dx);
    double q = 1.0 / (lam * lam) - fx * fx - fy * fy; double ph = 2.0 * M_PI * z * sqrt(q > 0 ? q : 0);
    G[2 * (a * N + b)] = cos(ph) / ((double)N * N); G[2 * (a * N + b) + 1] = sin(ph) / ((double)N * N);
  }
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: prop N check|step|fft2|parts\n"); return 1; }
  bench_init(); long N = atol(argv[1]); col_pfd = getenv("PFD") ? atoi(getenv("PFD")) : 4; if (getenv("TRNAIVE")) tr_naive = atoi(getenv("TRNAIVE")); if (getenv("NOFUSE")) use_fused = 0; NN = N; const char *mode = argv[2];
  size_t bytes = 16 * N * N;
  U = col_alloc(bytes); W = col_alloc(bytes); T = col_alloc(bytes); Gt = col_alloc(bytes);
  make_G(Gt, N); GT = col_alloc(bytes); transpose(Gt, GT, N); bench_fill_d(U, 2 * N * N);
  PP = prop_new(N, Gt, getenv("BROKEN") != NULL);
  printf("# prop N=%ld mode=%s SW=%ld PFD=%d col radices:", N, mode, PP->SW, col_pfd); for (int s = 0; s < PP->cp->S; s++) printf(" %d", PP->cp->R[s]); printf("\n");
  KD = kd_new(N);
  PF.N = N; PF.s = pffftd_new_setup((int)N, PFFFT_COMPLEX); PF.work = pffftd_aligned_malloc(16 * N);
  int want_fftw = !getenv("NOFFTW") && strcmp(mode, "parts") && strcmp(mode, "check") && !(getenv("NOFFTW"));
  if (want_fftw) {
    fftw_set_timelimit(getenv("FFTW_TL") ? atof(getenv("FFTW_TL")) : 60.0);
    const char *wf = getenv("FFTW_WF") ? getenv("FFTW_WF") : "fftw-d-wisdom.txt"; fftw_import_wisdom_from_filename(wf);
    uint64_t t0 = bench_now_ns();
    unsigned fl = getenv("FFTW_EST") ? FFTW_ESTIMATE : FFTW_MEASURE;
    FWf = fftw_plan_dft_2d((int)N, (int)N, (fftw_complex *)U, (fftw_complex *)W, FFTW_FORWARD, fl);
    FWb = fftw_plan_dft_2d((int)N, (int)N, (fftw_complex *)W, (fftw_complex *)U, FFTW_BACKWARD, fl);
    fftw_export_wisdom_to_filename(wf);
    printf("# FFTW plans (%s) took %.1f s\n", fl == FFTW_ESTIMATE ? "ESTIMATE" : "MEASURE, time limit", (bench_now_ns() - t0) * 1e-9);
    bench_fill_d(U, 2 * N * N);
  }
  if (!strcmp(mode, "plan")) return 0;
  if (!strcmp(mode, "mt") || !strcmp(mode, "mtcheck")) {
    if (getenv("SMEUNITS")) mt_sme_share_units = atoi(getenv("SMEUNITS")); if (getenv("KFRW")) mt_kfr = 1; if (getenv("TAIL")) mt_tail = atoi(getenv("TAIL"));
    printf("# cluster: NEON threads use %s; SME takes %d row units (8 rows each) at a time; tail %d\n", mt_kfr ? "KFR-d" : "PFFFT-d", mt_sme_share_units, mt_tail);
    mt_init(3); double *R = col_alloc(bytes), *X = col_alloc(bytes);
    memcpy(R, U, bytes); ref_fft2(R, N, 1, 1.0); gmul(R, Gt, N * N); ref_fft2(R, N, 0, 1.0);
    for (int v = 0; v < 2; v++) {
      memcpy(X, U, bytes); double *Us = U; U = X; mt_step(U, W, v == 0); U = Us;
      double em = 0, rm = 0; for (long a = 0; a < N; a++) for (long b = 0; b < N; b++) for (int c = 0; c < 2; c++) {
        double d = fabs(X[2 * (((N - a) % N) * N + (N - b) % N) + c] - R[2 * (a * N + b) + c]); if (!(d <= em)) em = d; rm = fmax(rm, fabs(R[2 * (a * N + b) + c])); }
      printf("CHECK cluster step (%s) read at -n vs pocketfft: max rel %.3e  %s\n", v == 0 ? "1 SME + 3 NEON" : "4 NEON", em / rm, em / rm < 1e-14 ? "PASS" : "FAIL");
    }
    free(R); free(X);
    if (!strcmp(mode, "mtcheck")) { MT.quit = 1; return 0; }
    bench_opts o = BENCH_OPTS_DEFAULT; o.units = 1e6; o.unit_name = "call(ms)"; bench_stat st[3]; char title[128];
    bench_prog pr[3] = {{"cluster 1 SME + 3 NEON", r_mt_sme, NULL, NULL}, {"cluster 4 NEON (same plan)", r_mt_neon, NULL, NULL}, {"1 SME thread (others spin)", r_sme_step, NULL, NULL}};
    snprintf(title, sizeof title, "cluster propagation step %ld^2 c128", N); bench_compare(title, pr, 3, &o, st);
    MT.quit = 1; return 0;
  }
  if (!strcmp(mode, "check")) {
    double *X = col_alloc(bytes), *R = col_alloc(bytes), rms;
    memcpy(R, U, bytes); ref_fft2(R, N, 1, 1.0);
    sme_fft2(PP, U, X);
    double e = relerr(X, R, 2 * N * N, &rms);
    printf("CHECK fft2 SME vs pocketfft (all %ld bins): max rel %.3e  rel RMS %.3e  %s\n", N * N, e, rms, e < 1e-14 ? "PASS" : "FAIL");
    pf_fft2(&PF, U, W, T);
    e = relerr(T, R, 2 * N * N, &rms); printf("CHECK fft2 PFFFT-d rows+transpose vs pocketfft: max rel %.3e  rel RMS %.3e  %s\n", e, rms, e < 1e-14 ? "PASS" : "FAIL");
    kf_fft2(U, W, T, N); e = relerr(T, R, 2 * N * N, &rms); printf("CHECK fft2 KFR-d rows+transpose vs pocketfft: max rel %.3e  rel RMS %.3e  %s\n", e, rms, e < 1e-14 ? "PASS" : "FAIL");
    // double-double DFT on sampled bins
    long ks[8][2] = {{0, 0}, {1, 0}, {0, 1}, {N / 2, N / 3}, {N - 1, 7}, {N / 3, N - 1}, {5, N / 2 + 3}, {N - 3, N - 5}};
    int nb = N >= 4096 ? 4 : 8; double esme = 0, epk = 0, xm = 0, s2s = 0, s2p = 0, s2x = 0;
    for (int b = 0; b < nb; b++) {
      long k1 = ks[b][0], k2 = ks[b][1]; double Xr, Xi; dd_bin(U, N, k1, k2, &Xr, &Xi);
      double ds = hypot(X[2 * (k1 * N + k2)] - Xr, X[2 * (k1 * N + k2) + 1] - Xi), dp = hypot(R[2 * (k1 * N + k2)] - Xr, R[2 * (k1 * N + k2) + 1] - Xi);
      esme = fmax(esme, ds); epk = fmax(epk, dp); s2s += ds * ds; s2p += dp * dp; s2x += Xr * Xr + Xi * Xi; xm = fmax(xm, hypot(Xr, Xi));
    }
    printf("CHECK fft2 vs double-double DFT on %d bins: SME rel RMS %.3e (max %.3e)  pocketfft rel RMS %.3e (max %.3e)  ratio %.2f  %s\n", nb,
           sqrt(s2s / s2x), esme / xm, sqrt(s2p / s2x), epk / xm, sqrt(s2s / s2p), sqrt(s2s / s2x) < 1e-14 ? "PASS" : "FAIL");
    // propagation step: SME (forward-only, reversed output) vs pocketfft ifft2(fft2(U) * H)
    memcpy(R, U, bytes); ref_fft2(R, N, 1, 1.0); gmul(R, Gt, N * N); ref_fft2(R, N, 0, 1.0);
    memcpy(X, U, bytes); prop_step(PP, X, W);
    double em = 0, rm = 0, se = 0, sr = 0;
    for (long a = 0; a < N; a++) for (long b = 0; b < N; b++) for (int c = 0; c < 2; c++) {
      double v = X[2 * (((N - a) % N) * N + (N - b) % N) + c], r = R[2 * (a * N + b) + c], d = fabs(v - r);
      if (!(d <= em)) em = d; rm = fmax(rm, fabs(r)); se += d * d; sr += r * r; }
    printf("CHECK propagation step SME (output read at -n) vs pocketfft ifft2(fft2(U)*H): max rel %.3e  rel RMS %.3e  %s\n", em / rm, sqrt(se / sr), em / rm < 1e-14 ? "PASS" : "FAIL");
    memcpy(X, U, bytes); pf_step(&PF, GT, X, W, T);
    e = relerr(X, R, 2 * N * N, &rms); printf("CHECK propagation step PFFFT-d (natural order) vs pocketfft: max rel %.3e  rel RMS %.3e  %s\n", e, rms, e < 1e-14 ? "PASS" : "FAIL");
    return 0;
  }
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = 1e6; o.unit_name = "call(ms)"; bench_stat st[10]; char title[128];
  if (!strcmp(mode, "step")) {
    bench_prog pr[4] = {{"SME step (rev. order)", r_sme_step, NULL, NULL}, {"PFFFT-d rows+transp", r_pf_step, NULL, NULL}, {"KFR-d rows+transp (rev.)", r_kf_step, NULL, NULL}, {"FFTW-d 2D MEASURE", r_fw_step, NULL, NULL}};
    snprintf(title, sizeof title, "propagation step %ld^2 c128", N); bench_compare(title, pr, want_fftw ? 4 : 3, &o, st);
  } else if (!strcmp(mode, "fft2")) {
    bench_prog pr[4] = {{"SME fft2", r_sme_fft2, NULL, NULL}, {"PFFFT-d fft2 (2 transp)", r_pf_fft2, NULL, NULL}, {"KFR-d fft2 (2 transp)", r_kf_fft2, NULL, NULL}, {"FFTW-d fft2", r_fw_fft2, NULL, NULL}};
    snprintf(title, sizeof title, "fft2 %ld^2 c128", N); bench_compare(title, pr, want_fftw ? 4 : 3, &o, st);
  } else if (!strcmp(mode, "cols")) {  // column phase only (tuning)
    bench_prog pr[3] = {{"SME cols (strips)", r_sme_cols, NULL, NULL}, {"SME cols+G+cols", r_sme_colsH, NULL, NULL}, {"SME cols, strip held in L2", r_strip_l2, NULL, NULL}};
    snprintf(title, sizeof title, "cols %ld^2 c128 PFD=%d SW=%ld", N, col_pfd, PP->SW); bench_compare(title, pr, 3, &o, st);
  } else if (!strcmp(mode, "parts")) {
    bench_prog pr[9] = {{"KFR-d rows pass", r_kf_rows, NULL, NULL}, {"SME rows pass", r_sme_rows, NULL, NULL}, {"SME cols (strips)", r_sme_cols, NULL, NULL}, {"SME cols+G+cols", r_sme_colsH, NULL, NULL},
                        {"PFFFT-d rows pass", r_pf_rows, NULL, NULL}, {"transpose", r_transpose, NULL, NULL}, {"memcpy", r_copy, NULL, NULL},
                        {"SME cols, strip held in L2", r_strip_l2, NULL, NULL}, {"SME cols+G+cols FUSED mid", r_sme_colsF, NULL, NULL}};
    snprintf(title, sizeof title, "parts %ld^2 c128", N); bench_compare(title, pr, 9, &o, st);
  }
  return 0;
}
