/* fc64.c -- av-1 (wave 15): fp64 / complex128 direct convolution and correlation on the M4's SME2 matrix unit,
 * a port of aj-1's fp32 Toeplitz-tile FIR (wave 11) to 8 x 8 fp64 FMOPA tiles (all 8 ZA.D tiles).
 *
 * One primitive (correlation of a zero-padded signal, any output range of the "full" result):
 *     y[t] = sum_{k<L} g[k] * xp[o0 + t + k],  t = 0..nout-1,
 *     xp[i] = x[i - (L-1)] for 0 <= i-(L-1) < nx, else 0.
 *   np.convolve(a, v, full) : g = v reversed, o0 = 0, nout = n + m - 1 (same / valid = sub-ranges)
 *   np.correlate(a, v, full): g = conj(v)
 * Paths: 1 = SME (Toeplitz tiles), 2 = NEON direct, 3 = scalar (reference order), 0 = auto.
 *
 * SME kernel (fp64): polyphase input P[p][r] = xp[s0 + 8 r + p] (p = 0..7, transposed through ZA in streaming
 *   mode); Toeplitz table T[c][i] = g[c - i] (0 outside), c = 0..Cn-1, Cn = roundup4(L + 7).
 *   ZA.D tile t (t = 0..7), pass rows R0..R0+63:  ZA_t[i][j] += T[c][i] * P[c % 8][R0 + c/8 + 8t + j]
 *   gives ZA_t[i][j] = y[8 (R0 + 8t + j) + i]; vertical slice j = 8 consecutive outputs. 512 outputs per pass.
 *   Per 4 taps-columns: 1 x4 load of T, 8 x4 loads of P, 32 FMOPA (64 MACs each).
 * Complex128 kernel: tiles 0..3 hold Re y, 4..7 hold Im y (256 complex outputs per pass); per c: 16 FMOPA
 *   (Re += Tr x Pr - Ti x Pi, Im += Ti x Pr + Tr x Pi; -Ti is a precomputed table).
 * Build (dylib for ctypes):
 *   clang -O3 -mcpu=apple-m4 -shared -fPIC fc64.c -framework Accelerate -o libfc64.dylib
 */
#include <Accelerate/Accelerate.h>
#include <pthread.h>
#include <arm_neon.h>
#include <arm_sme.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int fc_broken = 0; /* 1 = deliberately shift one Toeplitz row (the check must FAIL) */
long fc_chunk = 8192; /* outputs per chunk (rounded to 512) */
long fc_neon_max = 0; /* auto path: use NEON when nout * L <= this (0 = never) */
long fc_neon_maxn = 0; /* auto path: use NEON when nout <= this */

/* ---------------- scratch (per thread) ---------------- */
/* Each thread's buffers live in one heap block owned by a pthread key, so they are freed when that thread exits,
   whichever thread it is (a caller's or a corr_mt helper's). They used to be static __thread buffers, which nothing
   freed: every multi-thread call leaked up to about 1 MB, and on 12 October 2026 that filled memory and panicked the
   Mac. (A key destructor cannot free __thread buffers instead: on macOS they are already gone when it runs.) */
typedef struct { double *p; size_t n; } buf_t;
typedef struct { buf_t P, T, E, Y, S; } scratch_t;
static pthread_key_t skey;
static int skey_ok;
static pthread_once_t sonce = PTHREAD_ONCE_INIT;
static void scratch_free(void *v) {
  scratch_t *s = v;
  free(s->P.p); free(s->T.p); free(s->E.p); free(s->Y.p); free(s->S.p);
  free(s);
}
static void scratch_key(void) { skey_ok = pthread_key_create(&skey, scratch_free) == 0; }
static scratch_t *scratch(void) { /* this thread's block, NULL if out of memory */
  pthread_once(&sonce, scratch_key);
  if (!skey_ok) return NULL;
  scratch_t *s = pthread_getspecific(skey);
  if (!s && (s = calloc(1, sizeof *s)) && pthread_setspecific(skey, s)) { free(s); s = NULL; }
  return s;
}
static double *getbuf(buf_t *b, size_t n) {
  if (b->n < n) {
    free(b->p);
    size_t nn = n + n / 4 + 64;
    if (posix_memalign((void **)&b->p, 128, nn * sizeof(double))) { b->p = NULL; b->n = 0; return NULL; } /* p was freed */
    b->n = nn;
  }
  return b->p;
}

/* ---------------- SME: transpose through ZA ---------------- */
#define TRB(t)                                                                       \
  {                                                                                  \
    svfloat64x4_t v0_ = svld1_f64_x4(pc, xr), v1_ = svld1_f64_x4(pc, xr + 32);      \
    if (chk) { mx = svmax_f64_x4(mx, svmax_f64_x4(v0_, v1_)); mn = svmin_f64_x4(mn, svmin_f64_x4(v0_, v1_)); } \
    svwrite_hor_za64_f64_vg4(t, 0, v0_);                                             \
    svwrite_hor_za64_f64_vg4(t, 4, v1_);                                             \
    svst1_ver_za64(t, 0, pg, P + 0 * rs + r0); svst1_ver_za64(t, 1, pg, P + 1 * rs + r0); \
    svst1_ver_za64(t, 2, pg, P + 2 * rs + r0); svst1_ver_za64(t, 3, pg, P + 3 * rs + r0); \
    svst1_ver_za64(t, 4, pg, P + 4 * rs + r0); svst1_ver_za64(t, 5, pg, P + 5 * rs + r0); \
    svst1_ver_za64(t, 6, pg, P + 6 * rs + r0); svst1_ver_za64(t, 7, pg, P + 7 * rs + r0); \
  }
/* P[p*rs + r] = x[8 r + p], r in [0, rows), rows multiple of 8 */
static void tr_body(const double *x, double *P, long rs, long rows, int *bad) __arm_streaming __arm_inout("za") {
  svbool_t pg = svptrue_b64();
  svcount_t pc = svptrue_c64();
  int chk = bad != NULL;
  svfloat64_t z0 = svdup_f64(0.0);
  svfloat64x4_t mx = svcreate4(z0, z0, z0, z0), mn = mx;
  long r0 = 0;
  for (; r0 + 32 <= rows;) {
    const double *xr = x + 8 * r0;
    TRB(0) r0 += 8; xr += 64;
    TRB(1) r0 += 8; xr += 64;
    TRB(2) r0 += 8; xr += 64;
    TRB(3) r0 += 8;
  }
  for (; r0 < rows; r0 += 8) { const double *xr = x + 8 * r0; TRB(0) }
  if (chk) { /* FMAX / FMIN propagate NaN; +-inf survive: the window is finite iff -inf < min and max < inf */
    svfloat64_t a = svmax_f64_x(pg, svmax_f64_x(pg, svget4(mx, 0), svget4(mx, 1)), svmax_f64_x(pg, svget4(mx, 2), svget4(mx, 3)));
    svfloat64_t b = svmin_f64_x(pg, svmin_f64_x(pg, svget4(mn, 0), svget4(mn, 1)), svmin_f64_x(pg, svget4(mn, 2), svget4(mn, 3)));
    double hi = svmaxv_f64(pg, a), lo = svminv_f64(pg, b);
    *bad = !(hi < __builtin_inf() && lo > -__builtin_inf());
  }
}

/* complex: Pr[p][r] = Re x[8r+p], Pi[p][r] = Im x[8r+p]; x interleaved (re, im) */
#define TRC(tr, ti)                                                                   \
  {                                                                                   \
    for (int s = 0; s < 8; s++) {                                                     \
      svfloat64x2_t v = svld2_f64(pg, xr + 16 * s);                                   \
      svwrite_hor_za64_f64_m(tr, s, pg, svget2(v, 0));                                \
      svwrite_hor_za64_f64_m(ti, s, pg, svget2(v, 1));                                \
    }                                                                                 \
    for (int p = 0; p < 8; p++) {                                                     \
      svst1_ver_za64(tr, p, pg, Pr + p * rs + r0);                                    \
      svst1_ver_za64(ti, p, pg, Pi + p * rs + r0);                                    \
    }                                                                                 \
  }
static void trc_body(const double *x, double *Pr, double *Pi, long rs, long rows) __arm_streaming __arm_inout("za") {
  svbool_t pg = svptrue_b64();
  long r0 = 0;
  for (; r0 + 16 <= rows;) {
    const double *xr = x + 16 * r0;
    TRC(0, 1) r0 += 8; xr += 128;
    TRC(2, 3) r0 += 8;
  }
  for (; r0 < rows; r0 += 8) { const double *xr = x + 16 * r0; TRC(0, 1) }
}

/* ---------------- SME: real kernel ---------------- */
#define MOPA8(u)                                                                 \
  {                                                                              \
    const double *pp = P + (p + u) * rs + R0 + q;                                \
    svfloat64x4_t B0 = svld1_f64_x4(pc, pp), B1 = svld1_f64_x4(pc, pp + 32);     \
    svfloat64_t a = svget4(A, u);                                                \
    svmopa_za64_f64_m(0, pg, pg, a, svget4(B0, 0));                              \
    svmopa_za64_f64_m(1, pg, pg, a, svget4(B0, 1));                              \
    svmopa_za64_f64_m(2, pg, pg, a, svget4(B0, 2));                              \
    svmopa_za64_f64_m(3, pg, pg, a, svget4(B0, 3));                              \
    svmopa_za64_f64_m(4, pg, pg, a, svget4(B1, 0));                              \
    svmopa_za64_f64_m(5, pg, pg, a, svget4(B1, 1));                              \
    svmopa_za64_f64_m(6, pg, pg, a, svget4(B1, 2));                              \
    svmopa_za64_f64_m(7, pg, pg, a, svget4(B1, 3));                              \
  }
#define OUT8(t)                                                                  \
  {                                                                              \
    svst1_f64_x4(pc, y + 8 * (R0 + 8 * (t) + 0), svread_ver_za64_f64_vg4(t, 0)); \
    svst1_f64_x4(pc, y + 8 * (R0 + 8 * (t) + 4), svread_ver_za64_f64_vg4(t, 4)); \
  }
/* rows: multiple of 64 */
static void k_body(const double *P, long rs, const double *T, int Cn, double *y, long rows) __arm_streaming __arm_inout("za") {
  svbool_t pg = svptrue_b64();
  svcount_t pc = svptrue_c64();
  for (long R0 = 0; R0 < rows; R0 += 64) {
    svzero_za();
    for (int c = 0; c < Cn; c += 4) {
      int q = c >> 3, p = c & 7;
      svfloat64x4_t A = svld1_f64_x4(pc, T + 8 * c);
      MOPA8(0) MOPA8(1) MOPA8(2) MOPA8(3)
    }
    OUT8(0) OUT8(1) OUT8(2) OUT8(3) OUT8(4) OUT8(5) OUT8(6) OUT8(7)
  }
}

/* Blocked (two-level) accumulation for accuracy: tiles 0..3 accumulate a block of fc_blk tap-columns for 256 outputs
   (rows R0..R0+31), then are folded EXACTLY into the long accumulators 4..7 by unit-vector outer products
   (ZA_{t+4}[i][.] += e_i (x) ZA_t[i][.], 8 FMOPA per tile). Error: chains of fc_blk taps + L/fc_blk block sums,
   instead of one chain of L. Cost: 32 extra FMOPA per block (8 / fc_blk relative). rows: multiple of 32. */
long fc_blk = 0; /* tap-columns per block (multiple of 4); 0 = plain one-chain kernel */
static const double EYE8[64] = {1,0,0,0,0,0,0,0, 0,1,0,0,0,0,0,0, 0,0,1,0,0,0,0,0, 0,0,0,1,0,0,0,0,
                                0,0,0,0,1,0,0,0, 0,0,0,0,0,1,0,0, 0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,1};
#define MOPA4b(u)                                                                \
  {                                                                              \
    svfloat64x4_t B0 = svld1_f64_x4(pc, P + (p + u) * rs + R0 + q);              \
    svfloat64_t a = svget4(A, u);                                                \
    svmopa_za64_f64_m(0, pg, pg, a, svget4(B0, 0));                              \
    svmopa_za64_f64_m(1, pg, pg, a, svget4(B0, 1));                              \
    svmopa_za64_f64_m(2, pg, pg, a, svget4(B0, 2));                              \
    svmopa_za64_f64_m(3, pg, pg, a, svget4(B0, 3));                              \
  }
#define FOLD(t)                                                                  \
  {                                                                              \
    svfloat64x4_t h0 = svread_hor_za64_f64_vg4(t, 0), h1 = svread_hor_za64_f64_vg4(t, 4); \
    svmopa_za64_f64_m(t + 4, pg, pg, svget4(E0, 0), svget4(h0, 0));             \
    svmopa_za64_f64_m(t + 4, pg, pg, svget4(E0, 1), svget4(h0, 1));             \
    svmopa_za64_f64_m(t + 4, pg, pg, svget4(E0, 2), svget4(h0, 2));             \
    svmopa_za64_f64_m(t + 4, pg, pg, svget4(E0, 3), svget4(h0, 3));             \
    svmopa_za64_f64_m(t + 4, pg, pg, svget4(E1, 0), svget4(h1, 0));             \
    svmopa_za64_f64_m(t + 4, pg, pg, svget4(E1, 1), svget4(h1, 1));             \
    svmopa_za64_f64_m(t + 4, pg, pg, svget4(E1, 2), svget4(h1, 2));             \
    svmopa_za64_f64_m(t + 4, pg, pg, svget4(E1, 3), svget4(h1, 3));             \
  }
#define OUTb(t)                                                                  \
  {                                                                              \
    svst1_f64_x4(pc, y + 8 * (R0 + 8 * (t) + 0), svread_ver_za64_f64_vg4(t + 4, 0)); \
    svst1_f64_x4(pc, y + 8 * (R0 + 8 * (t) + 4), svread_ver_za64_f64_vg4(t + 4, 4)); \
  }
static void kb_body(const double *P, long rs, const double *T, int Cn, double *y, long rows, int blk) __arm_streaming __arm_inout("za") {
  svbool_t pg = svptrue_b64();
  svcount_t pc = svptrue_c64();
  svfloat64x4_t E0 = svld1_f64_x4(pc, EYE8), E1 = svld1_f64_x4(pc, EYE8 + 32);
  for (long R0 = 0; R0 < rows; R0 += 32) {
    svzero_za();
    for (int cb = 0; cb < Cn; cb += blk) {
      int ce = cb + blk < Cn ? cb + blk : Cn;
      if (cb) svzero_mask_za(0x0F);
      for (int c = cb; c < ce; c += 4) {
        int q = c >> 3, p = c & 7;
        svfloat64x4_t A = svld1_f64_x4(pc, T + 8 * c);
        MOPA4b(0) MOPA4b(1) MOPA4b(2) MOPA4b(3)
      }
      FOLD(0) FOLD(1) FOLD(2) FOLD(3)
    }
    OUTb(0) OUTb(1) OUTb(2) OUTb(3)
  }
}
__arm_locally_streaming __arm_new("za") static void smeb_chunk(const double *src, double *P, long rs, long prow8,
                                                                const double *T, int Cn, double *y, long full_rows,
                                                                double *ys, int blk, int *bad) {
  tr_body(src, P, rs, prow8, bad);
  if (full_rows) kb_body(P, rs, T, Cn, y, full_rows, blk);
  if (ys) kb_body(P + full_rows, rs, T, Cn, ys, 32, blk);
}

/* Blocked accumulation keeping all 8 tiles busy (path 5): each block of fc_blk5 tap-columns is accumulated from zero
   in ZA (512 outputs), then added to a running sum S kept in a scratch buffer in ZA-array order:
   ZA[array group w] += S[w] (multi-vector FADD into ZA, one rounding per block), S[w] = ZA[w]; the last block's sum is
   read out as usual. Same summation tree as path 4 (chains of fc_blk5 taps, then block sums in order). */
long fc_blk5 = 64;
static void kb8_body(const double *P, long rs, const double *T, int Cn, double *y, long rows, int blk, double *S) __arm_streaming __arm_inout("za") {
  svbool_t pg = svptrue_b64();
  svcount_t pc = svptrue_c64();
  for (long R0 = 0; R0 < rows; R0 += 64) {
    for (int cb = 0; cb < Cn; cb += blk) {
      int ce = cb + blk < Cn ? cb + blk : Cn;
      svzero_za();
      for (int c = cb; c < ce; c += 4) {
        int q = c >> 3, p = c & 7;
        svfloat64x4_t A = svld1_f64_x4(pc, T + 8 * c);
        MOPA8(0) MOPA8(1) MOPA8(2) MOPA8(3)
      }
      int last = ce >= Cn;
      if (cb) for (int w = 0; w < 16; w++) svadd_za64_f64_vg1x4(w, svld1_f64_x4(pc, S + 32 * w));
      if (!last) for (int w = 0; w < 16; w++) svst1_f64_x4(pc, S + 32 * w, svread_za64_f64_vg1x4(w));
    }
    OUT8(0) OUT8(1) OUT8(2) OUT8(3) OUT8(4) OUT8(5) OUT8(6) OUT8(7)
  }
}
__arm_locally_streaming __arm_new("za") static void smeb8_chunk(const double *src, double *P, long rs, long prow8,
                                                                 const double *T, int Cn, double *y, long full_rows,
                                                                 double *ys, int blk, double *S, int *bad) {
  tr_body(src, P, rs, prow8, bad);
  if (full_rows) kb8_body(P, rs, T, Cn, y, full_rows, blk, S);
  if (ys) kb8_body(P + full_rows, rs, T, Cn, ys, 64, blk, S);
}

/* one streaming entry per chunk: transpose, full passes, optional partial pass into ys */
__arm_locally_streaming __arm_new("za") static void sme_chunk(const double *src, double *P, long rs, long prow8,
                                                               const double *T, int Cn, double *y, long full_rows,
                                                               double *ys, int *bad) {
  tr_body(src, P, rs, prow8, bad);
  if (full_rows) k_body(P, rs, T, Cn, y, full_rows);
  if (ys) k_body(P + full_rows, rs, T, Cn, ys, 64);
}

/* ---------------- SME: complex kernel ---------------- */
#define MOPAC(u)                                                                 \
  {                                                                              \
    svfloat64x4_t Br = svld1_f64_x4(pc, Pr + (p + u) * rs + R0 + q);             \
    svfloat64x4_t Bi = svld1_f64_x4(pc, Pi + (p + u) * rs + R0 + q);             \
    svfloat64_t ar = svget4(Ar, u), ai = svget4(Ai, u), an = svget4(An, u);      \
    svmopa_za64_f64_m(0, pg, pg, ar, svget4(Br, 0));                             \
    svmopa_za64_f64_m(1, pg, pg, ar, svget4(Br, 1));                             \
    svmopa_za64_f64_m(2, pg, pg, ar, svget4(Br, 2));                             \
    svmopa_za64_f64_m(3, pg, pg, ar, svget4(Br, 3));                             \
    svmopa_za64_f64_m(4, pg, pg, ai, svget4(Br, 0));                             \
    svmopa_za64_f64_m(5, pg, pg, ai, svget4(Br, 1));                             \
    svmopa_za64_f64_m(6, pg, pg, ai, svget4(Br, 2));                             \
    svmopa_za64_f64_m(7, pg, pg, ai, svget4(Br, 3));                             \
    svmopa_za64_f64_m(0, pg, pg, an, svget4(Bi, 0));                             \
    svmopa_za64_f64_m(1, pg, pg, an, svget4(Bi, 1));                             \
    svmopa_za64_f64_m(2, pg, pg, an, svget4(Bi, 2));                             \
    svmopa_za64_f64_m(3, pg, pg, an, svget4(Bi, 3));                             \
    svmopa_za64_f64_m(4, pg, pg, ar, svget4(Bi, 0));                             \
    svmopa_za64_f64_m(5, pg, pg, ar, svget4(Bi, 1));                             \
    svmopa_za64_f64_m(6, pg, pg, ar, svget4(Bi, 2));                             \
    svmopa_za64_f64_m(7, pg, pg, ar, svget4(Bi, 3));                             \
  }
#define OUTC(t)                                                                  \
  for (int j = 0; j < 8; j += 4) {                                               \
    svfloat64x4_t re = svread_ver_za64_f64_vg4(t, j), im = svread_ver_za64_f64_vg4(t + 4, j); \
    double *yo = y + 16 * (R0 + 8 * (t) + j);                                    \
    svst2_f64(pg, yo, svcreate2(svget4(re, 0), svget4(im, 0)));                  \
    svst2_f64(pg, yo + 16, svcreate2(svget4(re, 1), svget4(im, 1)));             \
    svst2_f64(pg, yo + 32, svcreate2(svget4(re, 2), svget4(im, 2)));             \
    svst2_f64(pg, yo + 48, svcreate2(svget4(re, 3), svget4(im, 3)));             \
  }
/* rows multiple of 32; y interleaved complex, 256 complex outputs per pass */
static void kc_body(const double *Pr, const double *Pi, long rs, const double *Tr, const double *Ti, const double *Tn,
                    int Cn, double *y, long rows) __arm_streaming __arm_inout("za") {
  svbool_t pg = svptrue_b64();
  svcount_t pc = svptrue_c64();
  for (long R0 = 0; R0 < rows; R0 += 32) {
    svzero_za();
    for (int c = 0; c < Cn; c += 4) {
      int q = c >> 3, p = c & 7;
      svfloat64x4_t Ar = svld1_f64_x4(pc, Tr + 8 * c), Ai = svld1_f64_x4(pc, Ti + 8 * c),
                    An = svld1_f64_x4(pc, Tn + 8 * c);
      MOPAC(0) MOPAC(1) MOPAC(2) MOPAC(3)
    }
    OUTC(0) OUTC(1) OUTC(2) OUTC(3)
  }
}
__arm_locally_streaming __arm_new("za") static void smec_chunk(const double *src, double *Pr, double *Pi, long rs,
                                                                long prow8, const double *Tr, const double *Ti,
                                                                const double *Tn, int Cn, double *y, long full_rows,
                                                                double *ys) {
  trc_body(src, Pr, Pi, rs, prow8);
  if (full_rows) kc_body(Pr, Pi, rs, Tr, Ti, Tn, Cn, y, full_rows);
  if (ys) kc_body(Pr + full_rows, Pi + full_rows, rs, Tr, Ti, Tn, Cn, ys, 32);
}

/* ---------------- Toeplitz tables ---------------- */
static void make_T(double *T, const double *g, long gs, int L, int Cn) {
  for (int c = 0; c < Cn; c++)
    for (int i = 0; i < 8; i++) {
      int k = c - i;
      if (fc_broken && i == 3) k = c - i + 1; /* the deliberately broken copy: row 3 shifted by one tap */
      T[8 * c + i] = (k >= 0 && k < L) ? g[(long)k * gs] : 0.0;
    }
}

/* fill E[0..cnt) with xp[s0 .. s0+cnt), xp[i] = x[i - lp] (w doubles per element) */
static void fill_edge(double *E, const double *x, long nx, long lp, long s0, long cnt, int w) {
  long a = s0 - lp; /* x index of E[0] */
  long i0 = a < 0 ? -a : 0; if (i0 > cnt) i0 = cnt;
  long i1 = nx - a; if (i1 > cnt) i1 = cnt; if (i1 < i0) i1 = i0;
  if (i0) memset(E, 0, sizeof(double) * w * i0);
  if (i1 > i0) memcpy(E + w * i0, x + w * (a + i0), sizeof(double) * w * (i1 - i0));
  if (cnt > i1) memset(E + w * i1, 0, sizeof(double) * w * (cnt - i1));
}

/* ---------------- real: SME path ---------------- */
/* Non-finite inputs: the Toeplitz table holds structural zeros, and 0 * NaN = NaN, 0 * inf = NaN, so a NaN or inf in x
   would spread a few outputs beyond its true support. Any output touched by a non-finite product is itself non-finite,
   so a scan of each finished chunk (hot in cache) finds every affected chunk, which is then recomputed with the NEON
   direct path (bit-identical order, numpy's support). fc_nancheck = 0 disables the scan. */
int fc_nancheck = 3; /* 0 off, 1 scan outputs on the core, 2 scan inputs on the core, 3 check inside the streaming transpose (default) */
static void neon_valid(const double *xp, const double *g, int L, double *y, long n);
static int neon_real(const double *x, long nx, const double *g, int L, long o0, long nout, double *y);
static int scalar_cplx(const double *x, long nx, const double *g, int L, long o0, long nout, double *y);
static int any_nonfinite(const double *y, long n) {
  float64x2_t a0 = vdupq_n_f64(0), a1 = a0, a2 = a0, a3 = a0;
  long i = 0;
  for (; i + 8 <= n; i += 8) {
    float64x2_t v0 = vld1q_f64(y + i), v1 = vld1q_f64(y + i + 2), v2 = vld1q_f64(y + i + 4), v3 = vld1q_f64(y + i + 6);
    a0 = vaddq_f64(a0, vsubq_f64(v0, v0)); a1 = vaddq_f64(a1, vsubq_f64(v1, v1));
    a2 = vaddq_f64(a2, vsubq_f64(v2, v2)); a3 = vaddq_f64(a3, vsubq_f64(v3, v3));
  }
  double s = vaddvq_f64(vaddq_f64(vaddq_f64(a0, a1), vaddq_f64(a2, a3)));
  for (; i < n; i++) s += y[i] - y[i];
  return s != s;
}
long fc_blk4 = 64; /* tap-columns per block for path 4 (blocked, accurate) */
static int sme_real(const double *x, long nx, const double *g, int L, long o0, long nout, double *y, long blkset, int v8) {
  int Cn = (L + 7 + 3) & ~3;
  int qmax = (Cn - 1) >> 3;
  long ch = fc_chunk;
  if (ch < 4L * L) ch = 4L * L;
  ch = (ch + 511) & ~511L;
  long nr = (nout + 511) & ~511L;
  if (ch > nr) ch = nr;
  long rs = ((ch / 8 + qmax + 1 + 7) & ~7L) + 8;
  scratch_t *sc = scratch();
  if (!sc) return -1;
  double *T = getbuf(&sc->T, 8 * (size_t)Cn), *P = getbuf(&sc->P, 8 * (size_t)rs), *E = getbuf(&sc->E, 8 * (size_t)rs),
         *ys = getbuf(&sc->Y, 512);
  if (!T || !P || !E || !ys) return -1;
  if (fc_nancheck && any_nonfinite(g, L)) return neon_real(x, nx, g, L, o0, nout, y); /* non-finite taps: NEON, numpy's support */
  make_T(T, g, 1, L, Cn);
  long lp = L - 1;
  int blk = (int)((blkset + 3) & ~3L), useb = blk > 0 && Cn > blk;
  double *S = getbuf(&sc->S, 512);
  if (!S) return -1;
  for (long m = 0; m < nout; m += ch) {
    long cn = ch < nout - m ? ch : nout - m;
    long full_rows = (cn / 512) * 64, rem = cn - full_rows * 8;
    long krows = full_rows + (rem ? 64 : 0);
    long prow8 = (krows + qmax + 1 + 7) & ~7L;
    long s0 = o0 + m;
    const double *src;
    if (s0 >= lp && s0 + 8 * prow8 <= lp + nx) src = x + (s0 - lp);
    else { fill_edge(E, x, nx, lp, s0, 8 * prow8, 1); src = E; }
    int bad = 0, *pbad = fc_nancheck == 3 ? &bad : NULL;
    if (fc_nancheck == 2 && any_nonfinite(src, 8 * prow8)) { /* input scan variant: core reads x, never SME output */
      double *Es = sc->E.p; size_t En = sc->E.n; sc->E.p = NULL; sc->E.n = 0;
      int r = neon_real(x, nx, g, L, o0 + m, cn, y + m);
      free(sc->E.p); sc->E.p = Es; sc->E.n = En;
      if (r) return r;
      continue;
    }
    if (useb && v8) {
      smeb8_chunk(src, P, rs, prow8, T, Cn, y + m, full_rows, rem ? ys : NULL, blk, S, pbad);
      if (rem) memcpy(y + m + full_rows * 8, ys, 8 * rem);
    } else if (useb) {
      long fr = (cn / 256) * 32, rm = cn - fr * 8;
      smeb_chunk(src, P, rs, prow8, T, Cn, y + m, fr, rm ? ys : NULL, blk, pbad);
      if (rm) memcpy(y + m + fr * 8, ys, 8 * rm);
    } else {
      sme_chunk(src, P, rs, prow8, T, Cn, y + m, full_rows, rem ? ys : NULL, pbad);
      if (rem) memcpy(y + m + full_rows * 8, ys, 8 * rem);
    }
    if ((fc_nancheck == 1 && any_nonfinite(y + m, cn)) || bad) {
      double *Es = sc->E.p; size_t En = sc->E.n; sc->E.p = NULL; sc->E.n = 0; /* keep E for the remaining chunks */
      int r = neon_real(x, nx, g, L, o0 + m, cn, y + m);
      free(sc->E.p); sc->E.p = Es; sc->E.n = En;
      if (r) return r;
    }
  }
  return 0;
}

/* ---------------- complex: SME path (x, g, y interleaved complex) ---------------- */
static int sme_cplx(const double *x, long nx, const double *g, int L, long o0, long nout, double *y) {
  int Cn = (L + 7 + 3) & ~3;
  int qmax = (Cn - 1) >> 3;
  long ch = fc_chunk / 2;
  if (ch < 4L * L) ch = 4L * L;
  ch = (ch + 255) & ~255L;
  long nr = (nout + 255) & ~255L;
  if (ch > nr) ch = nr;
  long rs = ((ch / 8 + qmax + 1 + 7) & ~7L) + 8;
  scratch_t *sc = scratch();
  if (!sc) return -1;
  double *T = getbuf(&sc->T, 8 * (size_t)Cn * 3), *P = getbuf(&sc->P, 16 * (size_t)rs),
         *E = getbuf(&sc->E, 16 * (size_t)rs), *ys = getbuf(&sc->Y, 512);
  if (!T || !P || !E || !ys) return -1;
  double *Tr = T, *Ti = T + 8 * Cn, *Tn = T + 16 * Cn;
  make_T(Tr, g, 2, L, Cn);
  make_T(Ti, g + 1, 2, L, Cn);
  for (long i = 0; i < 8L * Cn; i++) Tn[i] = -Ti[i];
  double *Pr = P, *Pi = P + 8 * rs;
  long lp = L - 1;
  for (long m = 0; m < nout; m += ch) {
    long cn = ch < nout - m ? ch : nout - m;
    long full_rows = (cn / 256) * 32, rem = cn - full_rows * 8;
    long krows = full_rows + (rem ? 32 : 0);
    long prow8 = (krows + qmax + 1 + 7) & ~7L;
    long s0 = o0 + m;
    const double *src;
    if (s0 >= lp && s0 + 8 * prow8 <= lp + nx) src = x + 2 * (s0 - lp);
    else { fill_edge(E, x, nx, lp, s0, 8 * prow8, 2); src = E; }
    smec_chunk(src, Pr, Pi, rs, prow8, Tr, Ti, Tn, Cn, y + 2 * m, full_rows, rem ? ys : NULL);
    if (rem) memcpy(y + 2 * (m + full_rows * 8), ys, 16 * rem);
    if (fc_nancheck && any_nonfinite(y + 2 * m, 2 * cn)) scalar_cplx(x, nx, g, L, o0 + m, cn, y + 2 * m);
  }
  return 0;
}

/* ---------------- NEON direct (valid correlation over a padded copy) ---------------- */
static void neon_valid(const double *xp, const double *g, int L, double *y, long n) {
  long t = 0;
  for (; t + 16 <= n; t += 16) {
    float64x2_t a0 = vdupq_n_f64(0), a1 = a0, a2 = a0, a3 = a0, a4 = a0, a5 = a0, a6 = a0, a7 = a0;
    const double *xx = xp + t;
    for (int k = 0; k < L; k++) {
      float64x2_t gk = vdupq_n_f64(g[k]);
      const double *z = xx + k;
      a0 = vfmaq_f64(a0, vld1q_f64(z), gk);
      a1 = vfmaq_f64(a1, vld1q_f64(z + 2), gk);
      a2 = vfmaq_f64(a2, vld1q_f64(z + 4), gk);
      a3 = vfmaq_f64(a3, vld1q_f64(z + 6), gk);
      a4 = vfmaq_f64(a4, vld1q_f64(z + 8), gk);
      a5 = vfmaq_f64(a5, vld1q_f64(z + 10), gk);
      a6 = vfmaq_f64(a6, vld1q_f64(z + 12), gk);
      a7 = vfmaq_f64(a7, vld1q_f64(z + 14), gk);
    }
    vst1q_f64(y + t, a0); vst1q_f64(y + t + 2, a1); vst1q_f64(y + t + 4, a2); vst1q_f64(y + t + 6, a3);
    vst1q_f64(y + t + 8, a4); vst1q_f64(y + t + 10, a5); vst1q_f64(y + t + 12, a6); vst1q_f64(y + t + 14, a7);
  }
  for (; t < n; t++) {
    double s = 0;
    for (int k = 0; k < L; k++) s = __builtin_fma(g[k], xp[t + k], s);
    y[t] = s;
  }
}
static int neon_real(const double *x, long nx, const double *g, int L, long o0, long nout, double *y) {
  long lp = L - 1;
  long need = nout + L - 1;
  const double *xp;
  if (o0 >= lp && o0 + need <= lp + nx) xp = x + (o0 - lp);
  else {
    scratch_t *sc = scratch();
    double *E = sc ? getbuf(&sc->E, need) : NULL;
    if (!E) return -1;
    fill_edge(E, x, nx, lp, o0, need, 1);
    xp = E;
  }
  neon_valid(xp, g, L, y, nout);
  return 0;
}

/* ---------------- scalar reference order (one fma chain per output, k ascending) ---------------- */
static int scalar_real(const double *x, long nx, const double *g, int L, long o0, long nout, double *y) {
  long lp = L - 1;
  for (long t = 0; t < nout; t++) {
    double s = 0;
    for (int k = 0; k < L; k++) {
      long i = o0 + t + k - lp;
      if (i >= 0 && i < nx) s = __builtin_fma(g[k], x[i], s);
    }
    y[t] = s;
  }
  return 0;
}
static int scalar_cplx(const double *x, long nx, const double *g, int L, long o0, long nout, double *y) {
  long lp = L - 1;
  for (long t = 0; t < nout; t++) {
    double sr = 0, si = 0;
    for (int k = 0; k < L; k++) {
      long i = o0 + t + k - lp;
      if (i >= 0 && i < nx) {
        double xr = x[2 * i], xi = x[2 * i + 1], gr = g[2 * k], gi = g[2 * k + 1];
        sr += gr * xr - gi * xi;
        si += gr * xi + gi * xr;
      }
    }
    y[2 * t] = sr; y[2 * t + 1] = si;
  }
  return 0;
}

/* ---------------- entry points ---------------- */
static int corr_mt(const double *x, long nx, const double *g, int L, long o0, long nout, double *y, int k);
int fc_corr_d(const double *x, long nx, const double *g, int L, long o0, long nout, double *y, int path) {
  if (nout <= 0) return 0;
  if (L <= 0) { memset(y, 0, 8 * nout); return 0; }
  if (path == 0) path = (nout <= fc_neon_maxn || nout * (long)L <= fc_neon_max) ? 2 : 1;
  if (path == 1) return sme_real(x, nx, g, L, o0, nout, y, fc_blk, 0);
  if (path == 4) return sme_real(x, nx, g, L, o0, nout, y, fc_blk4, 0);
  if (path == 5) return sme_real(x, nx, g, L, o0, nout, y, fc_blk5, 1);
  if (path > 10 && path <= 18) return corr_mt(x, nx, g, L, o0, nout, y, path - 10);
  if (path == 2) return neon_real(x, nx, g, L, o0, nout, y);
  return scalar_real(x, nx, g, L, o0, nout, y);
}
int fc_corr_z(const double *x, long nx, const double *g, int L, long o0, long nout, double *y, int path) {
  if (nout <= 0) return 0;
  if (L <= 0) { memset(y, 0, 16 * nout); return 0; }
  if (path == 3) return scalar_cplx(x, nx, g, L, o0, nout, y);
  return sme_cplx(x, nx, g, L, o0, nout, y);
}

/* ---------------- whole P-cluster: the caller drives SME on the first fc_frac of the outputs, k helper threads run
   the NEON direct path on the rest (path 10 + k). Threads are created per call (about 10-20 us each); when one cannot
   be created (out of memory or threads), the caller computes that share itself. ---------------- */
double fc_frac = 0.8;
typedef struct { const double *x; long nx; const double *g; int L; long o0, nout; double *y; int r; } job_t;
static void *neon_job(void *p) {
  job_t *j = p;
  pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
  j->r = neon_real(j->x, j->nx, j->g, j->L, j->o0, j->nout, j->y);
  return NULL;
}
static int corr_mt(const double *x, long nx, const double *g, int L, long o0, long nout, double *y, int k) {
  long ns = ((long)(nout * fc_frac)) & ~511L;
  if (k < 1 || nout - ns < 64 * k) return sme_real(x, nx, g, L, o0, nout, y, fc_blk, 0);
  pthread_t th[8];
  job_t jb[8];
  long rest = nout - ns, per = (rest + k - 1) / k;
  int nt = 0, r = 0;
  for (int i = 0; i < k && i < 8; i++) {
    long a = ns + i * per, b = a + per < nout ? a + per : nout;
    if (b <= a) break;
    jb[nt] = (job_t){x, nx, g, L, o0 + a, b - a, y + a, 0};
    if (pthread_create(&th[nt], NULL, neon_job, &jb[nt]) == 0) nt++;
    else if ((jb[nt].r = neon_real(x, nx, g, L, o0 + a, b - a, y + a))) r = jb[nt].r; /* no thread: do it here */
  }
  int rs = ns ? sme_real(x, nx, g, L, o0, ns, y, fc_blk, 0) : 0;
  if (rs) r = rs;
  for (int i = 0; i < nt; i++) { pthread_join(th[i], NULL); if (jb[i].r) r = jb[i].r; }
  return r;
}

/* vDSP_convD baseline in the same primitive form (correlation: vDSP_convD with stride +1 IS correlation).
   Needs a padded copy when the range touches the zero padding. */
int fc_vdsp_d(const double *x, long nx, const double *g, int L, long o0, long nout, double *y) {
  long lp = L - 1, need = nout + L - 1;
  const double *xp;
  if (o0 >= lp && o0 + need <= lp + nx) xp = x + (o0 - lp);
  else {
    scratch_t *sc = scratch();
    double *E = sc ? getbuf(&sc->E, need) : NULL;
    if (!E) return -1;
    fill_edge(E, x, nx, lp, o0, need, 1);
    xp = E;
  }
  vDSP_convD(xp, 1, g, 1, y, 1, (vDSP_Length)nout, (vDSP_Length)L);
  return 0;
}
