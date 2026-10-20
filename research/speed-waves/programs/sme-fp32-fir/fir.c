/* fir.c -- single-channel fp32 FIR on the M4: SME2 Toeplitz-tile kernel vs vDSP_conv vs vDSP overlap-save.
 * Semantics ("valid" convolution, same as vDSP_conv with reversed taps):
 *   x has n + L - 1 samples, h has L taps, y[m] = sum_k h[k] x[m + L - 1 - k], m = 0..n-1.
 *   Internally g[k] = h[L-1-k] (correlation form): y[m] = sum_k g[k] x[m + k].
 * SME kernel (agent aj-1, wave 11):
 *   polyphase input P[p][r] = x[16 r + p] (NEON transpose, per chunk);
 *   Toeplitz table T[c][i] = g[c - i] (0 outside [0,L)), c = 0..Cn-1, Cn = roundup4(L + 15);
 *   ZA tile t (t=0..3), pass rows R0..R0+63: ZA_t[i][j] += T[c][i] * P[c%16][R0 + c/16 + 16t + j]
 *   gives ZA_t[i][j] = y[16 (R0 + 16t + j) + i]; vertical slice j = 16 contiguous outputs.
 * Build: clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$B/pffft/include fir.c $B/pffft/lib/libpffft.a -framework Accelerate -o fir
 *   (B=$BASELINES)
 * Use:   fir check L LOGN          (all three programs vs a double reference)
 *        fir time L LOGN same|separate [progs]   (progs: bitmask 1=vdsp_conv 2=ols 4=sme; default 7)
 *        fir olsscan L LOGN        (time overlap-save block sizes, print best)
 */
#include <Accelerate/Accelerate.h>
#include <arm_neon.h>
#include <arm_sme.h>
#include "../timing-harness/bench.h"
#include "pffft.h"

static int OLS_N_FORCE = 0; /* env OLS_N */
static long SME_CH = 0;     /* env SME_CH: outputs per chunk */

/* ---------------- reference ---------------- */
static double ref_one(const float *x, const float *h, int L, long m) {
  double s = 0;
  for (int k = 0; k < L; k++) s += (double)h[k] * (double)x[m + L - 1 - k];
  return s;
}

/* ---------------- vDSP_conv ---------------- */
typedef struct {
  bench_pool *xp, *yp;
  const float *h;
  int L;
  long n;
  /* ols */
  int N, log2N;
  vDSP_DFT_Setup fwd, inv;
  float *Hr, *Hi; /* scaled filter spectrum, packed (N/2 each) */
  float *sr, *si, *tr, *ti, *tmp;
  /* sme */
  float *T; int Cn; float *P; long rs; long ch;
  float *er, *od; /* ols2: whole input pre-split into even / odd samples */
  float *ys;      /* sme: scratch for the last partial pass */
  PFFFT_Setup *pf; float *pX, *pY, *pH, *pw, *pt; int pN;
} ctx_t;
static int OLS_PF_N = 0; /* env OLS_PF_N */
static int default_pf_N(int L) { int N = 512; while (N < 4 * L) N *= 2; return N; }
static void pf_setup(ctx_t *c, int N) {
  c->pN = N;
  c->pf = pffft_new_setup(N, PFFFT_REAL);
  if (!c->pf) { fprintf(stderr, "pffft setup failed N=%d\n", N); exit(1); }
  c->pX = bench_alloc(4 * (size_t)N); c->pY = bench_alloc(4 * (size_t)N); c->pH = bench_alloc(4 * (size_t)N);
  c->pw = bench_alloc(8 * (size_t)N); c->pt = bench_alloc(8 * (size_t)N);
  float *hp = bench_alloc(4 * (size_t)N);
  for (int k = 0; k < c->L; k++) hp[k] = c->h[k] / (float)N;
  pffft_transform(c->pf, hp, c->pH, c->pw, PFFFT_FORWARD);
  free(hp);
}
/* overlap-save with PFFFT (NEON, z-domain order, no reordering); step multiple of 4 (16-byte alignment) */
static void run_ols_pf(void *vc, long i) {
  ctx_t *c = vc;
  const float *x = bench_pool_get(c->xp, i);
  float *y = bench_pool_get(c->yp, i);
  int N = c->pN, L = c->L;
  long n = c->n, xlen = n + L - 1, step = (N - L + 1) & ~3L, m = 0;
  for (; m + step <= n && m + N <= xlen; m += step) {
    pffft_transform(c->pf, x + m, c->pX, c->pw, PFFFT_FORWARD);
    pffft_zconvolve(c->pf, c->pX, c->pH, c->pY);
    pffft_transform(c->pf, c->pY, c->pt, c->pw, PFFFT_BACKWARD);
    memcpy(y + m, c->pt + L - 1, 4 * step);
  }
  float *seg = c->pt + N;
  while (m < n) {
    long avail = xlen - m, cp = avail < N ? avail : N;
    memcpy(seg, x + m, 4 * cp); memset(seg + cp, 0, 4 * (N - cp));
    long k = n - m < N - L + 1 ? n - m : N - L + 1;
    pffft_transform(c->pf, seg, c->pX, c->pw, PFFFT_FORWARD);
    pffft_zconvolve(c->pf, c->pX, c->pH, c->pY);
    pffft_transform(c->pf, c->pY, c->pt, c->pw, PFFFT_BACKWARD);
    memcpy(y + m, c->pt + L - 1, 4 * k);
    m += k;
  }
}

static void run_vdsp(void *vc, long i) {
  ctx_t *c = vc;
  const float *x = bench_pool_get(c->xp, i);
  float *y = bench_pool_get(c->yp, i);
  vDSP_conv(x, 1, c->h + c->L - 1, -1, y, 1, c->n, c->L);
}

/* ---------------- overlap-save with vDSP real DFT (zrop) ---------------- */
static void ols_setup(ctx_t *c, int N) {
  c->N = N;
  c->fwd = vDSP_DFT_zrop_CreateSetup(NULL, N, vDSP_DFT_FORWARD);
  c->inv = vDSP_DFT_zrop_CreateSetup(c->fwd, N, vDSP_DFT_INVERSE);
  if (!c->fwd || !c->inv) { fprintf(stderr, "DFT setup failed N=%d\n", N); exit(1); }
  int h2 = N / 2;
  c->Hr = bench_alloc(4 * h2); c->Hi = bench_alloc(4 * h2);
  c->sr = bench_alloc(4 * h2); c->si = bench_alloc(4 * h2);
  c->tr = bench_alloc(4 * h2); c->ti = bench_alloc(4 * h2);
  c->tmp = bench_alloc(4 * 2 * (size_t)N);
  c->er = bench_alloc(4 * (size_t)((c->n + c->L) / 2 + 2)); c->od = bench_alloc(4 * (size_t)((c->n + c->L) / 2 + 2));
  float *hp = bench_alloc(4 * (size_t)N);
  memcpy(hp, c->h, 4 * (size_t)c->L);
  DSPSplitComplex s = {c->sr, c->si};
  vDSP_ctoz((const DSPComplex *)hp, 2, &s, 1, h2);
  vDSP_DFT_Execute(c->fwd, c->sr, c->si, c->Hr, c->Hi);
  /* forward zrop = 2 x DFT; inverse zrop of packed Z = N x IDFT * ... ; total scale for fwd*fwd*inv = 4 * N/... */
  float sc = 1.0f / (4.0f * (float)N); /* fwd zrop = 2 DFT (x2 twice), inv zrop = N IDFT: calibrated by check */
  for (int k = 0; k < h2; k++) { c->Hr[k] *= sc; c->Hi[k] *= sc; }
  free(hp);
}

static void ols_block(ctx_t *c, const float *seg, float *out, long nout) {
  int N = c->N, h2 = N / 2, L = c->L;
  DSPSplitComplex s = {c->sr, c->si}, t = {c->tr, c->ti}, H = {c->Hr, c->Hi};
  vDSP_ctoz((const DSPComplex *)seg, 2, &s, 1, h2);
  vDSP_DFT_Execute(c->fwd, c->sr, c->si, c->tr, c->ti);
  float dc = c->tr[0] * c->Hr[0], ny = c->ti[0] * c->Hi[0];
  vDSP_zvmul(&t, 1, &H, 1, &t, 1, h2, 1);
  c->tr[0] = dc; c->ti[0] = ny;
  vDSP_DFT_Execute(c->inv, c->tr, c->ti, c->sr, c->si);
  if (((L - 1) & 1) == 0 && nout == N - L + 1) {
    /* output starts at an even position: write straight from split */
    vDSP_ztoc(&s, 1, (DSPComplex *)c->tmp, 2, h2);
    memcpy(out, c->tmp + L - 1, 4 * nout);
  } else {
    vDSP_ztoc(&s, 1, (DSPComplex *)c->tmp, 2, h2);
    memcpy(out, c->tmp + L - 1, 4 * nout);
  }
}

/* Stronger overlap-save: split the whole input into even/odd ONCE (vDSP_ctoz), even block step, forward DFT
   straight from the split arrays, inverse into split, interleave (vDSP_ztoc) straight into y. */
static void run_ols2(void *vc, long i) {
  ctx_t *c = vc;
  const float *x = bench_pool_get(c->xp, i);
  float *y = bench_pool_get(c->yp, i);
  int N = c->N, L = c->L, h2 = N / 2;
  long n = c->n, xlen = n + L - 1;
  long step = (N - L + 1) & ~1L;
  DSPSplitComplex xs = {c->er, c->od};
  vDSP_ctoz((const DSPComplex *)x, 2, &xs, 1, xlen / 2);
  if (xlen & 1) { c->er[xlen / 2] = x[xlen - 1]; c->od[xlen / 2] = 0.f; }
  DSPSplitComplex t = {c->tr, c->ti}, H = {c->Hr, c->Hi};
  long m = 0;
  for (; m + step <= n && m + N <= xlen + 1; m += step) {
    vDSP_DFT_Execute(c->fwd, c->er + m / 2, c->od + m / 2, c->tr, c->ti);
    float dc = c->tr[0] * c->Hr[0], ny = c->ti[0] * c->Hi[0];
    vDSP_zvmul(&t, 1, &H, 1, &t, 1, h2, 1);
    c->tr[0] = dc; c->ti[0] = ny;
    vDSP_DFT_Execute(c->inv, c->tr, c->ti, c->sr, c->si);
    if (L & 1) { /* first useful output t = L-1 is even */
      DSPSplitComplex o = {c->sr + (L - 1) / 2, c->si + (L - 1) / 2};
      vDSP_ztoc(&o, 1, (DSPComplex *)(y + m), 2, step / 2);
    } else {
      y[m] = c->si[(L - 1) / 2];
      DSPSplitComplex o = {c->sr + L / 2, c->si + L / 2};
      vDSP_ztoc(&o, 1, (DSPComplex *)(y + m + 1), 2, (step - 2) / 2);
      y[m + step - 1] = c->sr[(L + step - 2) / 2];
    }
  }
  while (m < n) { /* tail through the simple path (zero-padded segment) */
    long avail = xlen - m, cp = avail < N ? avail : N;
    float *seg = c->tmp + N;
    memcpy(seg, x + m, 4 * cp);
    memset(seg + cp, 0, 4 * (N - cp));
    long k = n - m < N - L + 1 ? n - m : N - L + 1;
    ols_block(c, seg, y + m, k);
    m += k;
  }
}

static void run_ols(void *vc, long i) {
  ctx_t *c = vc;
  const float *x = bench_pool_get(c->xp, i);
  float *y = bench_pool_get(c->yp, i);
  int N = c->N, L = c->L;
  long step = N - L + 1, n = c->n, m = 0;
  for (; m + step <= n; m += step) ols_block(c, x + m, y + m, step);
  if (m < n) { /* last partial block: copy remaining input into a zero-padded segment */
    long avail = n + L - 1 - m;
    float *seg = c->tmp + N; /* tmp has 2N */
    memcpy(seg, x + m, 4 * avail);
    memset(seg + avail, 0, 4 * (N - avail));
    ols_block(c, seg, y + m, n - m);
  }
}

/* ---------------- SME2 Toeplitz kernel ---------------- */
static void sme_setup(ctx_t *c) {
  int L = c->L;
  c->Cn = (L + 15 + 3) & ~3;
  c->T = bench_alloc(4 * 16 * (size_t)(c->Cn + 4)); /* + one zero group for the short-block pass */
  for (int cc = 0; cc < c->Cn; cc++)
    for (int i = 0; i < 16; i++) {
      int k = cc - i;
      c->T[16 * cc + i] = (k >= 0 && k < L) ? c->h[L - 1 - k] : 0.f; /* g[k] = h[L-1-k] */
    }
  long ch = SME_CH ? SME_CH : 16384;
  if (ch < 4L * L) ch = 4L * L;
  ch = (ch + 1023) & ~1023L;
  if (ch > c->n) ch = (c->n + 1023) & ~1023L;
  c->ys = bench_alloc(4 * 1024);
  c->ch = ch;
  int qmax = (c->Cn - 1) / 16;
  c->rs = ((ch / 16 + qmax + 1) + 31) & ~31L; /* row length, padded */
  c->P = bench_alloc(4 * 16 * (size_t)c->rs);
}

/* P[p*rs + r] = x[16 r + p] for r in [0, rows); x valid for indices < xvalid (zero beyond). */
static void transpose16(const float *x, long xvalid, float *P, long rs, long rows) {
  long r = 0;
  long rfull = xvalid >= 16 ? (xvalid - 16) / 16 + 1 : 0; /* rows fully inside */
  if (rfull > rows) rfull = rows;
  for (; r + 4 <= rfull; r += 4) {
    const float *xr = x + 16 * r;
    for (int g = 0; g < 4; g++) {
      float32x4_t a0 = vld1q_f32(xr + 4 * g), a1 = vld1q_f32(xr + 16 + 4 * g), a2 = vld1q_f32(xr + 32 + 4 * g),
                  a3 = vld1q_f32(xr + 48 + 4 * g);
      float32x4_t t0 = vtrn1q_f32(a0, a1), t1 = vtrn2q_f32(a0, a1), t2 = vtrn1q_f32(a2, a3), t3 = vtrn2q_f32(a2, a3);
      float32x4_t o0 = vreinterpretq_f32_f64(vtrn1q_f64(vreinterpretq_f64_f32(t0), vreinterpretq_f64_f32(t2)));
      float32x4_t o1 = vreinterpretq_f32_f64(vtrn1q_f64(vreinterpretq_f64_f32(t1), vreinterpretq_f64_f32(t3)));
      float32x4_t o2 = vreinterpretq_f32_f64(vtrn2q_f64(vreinterpretq_f64_f32(t0), vreinterpretq_f64_f32(t2)));
      float32x4_t o3 = vreinterpretq_f32_f64(vtrn2q_f64(vreinterpretq_f64_f32(t1), vreinterpretq_f64_f32(t3)));
      vst1q_f32(P + (4 * g + 0) * rs + r, o0);
      vst1q_f32(P + (4 * g + 1) * rs + r, o1);
      vst1q_f32(P + (4 * g + 2) * rs + r, o2);
      vst1q_f32(P + (4 * g + 3) * rs + r, o3);
    }
  }
  for (; r < rows; r++)
    for (int p = 0; p < 16; p++) {
      long ix = 16 * r + p;
      P[p * rs + r] = ix < xvalid ? x[ix] : 0.f;
    }
}

#define MOPA4(u)                                                              \
  {                                                                           \
    svfloat32x4_t B = svld1_f32_x4(pc, P + (p + u) * rs + R0 + q);            \
    svfloat32_t a = svget4(A, u);                                             \
    svmopa_za32_f32_m(0, pg, pg, a, svget4(B, 0));                            \
    svmopa_za32_f32_m(1, pg, pg, a, svget4(B, 1));                            \
    svmopa_za32_f32_m(2, pg, pg, a, svget4(B, 2));                            \
    svmopa_za32_f32_m(3, pg, pg, a, svget4(B, 3));                            \
  }
#define OUT4(t, j)                                                            \
  {                                                                           \
    svfloat32x4_t q4 = svread_ver_za32_f32_vg4(t, j);                         \
    svst1_f32_x4(pc, y + 16 * (R0 + 16 * t + j), q4);                         \
  }

/* Streaming-mode transpose through ZA: P[p*rs + r] = x[16 r + p], r in [0, rows), rows multiple of 16,
   all x reads in range (caller guarantees). Variant: x4 loads + horizontal vg4 writes, vertical single-slice stores. */
#define TLOAD(t, s) { svfloat32x4_t v = svld1_f32_x4(pc, xr + 64 * (s)); svwrite_hor_za32_f32_vg4(t, 4 * (s), v); }
#define TST(t, p) svst1_ver_za32(t, p, pg, P + (p) * rs + r0);
static long SME_PF = 0; /* env SME_PF: prefetch distance in floats for the ZA transpose (0 = none) */
#define TPF if (pfd) { const char *pp = (const char *)(xr + pfd); for (int l = 0; l < 1024; l += 128) __builtin_prefetch(pp + l, 0, 0); }
#define TBLK(t) TPF \
  TLOAD(t, 0) TLOAD(t, 1) TLOAD(t, 2) TLOAD(t, 3) \
  TST(t, 0) TST(t, 1) TST(t, 2) TST(t, 3) TST(t, 4) TST(t, 5) TST(t, 6) TST(t, 7) \
  TST(t, 8) TST(t, 9) TST(t, 10) TST(t, 11) TST(t, 12) TST(t, 13) TST(t, 14) TST(t, 15)
__arm_locally_streaming __arm_new("za") static void sme_transpose(const float *x, float *P, long rs, long rows) {
  svbool_t pg = svptrue_b32();
  svcount_t pc = svptrue_c32();
  long r0 = 0, pfd = SME_PF;
  for (; r0 + 64 <= rows; ) {
    const float *xr = x + 16 * r0;
    TBLK(0) r0 += 16; xr += 256;
    TBLK(1) r0 += 16; xr += 256;
    TBLK(2) r0 += 16; xr += 256;
    TBLK(3) r0 += 16;
  }
  for (; r0 < rows; r0 += 16) { const float *xr = x + 16 * r0; TBLK(0) }
}
static int SME_TR = 1; /* env SME_TR=0: use the NEON transpose instead of the streaming ZA transpose */

/* rows: number of output rows (16 outputs each), multiple of 64 */
__arm_locally_streaming __arm_new("za") static void sme_kernel(const float *P, long rs, const float *T, int Cn,
                                                                float *y, long rows) {
  svbool_t pg = svptrue_b32();
  svcount_t pc = svptrue_c32();
  for (long R0 = 0; R0 < rows; R0 += 64) {
    svzero_za();
    for (int c = 0; c < Cn; c += 4) {
      int q = c >> 4, p = c & 15;
      svfloat32x4_t A = svld1_f32_x4(pc, T + 16 * c);
      MOPA4(0) MOPA4(1) MOPA4(2) MOPA4(3)
    }
    OUT4(0, 0) OUT4(0, 4) OUT4(0, 8) OUT4(0, 12)
    OUT4(1, 0) OUT4(1, 4) OUT4(1, 8) OUT4(1, 12)
    OUT4(2, 0) OUT4(2, 4) OUT4(2, 8) OUT4(2, 12)
    OUT4(3, 0) OUT4(3, 4) OUT4(3, 8) OUT4(3, 12)
  }
}


/* Small-block pass (at most 16 rows = 256 outputs): the TAPS are split over the 4 ZA tiles (4 independent
   accumulation chains for the same 16 rows), then the 4 tiles are summed. y gets 256 outputs (16 columns). */
__arm_locally_streaming __arm_new("za") static void sme_kernel_small(const float *P, long rs, const float *T, int Cn,
                                                                      float *y) {
  svbool_t pg = svptrue_b32();
  svcount_t pc = svptrue_c32();
  int ng = Cn / 4, gq = (ng + 3) / 4;
  svzero_za();
  for (int kk = 0; kk < gq; kk++) {
#define SMALLT(t)                                                                  \
    {                                                                              \
      int g = (t) * gq + kk;                                                       \
      if (g < ng) {                                                                \
        int c = 4 * g, q = c >> 4, p = c & 15;                                     \
        svfloat32x4_t A = svld1_f32_x4(pc, T + 16 * c);                            \
        svmopa_za32_f32_m(t, pg, pg, svget4(A, 0), svld1_f32(pg, P + (p + 0) * rs + q)); \
        svmopa_za32_f32_m(t, pg, pg, svget4(A, 1), svld1_f32(pg, P + (p + 1) * rs + q)); \
        svmopa_za32_f32_m(t, pg, pg, svget4(A, 2), svld1_f32(pg, P + (p + 2) * rs + q)); \
        svmopa_za32_f32_m(t, pg, pg, svget4(A, 3), svld1_f32(pg, P + (p + 3) * rs + q)); \
      }                                                                            \
    }
    SMALLT(0) SMALLT(1) SMALLT(2) SMALLT(3)
  }
  for (int j = 0; j < 16; j += 4) {
    svfloat32x4_t a = svread_ver_za32_f32_vg4(0, j), b = svread_ver_za32_f32_vg4(1, j);
    svfloat32x4_t c2 = svread_ver_za32_f32_vg4(2, j), d = svread_ver_za32_f32_vg4(3, j);
    svfloat32_t o0 = svadd_f32_x(pg, svadd_f32_x(pg, svget4(a, 0), svget4(b, 0)), svadd_f32_x(pg, svget4(c2, 0), svget4(d, 0)));
    svfloat32_t o1 = svadd_f32_x(pg, svadd_f32_x(pg, svget4(a, 1), svget4(b, 1)), svadd_f32_x(pg, svget4(c2, 1), svget4(d, 1)));
    svfloat32_t o2 = svadd_f32_x(pg, svadd_f32_x(pg, svget4(a, 2), svget4(b, 2)), svadd_f32_x(pg, svget4(c2, 2), svget4(d, 2)));
    svfloat32_t o3 = svadd_f32_x(pg, svadd_f32_x(pg, svget4(a, 3), svget4(b, 3)), svadd_f32_x(pg, svget4(c2, 3), svget4(d, 3)));
    svst1_f32_x4(pc, y + 16 * j, svcreate4(o0, o1, o2, o3));
  }
}
static void tr_body(const float *x, float *P, long rs, long rows) __arm_streaming __arm_inout("za") {
  svbool_t pg = svptrue_b32();
  svcount_t pc = svptrue_c32();
  long r0 = 0, pfd = SME_PF;
  for (; r0 + 64 <= rows; ) {
    const float *xr = x + 16 * r0;
    TBLK(0) r0 += 16; xr += 256;
    TBLK(1) r0 += 16; xr += 256;
    TBLK(2) r0 += 16; xr += 256;
    TBLK(3) r0 += 16;
  }
  for (; r0 < rows; r0 += 16) { const float *xr = x + 16 * r0; TBLK(0) }
}
static void small_body(const float *P, long rs, const float *T, int Cn, float *y) __arm_streaming __arm_inout("za") {
  svbool_t pg = svptrue_b32();
  svcount_t pc = svptrue_c32();
  int ng = Cn / 4, gq = (ng + 3) / 4;
  svzero_za();
  /* tile t takes groups [t*gq, (t+1)*gq); groups past ng use a zero Toeplitz block (T is padded) */
  for (int kk = 0; kk < gq; kk++) {
    int g0 = kk, g1 = gq + kk, g2 = 2 * gq + kk, g3 = 3 * gq + kk;
    svfloat32x4_t A0 = svld1_f32_x4(pc, T + 64 * g0), A1 = svld1_f32_x4(pc, T + 64 * (g1 < ng ? g1 : ng));
    svfloat32x4_t A2 = svld1_f32_x4(pc, T + 64 * (g2 < ng ? g2 : ng)), A3 = svld1_f32_x4(pc, T + 64 * (g3 < ng ? g3 : ng));
#define SB(g, u) svld1_f32(pg, P + ((4 * (g)) & 15) * rs + (u) * rs + ((4 * (g)) >> 4))
#define SU(u) \
    svmopa_za32_f32_m(0, pg, pg, svget4(A0, u), SB(g0, u)); \
    svmopa_za32_f32_m(1, pg, pg, svget4(A1, u), SB(g1 < ng ? g1 : ng, u)); \
    svmopa_za32_f32_m(2, pg, pg, svget4(A2, u), SB(g2 < ng ? g2 : ng, u)); \
    svmopa_za32_f32_m(3, pg, pg, svget4(A3, u), SB(g3 < ng ? g3 : ng, u));
    SU(0) SU(1) SU(2) SU(3)
  }
  for (int j = 0; j < 16; j += 4) {
    svfloat32x4_t a = svread_ver_za32_f32_vg4(0, j), b = svread_ver_za32_f32_vg4(1, j);
    svfloat32x4_t c2 = svread_ver_za32_f32_vg4(2, j), d = svread_ver_za32_f32_vg4(3, j);
    svfloat32_t o0 = svadd_f32_x(pg, svadd_f32_x(pg, svget4(a, 0), svget4(b, 0)), svadd_f32_x(pg, svget4(c2, 0), svget4(d, 0)));
    svfloat32_t o1 = svadd_f32_x(pg, svadd_f32_x(pg, svget4(a, 1), svget4(b, 1)), svadd_f32_x(pg, svget4(c2, 1), svget4(d, 1)));
    svfloat32_t o2 = svadd_f32_x(pg, svadd_f32_x(pg, svget4(a, 2), svget4(b, 2)), svadd_f32_x(pg, svget4(c2, 2), svget4(d, 2)));
    svfloat32_t o3 = svadd_f32_x(pg, svadd_f32_x(pg, svget4(a, 3), svget4(b, 3)), svadd_f32_x(pg, svget4(c2, 3), svget4(d, 3)));
    svst1_f32_x4(pc, y + 16 * j, svcreate4(o0, o1, o2, o3));
  }
}

/* one streaming-mode entry per short block: ZA transpose of the full rows, then the taps-split pass */
__arm_locally_streaming __arm_new("za") static void sme_small_fused(const float *x, float *P, long rs, long rz,
                                                                     const float *T, int Cn, float *y) {
  if (rz) tr_body(x, P, rs, rz);
  small_body(P, rs, T, Cn, y);
}
static int SME_FUSE = 1; /* env SME_FUSE=0: separate transpose and kernel calls for short blocks */
static int SME_SMALL = 1; /* env SME_SMALL=0 disables the small-block pass */
static void run_sme_konly(void *vc, long i) { /* timing only: kernel on a fixed P (results wrong) */
  ctx_t *c = vc;
  float *y = bench_pool_get(c->yp, i);
  long n = c->n;
  for (long m = 0; m < n; m += c->ch) {
    long cn = c->ch < n - m ? c->ch : n - m;
    sme_kernel(c->P, c->rs, c->T, c->Cn, y + m, cn / 16);
  }
}
static void do_transpose(ctx_t *c, const float *x, long xvalid, long rows) {
  if (SME_TR) {
    long rfull = xvalid >= 16 ? (xvalid - 16) / 16 + 1 : 0; /* rows fully inside x */
    long rz = rows < rfull ? rows : rfull;
    rz &= ~15L;
    if (rz) sme_transpose(x, c->P, c->rs, rz);
    if (rz < rows) { /* tail rows: rows inside x and the one partial row element by element, then zeros */
      long rlast = (xvalid + 15) / 16; /* rows touching valid x */
      if (rlast > rows) rlast = rows;
      for (long r = rz; r < rlast; r++)
        for (int p = 0; p < 16; p++) { long ix = 16 * r + p; c->P[p * c->rs + r] = ix < xvalid ? x[ix] : 0.f; }
      long r0z = rlast > rz ? rlast : rz;
      if (r0z < rows)
        for (int p = 0; p < 16; p++) memset(c->P + p * c->rs + r0z, 0, 4 * (rows - r0z));
    }
  } else transpose16(x, xvalid, c->P, c->rs, rows);
}
static void run_sme_tonly(void *vc, long i) { /* timing only: the transposes alone */
  ctx_t *c = vc;
  const float *x = bench_pool_get(c->xp, i);
  long n = c->n, xlen = n + c->L - 1, qmax = (c->Cn - 1) / 16;
  for (long m = 0; m < n; m += c->ch) {
    long cn = c->ch < n - m ? c->ch : n - m;
    do_transpose(c, x + m, xlen - m, cn / 16 + qmax + 1);
  }
}
static void run_sme(void *vc, long i) {
  ctx_t *c = vc;
  const float *x = bench_pool_get(c->xp, i);
  float *y = bench_pool_get(c->yp, i);
  long n = c->n, xlen = n + c->L - 1, qmax = (c->Cn - 1) / 16;
  for (long m = 0; m < n; m += c->ch) {
    long cn = c->ch < n - m ? c->ch : n - m;
    long rows = (cn / 1024) * 64, rem = cn - rows * 16; /* full passes, then a partial pass into scratch */
    if (rows == 0 && rem <= 256 && SME_SMALL && SME_FUSE && SME_TR) { /* short block: one streaming entry */
      long need = 16 + qmax + 1, xv = xlen - m;
      long rfull = xv >= 16 ? (xv - 16) / 16 + 1 : 0;
      long rz = (need < rfull ? need : rfull) & ~15L;
      long rlast = (xv + 15) / 16; if (rlast > need) rlast = need;
      for (long r = rz; r < rlast; r++)
        for (int p = 0; p < 16; p++) { long ix = 16 * r + p; c->P[p * c->rs + r] = ix < xv ? x[m + ix] : 0.f; }
      long r0z = rlast > rz ? rlast : rz;
      if (r0z < need) for (int p = 0; p < 16; p++) memset(c->P + p * c->rs + r0z, 0, 4 * (need - r0z));
      sme_small_fused(x + m, c->P, c->rs, rz, c->T, c->Cn, c->ys);
      memcpy(y + m, c->ys, 4 * rem);
      continue;
    }
    do_transpose(c, x + m, xlen - m, rows + (rem ? (rem <= 256 && SME_SMALL ? 16 : 64) : 0) + qmax + 1);
    if (rows) sme_kernel(c->P, c->rs, c->T, c->Cn, y + m, rows);
    if (rem && rem <= 256 && SME_SMALL) {
      sme_kernel_small(c->P + rows, c->rs, c->T, c->Cn, c->ys);
      memcpy(y + m + rows * 16, c->ys, 4 * rem);
    } else if (rem) {
      sme_kernel(c->P + rows, c->rs, c->T, c->Cn, c->ys, 64);
      memcpy(y + m + rows * 16, c->ys, 4 * rem);
    }
  }
}


/* ---------------- SME kernel as a reusable filter object + fast FIR algorithm (Karatsuba, 2-parallel FFA) ---------------- */
typedef struct ffa_s {
  int L, level;
  /* level 0: direct SME */
  float *T; int Cn; long ch, rs; float *P;
  /* level >= 1 */
  struct ffa_s *sub[3];
  float *x0, *x1, *sx, *A, *B, *C, *g0;
} ffa_t;

static void za_transpose(float *P, long rs, const float *x, long xvalid, long rows) {
  long rfull = xvalid >= 16 ? (xvalid - 16) / 16 + 1 : 0;
  long rz = rows < rfull ? rows : rfull;
  rz &= ~15L;
  if (rz) sme_transpose(x, P, rs, rz);
  for (long r = rz; r < rows; r++)
    for (int p = 0; p < 16; p++) { long ix = 16 * r + p; P[p * rs + r] = ix < xvalid ? x[ix] : 0.f; }
}

/* g: correlation taps (y[m] = sum_k g[k] x[m+k]); nmax: max outputs per call */
static ffa_t *ffa_new(const float *g, int L, long nmax, int level) {
  ffa_t *f = calloc(1, sizeof *f);
  f->L = L; f->level = level;
  if (level == 0) {
    f->Cn = (L + 15 + 3) & ~3;
    f->T = bench_alloc(4 * 16 * (size_t)f->Cn);
    for (int cc = 0; cc < f->Cn; cc++)
      for (int i = 0; i < 16; i++) { int k = cc - i; f->T[16 * cc + i] = (k >= 0 && k < L) ? g[k] : 0.f; }
    long ch = SME_CH ? SME_CH : 16384;
    if (ch < 4L * L) ch = 4L * L;
    ch = (ch + 1023) & ~1023L;
    if (ch > nmax) ch = nmax;
    f->ch = ch;
    int qmax = (f->Cn - 1) / 16;
    f->rs = ((ch / 16 + qmax + 1) + 31) & ~31L;
    f->P = bench_alloc(4 * 16 * (size_t)f->rs);
    return f;
  }
  int L2 = L / 2;
  float *g0 = bench_alloc(4 * L2), *g1 = bench_alloc(4 * L2), *g01 = bench_alloc(4 * L2);
  for (int b = 0; b < L2; b++) { g0[b] = g[2 * b]; g1[b] = g[2 * b + 1]; g01[b] = g0[b] + g1[b]; }
  f->sub[0] = ffa_new(g0, L2, nmax / 2, level - 1);
  f->sub[1] = ffa_new(g1, L2, nmax / 2, level - 1);
  f->sub[2] = ffa_new(g01, L2, nmax / 2, level - 1);
  f->g0 = g0;
  size_t hl = nmax / 2 + L2 + 64;
  f->x0 = bench_alloc(4 * hl); f->x1 = bench_alloc(4 * hl); f->sx = bench_alloc(4 * hl);
  f->A = bench_alloc(4 * hl); f->B = bench_alloc(4 * hl); f->C = bench_alloc(4 * hl);
  return f;
}

/* y[m] = sum_k g[k] x[m+k], m < n; x has n + L - 1 samples; n multiple of 1024 * 2^level */
static void ffa_run(ffa_t *f, const float *x, float *y, long n) {
  int L = f->L;
  long xlen = n + L - 1;
  if (f->level == 0) {
    int qmax = (f->Cn - 1) / 16;
    for (long m = 0; m < n; m += f->ch) {
      long cn = f->ch < n - m ? f->ch : n - m;
      za_transpose(f->P, f->rs, x + m, xlen - m, cn / 16 + qmax + 1);
      sme_kernel(f->P, f->rs, f->T, f->Cn, y + m, cn / 16);
    }
    return;
  }
  int L2 = L / 2;
  long h = n / 2;
  /* x0[j] = x[2j] (h + L2 samples), x1[j] = x[2j+1] (h + L2 - 1 samples); xlen is odd (L even) */
  DSPSplitComplex xs = {f->x0, f->x1};
  vDSP_ctoz((const DSPComplex *)x, 2, &xs, 1, xlen / 2);
  f->x0[xlen / 2] = x[xlen - 1];
  /* sx[j] = x1[j] + x0[j+1] */
  vDSP_vadd(f->x1, 1, f->x0 + 1, 1, f->sx, 1, h + L2 - 1);
  ffa_run(f->sub[0], f->x0, f->A, h);
  { double acc = 0; for (int b = 0; b < L2; b++) acc += (double)f->g0[b] * f->x0[h + b]; f->A[h] = (float)acc; }
  ffa_run(f->sub[1], f->x1, f->B, h);
  ffa_run(f->sub[2], f->sx, f->C, h);
  /* y[2a] = A[a] + B[a]; y[2a+1] = C[a] - A[a+1] - B[a] */
  const float *A = f->A, *B = f->B, *C = f->C;
  for (long a = 0; a < h; a += 4) {
    float32x4_t va = vld1q_f32(A + a), va1 = vld1q_f32(A + a + 1), vb = vld1q_f32(B + a), vc = vld1q_f32(C + a);
    float32x4x2_t o;
    o.val[0] = vaddq_f32(va, vb);
    o.val[1] = vsubq_f32(vsubq_f32(vc, va1), vb);
    vst2q_f32(y + 2 * a, o);
  }
}
static int FFA_LEVEL = 1; /* env FFA */
static ffa_t *FFA_OBJ = NULL;
static void run_sme_ffa(void *vc, long i) {
  ctx_t *c = vc;
  const float *x = bench_pool_get(c->xp, i);
  float *y = bench_pool_get(c->yp, i);
  ffa_run(FFA_OBJ, x, y, c->n);
}

/* ---------------- main ---------------- */
static void make_ctx(ctx_t *c, int L, long n, bench_pool *xp, bench_pool *yp, const float *h, int ols_N) {
  memset(c, 0, sizeof *c);
  c->L = L; c->n = n; c->xp = xp; c->yp = yp; c->h = h;
  sme_setup(c);
  ols_setup(c, ols_N);
  pf_setup(c, OLS_PF_N ? OLS_PF_N : default_pf_N(L));
}

static int default_ols_N(int L) {
  /* tuned by 'olsscan' (see ols_scan.txt); fallback 8L, at least 1024 */
  int N = 1024;
  while (N < 8 * L) N *= 2;
  return N;
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: fir check|time|olsscan L LOGN [same|separate] [mask]\n"); return 1; }
  const char *mode = argv[1];
  int L = atoi(argv[2]), logn = atoi(argv[3]);
  long n = (1L << logn) + (getenv("N_EXTRA") ? atol(getenv("N_EXTRA")) : 0);
  if (getenv("OLS_N")) OLS_N_FORCE = atoi(getenv("OLS_N"));
  if (getenv("SME_CH")) SME_CH = atol(getenv("SME_CH"));
  if (getenv("FFA")) FFA_LEVEL = atoi(getenv("FFA"));
  if (getenv("SME_SMALL")) SME_SMALL = atoi(getenv("SME_SMALL"));
  if (getenv("SME_FUSE")) SME_FUSE = atoi(getenv("SME_FUSE"));
  if (getenv("SME_PF")) SME_PF = atol(getenv("SME_PF"));
  if (getenv("OLS_PF_N")) OLS_PF_N = atoi(getenv("OLS_PF_N"));
  if (getenv("SME_TR")) SME_TR = atoi(getenv("SME_TR"));
  bench_init();
  float *h = bench_alloc(4 * (size_t)L);
  bench_fill_f(h, L);
  size_t xbytes = 4 * (size_t)(n + L - 1), ybytes = 4 * (size_t)n;
  int sep = argc > 4 && !strcmp(argv[4], "separate");
  bench_pool xp, yp;
  bench_pool_init(&xp, xbytes, (!strcmp(mode, "time") && sep) ? BENCH_SEPARATE : BENCH_SAME);
  bench_pool_init(&yp, ybytes, (!strcmp(mode, "time") && sep) ? BENCH_SEPARATE : BENCH_SAME);
  for (int b = 0; b < xp.count; b++) bench_fill_f(xp.buf[b], n + L - 1);
  int N = OLS_N_FORCE ? OLS_N_FORCE : default_ols_N(L);
  ctx_t c;
  make_ctx(&c, L, n, &xp, &yp, h, N);
  {
    float *g = bench_alloc(4 * (size_t)L);
    for (int k = 0; k < L; k++) g[k] = h[L - 1 - k];
    FFA_OBJ = ffa_new(g, L, n, FFA_LEVEL);
  }
  printf("# L=%d n=2^%d ols_N=%d sme_chunk=%ld Cn=%d useful MAC fraction of SME %.4f\n", L, logn, N, c.ch, c.Cn,
         (double)L / c.Cn);

  if (!strcmp(mode, "check")) {
    const float *x = xp.buf[0];
    /* reference: all outputs if n*L <= 2^28, else first 4096, last 4096 and 8192 random positions */
    long nchk = 0; long *idx = malloc(sizeof(long) * (n < 16384 ? n : 16384) + 8);
    if ((double)n * L <= 268435456.0) { free(idx); idx = malloc(sizeof(long) * n); for (long m = 0; m < n; m++) idx[nchk++] = m; }
    else {
      for (long m = 0; m < 4096; m++) idx[nchk++] = m;
      for (long m = n - 4096; m < n; m++) idx[nchk++] = m;
      for (int k = 0; k < 8192; k++) idx[nchk++] = (long)(bench_rand_u64() % (uint64_t)n);
    }
    double *ref = malloc(8 * nchk);
    for (long k = 0; k < nchk; k++) ref[k] = ref_one(x, h, L, idx[k]);
    float *got = malloc(4 * nchk);
    void (*fn[6])(void *, long) = {run_vdsp, run_ols, run_sme, run_ols2, run_ols_pf, run_sme_ffa};
    const char *nm[6] = {"vdsp_conv", "ols_vdsp", "sme_toeplitz", "ols2_vdsp", "ols_pffft", "sme_ffa"};
    char lab[96];
    for (int f = 0; f < 6; f++) {
      if (f == 5 && n % (1024L << FFA_LEVEL)) continue; /* Karatsuba variant needs n multiple of 1024*2^level */
      memset(yp.buf[0], 0xff, ybytes); /* NaN fill: catches unwritten outputs */
      fn[f](&c, 0);
      const float *y = yp.buf[0];
      for (long k = 0; k < nchk; k++) got[k] = y[idx[k]];
      snprintf(lab, sizeof lab, "%s L=%d n=2^%d (%ld outs)", nm[f], L, logn, nchk);
      bench_check_f(lab, got, ref, nchk, 1e-5);
      /* edges explicitly: first L and last L outputs */
      long ne = L < n ? L : n;
      double *re = malloc(8 * ne); float *ge = malloc(4 * ne);
      for (long k = 0; k < ne; k++) { re[k] = ref_one(x, h, L, k); ge[k] = y[k]; }
      snprintf(lab, sizeof lab, "%s first %ld outputs", nm[f], ne);
      bench_check_f(lab, ge, re, ne, 1e-5);
      for (long k = 0; k < ne; k++) { re[k] = ref_one(x, h, L, n - ne + k); ge[k] = y[n - ne + k]; }
      snprintf(lab, sizeof lab, "%s last %ld outputs", nm[f], ne);
      bench_check_f(lab, ge, re, ne, 1e-5);
      free(re); free(ge);
    }
    return 0;
  }
  if (!strcmp(mode, "olsscan")) {
    bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)n; o.unit_name = "output"; o.min_reps = 9; o.print_reps = 0;
    for (int NN = 256; NN <= 32768; NN *= 2) {
      if (NN < 2 * L) continue;
      ctx_t cc; make_ctx(&cc, L, n, &xp, &yp, h, NN);
      OLS_PF_N = NN;
      ctx_t cc2; make_ctx(&cc2, L, n, &xp, &yp, h, NN);
      bench_prog pr[2] = {{"ols2_vdsp", run_ols2, NULL, &cc}, {"ols_pffft", run_ols_pf, NULL, &cc2}};
      bench_stat st[2]; char t[64]; snprintf(t, sizeof t, "olsscan L=%d N=%d", L, NN);
      bench_compare(t, pr, 2, &o, st);
    }
    return 0;
  }
  int mask = argc > 5 ? atoi(argv[5]) : 75; /* 1 ols2, 2 vdsp_conv, 4 ols(simple), 8 sme, 16 sme kernel only, 32 transpose only, 64 ols_pffft, 128 sme_ffa (env FFA=levels) */
  bench_prog all[8] = {{"ols2_vdsp", run_ols2, NULL, &c}, {"vdsp_conv", run_vdsp, NULL, &c}, {"ols_vdsp", run_ols, NULL, &c}, {"sme_toeplitz", run_sme, NULL, &c},
                       {"sme_kernel_only", run_sme_konly, NULL, &c}, {"sme_transpose_only", run_sme_tonly, NULL, &c}, {"ols_pffft", run_ols_pf, NULL, &c}, {"sme_ffa", run_sme_ffa, NULL, &c}};
  bench_prog pr[8]; int np = 0;
  for (int f = 0; f < 8; f++) if (mask & (1 << f)) pr[np++] = all[f];
  bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)n; o.unit_name = "output";
  bench_stat st[8];
  char t[96]; snprintf(t, sizeof t, "fir L=%d n=2^%d %s olsN=%d", L, logn, sep ? "separate" : "same", N);
  bench_compare(t, pr, np, &o, st);
  for (int f = 0; f < np; f++)
    printf("GMACS %-40s %-14s best %8.4f ns/out = %7.1f GMAC/s (median %7.1f)\n", t, st[f].name, st[f].best,
           L / st[f].best, L / st[f].median);
  return 0;
}
