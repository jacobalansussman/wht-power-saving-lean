/* kern.h -- am-2 (wave 12): aj-2's M1 x 16 SME2 complex FFT (wave-11/aj-2/fftM16.c, copied and changed)
 * turned into a two-real-frames STFT kernel. Include with M1 (17..32) and SFX defined.
 * Change from aj-2: transform p reads its real part from sig + 2*hop*p and its imaginary part from
 * sig + 2*hop*p + hop (two consecutive STFT frames straight from the signal: no framing copy), and writes
 * Z = [re(N), im(N)] per pair to out + 2*N*p. Prefetch into L2, PFD pairs ahead.
 * N = 16*M1, n = 16 n1 + n2, k = k1 + M1 k2; tables and passes exactly as aj-2's fftM16.c. */
#define KN (16 * M1)
#define CAT_(a, b) a##b
#define CAT(a, b) CAT_(a, b)
#define FAt CAT(FAtab_, SFX)
#define FBt CAT(FBtab_, SFX)
#define W4t CAT(W4tab_, SFX)
#define S4t CAT(S4_, SFX)
#ifndef PFD
#define PFD 8
#endif
static float FAt[32 * 64] __attribute__((aligned(128)));
static float FBt[16 * 32] __attribute__((aligned(128)));
static float W4t[2 * 4 * 128] __attribute__((aligned(128)));
static _Thread_local float S4t[4 * 256] __attribute__((aligned(128)));

static void CAT(make_tables_, SFX)(void) {
  for (int n1 = 0; n1 < M1; n1++)
    for (int k = 0; k < 32; k++) {
      double th = 2 * M_PI * k * n1 / (double)M1;
      FAt[64 * n1 + k] = k < M1 ? (float)cos(th) : 0.f;
      FAt[64 * n1 + 32 + k] = k < M1 ? (float)(-sin(th)) : 0.f;
    }
  for (int n2 = 0; n2 < 16; n2++)
    for (int k = 0; k < 16; k++) {
      double th = 2 * M_PI * k * n2 / 16.0;
      FBt[32 * n2 + k] = (float)cos(th);
      FBt[32 * n2 + 16 + k] = (float)(-sin(th));
    }
  for (int b = 0; b < 2; b++)
    for (int j = 0; j < 4; j++)
      for (int m = 0; m < 4; m++) {
        int n2 = j + 4 * m;
        for (int c = 0; c < 16; c++) {
          int k1 = 16 * b + c;
          double th = 2 * M_PI * k1 * n2 / (double)KN;
          W4t[512 * b + 128 * j + 16 * m + c] = k1 < M1 ? (float)(-2.0 * sin(th / 2) * sin(th / 2)) : 0.f;
          W4t[512 * b + 128 * j + 64 + 16 * m + c] = k1 < M1 ? (float)(-sin(th)) : 0.f;
        }
      }
}

#ifndef TWB_DEFINED
#define TWB_DEFINED
#define TWB(j, t0, wbase)                                              \
  {                                                                    \
    svfloat32x4_t a = svread_za32_f32_vg1x4(4 * (j) + (t0));           \
    svfloat32x4_t c = svread_za32_f32_vg1x4(4 * (j) + (t0) + 1);       \
    svfloat32x4_t wm = svld1_f32_x4(pc, (wbase) + 128 * (j));          \
    svfloat32x4_t ws = svld1_f32_x4(pc, (wbase) + 128 * (j) + 64);     \
    svmla_za32_f32_vg1x4(4 * (j) + (t0), a, wm);                       \
    svmls_za32_f32_vg1x4(4 * (j) + (t0), c, ws);                       \
    svmla_za32_f32_vg1x4(4 * (j) + (t0) + 1, a, ws);                   \
    svmla_za32_f32_vg1x4(4 * (j) + (t0) + 1, c, wm);                   \
  }
#endif

__arm_locally_streaming __arm_new("za") static void CAT(sme_pairs_, SFX)(const float *sig, long hop, float *out, long B) {
  svbool_t pg = svptrue_b32(), p9 = svwhilelt_b32(0, M1 - 16);
  svcount_t pc = svptrue_c32();
  const long span = KN + hop; /* floats touched by one pair */
  for (int q = 1; q < PFD && q < B; q++)
    for (long l = 0; l < span; l += 32) __builtin_prefetch(sig + 2 * hop * q + l, 0, 2);
  for (long b = 0; b < B; b++) {
    const float *xr = sig + 2 * hop * b, *xi = xr + hop;
    float *yr = out + 2 * KN * b, *yi = yr + KN;
    if (PFD && b + PFD < B)
      for (long l = 0; l < span; l += 32) __builtin_prefetch(sig + 2 * hop * (b + PFD) + l, 0, 2);
    svzero_mask_za(0xff);
    {
      svfloat32x4_t f = svld1_f32_x4(pc, FAt);
      svfloat32_t ar = svld1_f32(pg, xr), ai = svld1_f32(pg, xi);
      svmopa_za32_f32_m(0, pg, pg, ar, svget4(f, 0)); svmopa_za32_f32_m(1, pg, pg, ai, svget4(f, 0));
      svmopa_za32_f32_m(2, pg, pg, ar, svget4(f, 1)); svmopa_za32_f32_m(3, pg, pg, ai, svget4(f, 1));
    }
    for (int n1 = 1; n1 < M1; n1++) {
      svfloat32x4_t f = svld1_f32_x4(pc, FAt + 64 * n1);
      svfloat32_t ar = svld1_f32(pg, xr + 16 * n1), ai = svld1_f32(pg, xi + 16 * n1);
      svfloat32_t fra = svget4(f, 0), frb = svget4(f, 1), fia = svget4(f, 2), fib = svget4(f, 3);
      svmopa_za32_f32_m(0, pg, pg, ar, fra); svmopa_za32_f32_m(1, pg, pg, ar, fia);
      svmopa_za32_f32_m(2, pg, pg, ar, frb); svmopa_za32_f32_m(3, pg, pg, ar, fib);
      svmops_za32_f32_m(0, pg, pg, ai, fia); svmopa_za32_f32_m(1, pg, pg, ai, fra);
      svmops_za32_f32_m(2, pg, pg, ai, fib); svmopa_za32_f32_m(3, pg, pg, ai, frb);
    }
    TWB(0, 0, W4t) TWB(1, 0, W4t) TWB(2, 0, W4t) TWB(3, 0, W4t)
    TWB(0, 2, W4t + 512) TWB(1, 2, W4t + 512) TWB(2, 2, W4t + 512) TWB(3, 2, W4t + 512)
    for (int n2 = 0; n2 < 16; n2++) {
      svstr_za(4 * n2 + 0, S4t + 16 * n2); svstr_za(4 * n2 + 1, S4t + 256 + 16 * n2);
      svstr_za(4 * n2 + 2, S4t + 512 + 16 * n2); svstr_za(4 * n2 + 3, S4t + 768 + 16 * n2);
    }
    svzero_mask_za(0xff);
    for (int n2 = 0; n2 < 16; n2++) {
      svfloat32x2_t f = svld1_f32_x2(pc, FBt + 32 * n2);
      svfloat32_t fr = svget2(f, 0), fi = svget2(f, 1);
      svfloat32_t yra = svld1_f32(pg, S4t + 16 * n2), yia = svld1_f32(pg, S4t + 256 + 16 * n2);
      svfloat32_t yrb = svld1_f32(pg, S4t + 512 + 16 * n2), yib = svld1_f32(pg, S4t + 768 + 16 * n2);
      svmopa_za32_f32_m(0, pg, pg, fr, yra); svmopa_za32_f32_m(1, pg, pg, fr, yia);
      svmopa_za32_f32_m(2, pg, pg, fr, yrb); svmopa_za32_f32_m(3, pg, pg, fr, yib);
      if (n2 != 0 && n2 != 8) {
        svmops_za32_f32_m(0, pg, pg, fi, yia); svmopa_za32_f32_m(1, pg, pg, fi, yra);
        svmops_za32_f32_m(2, pg, pg, fi, yib); svmopa_za32_f32_m(3, pg, pg, fi, yrb);
      }
    }
    for (int k2 = 0; k2 < 16; k2++) {
      svst1_hor_za32(0, k2, pg, yr + M1 * k2); svst1_hor_za32(1, k2, pg, yi + M1 * k2);
      svst1_hor_za32(2, k2, p9, yr + M1 * k2 + 16); svst1_hor_za32(3, k2, p9, yi + M1 * k2 + 16);
    }
  }
}
__arm_locally_streaming __arm_new("za") static void CAT(sme_pairs_w_, SFX)(const float *sig, long hop, const float *win, float *out, long B) {
  svbool_t pg = svptrue_b32(), p9 = svwhilelt_b32(0, M1 - 16);
  svcount_t pc = svptrue_c32();
  const long span = KN + hop; /* floats touched by one pair */
  for (int q = 1; q < PFD && q < B; q++)
    for (long l = 0; l < span; l += 32) __builtin_prefetch(sig + 2 * hop * q + l, 0, 2);
  for (long b = 0; b < B; b++) {
    const float *xr = sig + 2 * hop * b, *xi = xr + hop;
    float *yr = out + 2 * KN * b, *yi = yr + KN;
    if (PFD && b + PFD < B)
      for (long l = 0; l < span; l += 32) __builtin_prefetch(sig + 2 * hop * (b + PFD) + l, 0, 2);
    svzero_mask_za(0xff);
    {
      svfloat32x4_t f = svld1_f32_x4(pc, FAt);
      svfloat32_t w0 = svld1_f32(pg, win); svfloat32_t ar = svmul_f32_x(pg, svld1_f32(pg, xr), w0), ai = svmul_f32_x(pg, svld1_f32(pg, xi), w0);
      svmopa_za32_f32_m(0, pg, pg, ar, svget4(f, 0)); svmopa_za32_f32_m(1, pg, pg, ai, svget4(f, 0));
      svmopa_za32_f32_m(2, pg, pg, ar, svget4(f, 1)); svmopa_za32_f32_m(3, pg, pg, ai, svget4(f, 1));
    }
    for (int n1 = 1; n1 < M1; n1++) {
      svfloat32x4_t f = svld1_f32_x4(pc, FAt + 64 * n1);
      svfloat32_t wn = svld1_f32(pg, win + 16 * n1); svfloat32_t ar = svmul_f32_x(pg, svld1_f32(pg, xr + 16 * n1), wn), ai = svmul_f32_x(pg, svld1_f32(pg, xi + 16 * n1), wn);
      svfloat32_t fra = svget4(f, 0), frb = svget4(f, 1), fia = svget4(f, 2), fib = svget4(f, 3);
      svmopa_za32_f32_m(0, pg, pg, ar, fra); svmopa_za32_f32_m(1, pg, pg, ar, fia);
      svmopa_za32_f32_m(2, pg, pg, ar, frb); svmopa_za32_f32_m(3, pg, pg, ar, fib);
      svmops_za32_f32_m(0, pg, pg, ai, fia); svmopa_za32_f32_m(1, pg, pg, ai, fra);
      svmops_za32_f32_m(2, pg, pg, ai, fib); svmopa_za32_f32_m(3, pg, pg, ai, frb);
    }
    TWB(0, 0, W4t) TWB(1, 0, W4t) TWB(2, 0, W4t) TWB(3, 0, W4t)
    TWB(0, 2, W4t + 512) TWB(1, 2, W4t + 512) TWB(2, 2, W4t + 512) TWB(3, 2, W4t + 512)
    for (int n2 = 0; n2 < 16; n2++) {
      svstr_za(4 * n2 + 0, S4t + 16 * n2); svstr_za(4 * n2 + 1, S4t + 256 + 16 * n2);
      svstr_za(4 * n2 + 2, S4t + 512 + 16 * n2); svstr_za(4 * n2 + 3, S4t + 768 + 16 * n2);
    }
    svzero_mask_za(0xff);
    for (int n2 = 0; n2 < 16; n2++) {
      svfloat32x2_t f = svld1_f32_x2(pc, FBt + 32 * n2);
      svfloat32_t fr = svget2(f, 0), fi = svget2(f, 1);
      svfloat32_t yra = svld1_f32(pg, S4t + 16 * n2), yia = svld1_f32(pg, S4t + 256 + 16 * n2);
      svfloat32_t yrb = svld1_f32(pg, S4t + 512 + 16 * n2), yib = svld1_f32(pg, S4t + 768 + 16 * n2);
      svmopa_za32_f32_m(0, pg, pg, fr, yra); svmopa_za32_f32_m(1, pg, pg, fr, yia);
      svmopa_za32_f32_m(2, pg, pg, fr, yrb); svmopa_za32_f32_m(3, pg, pg, fr, yib);
      if (n2 != 0 && n2 != 8) {
        svmops_za32_f32_m(0, pg, pg, fi, yia); svmopa_za32_f32_m(1, pg, pg, fi, yra);
        svmops_za32_f32_m(2, pg, pg, fi, yib); svmopa_za32_f32_m(3, pg, pg, fi, yrb);
      }
    }
    for (int k2 = 0; k2 < 16; k2++) {
      svst1_hor_za32(0, k2, pg, yr + M1 * k2); svst1_hor_za32(1, k2, pg, yi + M1 * k2);
      svst1_hor_za32(2, k2, p9, yr + M1 * k2 + 16); svst1_hor_za32(3, k2, p9, yi + M1 * k2 + 16);
    }
  }
}
#undef KN
#undef FAt
#undef FBt
#undef W4t
#undef S4t
