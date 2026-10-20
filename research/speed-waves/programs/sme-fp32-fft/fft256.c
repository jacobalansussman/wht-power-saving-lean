/* fft256.c -- aj-2 (wave 11): batched 256-point complex FFT (fp32, forward, unscaled, e^{-i}) on the
 * M4's SME2 matrix unit, against vDSP (DFT_Execute, fft_zop, fftm_zop) and PFFFT z-order.
 *
 * Layout (same as vDSP split complex): a batch of B transforms, re[B*256] and im[B*256], transform b at
 * offset 256*b, natural order in and out. Out of place.
 *
 * Kernel v1 ("sme1"): four-step 256 = 16 x 16, n = 16 n1 + n2, k = k1 + 16 k2, per transform:
 *   pass A  Y = F16 * X        (X[n1][n2] = x[16 n1 + n2]; rows of X are contiguous 16-float rows)
 *           tiles ZA0 = Re Y, ZA1 = Im Y (rows k1, cols n2), outer products F[:,n1] (x) X[n1,:]
 *   twiddle Y'[k1][n2] = Y[k1][n2] * w256^(k1 n2), IN PLACE in ZA with multi-vector FMLA/FMLS on ZA
 *           array vectors: Y' = Y + Y*(W-1) (re/im mixed), read by 4-row groups (MOVA vg1x4)
 *   pass B  Z = Y' * F16       (tiles ZA2 = Re Z, ZA3 = Im Z, rows k1, cols k2), left operands are
 *           vertical slices of ZA0/ZA1 (columns n2 of Y'), right operands rows of F16
 *   output  out[k1 + 16 k2] = Z[k1][k2]: vertical slice k2 of ZA2/ZA3 is out[16 k2 .. 16 k2 + 15],
 *           stored straight from ZA (ST1W za vertical slice): NATURAL order.
 * F16 is symmetric; columns n and 16-n share Re and negate Im, so 16 distinct vectors are kept.
 *
 * Usage: fft256 check | time same|separate B...   Run timings under the timing lock.
 */
#include "../timing-harness/bench.h"
#include <Accelerate/Accelerate.h>
#include <arm_sme.h>
#include "pffft.h"

#define N 256
static float Ftab[32 * 16];    /* Fr[n] (n=0..8) at 16n; Fi[n] (n=1..7) at 16*(8+n) ; padded */
static float Wtab[2 * 256];    /* twiddle groups: for j=0..3: Re(W)-1 rows k1=j,j+4,j+8,j+12 (64 floats), then Im rows */

static void make_tables(void) {
  for (int n = 0; n <= 8; n++)
    for (int k = 0; k < 16; k++) Ftab[16 * n + k] = (float)cos(2 * M_PI * k * n / 16.0);
  for (int n = 1; n <= 7; n++)
    for (int k = 0; k < 16; k++) Ftab[16 * (8 + n) + k] = (float)(-sin(2 * M_PI * k * n / 16.0));
  for (int j = 0; j < 4; j++)
    for (int m = 0; m < 4; m++) {
      int k1 = j + 4 * m;
      for (int n2 = 0; n2 < 16; n2++) {
        double th = 2 * M_PI * k1 * n2 / 256.0;
        Wtab[128 * j + 16 * m + n2] = (float)(-2.0 * sin(th / 2) * sin(th / 2)); /* cos - 1, accurately */
        Wtab[128 * j + 64 + 16 * m + n2] = (float)(-sin(th));
      }
    }
}

#include "fft256_gen.h"
#define LOADF \
  svfloat32_t fr0 = svld1_f32(pg, Ftab + 0), fr1 = svld1_f32(pg, Ftab + 16), fr2 = svld1_f32(pg, Ftab + 32), \
    fr3 = svld1_f32(pg, Ftab + 48), fr4 = svld1_f32(pg, Ftab + 64), fr5 = svld1_f32(pg, Ftab + 80), \
    fr6 = svld1_f32(pg, Ftab + 96), fr7 = svld1_f32(pg, Ftab + 112), fr8 = svld1_f32(pg, Ftab + 128); \
  svfloat32_t fi1 = svld1_f32(pg, Ftab + 144), fi2 = svld1_f32(pg, Ftab + 160), fi3 = svld1_f32(pg, Ftab + 176), \
    fi4 = svld1_f32(pg, Ftab + 192), fi5 = svld1_f32(pg, Ftab + 208), fi6 = svld1_f32(pg, Ftab + 224), \
    fi7 = svld1_f32(pg, Ftab + 240);
/* twiddle in place: group j = array vectors 4j (+16,+32,+48) = tile0 rows j,j+4,j+8,j+12; 4j+1 = same rows of tile1 */
#define TW(j)                                                          \
    {                                                                  \
      svfloat32x4_t a = svread_za32_f32_vg1x4(4 * (j));                \
      svfloat32x4_t c = svread_za32_f32_vg1x4(4 * (j) + 1);            \
      svfloat32x4_t wm = svld1_f32_x4(pc, Wtab + 128 * (j));           \
      svfloat32x4_t ws = svld1_f32_x4(pc, Wtab + 128 * (j) + 64);      \
      svmla_za32_f32_vg1x4(4 * (j), a, wm);                            \
      svmls_za32_f32_vg1x4(4 * (j), c, ws);                            \
      svmla_za32_f32_vg1x4(4 * (j) + 1, a, ws);                        \
      svmla_za32_f32_vg1x4(4 * (j) + 1, c, wm);                        \
    }

__arm_locally_streaming __arm_new("za") void fft256_sme1(const float *re, const float *im, float *ore, float *oim, long B) {
  svbool_t pg = svptrue_b32();
  svcount_t pc = svptrue_c32();
  LOADF
  for (long b = 0; b < B; b++) {
    const float *xr = re + N * b, *xi = im + N * b;
    float *yr = ore + N * b, *yi = oim + N * b;
    svzero_mask_za(0xff);
    PASSA_G0 PASSA_G1 PASSA_G2 PASSA_G3
    TW(0) TW(1) TW(2) TW(3)
    PASSB_G0 PASSB_G1 PASSB_G2 PASSB_G3
    for (int k2 = 0; k2 < 16; k2++) {
      svst1_ver_za32(2, k2, pg, yr + 16 * k2);
      svst1_ver_za32(3, k2, pg, yi + 16 * k2);
    }
  }
}


/* Kernel v2 ("sme2"): software pipeline. Pass B of transform t (tiles 2,3; left operands = columns of
 * Y'(t) loaded from a 4 KB scratch with 4-vector loads) is interleaved FMOPA by FMOPA with pass A of
 * transform t+1 (tiles 0,1), so all four tiles accumulate at once. Y'(t+1) (twiddled in place) is
 * then stored to the scratch by ZA vertical-slice stores; Z(t) is stored to the output the same way. */
#include "fft256_gen2.h"
static float Sbuf[2 * 512] __attribute__((aligned(128)));
#define STORE_YP(S) for (int n2 = 0; n2 < 16; n2++) { svst1_ver_za32(0, n2, pg, (S) + 16 * n2); svst1_ver_za32(1, n2, pg, (S) + 256 + 16 * n2); }
#define STORE_Z(yr, yi) for (int k2 = 0; k2 < 16; k2++) { svst1_ver_za32(2, k2, pg, (yr) + 16 * k2); svst1_ver_za32(3, k2, pg, (yi) + 16 * k2); }
__arm_locally_streaming __arm_new("za") void fft256_sme2(const float *re, const float *im, float *ore, float *oim, long B) {
  svbool_t pg = svptrue_b32();
  svcount_t pc = svptrue_c32();
  LOADF
  {
    const float *xr = re, *xi = im;
    svzero_mask_za(0x33);
    PASSA_G0 PASSA_G1 PASSA_G2 PASSA_G3
    TW(0) TW(1) TW(2) TW(3)
    STORE_YP(Sbuf)
  }
  for (long t = 0; t < B; t++) {
    const float *sr = Sbuf + 512 * (t & 1), *si = sr + 256;
    svzero_mask_za(0xcc);
    if (t + 1 < B) {
      const float *xr = re + N * (t + 1), *xi = im + N * (t + 1);
      svzero_mask_za(0x33);
      AB_G0 AB_G1 AB_G2 AB_G3
      TW(0) TW(1) TW(2) TW(3)
      float *S2 = Sbuf + 512 * ((t + 1) & 1);
      STORE_YP(S2)
    } else {
      BM_G0 BM_G1 BM_G2 BM_G3
    }
    STORE_Z(ore + N * t, oim + N * t)
  }
}


/* Kernel v3 ("sme3"): v2 with transposed tiles. Pass A computes Y^T (rows n2, cols k1), pass B computes
 * Z^T (rows k2, cols k1). Then the scratch rows (Y' columns) and the output rows out[16k2..16k2+15] are
 * HORIZONTAL slices, i.e. whole ZA array vectors, stored with STR ZA (array vector 4i+tile). The twiddle
 * table is the same because w^(k1 n2) is symmetric. */
#include "fft256_gen3.h"
#define STORE_YP3(S) for (int n2 = 0; n2 < 16; n2++) { svstr_za(4 * n2 + 0, (S) + 16 * n2); svstr_za(4 * n2 + 1, (S) + 256 + 16 * n2); }
#define STORE_Z3(yr, yi) for (int k2 = 0; k2 < 16; k2++) { svstr_za(4 * k2 + 2, (yr) + 16 * k2); svstr_za(4 * k2 + 3, (yi) + 16 * k2); }
__arm_locally_streaming __arm_new("za") void fft256_sme3(const float *re, const float *im, float *ore, float *oim, long B) {
  svbool_t pg = svptrue_b32();
  svcount_t pc = svptrue_c32();
  LOADF
  {
    const float *xr = re, *xi = im;
    svzero_mask_za(0x33);
    T_A_G0 T_A_G1 T_A_G2 T_A_G3
    TW(0) TW(1) TW(2) TW(3)
    STORE_YP3(Sbuf)
  }
  for (long t = 0; t < B; t++) {
    const float *sr = Sbuf + 512 * (t & 1), *si = sr + 256;
    svzero_mask_za(0xcc);
    if (t + 1 < B) {
      const float *xr = re + N * (t + 1), *xi = im + N * (t + 1);
      svzero_mask_za(0x33);
      T_AB_G0 T_AB_G1 T_AB_G2 T_AB_G3
      TW(0) TW(1) TW(2) TW(3)
      float *S2 = Sbuf + 512 * ((t + 1) & 1);
      STORE_YP3(S2)
    } else {
      T_B_G0 T_B_G1 T_B_G2 T_B_G3
    }
    STORE_Z3(ore + N * t, oim + N * t)
  }
}


/* Kernel v4/v5: v3 with a different store path. v4: both the Y' scratch and the output are stored as
 * MOVA 4-row reads (svread_hor_za32_f32_vg4) + 4-vector stores (ST1W x4). v5: output that way, Y' by STR. */
#define STX4(tile, ptr) { \
  svst1_f32_x4(pc, (ptr) + 0,   svread_hor_za32_f32_vg4(tile, 0)); svst1_f32_x4(pc, (ptr) + 64,  svread_hor_za32_f32_vg4(tile, 4)); \
  svst1_f32_x4(pc, (ptr) + 128, svread_hor_za32_f32_vg4(tile, 8)); svst1_f32_x4(pc, (ptr) + 192, svread_hor_za32_f32_vg4(tile, 12)); }
#define KER45(name, STYP, STZ) \
__arm_locally_streaming __arm_new("za") void name(const float *re, const float *im, float *ore, float *oim, long B) { \
  svbool_t pg = svptrue_b32(); svcount_t pc = svptrue_c32(); LOADF \
  { const float *xr = re, *xi = im; svzero_mask_za(0x33); T_A_G0 T_A_G1 T_A_G2 T_A_G3 TW(0) TW(1) TW(2) TW(3) STYP(Sbuf) } \
  for (long t = 0; t < B; t++) { \
    const float *sr = Sbuf + 512 * (t & 1), *si = sr + 256; \
    svzero_mask_za(0xcc); \
    if (t + 1 < B) { \
      const float *xr = re + N * (t + 1), *xi = im + N * (t + 1); \
      svzero_mask_za(0x33); \
      T_AB_G0 T_AB_G1 T_AB_G2 T_AB_G3 \
      TW(0) TW(1) TW(2) TW(3) \
      float *S2 = Sbuf + 512 * ((t + 1) & 1); \
      STYP(S2) \
    } else { T_B_G0 T_B_G1 T_B_G2 T_B_G3 } \
    STZ(ore + N * t, oim + N * t) \
  } }
#define STYP_X4(S) { STX4(0, (S)) STX4(1, (S) + 256) }
#define STZ_X4(yr, yi) { STX4(2, (yr)) STX4(3, (yi)) }
KER45(fft256_sme4, STYP_X4, STZ_X4)
KER45(fft256_sme5, STORE_YP3, STZ_X4)


/* Kernel v6: v5 plus software prefetch (PRFM) of the input of transform t+PFD (re and im, 32 lines of
 * 128 bytes... 16 per array) and a store-prefetch of output t+PFD, for the separate-arrays (DRAM) case. */
#ifndef PFD
#define PFD 4
#endif
__arm_locally_streaming __arm_new("za") void fft256_sme6(const float *re, const float *im, float *ore, float *oim, long B) {
  svbool_t pg = svptrue_b32(); svcount_t pc = svptrue_c32(); LOADF
  for (int q = 1; q < PFD && q < B; q++) for (int l = 0; l < 8; l++) { __builtin_prefetch(re + N * q + 32 * l, 0, 3); __builtin_prefetch(im + N * q + 32 * l, 0, 3); }
  { const float *xr = re, *xi = im; svzero_mask_za(0x33); T_A_G0 T_A_G1 T_A_G2 T_A_G3 TW(0) TW(1) TW(2) TW(3) STORE_YP3(Sbuf) }
  for (long t = 0; t < B; t++) {
    const float *sr = Sbuf + 512 * (t & 1), *si = sr + 256;
    if (t + PFD < B) for (int l = 0; l < 8; l++) {
      __builtin_prefetch(re + N * (t + PFD) + 32 * l, 0, 3); __builtin_prefetch(im + N * (t + PFD) + 32 * l, 0, 3);
      __builtin_prefetch(ore + N * (t + PFD) + 32 * l, 1, 3); __builtin_prefetch(oim + N * (t + PFD) + 32 * l, 1, 3); }
    svzero_mask_za(0xcc);
    if (t + 1 < B) {
      const float *xr = re + N * (t + 1), *xi = im + N * (t + 1);
      svzero_mask_za(0x33);
      T_AB_G0 T_AB_G1 T_AB_G2 T_AB_G3
      TW(0) TW(1) TW(2) TW(3)
      float *S2 = Sbuf + 512 * ((t + 1) & 1);
      STORE_YP3(S2)
    } else { T_B_G0 T_B_G1 T_B_G2 T_B_G3 }
    STZ_X4(ore + N * t, oim + N * t)
  }
}


/* Kernel v7/v8: v5 plus input-only prefetch into L2 (PRFM PLDL2KEEP) at distance 4 (v7) or 8 (v8)
 * transforms; no output prefetch. */
#define KER78(name, DIST) \
__arm_locally_streaming __arm_new("za") void name(const float *re, const float *im, float *ore, float *oim, long B) { \
  svbool_t pg = svptrue_b32(); svcount_t pc = svptrue_c32(); LOADF \
  for (int q = 1; q < DIST && q < B; q++) for (int l = 0; l < 8; l++) { __builtin_prefetch(re + N * q + 32 * l, 0, 2); __builtin_prefetch(im + N * q + 32 * l, 0, 2); } \
  { const float *xr = re, *xi = im; svzero_mask_za(0x33); T_A_G0 T_A_G1 T_A_G2 T_A_G3 TW(0) TW(1) TW(2) TW(3) STORE_YP3(Sbuf) } \
  for (long t = 0; t < B; t++) { \
    const float *sr = Sbuf + 512 * (t & 1), *si = sr + 256; \
    if (t + DIST < B) for (int l = 0; l < 8; l++) { __builtin_prefetch(re + N * (t + DIST) + 32 * l, 0, 2); __builtin_prefetch(im + N * (t + DIST) + 32 * l, 0, 2); } \
    svzero_mask_za(0xcc); \
    if (t + 1 < B) { \
      const float *xr = re + N * (t + 1), *xi = im + N * (t + 1); \
      svzero_mask_za(0x33); \
      T_AB_G0 T_AB_G1 T_AB_G2 T_AB_G3 \
      TW(0) TW(1) TW(2) TW(3) \
      float *S2 = Sbuf + 512 * ((t + 1) & 1); \
      STORE_YP3(S2) \
    } else { T_B_G0 T_B_G1 T_B_G2 T_B_G3 } \
    STZ_X4(ore + N * t, oim + N * t) \
  } }
KER78(fft256_sme7, 4)
KER78(fft256_sme8, 8)


/* Kernel v9: drain-overlapped schedule (gen9.py). Loop invariant: tiles 0,1 hold the raw pass-A result
 * Y^T(t); tiles 2,3 are zero. Per transform t:
 *   1. twiddle + STR of Y'(t) row groups j=0,1 to the scratch (slot 4j+m holds row n2=j+4m)
 *   2. pass B(t) on row groups 0,1 (tiles 2,3) interleaved with twiddle + STR of groups 2,3
 *   3. zero tiles 0,1
 *   4. pass B(t) on groups 2,3 interleaved FMOPA by FMOPA with pass A(t+1) rows 0..7 (tiles 0,1)
 *   5. pass A(t+1) rows 8..15 interleaved with the output stores of Z(t) (MOVA x4 + ST1W x4); zero 2,3
 * so the twiddle/staging drain overlaps pass B and the output drain overlaps pass A. Input prefetch to L2
 * 8 transforms ahead as in v8. */
#include "fft256_gen9.h"
static float S9[512] __attribute__((aligned(128)));
__arm_locally_streaming __arm_new("za") void fft256_sme9(const float *re, const float *im, float *ore, float *oim, long B) {
  svbool_t pg = svptrue_b32(); svcount_t pc = svptrue_c32(); LOADF
  float *S = S9;
  for (int q = 1; q < 8 && q < B; q++) for (int l = 0; l < 8; l++) { __builtin_prefetch(re + N * q + 32 * l, 0, 2); __builtin_prefetch(im + N * q + 32 * l, 0, 2); }
  svzero_mask_za(0xff);
  { const float *xr = re, *xi = im; T_A_G0 T_A_G1 T_A_G2 T_A_G3 }
  for (long t = 0; t < B; t++) {
    float *yo_r = ore + N * t, *yo_i = oim + N * t;
    if (t + 8 < B) for (int l = 0; l < 8; l++) { __builtin_prefetch(re + N * (t + 8) + 32 * l, 0, 2); __builtin_prefetch(im + N * (t + 8) + 32 * l, 0, 2); }
    V9_STEP1
    if (t + 1 < B) {
      const float *xr = re + N * (t + 1), *xi = im + N * (t + 1);
      V9_STEP2
      svzero_mask_za(0x33);
      V9_STEP4
      V9_STEP5
      svzero_mask_za(0xcc);
    } else {
      V9_STEP2_LAST
      V9_STEP4_LAST
      V9_STEP5_LAST
    }
  }
}

/* ---------------- reference ---------------- */
static void dft_ref(const float *re, const float *im, double *Xr, double *Xi) {
  for (int k = 0; k < N; k++) {
    double sr = 0, si = 0;
    for (int n = 0; n < N; n++) {
      double th = -2 * M_PI * (double)((long)k * n % N) / N, c = cos(th), s = sin(th);
      sr += re[n] * c - im[n] * s;
      si += re[n] * s + im[n] * c;
    }
    Xr[k] = sr; Xi[k] = si;
  }
}

/* ---------------- programs ---------------- */
enum { P_SME1, P_SME2, P_SME3, P_SME4, P_SME5, P_SME6, P_SME7, P_SME8, P_SME9, P_VDFT, P_VFFT, P_VFFTM, P_PFFFTZ, NP };
static const char *pname[NP] = {"SME2 four-step v1", "SME2 four-step v2 pipelined", "SME2 v3 pipelined+STR", "SME2 v4 (x4 stores)", "SME2 v5 (STR Y', x4 out)", "SME2 v6 (v5 + prefetch)", "SME2 v7 (v5 + L2 pf 4)", "SME2 v8 (v5 + L2 pf 8)", "SME2 v9 (drain overlapped)", "vDSP_DFT_Execute xB", "vDSP_fft_zop xB", "vDSP_fftm_zop (1 call)", "PFFFT z-order xB (interl.)"};
typedef struct {
  int kind;
  long B;
  bench_pool in, out; /* split: re[B*N] then im[B*N]; pffft: interleaved 2*B*N */
  vDSP_DFT_Setup ds;
  FFTSetup fs;
  PFFFT_Setup *ps;
  float *work;
} ctx;

static void prun(void *v, long i) {
  ctx *c = v;
  float *in = bench_pool_get(&c->in, i), *out = bench_pool_get(&c->out, i);
  long B = c->B, M = B * N;
  switch (c->kind) {
    case P_SME1: fft256_sme1(in, in + M, out, out + M, B); break;
    case P_SME2: fft256_sme2(in, in + M, out, out + M, B); break;
    case P_SME3: fft256_sme3(in, in + M, out, out + M, B); break;
    case P_SME4: fft256_sme4(in, in + M, out, out + M, B); break;
    case P_SME5: fft256_sme5(in, in + M, out, out + M, B); break;
    case P_SME6: fft256_sme6(in, in + M, out, out + M, B); break;
    case P_SME7: fft256_sme7(in, in + M, out, out + M, B); break;
    case P_SME8: fft256_sme8(in, in + M, out, out + M, B); break;
    case P_SME9: fft256_sme9(in, in + M, out, out + M, B); break;
    case P_VDFT:
      for (long b = 0; b < B; b++) vDSP_DFT_Execute(c->ds, in + N * b, in + M + N * b, out + N * b, out + M + N * b);
      break;
    case P_VFFT:
      for (long b = 0; b < B; b++) {
        DSPSplitComplex a = {in + N * b, in + M + N * b}, o = {out + N * b, out + M + N * b};
        vDSP_fft_zop(c->fs, &a, 1, &o, 1, 8, kFFTDirection_Forward);
      }
      break;
    case P_VFFTM: {
      DSPSplitComplex a = {in, in + M}, o = {out, out + M};
      vDSP_fftm_zop(c->fs, &a, 1, N, &o, 1, N, 8, B, kFFTDirection_Forward);
    } break;
    case P_PFFFTZ:
      for (long b = 0; b < B; b++) pffft_transform(c->ps, in + 2 * N * b, out + 2 * N * b, c->work, PFFFT_FORWARD);
      break;
  }
}

static void setup(ctx *c, int kind, long B, int mode) {
  memset(c, 0, sizeof *c);
  c->kind = kind; c->B = B;
  bench_pool_init(&c->in, sizeof(float) * 2 * N * B, mode);
  bench_pool_init(&c->out, sizeof(float) * 2 * N * B, mode);
  for (int p = 0; p < c->in.count; p++) bench_fill_f(c->in.buf[p], 2 * N * B);
  if (kind == P_VDFT) c->ds = vDSP_DFT_zop_CreateSetup(NULL, N, vDSP_DFT_FORWARD);
  if (kind == P_VFFT || kind == P_VFFTM) c->fs = vDSP_create_fftsetup(8, kFFTRadix2);
  if (kind == P_PFFFTZ) { c->ps = pffft_new_setup(N, PFFFT_COMPLEX); c->work = pffft_aligned_malloc(sizeof(float) * 2 * N); }
}

/* check one program on a batch: compare every transform with the fp64 DFT */
static int check(int kind, long B) {
  ctx c; setup(&c, kind, B, BENCH_SAME);
  float *in = c.in.buf[0], *out = c.out.buf[0];
  long M = B * N;
  if (kind == P_PFFFTZ) { /* interleave the same split data */
    float *t = bench_alloc(sizeof(float) * 2 * M);
    memcpy(t, in, sizeof(float) * 2 * M);
    for (long j = 0; j < M; j++) { in[2 * j] = t[j]; in[2 * j + 1] = t[M + j]; }
    free(t);
  }
  prun(&c, 0);
  double emax = 0, rmax = 0, Xr[N], Xi[N];
  float ir[N], ii[N];
  float *tmp = bench_alloc(sizeof(float) * 2 * N);
  for (long b = 0; b < B; b++) {
    for (int n = 0; n < N; n++) {
      if (kind == P_PFFFTZ) { ir[n] = in[2 * (N * b + n)]; ii[n] = in[2 * (N * b + n) + 1]; }
      else { ir[n] = in[N * b + n]; ii[n] = in[M + N * b + n]; }
    }
    dft_ref(ir, ii, Xr, Xi);
    const float *gr, *gi; int stride = 1;
    if (kind == P_PFFFTZ) { pffft_zreorder(c.ps, out + 2 * N * b, tmp, PFFFT_FORWARD); gr = tmp; gi = tmp + 1; stride = 2; }
    else { gr = out + N * b; gi = out + M + N * b; }
    for (int k = 0; k < N; k++) {
      double e = hypot(gr[stride * k] - Xr[k], gi[stride * k] - Xi[k]);
      if (!(e <= emax)) emax = e;
      double m = hypot(Xr[k], Xi[k]);
      if (m > rmax) rmax = m;
    }
  }
  int ok = emax / rmax < 1e-5;
  printf("CHECK %-28s B=%-4ld every transform vs fp64 DFT, natural order: max abs err %.3e  rel (to max |X|) %.3e  %s\n",
         pname[kind], B, emax, emax / rmax, ok ? "PASS" : "FAIL");
  free(tmp);
  return ok;
}

int main(int argc, char **argv) {
  bench_init();
  make_tables();
  if (argc < 2) { fprintf(stderr, "usage: fft256 check | time same|separate B...\n"); return 2; }
  if (!strcmp(argv[1], "check")) {
    int ok = 1;
    for (int k = 0; k < NP; k++) { ok &= check(k, 64); }
    ok &= check(P_SME1, 1); ok &= check(P_SME1, 3); ok &= check(P_SME2, 1); ok &= check(P_SME2, 2); ok &= check(P_SME2, 3); ok &= check(P_SME3, 1); ok &= check(P_SME3, 2); ok &= check(P_SME3, 5); ok &= check(P_SME4, 1); ok &= check(P_SME4, 5); ok &= check(P_SME5, 1); ok &= check(P_SME5, 5); ok &= check(P_SME6, 1); ok &= check(P_SME6, 7); ok &= check(P_SME7, 9); ok &= check(P_SME8, 9); ok &= check(P_SME9, 1); ok &= check(P_SME9, 2); ok &= check(P_SME9, 11);
    return ok ? 0 : 1;
  }
  int mode = !strcmp(argv[2], "separate") ? BENCH_SEPARATE : BENCH_SAME;
  for (int a = 3; a < argc; a++) {
    long B = atol(argv[a]);
    ctx c[NP]; bench_prog p[NP]; bench_stat st[NP];
    int np = 0; const char *sel = getenv("PROGS");
    for (int k = 0; k < NP; k++) {
      if (sel) { char key[8]; snprintf(key, sizeof key, ",%d,", k); char buf[128]; snprintf(buf, sizeof buf, ",%s,", sel); if (!strstr(buf, key)) continue; }
      setup(&c[np], k, B, mode); p[np] = (bench_prog){pname[k], prun, NULL, &c[np]}; np++; }
    char title[128];
    snprintf(title, sizeof title, "cfft256 B=%ld %s", B, argv[2]);
    bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)N * B;
    bench_compare(title, p, np, &o, st);
    for (int k = 0; k < np; k++) { bench_pool_free(&c[k].in); bench_pool_free(&c[k].out); }
  }
  return 0;
}
