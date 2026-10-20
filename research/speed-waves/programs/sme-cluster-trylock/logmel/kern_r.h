/* kern_r.h -- am-2 (wave 12): REAL-input two-frame kernel for N = 400 = 25 x 16 on the SME2 unit.
 * For a real frame, pass A's output is Hermitian in k1 (Y[25-k1] = conj Y[k1]), so only k1 = 0..12 is
 * computed: ONE block of 16 k1 lanes per frame (lanes 13..15 are valid but unused). Two frames then fill
 * the four tiles: tiles 0/1 = frame a re/im, tiles 2/3 = frame b re/im.
 *   pass A^T: 4 FMOPA per n1 (100 per pair) instead of 8 (200) for two frames packed as one complex input.
 *   twiddle:  block-0 twiddles applied to both tile pairs, in place (aj-2's TWB).
 *   pass B^T: aj-2's pass B unchanged (its tiles 2/3 now hold frame b instead of k1 block 16..24).
 *   output:   every tile row (k2) stored whole with STR ZA: Z[t][k2][k1], 1024 floats per pair.
 * The power |X[k]|^2 = |Z[k1 + 25 k2]|^2 for k1 <= 12; for k1 >= 13 it is |Z[(25-k1) + 25 (15-k2)]|^2
 * (|X[N-k]| = |X[k]| for real input), so the NEON post-pass is a square-sum plus a fixed permutation.
 * WIN: multiply each loaded row by the window inside streaming mode (2 FMUL per n1). */
#ifndef PFD
#define PFD 8
#endif
#define RKN 400
#define RFUN(name, WIN)                                                                                         \
  __arm_locally_streaming __arm_new("za") static void name(const float *sig, long hop, const float *win, float *out, long B) { \
    svbool_t pg = svptrue_b32();                                                                                \
    svcount_t pc = svptrue_c32();                                                                               \
    const float *W4 = W4tab_400, *FA = FAtab_400, *FB = FBtab_400; float *S = S4_400;                         \
    const long span = RKN + hop;                                                                                \
    for (int q = 1; q < PFD && q < B; q++)                                                                      \
      for (long l = 0; l < span; l += 32) __builtin_prefetch(sig + 2 * hop * q + l, 0, 2);                      \
    for (long b = 0; b < B; b++) {                                                                              \
      const float *xa = sig + 2 * hop * b, *xb = xa + hop;                                                      \
      float *o = out + 1024 * b;                                                                                \
      if (PFD && b + PFD < B)                                                                                   \
        for (long l = 0; l < span; l += 32) __builtin_prefetch(sig + 2 * hop * (b + PFD) + l, 0, 2);            \
      svzero_mask_za(0xff);                                                                                     \
      for (int n1 = 0; n1 < 25; n1++) {                                                                         \
        svfloat32_t fr = svld1_f32(pg, FA + 64 * n1), fi = svld1_f32(pg, FA + 64 * n1 + 32);                   \
        svfloat32_t a = svld1_f32(pg, xa + 16 * n1), c = svld1_f32(pg, xb + 16 * n1);                          \
        if (WIN) { svfloat32_t w = svld1_f32(pg, win + 16 * n1); a = svmul_f32_x(pg, a, w); c = svmul_f32_x(pg, c, w); } \
        svmopa_za32_f32_m(0, pg, pg, a, fr); svmopa_za32_f32_m(1, pg, pg, a, fi);                              \
        svmopa_za32_f32_m(2, pg, pg, c, fr); svmopa_za32_f32_m(3, pg, pg, c, fi);                              \
      }                                                                                                         \
      TWB(0, 0, W4) TWB(1, 0, W4) TWB(2, 0, W4) TWB(3, 0, W4)                                                   \
      TWB(0, 2, W4) TWB(1, 2, W4) TWB(2, 2, W4) TWB(3, 2, W4)                                                   \
      for (int n2 = 0; n2 < 16; n2++) {                                                                         \
        svstr_za(4 * n2 + 0, S + 16 * n2); svstr_za(4 * n2 + 1, S + 256 + 16 * n2);                            \
        svstr_za(4 * n2 + 2, S + 512 + 16 * n2); svstr_za(4 * n2 + 3, S + 768 + 16 * n2);                      \
      }                                                                                                         \
      svzero_mask_za(0xff);                                                                                     \
      for (int n2 = 0; n2 < 16; n2++) {                                                                         \
        svfloat32x2_t f = svld1_f32_x2(pc, FB + 32 * n2);                                                       \
        svfloat32_t fr = svget2(f, 0), fi = svget2(f, 1);                                                       \
        svfloat32_t yra = svld1_f32(pg, S + 16 * n2), yia = svld1_f32(pg, S + 256 + 16 * n2);                  \
        svfloat32_t yrb = svld1_f32(pg, S + 512 + 16 * n2), yib = svld1_f32(pg, S + 768 + 16 * n2);            \
        svmopa_za32_f32_m(0, pg, pg, fr, yra); svmopa_za32_f32_m(1, pg, pg, fr, yia);                          \
        svmopa_za32_f32_m(2, pg, pg, fr, yrb); svmopa_za32_f32_m(3, pg, pg, fr, yib);                          \
        if (n2 != 0 && n2 != 8) {                                                                               \
          svmops_za32_f32_m(0, pg, pg, fi, yia); svmopa_za32_f32_m(1, pg, pg, fi, yra);                        \
          svmops_za32_f32_m(2, pg, pg, fi, yib); svmopa_za32_f32_m(3, pg, pg, fi, yrb);                        \
        }                                                                                                       \
      }                                                                                                         \
      for (int k2 = 0; k2 < 16; k2++) {                                                                         \
        svstr_za(4 * k2 + 0, o + 16 * k2); svstr_za(4 * k2 + 1, o + 256 + 16 * k2);                             \
        svstr_za(4 * k2 + 2, o + 512 + 16 * k2); svstr_za(4 * k2 + 3, o + 768 + 16 * k2);                       \
      }                                                                                                         \
    }                                                                                                           \
  }
RFUN(sme_real400_w, 1)
RFUN(sme_real400, 0)

/* power of one frame from its two 16x16 tiles (R, I: [k2][k1]) into p[0..200] */
static inline void power_r400(const float *restrict R, const float *restrict I, float *restrict p) {
  for (int k2 = 0; k2 < 8; k2++) {
    const float *r = R + 16 * k2, *i = I + 16 * k2;
    for (int k1 = 0; k1 < 13; k1++) p[25 * k2 + k1] = r[k1] * r[k1] + i[k1] * i[k1];
    const float *rr = R + 16 * (15 - k2), *ii = I + 16 * (15 - k2);
    for (int k1 = 13; k1 < 25; k1++) p[25 * k2 + k1] = rr[25 - k1] * rr[25 - k1] + ii[25 - k1] * ii[25 - k1];
  }
  p[200] = R[16 * 8] * R[16 * 8] + I[16 * 8] * I[16 * 8];
}

#include <arm_neon.h>
/* NEON version of power_r400 (same result) */
static inline float32x4_t sq4(float32x4_t r, float32x4_t i) { return vfmaq_f32(vmulq_f32(r, r), i, i); }
static inline float32x4_t rev4(float32x4_t v) { v = vrev64q_f32(v); return vextq_f32(v, v, 2); }
static inline void power_r400_neon(const float *restrict R, const float *restrict I, float *restrict p) {
  for (int k2 = 0; k2 < 8; k2++) {
    const float *r = R + 16 * k2, *i = I + 16 * k2;
    float *q = p + 25 * k2;
    vst1q_f32(q + 0, sq4(vld1q_f32(r + 0), vld1q_f32(i + 0)));
    vst1q_f32(q + 4, sq4(vld1q_f32(r + 4), vld1q_f32(i + 4)));
    vst1q_f32(q + 8, sq4(vld1q_f32(r + 8), vld1q_f32(i + 8)));
    q[12] = r[12] * r[12] + i[12] * i[12];
    const float *rr = R + 16 * (15 - k2), *ii = I + 16 * (15 - k2);
    vst1q_f32(q + 13, rev4(sq4(vld1q_f32(rr + 9), vld1q_f32(ii + 9))));
    vst1q_f32(q + 17, rev4(sq4(vld1q_f32(rr + 5), vld1q_f32(ii + 5))));
    vst1q_f32(q + 21, rev4(sq4(vld1q_f32(rr + 1), vld1q_f32(ii + 1))));
  }
  p[200] = R[16 * 8] * R[16 * 8] + I[16 * 8] * I[16 * 8];
}

/* |Z|^2 of the frame in tiles (TR, TI) into tile TI, then the 201 power bins written straight to p:
 * rows k2 = 0..7 give bins 25k2 + 0..15 directly; bins 25k2 + 16..24 come from row 15-k2, lanes 9..1,
 * reversed (|X[N-k]| = |X[k]|); bin 200 = row 8 lane 0. */
#define RP_POWER(TR, TI, MASK, p)                                                                       \
  {                                                                                                     \
    svfloat32x4_t i0 = svread_za32_f32_vg1x4(0 + TI), i1 = svread_za32_f32_vg1x4(4 + TI);               \
    svfloat32x4_t i2 = svread_za32_f32_vg1x4(8 + TI), i3 = svread_za32_f32_vg1x4(12 + TI);              \
    svzero_mask_za(MASK);                                                                               \
    svmla_za32_f32_vg1x4(0 + TI, i0, i0); svmla_za32_f32_vg1x4(4 + TI, i1, i1);                         \
    svmla_za32_f32_vg1x4(8 + TI, i2, i2); svmla_za32_f32_vg1x4(12 + TI, i3, i3);                        \
    svfloat32x4_t r0 = svread_za32_f32_vg1x4(0 + TR), r1 = svread_za32_f32_vg1x4(4 + TR);               \
    svfloat32x4_t r2 = svread_za32_f32_vg1x4(8 + TR), r3 = svread_za32_f32_vg1x4(12 + TR);              \
    svmla_za32_f32_vg1x4(0 + TI, r0, r0); svmla_za32_f32_vg1x4(4 + TI, r1, r1);                         \
    svmla_za32_f32_vg1x4(8 + TI, r2, r2); svmla_za32_f32_vg1x4(12 + TI, r3, r3);                        \
    for (int k2 = 0; k2 < 8; k2++) svst1_hor_za32(TI, k2, pg, (p) + 25 * k2);                         \
    svst1_hor_za32(TI, 8, p1, (p) + 200);                                                               \
    for (int k2 = 0; k2 < 8; k2++) {                                                                    \
      svfloat32_t v = svread_hor_za32_f32_m(svundef_f32(), pg, TI, 15 - k2);                            \
      svst1_f32(prv, (p) + 25 * k2 + 10, svrev_f32(v));                                                 \
    }                                                                                                   \
  }
#define RPFUN(name, WIN)                                                                                         \
  __arm_locally_streaming __arm_new("za") static void name(const float *sig, long hop, const float *win, float *P, long B) { \
    svbool_t pg = svptrue_b32(), p1 = svwhilelt_b32(0, 1);                                                      \
    svbool_t prv = svbic_b_z(pg, svwhilelt_b32(0, 15), svwhilelt_b32(0, 6));                                                                                \
    svcount_t pc = svptrue_c32();                                                                               \
    const float *W4 = W4tab_400, *FA = FAtab_400, *FB = FBtab_400; float *S = S4_400;                         \
    const long span = RKN + hop;                                                                                \
    for (int q = 1; q < PFD && q < B; q++)                                                                      \
      for (long l = 0; l < span; l += 32) __builtin_prefetch(sig + 2 * hop * q + l, 0, 2);                      \
    for (long b = 0; b < B; b++) {                                                                              \
      const float *xa = sig + 2 * hop * b, *xb = xa + hop;                                                      \
      float *pa = P + 402 * b, *pb = pa + 201;                                                                                \
      if (PFD && b + PFD < B)                                                                                   \
        for (long l = 0; l < span; l += 32) __builtin_prefetch(sig + 2 * hop * (b + PFD) + l, 0, 2);            \
      svzero_mask_za(0xff);                                                                                     \
      for (int n1 = 0; n1 < 25; n1++) {                                                                         \
        svfloat32_t fr = svld1_f32(pg, FA + 64 * n1), fi = svld1_f32(pg, FA + 64 * n1 + 32);                   \
        svfloat32_t a = svld1_f32(pg, xa + 16 * n1), c = svld1_f32(pg, xb + 16 * n1);                          \
        if (WIN) { svfloat32_t w = svld1_f32(pg, win + 16 * n1); a = svmul_f32_x(pg, a, w); c = svmul_f32_x(pg, c, w); } \
        svmopa_za32_f32_m(0, pg, pg, a, fr); svmopa_za32_f32_m(1, pg, pg, a, fi);                              \
        svmopa_za32_f32_m(2, pg, pg, c, fr); svmopa_za32_f32_m(3, pg, pg, c, fi);                              \
      }                                                                                                         \
      TWB(0, 0, W4) TWB(1, 0, W4) TWB(2, 0, W4) TWB(3, 0, W4)                                                   \
      TWB(0, 2, W4) TWB(1, 2, W4) TWB(2, 2, W4) TWB(3, 2, W4)                                                   \
      for (int n2 = 0; n2 < 16; n2++) {                                                                         \
        svstr_za(4 * n2 + 0, S + 16 * n2); svstr_za(4 * n2 + 1, S + 256 + 16 * n2);                            \
        svstr_za(4 * n2 + 2, S + 512 + 16 * n2); svstr_za(4 * n2 + 3, S + 768 + 16 * n2);                      \
      }                                                                                                         \
      svzero_mask_za(0xff);                                                                                     \
      for (int n2 = 0; n2 < 16; n2++) {                                                                         \
        svfloat32x2_t f = svld1_f32_x2(pc, FB + 32 * n2);                                                       \
        svfloat32_t fr = svget2(f, 0), fi = svget2(f, 1);                                                       \
        svfloat32_t yra = svld1_f32(pg, S + 16 * n2), yia = svld1_f32(pg, S + 256 + 16 * n2);                  \
        svfloat32_t yrb = svld1_f32(pg, S + 512 + 16 * n2), yib = svld1_f32(pg, S + 768 + 16 * n2);            \
        svmopa_za32_f32_m(0, pg, pg, fr, yra); svmopa_za32_f32_m(1, pg, pg, fr, yia);                          \
        svmopa_za32_f32_m(2, pg, pg, fr, yrb); svmopa_za32_f32_m(3, pg, pg, fr, yib);                          \
        if (n2 != 0 && n2 != 8) {                                                                               \
          svmops_za32_f32_m(0, pg, pg, fi, yia); svmopa_za32_f32_m(1, pg, pg, fi, yra);                        \
          svmops_za32_f32_m(2, pg, pg, fi, yib); svmopa_za32_f32_m(3, pg, pg, fi, yrb);                        \
        }                                                                                                       \
      }                                                                                                         \
      RP_POWER(0, 1, 0x22, pa) RP_POWER(2, 3, 0x88, pb)                                                       \
    }                                                                                                           \
  }
RPFUN(sme_real400_pw, 1)
