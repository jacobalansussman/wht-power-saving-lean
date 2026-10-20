/* wht_neon.h -- a careful NEON textbook Walsh-Hadamard transform (fp32, in place, unnormalised,
 * natural (Hadamard) order), header-only. speed wave 11, agent aj-base-cpu. Correct for 2^1..2^24
 * (checked by wht_check.c against a plain double-precision reference).
 *
 * Plan:
 *   1. "base64" kernel on every block of 64 floats (16 q-registers): the 4 cross-vector stages
 *      (strides 4,8,16,32 floats) as register butterflies, then the 2 in-vector stages (strides 1,2)
 *      on PAIRS of vectors with trn1/trn2 (64-bit) zip and uzp (3 shuffles per vector).
 *   2. remaining stages in radix-16 / 8 / 4 / 2 register passes: 16 vectors at stride S floats,
 *      4 butterfly stages in registers, store.
 *   3. cache blocking: for n > WHT_NEON_BLOCK_LOG (default 14 = 64 KiB of floats, inside the M4
 *      P-core 128 KiB L1D) each block of 2^B floats is finished first (all its low stages), then the
 *      high stages run as radix-16 passes over the whole array.
 * Usage: wht_neon_f32(float *buf, int log_n);  log_n in 0..30 (tested to 24).
 */
#ifndef WHT_NEON_H
#define WHT_NEON_H
#include <arm_neon.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifndef WHT_NEON_BLOCK_LOG
#define WHT_NEON_BLOCK_LOG 14
#endif
static int wht_neon_block_log = WHT_NEON_BLOCK_LOG;

#define WB(a, b)                      \
  do {                                \
    float32x4_t _t = (a);             \
    (a) = vaddq_f32(_t, (b));         \
    (b) = vsubq_f32(_t, (b));         \
  } while (0)

/* in-vector stages (strides 1 and 2) on two vectors x, y at once */
static inline void wht_inlane_pair(float32x4_t *x, float32x4_t *y) {
  float32x4_t t1 = vtrn1q_f32(*x, *y); /* x0 y0 x2 y2 */
  float32x4_t t2 = vtrn2q_f32(*x, *y); /* x1 y1 x3 y3 */
  float32x4_t s = vaddq_f32(t1, t2);   /* x01+ y01+ x23+ y23+ */
  float32x4_t d = vsubq_f32(t1, t2);   /* x01- y01- x23- y23- */
  float32x4_t u1 = vreinterpretq_f32_f64(vzip1q_f64(vreinterpretq_f64_f32(s), vreinterpretq_f64_f32(d)));
  float32x4_t u2 = vreinterpretq_f32_f64(vzip2q_f64(vreinterpretq_f64_f32(s), vreinterpretq_f64_f32(d)));
  float32x4_t p = vaddq_f32(u1, u2); /* X0 Y0 X1 Y1 */
  float32x4_t m = vsubq_f32(u1, u2); /* X2 Y2 X3 Y3 */
  *x = vuzp1q_f32(p, m);
  *y = vuzp2q_f32(p, m);
}

static inline void wht_base64(float *p) {
  float32x4_t v[16];
  for (int j = 0; j < 16; j++) v[j] = vld1q_f32(p + 4 * j);
  for (int h = 1; h < 16; h <<= 1)
    for (int i = 0; i < 16; i += 2 * h)
      for (int j = i; j < i + h; j++) WB(v[j], v[j + h]);
  for (int j = 0; j < 16; j += 2) wht_inlane_pair(&v[j], &v[j + 1]);
  for (int j = 0; j < 16; j++) vst1q_f32(p + 4 * j, v[j]);
}

/* radix-R cross pass (R = 2,4,8,16): stages with strides S, 2S, ... (R-1)S... over [p, p+R*S) */
#define WHT_DEF_RADIX(R)                                                          \
  static inline void wht_r##R(float *p, size_t S) {                               \
    for (size_t c = 0; c < S; c += 4) {                                           \
      float32x4_t v[R];                                                           \
      for (int j = 0; j < R; j++) v[j] = vld1q_f32(p + c + j * S);                \
      for (int h = 1; h < R; h <<= 1)                                             \
        for (int i = 0; i < R; i += 2 * h)                                        \
          for (int j = i; j < i + h; j++) WB(v[j], v[j + h]);                     \
      for (int j = 0; j < R; j++) vst1q_f32(p + c + j * S, v[j]);                 \
    }                                                                             \
  }
WHT_DEF_RADIX(2)
WHT_DEF_RADIX(4)
WHT_DEF_RADIX(8)
WHT_DEF_RADIX(16)

/* stages lo..hi-1 (stride 2^lo .. 2^(hi-1), lo >= 2) over the array buf of 2^log_n floats */
static inline void wht_cross_stages(float *buf, int log_n, int lo, int hi, int maxr) {
  size_t n = (size_t)1 << log_n;
  int s = lo;
  while (s < hi) {
    int left = hi - s, passes = (left + maxr - 1) / maxr;
    int r = (left + passes - 1) / passes; /* spread evenly, max radix 2^maxr */
    size_t S = (size_t)1 << s, span = S << r;
    for (size_t b = 0; b < n; b += span) {
      float *q = buf + b;
      switch (r) {
        case 1: wht_r2(q, S); break;
        case 2: wht_r4(q, S); break;
        case 3: wht_r8(q, S); break;
        default: wht_r16(q, S); break;
      }
    }
    s += r;
  }
}

static void wht_small_scalar(float *x, int log_n) {
  size_t n = (size_t)1 << log_n;
  for (size_t h = 1; h < n; h <<= 1)
    for (size_t i = 0; i < n; i += 2 * h)
      for (size_t j = i; j < i + h; j++) {
        float a = x[j], b = x[j + h];
        x[j] = a + b;
        x[j + h] = a - b;
      }
}

/* all stages of a block of 2^log_n floats, log_n >= 6 */
static inline void wht_block(float *buf, int log_n) {
  size_t n = (size_t)1 << log_n;
  for (size_t b = 0; b < n; b += 64) wht_base64(buf + b);
  wht_cross_stages(buf, log_n, 6, log_n, 4);
}

/* high stages (strides >= 2^B) by column strips: gather a strip of W floats from each of the 2^R
   rows into a small contiguous buffer (<= 64 KiB, L1-resident), do all R stages there, scatter back.
   This avoids the L1 set conflicts of 16 rows at a power-of-two stride (the cause of a 3x slowdown
   at 2^22 in the direct version). */
static int wht_neon_strips = 0;
/* max radix (log2) of the direct high-stage passes: 3 = radix 8, so the 8 rows at a power-of-two
   stride fit the 8 ways of an L1 set; radix 16 thrashes (measured 2-3x slower at 2^18, 2^22) */
static int wht_neon_high_maxr = 3;
static void wht_high_strips(float *buf, int log_n, int B) {
  int R = log_n - B;
  size_t rows = (size_t)1 << R, cols = (size_t)1 << B;
  int LW = 14 - R;
  if (LW > 6) LW = 6;
  if (LW < 2) LW = 2;
  size_t W = (size_t)1 << LW;
  static float *tmp = NULL;
  static size_t tmpn = 0;
  if (tmpn < rows * W) {
    free(tmp);
    if (posix_memalign((void **)&tmp, 128, rows * W * sizeof(float))) abort();
    tmpn = rows * W;
  }
  for (size_t c = 0; c < cols; c += W) {
    for (size_t r = 0; r < rows; r++) memcpy(tmp + r * W, buf + r * cols + c, W * sizeof(float));
    wht_cross_stages(tmp, LW + R, LW, LW + R, 4);
    for (size_t r = 0; r < rows; r++) memcpy(buf + r * cols + c, tmp + r * W, W * sizeof(float));
  }
}

static void wht_neon_f32(float *buf, int log_n) {
  if (log_n <= 0) return;
  if (log_n < 6) {
    wht_small_scalar(buf, log_n);
    return;
  }
  int B = wht_neon_block_log;
  if (log_n <= B) {
    wht_block(buf, log_n);
    return;
  }
  size_t n = (size_t)1 << log_n, bs = (size_t)1 << B;
  for (size_t b = 0; b < n; b += bs) wht_block(buf + b, B);
  if (wht_neon_strips) wht_high_strips(buf, log_n, B);
  else wht_cross_stages(buf, log_n, B, log_n, wht_neon_high_maxr);
}
#endif
