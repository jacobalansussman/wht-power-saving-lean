// ap0lib.c -- ap-0: numpy-layout front end for the SME fp64 FFT. complex128 rows (interleaved re,im), forward,
// unnormalised, out of place, natural order, N = R * 8^P with R in {1, 2, 4} (512, 1024, 2048, 4096, 8192, 16384, 32768).
// Per chunk of rows: NEON pass (deinterleave + one radix-R DIF step with twiddles) -> split scratch; SME pass
// (smeg_fft on R*rows transforms of 8^P points, streaming mode entered once per chunk); NEON pass (R-way output
// interleave + re/im interleave) -> y.
#include "smegi.h"
#include <arm_neon.h>
#include <stdio.h>
typedef struct { smeg_plan *pl; long N, M; int R; double *wr, *wi; double *xr, *xi, *yr, *yi; long CB; } ap0_t;

ap0_t *ap0_new(long N) {
  ap0_t *h = calloc(1, sizeof *h); h->N = N;
  int P = 0; long m = N; while (m % 8 == 0 && m > 1) { m /= 8; P++; }
  // N = m * 8^P; need m in {1,2,4}; if m == 1 use R = 1; else R = m
  if (!(m == 1 || m == 2 || m == 4) || P < 3) { if (m == 1 && P >= 4) {} else { fprintf(stderr, "ap0: unsupported N=%ld\n", N); return NULL; } }
  h->R = (int)m; h->M = N / h->R; h->pl = smeg_new(P, 0);
  h->wr = smeg_alloc(8 * N); h->wi = smeg_alloc(8 * N);  // w_N^{n r} for r = 1..R-1, n < M: index (r-1)*M + n
  for (int r = 1; r < h->R; r++) for (long n = 0; n < h->M; n++) smeg_w((long)r * n, N, &h->wr[(r - 1) * h->M + n], &h->wi[(r - 1) * h->M + n]);
  long chunk = getenv("AP0_CHUNK") ? atol(getenv("AP0_CHUNK")) : 131072; h->CB = chunk / N; if (h->CB < 2) h->CB = 2;
  long tot = h->CB * N; h->xr = smeg_alloc(8 * tot); h->xi = smeg_alloc(8 * tot); h->yr = smeg_alloc(8 * tot); h->yi = smeg_alloc(8 * tot);
  return h;
}

static void pre(const ap0_t *h, const double *x, long cb) {
  long N = h->N, M = h->M;
  for (long t = 0; t < cb; t++) {
    const double *xt = x + 2 * t * N;
    if (h->R == 1) {
      double *dr = h->xr + t * N, *di = h->xi + t * N;
      for (long n = 0; n < N; n += 2) { float64x2x2_t v = vld2q_f64(xt + 2 * n); vst1q_f64(dr + n, v.val[0]); vst1q_f64(di + n, v.val[1]); }
    } else if (h->R == 2) {
      double *d0r = h->xr + (2 * t) * M, *d0i = h->xi + (2 * t) * M, *d1r = d0r + M, *d1i = d0i + M;
      for (long n = 0; n < M; n += 2) {
        float64x2x2_t a = vld2q_f64(xt + 2 * n), b = vld2q_f64(xt + 2 * (n + M));
        float64x2_t sr = vaddq_f64(a.val[0], b.val[0]), si = vaddq_f64(a.val[1], b.val[1]);
        float64x2_t dr = vsubq_f64(a.val[0], b.val[0]), di = vsubq_f64(a.val[1], b.val[1]);
        float64x2_t wr = vld1q_f64(h->wr + n), wi = vld1q_f64(h->wi + n);
        vst1q_f64(d0r + n, sr); vst1q_f64(d0i + n, si);
        vst1q_f64(d1r + n, vfmsq_f64(vmulq_f64(dr, wr), di, wi)); vst1q_f64(d1i + n, vfmaq_f64(vmulq_f64(dr, wi), di, wr));
      }
    } else {  // R == 4
      double *o[4][2]; for (int r = 0; r < 4; r++) { o[r][0] = h->xr + (4 * t + r) * M; o[r][1] = h->xi + (4 * t + r) * M; }
      for (long n = 0; n < M; n += 2) {
        float64x2x2_t a0 = vld2q_f64(xt + 2 * n), a1 = vld2q_f64(xt + 2 * (n + M)), a2 = vld2q_f64(xt + 2 * (n + 2 * M)), a3 = vld2q_f64(xt + 2 * (n + 3 * M));
        float64x2_t s02r = vaddq_f64(a0.val[0], a2.val[0]), s02i = vaddq_f64(a0.val[1], a2.val[1]);
        float64x2_t d02r = vsubq_f64(a0.val[0], a2.val[0]), d02i = vsubq_f64(a0.val[1], a2.val[1]);
        float64x2_t s13r = vaddq_f64(a1.val[0], a3.val[0]), s13i = vaddq_f64(a1.val[1], a3.val[1]);
        float64x2_t d13r = vsubq_f64(a1.val[0], a3.val[0]), d13i = vsubq_f64(a1.val[1], a3.val[1]);
        // r=0: s02 + s13; r=2: s02 - s13; r=1: d02 - i d13; r=3: d02 + i d13   (-i (a+bi) = b - ai)
        float64x2_t y0r = vaddq_f64(s02r, s13r), y0i = vaddq_f64(s02i, s13i);
        float64x2_t y2r = vsubq_f64(s02r, s13r), y2i = vsubq_f64(s02i, s13i);
        float64x2_t y1r = vaddq_f64(d02r, d13i), y1i = vsubq_f64(d02i, d13r);
        float64x2_t y3r = vsubq_f64(d02r, d13i), y3i = vaddq_f64(d02i, d13r);
        vst1q_f64(o[0][0] + n, y0r); vst1q_f64(o[0][1] + n, y0i);
        float64x2_t w1r = vld1q_f64(h->wr + n), w1i = vld1q_f64(h->wi + n), w2r = vld1q_f64(h->wr + M + n), w2i = vld1q_f64(h->wi + M + n);
        float64x2_t w3r = vld1q_f64(h->wr + 2 * M + n), w3i = vld1q_f64(h->wi + 2 * M + n);
        vst1q_f64(o[1][0] + n, vfmsq_f64(vmulq_f64(y1r, w1r), y1i, w1i)); vst1q_f64(o[1][1] + n, vfmaq_f64(vmulq_f64(y1r, w1i), y1i, w1r));
        vst1q_f64(o[2][0] + n, vfmsq_f64(vmulq_f64(y2r, w2r), y2i, w2i)); vst1q_f64(o[2][1] + n, vfmaq_f64(vmulq_f64(y2r, w2i), y2i, w2r));
        vst1q_f64(o[3][0] + n, vfmsq_f64(vmulq_f64(y3r, w3r), y3i, w3i)); vst1q_f64(o[3][1] + n, vfmaq_f64(vmulq_f64(y3r, w3i), y3i, w3r));
      }
    }
  }
}
static void post(const ap0_t *h, double *y, long cb) {
  long N = h->N, M = h->M; int R = h->R;
  for (long t = 0; t < cb; t++) {
    double *yt = y + 2 * t * N;
    if (R == 1) {
      const double *sr = h->yr + t * N, *si = h->yi + t * N;
      for (long n = 0; n < N; n += 2) { float64x2x2_t v; v.val[0] = vld1q_f64(sr + n); v.val[1] = vld1q_f64(si + n); vst2q_f64(yt + 2 * n, v); }
    } else if (R == 2) {  // y[2k + r] = sub_r[k]
      const double *s0r = h->yr + 2 * t * M, *s0i = h->yi + 2 * t * M, *s1r = s0r + M, *s1i = s0i + M;
      for (long k = 0; k < M; k += 2) {
        float64x2_t ar = vld1q_f64(s0r + k), ai = vld1q_f64(s0i + k), br = vld1q_f64(s1r + k), bi = vld1q_f64(s1i + k);
        // order: (a0r a0i b0r b0i a1r a1i b1r b1i)
        float64x2_t a0 = vzip1q_f64(ar, ai), a1 = vzip2q_f64(ar, ai), b0 = vzip1q_f64(br, bi), b1 = vzip2q_f64(br, bi);
        vst1q_f64(yt + 4 * k, a0); vst1q_f64(yt + 4 * k + 2, b0); vst1q_f64(yt + 4 * k + 4, a1); vst1q_f64(yt + 4 * k + 6, b1);
      }
    } else {
      const double *sr[4], *si[4]; for (int r = 0; r < 4; r++) { sr[r] = h->yr + (4 * t + r) * M; si[r] = h->yi + (4 * t + r) * M; }
      for (long k = 0; k < M; k += 2) {
        for (int r = 0; r < 4; r++) {
          float64x2_t ar = vld1q_f64(sr[r] + k), ai = vld1q_f64(si[r] + k);
          vst1q_f64(yt + 8 * k + 2 * r, vzip1q_f64(ar, ai)); vst1q_f64(yt + 8 * k + 8 + 2 * r, vzip2q_f64(ar, ai));
        }
      }
    }
  }
}
int ap0_skip = 0;  // diagnostics only: 1 = skip pre pass, 2 = skip SME pass, 4 = skip post pass
// forward FFT of B rows of N complex128 values: x, y interleaved, row stride N
void ap0_fft(ap0_t *h, const double *x, double *y, long B) {
  if (h->R == 1 && !(ap0_skip & 8)) { smeg_fft2(h->pl, x, NULL, y, NULL, B, 1); return; }  // 8 = force the conversion path
  for (long t0 = 0; t0 < B; t0 += h->CB) {
    long cb = B - t0 < h->CB ? B - t0 : h->CB;
    if (!(ap0_skip & 1)) pre(h, x + 2 * t0 * h->N, cb);
    if (!(ap0_skip & 2)) smeg_fft2(h->pl, h->xr, h->xi, h->yr, h->yi, cb * h->R, 0);
    if (!(ap0_skip & 4)) post(h, y + 2 * t0 * h->N, cb);
  }
}
void *ap0_new_v(long N) { return ap0_new(N); }
void ap0_fft_v(void *h, const double *x, double *y, long B) { ap0_fft((ap0_t *)h, x, y, B); }
