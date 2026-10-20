// ap0p.c -- ap-0: numpy-layout front end for N = 2 * 8^P (1,024, 8,192) by decimation in time: SME does the two
// half-length transforms of every row straight from the interleaved input (smegp.h, LD4D), then one NEON pass does
// the radix-2 butterfly with twiddles and writes interleaved output. No pre-pass.
#include "smegp.h"
#include <arm_neon.h>
typedef struct { smegp_plan *pl; long N, M; double *wr, *wi, *yr, *yi; long CB; } ap0p_t;
ap0p_t *ap0p_new(long N) {
  ap0p_t *h = calloc(1, sizeof *h); h->N = N; h->M = N / 2; int P = 0; long m = h->M; while (m > 1) { m /= 8; P++; }
  h->pl = smegp_new(P, 0); h->wr = smegp_alloc(8 * h->M); h->wi = smegp_alloc(8 * h->M);
  for (long k = 0; k < h->M; k++) smegp_w(k, N, &h->wr[k], &h->wi[k]);
  long chunk = getenv("AP0_CHUNK") ? atol(getenv("AP0_CHUNK")) : 65536; h->CB = chunk / N; if (h->CB < 2) h->CB = 2;
  h->yr = smegp_alloc(8 * h->CB * N); h->yi = smegp_alloc(8 * h->CB * N);
  return h;
}
void ap0p_fft(ap0p_t *h, const double *x, double *y, long B) {
  long N = h->N, M = h->M;
  for (long t0 = 0; t0 < B; t0 += h->CB) {
    long cb = B - t0 < h->CB ? B - t0 : h->CB;
    smegp_fft2(h->pl, x + 2 * t0 * N, h->yr, h->yi, cb);
    for (long t = 0; t < cb; t++) {
      const double *a_r = h->yr + t * N, *a_i = h->yi + t * N, *b_r = a_r + M, *b_i = a_i + M; double *yt = y + 2 * (t0 + t) * N;
      for (long k = 0; k < M; k += 2) {
        float64x2_t y0r = vld1q_f64(a_r + k), y0i = vld1q_f64(a_i + k), y1r = vld1q_f64(b_r + k), y1i = vld1q_f64(b_i + k);
        float64x2_t wr = vld1q_f64(h->wr + k), wi = vld1q_f64(h->wi + k);
        float64x2_t tr = vfmsq_f64(vmulq_f64(y1r, wr), y1i, wi), ti = vfmaq_f64(vmulq_f64(y1r, wi), y1i, wr);
        float64x2_t lr = vaddq_f64(y0r, tr), li = vaddq_f64(y0i, ti), hr = vsubq_f64(y0r, tr), hi = vsubq_f64(y0i, ti);
        vst1q_f64(yt + 2 * k, vzip1q_f64(lr, li)); vst1q_f64(yt + 2 * k + 2, vzip2q_f64(lr, li));
        vst1q_f64(yt + 2 * (k + M), vzip1q_f64(hr, hi)); vst1q_f64(yt + 2 * (k + M) + 2, vzip2q_f64(hr, hi));
      }
    }
  }
}
