// as0.c -- wave 14, agent as-0: one-pass float64 Welch PSD (scipy.signal.welch semantics: given window, constant
// detrend, density scaling, one-sided, average='mean', nfft = nperseg even, any noverlap) with three FFT engines:
//   SME  : ap-0's fp64 radix-8 SME kernel (smegi.h / ap0lib2.c, copied from wave-13/ap-0); two real frames packed per
//          complex row; chunk of rows = NEON prep (mean, detrend, window, pack) -> SME FFT -> |Z|^2 per bin
//          mode 0 ("SME-A"): SME writes Z (L2 scratch), NEON pass accumulates |Z|^2
//          mode 1 ("SME-F", N = 512 / 4096 only): last SME pass accumulates |Z|^2 straight from ZA (no Z written)
//   VDSP : per frame pair: same NEON prep into split re/im (L1) -> vDSP_DFT_ExecuteD (split, fp64) -> NEON |Z|^2 acc
//   PFD  : per frame pair: NEON prep (interleaved) -> pffftd_transform (z-order, NEON double) -> |Z|^2 acc in z-order,
//          one pffftd_zreorder of the accumulator at the end
// Fold at the end: |X1_k|^2 + |X2_k|^2 = (|Z_k|^2 + |Z_{N-k}|^2) / 2. An odd segment count pads the last row with a
// zero imaginary frame (then the fold gives |X_k|^2 exactly).
#include "ap0lib2.c"
#include <Accelerate/Accelerate.h>
#include <arm_neon.h>
#include <stdio.h>
#include "pffft_double.h"
// KFR (kfrlib/kfr master 7c1e55c, GPL-2 / commercial, private measurement only): ar-0's wrapper kfrd.cc / kfrd.o
void *kd_new(long n); void kd_exec(void *v, const double *in, double *out);
// FFTW 3.3.10 double built WITH NEON by ar-f (wave 13): $FFTW_DOUBLE_NEON (GPL, private measurement only)
#include "fftw3.h"

typedef struct {
  long N, step;
  double *win, scale;  // scale = 1 / (fs * sum w^2)
  void *sme; long CB;  // pairs per SME chunk
  double *zin, *zout, *acc, *acc2;
  vDSP_DFT_SetupD vs; double *vr, *vi, *vor, *voi;
  PFFFTD_Setup *pf; double *pin, *pout, *pwork; int simd;
  int broken;  // 1 = skip the k <-> N-k fold (must FAIL)
  // SME-D (N = R * 8^P, R in {2, 4}): radix-R DIF step fused into the NEON prep, R sub-transforms of M = N / R on the
  // SME kernel with the |Z|^2 accumulation in its last pass (smegi.h inter = 3), accumulator in [r][k] order
  int R; long M; smeg_plan *plM; double *dxr, *dxi, *accp, *twr, *twi, *zeros;
  void *kd;    // KFR complex<double> plan
  double *mu;  // SME-W: frame means of a chunk
  fftw_plan fw; double *fin, *fout;
  int skip;    // diagnostics only: 1 = skip prep, 2 = skip SME
} as0_t;

static double *as0_alloc(size_t n) { void *p = NULL; posix_memalign(&p, 128, n * sizeof(double) + 128); memset(p, 0, n * sizeof(double)); return p; }

as0_t *as0_new(long N, long noverlap, const double *win, double fs) {
  as0_t *h = calloc(1, sizeof *h); h->N = N; h->step = N - noverlap;
  h->win = as0_alloc(N); memcpy(h->win, win, N * sizeof(double));
  double s2 = 0; for (long n = 0; n < N; n++) s2 += win[n] * win[n]; h->scale = 1.0 / (fs * s2);
  long m = N; int P = 0; while (m % 8 == 0 && m > 1) { m /= 8; P++; }
  int ok = (m == 1 && P >= 3) || ((m == 2 || m == 4) && P >= 3);
  if (ok) { h->sme = ap0v2_new(N); }
  long chunk = getenv("AS0_CHUNK") ? atol(getenv("AS0_CHUNK")) : 131072; h->CB = chunk / N; if (h->CB < 2) h->CB = 2;
  h->zin = as0_alloc(2 * h->CB * N); h->zout = as0_alloc(2 * h->CB * N); h->acc = as0_alloc(N); h->acc2 = as0_alloc(2 * N);
  h->vs = vDSP_DFT_zop_CreateSetupD(NULL, N, vDSP_DFT_FORWARD);
  h->vr = as0_alloc(N); h->vi = as0_alloc(N); h->vor = as0_alloc(N); h->voi = as0_alloc(N);
  h->pf = pffftd_new_setup((int)N, PFFFT_COMPLEX); h->pin = pffftd_aligned_malloc(16 * N); h->pout = pffftd_aligned_malloc(16 * N);
  h->pwork = pffftd_aligned_malloc(16 * N); h->simd = pffftd_simd_size();
  h->broken = getenv("AS0_BROKEN") ? atoi(getenv("AS0_BROKEN")) : 0;
  h->skip = getenv("AS0_SKIP") ? atoi(getenv("AS0_SKIP")) : 0;
  h->zeros = as0_alloc(N); h->kd = kd_new(N); h->mu = as0_alloc(2 * h->CB + 2);
  h->fin = fftw_malloc(16 * N); h->fout = fftw_malloc(16 * N);
  h->fw = fftw_plan_dft_1d((int)N, (fftw_complex *)h->fin, (fftw_complex *)h->fout, FFTW_FORWARD, getenv("AS0_FFTW_PATIENT") ? FFTW_PATIENT : FFTW_MEASURE);
  if ((m == 2 || m == 4) && P >= 3) {
    h->R = (int)m; h->M = N / m; h->plM = smeg_new(P, 0);
    h->dxr = as0_alloc(h->CB * N); h->dxi = as0_alloc(h->CB * N); h->accp = as0_alloc(N);
    h->twr = as0_alloc(N); h->twi = as0_alloc(N);
    for (int r = 1; r < h->R; r++) for (long n = 0; n < h->M; n++) smeg_w((long)r * n, N, &h->twr[(r - 1) * h->M + n], &h->twi[(r - 1) * h->M + n]);
  }
  return h;
}

static inline double fmean(const double *f, long N) {
  float64x2_t s0 = vdupq_n_f64(0), s1 = s0, s2 = s0, s3 = s0;
  for (long n = 0; n < N; n += 8) {
    s0 = vaddq_f64(s0, vld1q_f64(f + n)); s1 = vaddq_f64(s1, vld1q_f64(f + n + 2));
    s2 = vaddq_f64(s2, vld1q_f64(f + n + 4)); s3 = vaddq_f64(s3, vld1q_f64(f + n + 6));
  }
  return vaddvq_f64(vaddq_f64(vaddq_f64(s0, s1), vaddq_f64(s2, s3))) / (double)N;
}
// z[n] = ((a[n]-ma) w[n], (b[n]-mb) w[n]) interleaved; b == NULL -> imaginary part 0
static inline void prep_int(const double *a, const double *b, const double *w, double *z, long N) {
  float64x2_t ma = vdupq_n_f64(fmean(a, N)), mb = vdupq_n_f64(b ? fmean(b, N) : 0.0);
  if (b) for (long n = 0; n < N; n += 2) {
    float64x2_t ww = vld1q_f64(w + n); float64x2x2_t v;
    v.val[0] = vmulq_f64(vsubq_f64(vld1q_f64(a + n), ma), ww); v.val[1] = vmulq_f64(vsubq_f64(vld1q_f64(b + n), mb), ww);
    vst2q_f64(z + 2 * n, v);
  } else for (long n = 0; n < N; n += 2) {
    float64x2_t ww = vld1q_f64(w + n); float64x2x2_t v;
    v.val[0] = vmulq_f64(vsubq_f64(vld1q_f64(a + n), ma), ww); v.val[1] = vdupq_n_f64(0); vst2q_f64(z + 2 * n, v);
  }
}
static inline void prep_split(const double *a, const double *b, const double *w, double *zr, double *zi, long N) {
  float64x2_t ma = vdupq_n_f64(fmean(a, N)), mb = vdupq_n_f64(b ? fmean(b, N) : 0.0);
  for (long n = 0; n < N; n += 2) {
    float64x2_t ww = vld1q_f64(w + n);
    vst1q_f64(zr + n, vmulq_f64(vsubq_f64(vld1q_f64(a + n), ma), ww));
    vst1q_f64(zi + n, b ? vmulq_f64(vsubq_f64(vld1q_f64(b + n), mb), ww) : vdupq_n_f64(0));
  }
}
static inline void acc_int(double *acc, const double *z, long N) {  // acc[k] += |z_k|^2, z interleaved
  for (long k = 0; k < N; k += 4) {
    float64x2x2_t u = vld2q_f64(z + 2 * k), v = vld2q_f64(z + 2 * k + 4);
    vst1q_f64(acc + k, vfmaq_f64(vfmaq_f64(vld1q_f64(acc + k), u.val[0], u.val[0]), u.val[1], u.val[1]));
    vst1q_f64(acc + k + 2, vfmaq_f64(vfmaq_f64(vld1q_f64(acc + k + 2), v.val[0], v.val[0]), v.val[1], v.val[1]));
  }
}
static inline void acc_split(double *acc, const double *zr, const double *zi, long N) {
  for (long k = 0; k < N; k += 4) {
    float64x2_t r0 = vld1q_f64(zr + k), i0 = vld1q_f64(zi + k), r1 = vld1q_f64(zr + k + 2), i1 = vld1q_f64(zi + k + 2);
    vst1q_f64(acc + k, vfmaq_f64(vfmaq_f64(vld1q_f64(acc + k), r0, r0), i0, i0));
    vst1q_f64(acc + k + 2, vfmaq_f64(vfmaq_f64(vld1q_f64(acc + k + 2), r1, r1), i1, i1));
  }
}
static void finish(const as0_t *h, const double *acc, long nseg, double *psd) {
  long N = h->N; double sc = h->scale / (double)nseg;
  for (long k = 0; k <= N / 2; k++) {
    double s = h->broken ? acc[k] : 0.5 * (acc[k] + acc[(N - k) % N]);
    psd[k] = s * sc * ((k == 0 || k == N / 2) ? 1.0 : 2.0);
  }
}
long as0_nseg(const as0_t *h, long L) { return L < h->N ? 0 : (L - h->N) / h->step + 1; }

// mode 0: SME writes Z, NEON accumulates; mode 1: fused accumulation in the last SME pass (N = 8^P only)
int as0_welch_sme(as0_t *h, const double *x, long L, double *psd, int mode) {
  if (!h->sme) return -1;
  ap0v2_t *hv = h->sme; if (mode == 1 && !(hv->a && hv->a->R == 1)) return -2;
  long N = h->N, st = h->step, nseg = as0_nseg(h, L), nrow = (nseg + 1) / 2; double *acc = h->acc;
  memset(acc, 0, N * sizeof(double));
  for (long r0 = 0; r0 < nrow; r0 += h->CB) {
    long cb = nrow - r0 < h->CB ? nrow - r0 : h->CB;
    for (long r = 0; r < cb; r++) {
      long s = 2 * (r0 + r); const double *a = x + s * st, *b = (s + 1 < nseg) ? a + st : NULL;
      if (!(h->skip & 1)) prep_int(a, b, h->win, h->zin + 2 * r * N, N);
    }
    if (mode == 1) { if (!(h->skip & 2)) smeg_fft2(hv->a->pl, h->zin, NULL, acc, NULL, cb, 2); }
    else { ap0v2_fft(hv, h->zin, h->zout, cb); for (long r = 0; r < cb; r++) acc_int(acc, h->zout + 2 * r * N, N); }
  }
  finish(h, acc, nseg, psd); return 0;
}
int as0_welch_vdsp(as0_t *h, const double *x, long L, double *psd) {
  long N = h->N, st = h->step, nseg = as0_nseg(h, L); double *acc = h->acc;
  memset(acc, 0, N * sizeof(double));
  for (long s = 0; s < nseg; s += 2) {
    const double *a = x + s * st, *b = (s + 1 < nseg) ? a + st : NULL;
    prep_split(a, b, h->win, h->vr, h->vi, N);
    vDSP_DFT_ExecuteD(h->vs, h->vr, h->vi, h->vor, h->voi);
    acc_split(acc, h->vor, h->voi, N);
  }
  finish(h, acc, nseg, psd); return 0;
}
int as0_welch_pfd(as0_t *h, const double *x, long L, double *psd) {
  long N = h->N, st = h->step, nseg = as0_nseg(h, L); double *acc = h->acc; int S = h->simd;
  memset(acc, 0, N * sizeof(double));
  for (long s = 0; s < nseg; s += 2) {
    const double *a = x + s * st, *b = (s + 1 < nseg) ? a + st : NULL;
    prep_int(a, b, h->win, h->pin, N);
    pffftd_transform(h->pf, h->pin, h->pout, h->pwork, PFFFT_FORWARD);
    const double *z = h->pout;  // z-order: blocks of S re then S im
    if (S == 2) for (long k = 0; k < N; k += 2) {
      float64x2_t r = vld1q_f64(z + 2 * k), i = vld1q_f64(z + 2 * k + 2);
      vst1q_f64(acc + k, vfmaq_f64(vfmaq_f64(vld1q_f64(acc + k), r, r), i, i));
    } else if (S == 4) for (long k = 0; k < N; k += 4) {
      float64x2_t r0 = vld1q_f64(z + 2 * k), r1 = vld1q_f64(z + 2 * k + 2), i0 = vld1q_f64(z + 2 * k + 4), i1 = vld1q_f64(z + 2 * k + 6);
      vst1q_f64(acc + k, vfmaq_f64(vfmaq_f64(vld1q_f64(acc + k), r0, r0), i0, i0));
      vst1q_f64(acc + k + 2, vfmaq_f64(vfmaq_f64(vld1q_f64(acc + k + 2), r1, r1), i1, i1));
    } else for (long k = 0; k < N; k++) { long blk = k / S, l = k % S; double re = z[2 * blk * S + l], im = z[2 * blk * S + S + l]; acc[k] += re * re + im * im; }
  }
  // back to natural order: put acc in z layout (re = acc, im = 0) and reorder
  double *zz = h->acc2; for (long k = 0; k < N; k++) { long blk = k / S, l = k % S; zz[2 * blk * S + l] = acc[k]; zz[2 * blk * S + S + l] = 0; }
  pffftd_zreorder(h->pf, zz, h->pout, PFFFT_FORWARD);
  double *nat = h->zout; for (long k = 0; k < N; k++) nat[k] = h->pout[2 * k];
  finish(h, nat, nseg, psd); return 0;
}
// SME-D prep: frames a, b (b may be h->zeros) -> mean, detrend, window, pack, one radix-R DIF step with twiddles ->
// split rows dxr/dxi[(r * CB + t) * M + n]
#define CMULW(yr, yi, wr_, wi_, dr, di) do { float64x2_t wr__ = (wr_), wi__ = (wi_); \
  dr = vfmsq_f64(vmulq_f64(yr, wr__), yi, wi__); di = vfmaq_f64(vmulq_f64(yr, wi__), yi, wr__); } while (0)
static inline void prep_dif(const as0_t *h, const double *a, const double *b, long t) {
  long N = h->N, M = h->M, CB = h->CB; const double *w = h->win;
  float64x2_t ma = vdupq_n_f64(fmean(a, N)), mb = vdupq_n_f64(b == h->zeros ? 0.0 : fmean(b, N));
  double *xr[4], *xi[4]; for (int r = 0; r < h->R; r++) { xr[r] = h->dxr + (r * CB + t) * M; xi[r] = h->dxi + (r * CB + t) * M; }
#define ZQ(q, zr, zi) float64x2_t zr, zi; { float64x2_t ww = vld1q_f64(w + n + (q) * M); \
    zr = vmulq_f64(vsubq_f64(vld1q_f64(a + n + (q) * M), ma), ww); zi = vmulq_f64(vsubq_f64(vld1q_f64(b + n + (q) * M), mb), ww); }
  if (h->R == 2) for (long n = 0; n < M; n += 2) {
    ZQ(0, z0r, z0i) ZQ(1, z1r, z1i)
    vst1q_f64(xr[0] + n, vaddq_f64(z0r, z1r)); vst1q_f64(xi[0] + n, vaddq_f64(z0i, z1i));
    float64x2_t dr = vsubq_f64(z0r, z1r), di = vsubq_f64(z0i, z1i), o_r, o_i;
    CMULW(dr, di, vld1q_f64(h->twr + n), vld1q_f64(h->twi + n), o_r, o_i); vst1q_f64(xr[1] + n, o_r); vst1q_f64(xi[1] + n, o_i);
  } else for (long n = 0; n < M; n += 2) {
    ZQ(0, a0r, a0i) ZQ(1, a1r, a1i) ZQ(2, a2r, a2i) ZQ(3, a3r, a3i)
    float64x2_t s02r = vaddq_f64(a0r, a2r), s02i = vaddq_f64(a0i, a2i), d02r = vsubq_f64(a0r, a2r), d02i = vsubq_f64(a0i, a2i);
    float64x2_t s13r = vaddq_f64(a1r, a3r), s13i = vaddq_f64(a1i, a3i), d13r = vsubq_f64(a1r, a3r), d13i = vsubq_f64(a1i, a3i);
    float64x2_t y1r = vaddq_f64(d02r, d13i), y1i = vsubq_f64(d02i, d13r), y2r = vsubq_f64(s02r, s13r), y2i = vsubq_f64(s02i, s13i);
    float64x2_t y3r = vsubq_f64(d02r, d13i), y3i = vaddq_f64(d02i, d13r), o_r, o_i;
    vst1q_f64(xr[0] + n, vaddq_f64(s02r, s13r)); vst1q_f64(xi[0] + n, vaddq_f64(s02i, s13i));
    CMULW(y1r, y1i, vld1q_f64(h->twr + n), vld1q_f64(h->twi + n), o_r, o_i); vst1q_f64(xr[1] + n, o_r); vst1q_f64(xi[1] + n, o_i);
    CMULW(y2r, y2i, vld1q_f64(h->twr + M + n), vld1q_f64(h->twi + M + n), o_r, o_i); vst1q_f64(xr[2] + n, o_r); vst1q_f64(xi[2] + n, o_i);
    CMULW(y3r, y3i, vld1q_f64(h->twr + 2 * M + n), vld1q_f64(h->twi + 2 * M + n), o_r, o_i); vst1q_f64(xr[3] + n, o_r); vst1q_f64(xi[3] + n, o_i);
  }
#undef ZQ
}
int as0_welch_smed(as0_t *h, const double *x, long L, double *psd) {
  if (!h->plM) return -3;
  long N = h->N, M = h->M, st = h->step, nseg = as0_nseg(h, L), nrow = (nseg + 1) / 2; int R = h->R;
  memset(h->accp, 0, N * sizeof(double));
  for (long r0 = 0; r0 < nrow; r0 += h->CB) {
    long cb = nrow - r0 < h->CB ? nrow - r0 : h->CB;
    if (!(h->skip & 1)) for (long r = 0; r < cb; r++) {
      long s = 2 * (r0 + r); const double *a = x + s * st, *b = (s + 1 < nseg) ? a + st : h->zeros;
      prep_dif(h, a, b, r);
    }
    if (!(h->skip & 2)) for (int q = 0; q < R; q++) smeg_fft2(h->plM, h->dxr + q * h->CB * M, h->dxi + q * h->CB * M, h->accp + q * M, NULL, cb, 3);
  }
  double *acc = h->acc; for (int q = 0; q < R; q++) for (long k = 0; k < M; k++) acc[R * k + q] = h->accp[q * M + k];
  finish(h, acc, nseg, psd); return 0;
}
int as0_welch_kfr(as0_t *h, const double *x, long L, double *psd) {
  long N = h->N, st = h->step, nseg = as0_nseg(h, L); double *acc = h->acc;
  memset(acc, 0, N * sizeof(double));
  for (long s = 0; s < nseg; s += 2) {
    const double *a = x + s * st, *b = (s + 1 < nseg) ? a + st : NULL;
    prep_int(a, b, h->win, h->pin, N); kd_exec(h->kd, h->pin, h->pout); acc_int(acc, h->pout, N);
  }
  finish(h, acc, nseg, psd); return 0;
}
// SME-W (N = 8^P): mean in a NEON pass, detrend + window inside SME pass 1 (raw frames read by the kernel), |Z|^2 in the
// last pass; an odd last frame goes through the SME-F path with a zero imaginary frame.
int as0_welch_smew(as0_t *h, const double *x, long L, double *psd) {
  if (!h->sme) return -1;
  ap0v2_t *hv = h->sme; if (!(hv->a && hv->a->R == 1)) return -2;
  smeg_plan *pl = hv->a->pl; long N = h->N, st = h->step, nseg = as0_nseg(h, L), npair = nseg / 2; double *acc = h->acc;
  memset(acc, 0, N * sizeof(double));
  pl->rs = 2 * st; pl->st = st; pl->win = h->win; pl->mu = h->mu;
  for (long r0 = 0; r0 < npair; r0 += h->CB) {
    long cb = npair - r0 < h->CB ? npair - r0 : h->CB; const double *x0 = x + 2 * r0 * st;
    if (!(h->skip & 1)) for (long r = 0; r < 2 * cb; r++) h->mu[r] = fmean(x0 + r * st, N);
    if (!(h->skip & 2)) smeg_fft2(pl, x0, NULL, acc, NULL, cb, 4);
  }
  if (nseg & 1) { prep_int(x + (nseg - 1) * st, NULL, h->win, h->zin, N); smeg_fft2(pl, h->zin, NULL, acc, NULL, 1, 2); }
  finish(h, acc, nseg, psd); return 0;
}
int as0_welch_fftw(as0_t *h, const double *x, long L, double *psd) {
  long N = h->N, st = h->step, nseg = as0_nseg(h, L); double *acc = h->acc;
  memset(acc, 0, N * sizeof(double));
  for (long s = 0; s < nseg; s += 2) {
    const double *a = x + s * st, *b = (s + 1 < nseg) ? a + st : NULL;
    prep_int(a, b, h->win, h->fin, N); fftw_execute(h->fw); acc_int(acc, h->fout, N);
  }
  finish(h, acc, nseg, psd); return 0;
}
// multi-channel: C rows of length L (row stride L), psd rows of N/2+1; engine 0 = SME-A, 1 = SME-F, 2 = vDSP, 3 = PFFFT-d, 4 = SME-D, 5 = KFR-d, 6 = SME-W, 7 = FFTW-d NEON
int as0_welch_rows(as0_t *h, const double *x, long C, long L, double *psd, int engine) {
  long F = h->N / 2 + 1; int rc = 0;
  for (long c = 0; c < C; c++) {
    const double *xc = x + c * L; double *pc = psd + c * F;
    if (engine <= 1) rc |= as0_welch_sme(h, xc, L, pc, engine);
    else if (engine == 4) rc |= as0_welch_smed(h, xc, L, pc);
    else if (engine == 5) rc |= as0_welch_kfr(h, xc, L, pc);
    else if (engine == 6) rc |= as0_welch_smew(h, xc, L, pc);
    else if (engine == 7) rc |= as0_welch_fftw(h, xc, L, pc);
    else if (engine == 2) rc |= as0_welch_vdsp(h, xc, L, pc);
    else rc |= as0_welch_pfd(h, xc, L, pc);
  }
  return rc;
}
