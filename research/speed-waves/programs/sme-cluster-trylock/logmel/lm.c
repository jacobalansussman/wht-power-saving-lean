/* lm.c -- ap-2 (wave 13): the WHOLE Whisper log-mel front end (n_fft 400, hop 160, periodic Hann, 80 Slaney mels,
 * log10(max(mel,1e-10)), max(x, max-8), (x+4)/4) on the M4, one thread and whole chip.
 * Built on am-2's real-input SME2 STFT kernel (kern_r.h, copied unchanged from wave-12/am-2).
 *
 * Programs (all produce out[F][80] = the Whisper-normalised log-mel; checked against fp64 before timing):
 *   B-sgemm-vlog  FFTW-c2-TW2 (FFTW 3.3.10 MEASURE, two windowed frames per complex DFT, fast post-pass) for all
 *                 frames, Accelerate cblas_sgemm mel, vvlog10f + vDSP clamps (= am-2's WM-FFTW-c2-TW2)
 *   B-sgemm-flog  the same with the NEON fast log (below) instead of vvlog10f
 *   B-nmel-flog   FFTW-c2-TW2 pair by pair; NEON sparse mel on the pair while it is in L1; NEON fast log
 *   S-sgemm-vlog  am-2's SMEr-P for all frames, then sgemm + vvlog10f (= am-2's WM-SMEr-P)
 *   S-sgemm-flog  SMEr-P + sgemm + fast log
 *   S-nmel-flog   SMEr-P in chunks of C pairs, NEON sparse mel + fast log per chunk
 *   S-zamel-flog  NEW: one streaming call does STFT + power + mel (outer products over bins in ZA, 32 frames per
 *                 group); only 80 values per frame leave streaming mode; NEON fast log + clamp pass after.
 * Fast log: log2 by exponent split + degree-8 polynomial (fitlog.py; 4.8e-8 abs in log2), folded with the
 *   (x+4)/4 scaling; the clamp max-8 becomes max(y, ymax-2) on the scaled value (identical in exact arithmetic).
 * Usage: lm check | time | mt SPEC   (env F frames=3000, C pairs per chunk=32, PROGS=comma list, BREAK=1 broken copy)
 */
#include "../../timing-harness/bench.h"
#include <Accelerate/Accelerate.h>
#include <arm_sme.h>
#include <arm_neon.h>
#include <fftw3.h>
#include <pthread.h>
#include <stdatomic.h>

#define M1 25
#define SFX 400
#include "kern.h"
#undef M1
#undef SFX
#include "kern_r.h"

#define NB 201  /* bins */
#ifndef NM
#define NM 80   /* mels (Whisper large-v3: 128, build with -DNM=128) */
#endif
#define PR 208  /* padded power row in the ZA-mel kernel */

/* ---------------- mel filterbank (librosa.filters.mel sr=16000 n_fft=400 n_mels=80, Slaney scale + norm) ------ */
static double hz2mel(double f) { double fsp = 200.0 / 3, mhz = 1000.0, mm = mhz / fsp, ls = log(6.4) / 27.0; return f < mhz ? f / fsp : mm + log(f / mhz) / ls; }
static double mel2hz(double m) { double fsp = 200.0 / 3, mhz = 1000.0, mm = mhz / fsp, ls = log(6.4) / 27.0; return m < mm ? fsp * m : mhz * exp(ls * (m - mm)); }
static double Wd[NM * NB]; /* fp64 filterbank [m][k] */
static void mel_filters(int brk) {
  double melf[NM + 2], fft_f[NB];
  for (int k = 0; k < NB; k++) fft_f[k] = 16000.0 / 400 * k;
  double m0 = hz2mel(0), m1 = hz2mel(8000.0);
  for (int i = 0; i < NM + 2; i++) melf[i] = mel2hz(m0 + (m1 - m0) * i / (NM + 1));
  for (int i = 0; i < NM; i++) {
    double fd0 = melf[i + 1] - melf[i], fd1 = melf[i + 2] - melf[i + 1], en = 2.0 / (melf[i + 2] - melf[i]);
    for (int k = 0; k < NB; k++) {
      double lo = -(melf[i] - fft_f[k]) / fd0, up = (melf[i + 2] - fft_f[k]) / fd1;
      double v = lo < up ? lo : up; Wd[i * NB + k] = (v > 0 ? v : 0) * en;
    }
  }
  if (brk) { /* broken copy: mel band 40 shifted up by one bin */
    for (int k = NB - 1; k > 0; k--) Wd[40 * NB + k] = Wd[40 * NB + k - 1];
    Wd[40 * NB] = 0;
  }
}
/* tables derived from Wd */
static float Wf[NM * NB] __attribute__((aligned(128)));        /* fp32 dense [m][k] for sgemm */
static float WT[(PR + 16) * NM] __attribute__((aligned(128))); /* fp32 transposed [k][m], zero for k >= 201 */
static int MT_LO[NM / 16], MT_LEN[NM / 16];                                /* bin range per tile of 16 mels (len multiple of 4) */
static float WC[NM * 48] __attribute__((aligned(128)));        /* compact NEON weights */
static int MLO[NM], MLEN4[NM], MOFF[NM];
static void make_mel_tables(int brk) {
  mel_filters(brk);
  memset(WT, 0, sizeof WT); memset(WC, 0, sizeof WC);
  for (int m = 0; m < NM; m++) for (int k = 0; k < NB; k++) { Wf[m * NB + k] = (float)Wd[m * NB + k]; WT[k * NM + m] = (float)Wd[m * NB + k]; }
  int off = 0;
  for (int m = 0; m < NM; m++) {
    int lo = NB, hi = 0;
    for (int k = 0; k < NB; k++) if (Wd[m * NB + k] != 0) { if (k < lo) lo = k; hi = k + 1; }
    if (hi == 0) { lo = 0; hi = 1; }
    int l4 = (hi - lo + 3) & ~3;
    MLO[m] = lo; MLEN4[m] = l4; MOFF[m] = off;
    for (int j = 0; j < l4; j++) WC[off + j] = (lo + j < NB) ? (float)Wd[m * NB + lo + j] : 0.f;
    off += l4;
  }
  for (int t = 0; t < NM / 16; t++) {
    int lo = NB, hi = 0;
    for (int m = 16 * t; m < 16 * t + 16; m++) { if (MLO[m] < lo) lo = MLO[m]; if (MLO[m] + MLEN4[m] > hi) hi = MLO[m] + MLEN4[m]; }
    int len = (hi - lo + 3) & ~3;
    if (lo + len > PR) lo = PR - len;
    MT_LO[t] = lo; MT_LEN[t] = len;
  }
}

/* ---------------- NEON fast log, Whisper-normalised: y = (log10(max(x,1e-10)) + 4) / 4 ---------------- */
static const float LOGC[8] = {1.442694772e+00f, -7.213571490e-01f, 4.809394436e-01f, -3.600872104e-01f, 2.867074859e-01f, -2.500691317e-01f, 2.368901786e-01f, -1.457439579e-01f};
static inline float32x4_t wlog4(float32x4_t x) {
  x = vmaxq_f32(x, vdupq_n_f32(1e-10f));
  int32x4_t xi = vreinterpretq_s32_f32(x);
  int32x4_t e = vshrq_n_s32(vsubq_s32(xi, vdupq_n_s32(0x3f3504f3)), 23);
  float32x4_t t = vsubq_f32(vreinterpretq_f32_s32(vsubq_s32(xi, vshlq_n_s32(e, 23))), vdupq_n_f32(1.f));
  float32x4_t p = vfmaq_f32(vdupq_n_f32(LOGC[6]), vdupq_n_f32(LOGC[7]), t);
  p = vfmaq_f32(vdupq_n_f32(LOGC[5]), p, t);
  p = vfmaq_f32(vdupq_n_f32(LOGC[4]), p, t);
  p = vfmaq_f32(vdupq_n_f32(LOGC[3]), p, t);
  p = vfmaq_f32(vdupq_n_f32(LOGC[2]), p, t);
  p = vfmaq_f32(vdupq_n_f32(LOGC[1]), p, t);
  p = vfmaq_f32(vdupq_n_f32(LOGC[0]), p, t);
  float32x4_t l2 = vfmaq_f32(vcvtq_f32_s32(e), p, t);       /* log2 x */
  return vfmaq_f32(vdupq_n_f32(1.f), l2, vdupq_n_f32(0.075257498915995f)); /* log10(2)/4 */
}
/* in-place over n (multiple of 16) values, returns running max */
static inline float32x4_t wlog_block(float *restrict v, long n, float32x4_t mx) {
  float32x4_t m0 = mx, m1 = mx;
  for (long i = 0; i < n; i += 8) {
    float32x4_t a = wlog4(vld1q_f32(v + i)), b = wlog4(vld1q_f32(v + i + 4));
    vst1q_f32(v + i, a); vst1q_f32(v + i + 4, b);
    m0 = vmaxq_f32(m0, a); m1 = vmaxq_f32(m1, b);
  }
  return vmaxq_f32(m0, m1);
}
static void clamp_pass(float *restrict v, long n, float mx) {
  float32x4_t th = vdupq_n_f32(mx - 2.f);
  for (long i = 0; i < n; i += 8) { vst1q_f32(v + i, vmaxq_f32(vld1q_f32(v + i), th)); vst1q_f32(v + i + 4, vmaxq_f32(vld1q_f32(v + i + 4), th)); }
}

/* ---------------- NEON sparse mel for two frames ---------------- */
static inline void mel_neon2_v1(const float *restrict pa, const float *restrict pb, float *restrict ma, float *restrict mb) {
  for (int m = 0; m < NM; m++) {
    const float *restrict w = WC + MOFF[m];
    const float *restrict a = pa + MLO[m], *restrict b = pb + MLO[m];
    float32x4_t sa = vdupq_n_f32(0), sb = sa;
    for (int j = 0; j < MLEN4[m]; j += 4) { float32x4_t wv = vld1q_f32(w + j); sa = vfmaq_f32(sa, wv, vld1q_f32(a + j)); sb = vfmaq_f32(sb, wv, vld1q_f32(b + j)); }
    ma[m] = vaddvq_f32(sa); mb[m] = vaddvq_f32(sb);
  }
}
/* v2: groups of 4 mels, the 4 dot products reduced together with 2 levels of pairwise adds, one vector store */
static inline void mel_neon2(const float *restrict pa, const float *restrict pb, float *restrict ma, float *restrict mb) {
  for (int g = 0; g < NM / 4; g++) {
    float32x4_t A0, A1, A2, A3, B0, B1, B2, B3;
#define MELQ(q, SA, SB)                                                                                          \
    {                                                                                                            \
      const int m = 4 * g + q; const float *restrict w = WC + MOFF[m];                                           \
      const float *restrict a = pa + MLO[m], *restrict b = pb + MLO[m];                                          \
      float32x4_t wv = vld1q_f32(w); SA = vmulq_f32(wv, vld1q_f32(a)); SB = vmulq_f32(wv, vld1q_f32(b));        \
      for (int j = 4; j < MLEN4[m]; j += 4) { wv = vld1q_f32(w + j); SA = vfmaq_f32(SA, wv, vld1q_f32(a + j)); SB = vfmaq_f32(SB, wv, vld1q_f32(b + j)); } \
    }
    MELQ(0, A0, B0) MELQ(1, A1, B1) MELQ(2, A2, B2) MELQ(3, A3, B3)
#undef MELQ
    vst1q_f32(ma + 4 * g, vpaddq_f32(vpaddq_f32(A0, A1), vpaddq_f32(A2, A3)));
    vst1q_f32(mb + 4 * g, vpaddq_f32(vpaddq_f32(B0, B1), vpaddq_f32(B2, B3)));
  }
}

#if NM == 80
#include "mel_gen80.h"
#else
#define mel_gen2 mel_neon2
#endif
/* ---------------- FFTW two-frame post-pass (am-2's post_pair2i) ---------------- */
static void post_pair2i(const float *restrict o, int N, float *restrict pa, float *restrict pb) {
  int h = N / 2;
  pa[0] = o[0] * o[0]; pb[0] = o[1] * o[1];
  pa[h] = o[2 * h] * o[2 * h]; pb[h] = o[2 * h + 1] * o[2 * h + 1];
  for (int k = 1; k < h; k++) {
    float xr = o[2 * k], xi = o[2 * k + 1], yr = o[2 * (N - k)], yi = o[2 * (N - k) + 1];
    float ar = xr + yr, ai = xi - yi, br = xi + yi, bi = yr - xr;
    pa[k] = 0.25f * (ar * ar + ai * ai);
    pb[k] = 0.25f * (br * br + bi * bi);
  }
}

/* ---------------- NEW: STFT + power + mel in one streaming call ---------------- */
#define GPF 16 /* pairs per mel group = 32 frames */
static _Thread_local float PG[32 * PR + 64] __attribute__((aligned(128)));
__arm_locally_streaming __arm_new("za") static void sme_lm400(const float *sig, long hop, const float *win, float *mel, long B) {
  svbool_t pg = svptrue_b32(), p1 = svwhilelt_b32(0, 1);
  svbool_t prv = svbic_b_z(pg, svwhilelt_b32(0, 15), svwhilelt_b32(0, 6));
  svcount_t pc = svptrue_c32();
  const float *W4 = W4tab_400, *FA = FAtab_400, *FB = FBtab_400; float *S = S4_400;
  float *Pg = PG;
  const long span = RKN + hop;
  for (int q = 1; q < PFD && q < B; q++)
    for (long l = 0; l < span; l += 32) __builtin_prefetch(sig + 2 * hop * q + l, 0, 2);
  for (long b = 0; b < B; b++) {
    const float *xa = sig + 2 * hop * b, *xb = xa + hop;
    float *pa = Pg + 2 * PR * (b % GPF), *pb = pa + PR;
    if (PFD && b + PFD < B)
      for (long l = 0; l < span; l += 32) __builtin_prefetch(sig + 2 * hop * (b + PFD) + l, 0, 2);
    svzero_mask_za(0xff);
    for (int n1 = 0; n1 < 25; n1++) {
      svfloat32_t fr = svld1_f32(pg, FA + 64 * n1), fi = svld1_f32(pg, FA + 64 * n1 + 32);
      svfloat32_t a = svld1_f32(pg, xa + 16 * n1), c = svld1_f32(pg, xb + 16 * n1);
      svfloat32_t w = svld1_f32(pg, win + 16 * n1); a = svmul_f32_x(pg, a, w); c = svmul_f32_x(pg, c, w);
      svmopa_za32_f32_m(0, pg, pg, a, fr); svmopa_za32_f32_m(1, pg, pg, a, fi);
      svmopa_za32_f32_m(2, pg, pg, c, fr); svmopa_za32_f32_m(3, pg, pg, c, fi);
    }
    TWB(0, 0, W4) TWB(1, 0, W4) TWB(2, 0, W4) TWB(3, 0, W4)
    TWB(0, 2, W4) TWB(1, 2, W4) TWB(2, 2, W4) TWB(3, 2, W4)
    for (int n2 = 0; n2 < 16; n2++) {
      svstr_za(4 * n2 + 0, S + 16 * n2); svstr_za(4 * n2 + 1, S + 256 + 16 * n2);
      svstr_za(4 * n2 + 2, S + 512 + 16 * n2); svstr_za(4 * n2 + 3, S + 768 + 16 * n2);
    }
    svzero_mask_za(0xff);
    for (int n2 = 0; n2 < 16; n2++) {
      svfloat32x2_t f = svld1_f32_x2(pc, FB + 32 * n2);
      svfloat32_t fr = svget2(f, 0), fi = svget2(f, 1);
      svfloat32_t yra = svld1_f32(pg, S + 16 * n2), yia = svld1_f32(pg, S + 256 + 16 * n2);
      svfloat32_t yrb = svld1_f32(pg, S + 512 + 16 * n2), yib = svld1_f32(pg, S + 768 + 16 * n2);
      svmopa_za32_f32_m(0, pg, pg, fr, yra); svmopa_za32_f32_m(1, pg, pg, fr, yia);
      svmopa_za32_f32_m(2, pg, pg, fr, yrb); svmopa_za32_f32_m(3, pg, pg, fr, yib);
      if (n2 != 0 && n2 != 8) {
        svmops_za32_f32_m(0, pg, pg, fi, yia); svmopa_za32_f32_m(1, pg, pg, fi, yra);
        svmops_za32_f32_m(2, pg, pg, fi, yib); svmopa_za32_f32_m(3, pg, pg, fi, yrb);
      }
    }
    RP_POWER(0, 1, 0x22, pa) RP_POWER(2, 3, 0x88, pb)
    if (b % GPF == GPF - 1 || b == B - 1) {
      /* mel phase: 32 frames' power rows (Pg, stride PR) -> mel[32][80]. Tiles 2/3: transposition of frames
       * 0..15 / 16..31 (LDR rows, read vertical slices = one bin over 16 frames); tiles 0/1: accumulators
       * (frames x 16 mels) for one tile of 16 mels at a time, over that tile's bin range only. */
      float *mo = mel + NM * 2 * (b - b % GPF);
      for (int t = 0; t < NM / 16; t++) {
        svzero_mask_za(0x33);
        const int lo = MT_LO[t], end = MT_LO[t] + MT_LEN[t];
        for (int k0 = lo; k0 < end; k0 += 16) {
          int nb = end - k0; if (nb > 16) nb = 16;
          for (int f = 0; f < 16; f++) { svldr_za(4 * f + 2, Pg + PR * f + k0); svldr_za(4 * f + 3, Pg + PR * (16 + f) + k0); }
          for (int j = 0; j < nb; j += 4) {
            svfloat32x4_t va = svread_ver_za32_f32_vg4(2, j), vb = svread_ver_za32_f32_vg4(3, j);
            const float *w = WT + NM * (k0 + j) + 16 * t;
            svfloat32_t w0 = svld1_f32(pg, w), w1 = svld1_f32(pg, w + NM), w2 = svld1_f32(pg, w + 2 * NM), w3 = svld1_f32(pg, w + 3 * NM);
            svmopa_za32_f32_m(0, pg, pg, svget4(va, 0), w0); svmopa_za32_f32_m(1, pg, pg, svget4(vb, 0), w0);
            svmopa_za32_f32_m(0, pg, pg, svget4(va, 1), w1); svmopa_za32_f32_m(1, pg, pg, svget4(vb, 1), w1);
            svmopa_za32_f32_m(0, pg, pg, svget4(va, 2), w2); svmopa_za32_f32_m(1, pg, pg, svget4(vb, 2), w2);
            svmopa_za32_f32_m(0, pg, pg, svget4(va, 3), w3); svmopa_za32_f32_m(1, pg, pg, svget4(vb, 3), w3);
          }
        }
        for (int f = 0; f < 16; f++) { svst1_hor_za32(0, f, pg, mo + NM * f + 16 * t); svst1_hor_za32(1, f, pg, mo + NM * (16 + f) + 16 * t); }
      }
    }
  }
}

/* ---------------- programs ---------------- */
enum { K_B_SGV, K_B_SGF, K_B_NMF, K_S_SGV, K_S_SGF, K_S_NMF, K_S_ZAF, K_B_GMF, K_X_STFT, K_X_ZANL, K_X_FSTFT, K_X_FLOG, NK };
static const char *kname[NK] = {"B-sgemm-vlog", "B-sgemm-flog", "B-nmel-flog", "S-sgemm-vlog", "S-sgemm-flog", "S-nmel-flog", "S-zamel-flog", "B-gmel-flog", "XS-stft-only", "XS-zamel-nolog", "XB-stft-only", "X-flog-only"};
typedef struct {
  int kind, F, C, hop;
  const float *sig; float *win;
  float *P, *Pc, *out, *buf, *buf2;
  fftwf_plan fp;
  volatile long frames; int qos; /* mt */
} lctx;

static float *g_sig; static long g_need; static float *g_win;
static void sgemm_mel(lctx *c) { cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, c->F, NM, NB, 1.f, c->P, NB, Wf, NB, 0.f, c->out, NM); }
static void vlog_all(lctx *c) {
  int M = c->F * NM;
  float lo = 1e-10f; vDSP_vthr(c->out, 1, &lo, c->out, 1, M);
  vvlog10f(c->out, c->out, &M);
  float mx; vDSP_maxv(c->out, 1, &mx, M);
  float th = mx - 8.f; vDSP_vthr(c->out, 1, &th, c->out, 1, M);
  float a = 0.25f, b1 = 1.f; vDSP_vsmsa(c->out, 1, &a, &b1, c->out, 1, M);
}
static void flog_all(lctx *c) {
  long M = (long)c->F * NM;
  float32x4_t mx = wlog_block(c->out, M, vdupq_n_f32(-1e30f));
  clamp_pass(c->out, M, vmaxvq_f32(mx));
}
static void fftw_c2_all(lctx *c) { /* FFTW-c2-TW2 into P[F][201] */
  const int N = 400, hop = c->hop; float *restrict b = c->buf, *restrict o = c->buf2; const float *restrict win = c->win;
  for (long pp = 0; pp < c->F / 2; pp++) {
    const float *restrict xa = c->sig + 2L * hop * pp, *restrict xb = xa + hop;
    for (int n = 0; n < N; n++) { b[2 * n] = xa[n] * win[n]; b[2 * n + 1] = xb[n] * win[n]; }
    fftwf_execute_dft(c->fp, (fftwf_complex *)b, (fftwf_complex *)o);
    post_pair2i(o, N, c->P + NB * 2 * pp, c->P + NB * (2 * pp + 1));
  }
}
/* macOS allocates a thread's TLS block lazily on first access (malloc): that must not happen in streaming mode
 * (SIGILL), so touch the kernels' _Thread_local scratch from normal mode first. */
static void tls_touch(void) { volatile float *a = S4_400, *b = PG; a[0] = 0; b[32 * PR + 63] = 0; }
static void run_core(lctx *c) {
  const int hop = c->hop, F = c->F;
  tls_touch();
  switch (c->kind) {
    case K_B_SGV: fftw_c2_all(c); sgemm_mel(c); vlog_all(c); break;
    case K_B_SGF: fftw_c2_all(c); sgemm_mel(c); flog_all(c); break;
    case K_S_SGV: sme_real400_pw(c->sig, hop, c->win, c->P, F / 2); sgemm_mel(c); vlog_all(c); break;
    case K_S_SGF: sme_real400_pw(c->sig, hop, c->win, c->P, F / 2); sgemm_mel(c); flog_all(c); break;
    case K_B_NMF:
    case K_B_GMF: {
      /* software-pipelined: the mel of pair pp-1 runs after the FFT of pair pp, so its power rows were stored long
       * before they are loaded (an unaligned vector load right after scalar stores to the same bytes stalls) */
      const int N = 400; float *restrict b = c->buf, *restrict o = c->buf2; const float *restrict win = c->win;
      float *pr[2] = {c->Pc, c->Pc + 2 * PR};
      float32x4_t mx = vdupq_n_f32(-1e30f);
      for (long pp = 0; pp < F / 2; pp++) {
        const float *restrict xa = c->sig + 2L * hop * pp, *restrict xb = xa + hop;
        for (int n = 0; n < N; n++) { b[2 * n] = xa[n] * win[n]; b[2 * n + 1] = xb[n] * win[n]; }
        fftwf_execute_dft(c->fp, (fftwf_complex *)b, (fftwf_complex *)o);
        if (pp > 0) {
          float *ma = c->out + NM * 2 * (pp - 1), *q = pr[(pp - 1) & 1];
          if (c->kind == K_B_GMF) mel_gen2(q, q + PR, ma, ma + NM); else mel_neon2(q, q + PR, ma, ma + NM);
          mx = wlog_block(ma, 2 * NM, mx);
        }
        float *q = pr[pp & 1];
        post_pair2i(o, N, q, q + PR);
      }
      { long pp = F / 2; float *ma = c->out + NM * 2 * (pp - 1), *q = pr[(pp - 1) & 1]; if (c->kind == K_B_GMF) mel_gen2(q, q + PR, ma, ma + NM); else mel_neon2(q, q + PR, ma, ma + NM); mx = wlog_block(ma, 2 * NM, mx); }
      clamp_pass(c->out, (long)F * NM, vmaxvq_f32(mx));
    } break;
    case K_S_NMF: {
      float32x4_t mx = vdupq_n_f32(-1e30f);
      long np = F / 2;
      for (long p0 = 0; p0 < np; p0 += c->C) {
        long cc = np - p0 < c->C ? np - p0 : c->C;
        sme_real400_pw(c->sig + 2L * hop * p0, hop, c->win, c->Pc, cc);
        for (long i = 0; i < cc; i++) {
          float *ma = c->out + NM * 2 * (p0 + i);
          mel_neon2(c->Pc + 2 * NB * i, c->Pc + 2 * NB * i + NB, ma, ma + NM);
        }
        mx = wlog_block(c->out + NM * 2 * p0, 2 * NM * cc, mx);
      }
      clamp_pass(c->out, (long)F * NM, vmaxvq_f32(mx));
    } break;
    case K_S_ZAF: sme_lm400(c->sig, hop, c->win, c->out, F / 2); flog_all(c); break;
    /* ablations (outputs are not the log-mel; not checked) */
    case K_X_STFT: sme_real400_pw(c->sig, hop, c->win, c->P, F / 2); break;
    case K_X_ZANL: sme_lm400(c->sig, hop, c->win, c->out, F / 2); break;
    case K_X_FSTFT: fftw_c2_all(c); break;
    case K_X_FLOG: flog_all(c); break; /* in place on a buffer filled once in setup: values stay finite */
  }
}
static void run_prog(void *v, long it) { (void)it; run_core((lctx *)v); }

static int setup(lctx *c, int kind, int F) {
  memset(c, 0, sizeof *c);
  c->kind = kind; c->F = F; c->hop = 160;
  c->C = getenv("C") ? atoi(getenv("C")) : 32;
  long need = 160L * (F - 1) + 400 + 64;
  if (!g_sig || g_need < need) {
    g_sig = bench_alloc(sizeof(float) * need); g_need = need;
    bench_rng_state = 0x1234567ull; bench_fill_f(g_sig, need);
    g_win = bench_alloc(sizeof(float) * 400);
    for (int n = 0; n < 400; n++) g_win[n] = (float)(0.5 - 0.5 * cos(2 * M_PI * n / 400));
  }
  c->sig = g_sig; c->win = g_win;
  c->P = bench_alloc(sizeof(float) * ((long)NB * F + 64));
  c->Pc = bench_alloc(sizeof(float) * (2L * NB * c->C + 4 * PR + 64));
  c->out = bench_alloc(sizeof(float) * (long)NM * (F + 34));
  if (kind == K_X_FLOG) for (long i = 0; i < (long)F * NM; i++) c->out[i] = 1.f + (float)(i % 977);
  c->buf = bench_alloc(sizeof(float) * 2 * 400 + 64); c->buf2 = bench_alloc(sizeof(float) * 2 * 400 + 64);
  if (kind == K_B_SGV || kind == K_B_SGF || kind == K_B_NMF || kind == K_B_GMF || kind == K_X_FSTFT) {
    static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER; /* FFTW planner is not thread safe */
    pthread_mutex_lock(&mu);
    c->fp = fftwf_plan_dft_1d(400, (fftwf_complex *)c->buf, (fftwf_complex *)c->buf2, FFTW_FORWARD, FFTW_MEASURE);
    pthread_mutex_unlock(&mu);
    if (!c->fp) return 0;
  }
  return 1;
}

/* ---------------- fp64 reference ---------------- */
static double *g_ref; static int g_refF;
static void make_ref(int F) {
  if (g_ref && g_refF == F) return;
  free(g_ref); g_ref = malloc(sizeof(double) * NM * F); g_refF = F;
  static double ct[400], st[400];
  for (int i = 0; i < 400; i++) { ct[i] = cos(2 * M_PI * i / 400); st[i] = -sin(2 * M_PI * i / 400); }
  double mx = -1e300;
  for (long f = 0; f < F; f++) {
    double x[400], p[NB];
    for (int n = 0; n < 400; n++) x[n] = (double)g_sig[160 * f + n] * (0.5 - 0.5 * cos(2 * M_PI * n / 400));
    for (int k = 0; k < NB; k++) { double re = 0, im = 0; for (int n = 0; n < 400; n++) { int j = (k * n) % 400; re += x[n] * ct[j]; im += x[n] * st[j]; } p[k] = re * re + im * im; }
    for (int m = 0; m < NM; m++) { double s = 0; for (int k = 0; k < NB; k++) s += Wd[m * NB + k] * p[k]; double l = log10(s > 1e-10 ? s : 1e-10); g_ref[f * NM + m] = l; if (l > mx) mx = l; }
  }
  for (long t = 0; t < (long)NM * F; t++) { double a = g_ref[t] > mx - 8 ? g_ref[t] : mx - 8; g_ref[t] = (a + 4) / 4; }
}
static int check_out(const char *label, const float *out, int F) {
  double e = 0, s2 = 0; long at = -1;
  for (long t = 0; t < (long)NM * F; t++) { double d = fabs((double)out[t] - g_ref[t]); if (!(d <= e)) { e = d; at = t; } s2 += d * d; }
  int ok = e < 1e-4;
  printf("CHECK %-14s F=%d all %d log-mel values vs fp64 (DFT, Slaney mel, Whisper log): max abs %.3e (frame %ld mel %ld) rms %.3e %s\n",
         label, F, NM * F, e, at / NM, at % NM, sqrt(s2 / ((double)NM * F)), ok ? "PASS" : "FAIL");
  return ok;
}

static int want(int k) {
  const char *s = getenv("PROGS"); if (!s) return 1;
  char buf[512]; snprintf(buf, sizeof buf, ",%s,", s); char key[64]; snprintf(key, sizeof key, ",%s,", kname[k]);
  return strstr(buf, key) != NULL;
}
static int find_prog(const char *n) { for (int k = 0; k < NK; k++) if (!strcmp(kname[k], n)) return k; return -1; }

/* ---------------- whole chip: a spinning pool; each thread runs whole 30 s chunks of its own program ---------------- */
typedef struct { lctx c; atomic_long frames; int qos; int stop; } mtt;
static atomic_int g_go, g_stop;
static void *mt_thread(void *v) {
  mtt *a = v;
  pthread_set_qos_class_self_np(a->qos ? QOS_CLASS_BACKGROUND : QOS_CLASS_USER_INTERACTIVE, 0);
  while (!atomic_load(&g_go)) ;
  while (!atomic_load(&g_stop)) { run_core(&a->c); atomic_fetch_add(&a->frames, a->c.F); }
  return NULL;
}
/* SPEC: comma list of PROG@P or PROG@E with optional *n, e.g. "S-zamel-flog@P,B-nmel-flog@P*3" */
static void mt_mode(const char *spec, int F) {
  mtt *th = calloc(32, sizeof(mtt)); int T = 0;
  char s[512]; snprintf(s, sizeof s, "%s", spec);
  for (char *tok = strtok(s, ","); tok; tok = strtok(NULL, ",")) {
    char name[64]; int n = 1; char cl = 'P';
    char *at = strchr(tok, '@'); if (!at) { printf("bad spec\n"); return; }
    snprintf(name, sizeof name, "%.*s", (int)(at - tok), tok); cl = at[1];
    char *st = strchr(at, '*'); if (st) n = atoi(st + 1);
    int k = find_prog(name); if (k < 0) { printf("unknown %s\n", name); return; }
    for (int i = 0; i < n; i++) { if (!setup(&th[T].c, k, F)) { printf("setup failed\n"); return; } th[T].qos = (cl == 'E'); T++; }
  }
  pthread_t pt[32];
  atomic_store(&g_go, 0); atomic_store(&g_stop, 0);
  for (int t = 0; t < T; t++) { atomic_store(&th[t].frames, 0); pthread_create(&pt[t], NULL, mt_thread, &th[t]); }
  atomic_store(&g_go, 1);
  struct timespec warm = {0, 300000000}; nanosleep(&warm, NULL);
  int R = getenv("R") ? atoi(getenv("R")) : 9; double res[64];
  long f0[32], f1[32];
  for (int r = 0; r < R; r++) {
    for (int t = 0; t < T; t++) f0[t] = atomic_load(&th[t].frames);
    uint64_t t0 = bench_now_ns();
    struct timespec w = {0, 400000000}; nanosleep(&w, NULL);
    uint64_t t1 = bench_now_ns();
    long tot = 0; for (int t = 0; t < T; t++) { f1[t] = atomic_load(&th[t].frames); tot += f1[t] - f0[t]; }
    res[r] = (double)(t1 - t0) / tot;
    printf("rep %d  MT %s  %.4f ns/frame aggregate = %.3f Mframes/s; per thread (Mframes/s):", r, spec, res[r], 1e3 / res[r]);
    for (int t = 0; t < T; t++) printf(" %.2f", (f1[t] - f0[t]) * 1e3 / (double)(t1 - t0));
    printf("\n");
  }
  atomic_store(&g_stop, 1);
  for (int t = 0; t < T; t++) pthread_join(pt[t], NULL);
  qsort(res, R, sizeof(double), bench_cmp_d);
  printf("RESULT MT %-60s T=%d best %.4f median %.4f ns/frame aggregate (median %.3f Mframes/s)\n", spec, T, res[0], res[R / 2], 1e3 / res[R / 2]);
  int allok = 1; char lab[64];
  for (int t = 0; t < T; t++) { snprintf(lab, sizeof lab, "t%d:%s", t, kname[th[t].c.kind]); allok &= check_out(lab, th[t].c.out, F); }
  printf("MTCHECK %s: all %d threads' last outputs %s\n", spec, T, allok ? "PASS" : "FAIL");
  fflush(stdout);
}


/* ---------------- end to end: 16 kHz int16 WAV from disk -> fp32 -> log-mel, 30 s chunks ---------------- */
#include <fcntl.h>
#include <unistd.h>
static void e2e_mode(const char *path, int F) {
  const long CH = 160L * F;                /* 480,000 new samples per 30 s chunk */
  const long NEED = 160L * (F - 1) + 400;  /* samples the 3000 frames read */
  int16_t *raw = bench_alloc(sizeof(int16_t) * (NEED + 64));
  float *sig = bench_alloc(sizeof(float) * (NEED + 64));
  for (int k = 0; k < NK; k++) if (want(k)) {
    lctx c; setup(&c, k, F); c.sig = sig;
    for (int pass = 0; pass < 3; pass++) {
      int fd = open(path, O_RDONLY); if (fd < 0) { perror(path); return; }
      off_t fsz = lseek(fd, 0, SEEK_END); long ns = (fsz - 44) / 2;
      double tr = 0, tc = 0, tf = 0; long chunks = 0; float chk = 0;
      for (long s0 = 0; s0 + NEED <= ns; s0 += CH) {
        uint64_t t0 = bench_now_ns();
        if (pread(fd, raw, sizeof(int16_t) * NEED, 44 + 2 * s0) != (ssize_t)(sizeof(int16_t) * NEED)) { printf("short read\n"); return; }
        uint64_t t1 = bench_now_ns();
        vDSP_vflt16(raw, 1, sig, 1, NEED); float sc = 1.f / 32768; vDSP_vsmul(sig, 1, &sc, sig, 1, NEED);
        uint64_t t2 = bench_now_ns();
        run_core(&c);
        uint64_t t3 = bench_now_ns();
        tr += t1 - t0; tc += t2 - t1; tf += t3 - t2; chunks++; chk += c.out[12345];
      }
      close(fd);
      double tot = tr + tc + tf, sec = chunks * 30.0;
      printf("E2E %-13s pass %d: %ld chunks (%.0f s of audio): read %.1f ms, int16->fp32 %.1f ms, log-mel %.1f ms, total %.1f ms per hour of audio; log-mel share %.1f %%; %.3f ns/frame log-mel (checksum %.4f)\n",
             kname[k], pass, chunks, sec, tr / 1e6 * 3600 / sec, tc / 1e6 * 3600 / sec, tf / 1e6 * 3600 / sec, tot / 1e6 * 3600 / sec, 100 * tf / tot, tf / (chunks * (double)F), chk);
    }
  }
  /* 48 kHz -> 16 kHz resampling of 30 s (1.44 M input samples), vDSP_desamp, windowed-sinc FIR, cutoff 7.2 kHz */
  int taps[2] = {48, 96};
  for (int ti = 0; ti < 2; ti++) {
    int P = taps[ti]; long nout = 480000; long nin = 3 * nout + P;
    float *x = bench_alloc(sizeof(float) * nin), *y = bench_alloc(sizeof(float) * nout), *h = bench_alloc(sizeof(float) * P);
    bench_fill_f(x, nin);
    for (int i = 0; i < P; i++) { double t = i - (P - 1) / 2.0, fc = 7200.0 / 48000; double sinc = t == 0 ? 2 * fc : sin(2 * M_PI * fc * t) / (M_PI * t); h[i] = (float)(sinc * (0.54 - 0.46 * cos(2 * M_PI * i / (P - 1)))); }
    double best = 1e30;
    for (int r = 0; r < 9; r++) { uint64_t t0 = bench_now_ns(); vDSP_desamp(x, 3, h, y, nout, P); double d = bench_now_ns() - t0; if (d < best) best = d; }
    printf("RESAMPLE 48k->16k vDSP_desamp %d taps: best %.2f ms per 30 s = %.1f ms per hour of audio = %.3f ns per 16 kHz output sample\n", P, best / 1e6, best / 1e6 * 120, best / nout);
  }
}

int main(int argc, char **argv) {
  setvbuf(stdout, NULL, _IOLBF, 0); bench_init(); make_tables_400();
  BLASSetThreading(BLAS_THREADING_SINGLE_THREADED);
  fftwf_import_wisdom_from_filename("$BASELINES/fftw/wisdom-m4-fftwf-measure.txt");
  int F = getenv("F") ? atoi(getenv("F")) : 3000; if (F & 1) F++;
  int brk = getenv("BREAK") ? atoi(getenv("BREAK")) : 0;
  make_mel_tables(0);
  if (argc < 2) return 2;
  if (!strcmp(argv[1], "tables")) {
    for (int t = 0; t < NM / 16; t++) printf("mel tile %d: bins %d..%d (%d)\n", t, MT_LO[t], MT_LO[t] + MT_LEN[t] - 1, MT_LEN[t]);
    int nz = 0; for (int i = 0; i < NM * NB; i++) nz += Wd[i] != 0; printf("nonzeros %d; NEON compact weights %d\n", nz, MOFF[NM - 1] + MLEN4[NM - 1]);
    FILE *fp = fopen(argc > 2 ? argv[2] : "melW_c.bin", "wb"); fwrite(Wd, sizeof(double), NM * NB, fp); fclose(fp);
    return 0;
  }
  lctx c0; setup(&c0, K_B_NMF, F); /* makes the signal */
  make_ref(F);
  if (brk == 1) make_mel_tables(1); /* the programs use the broken filterbank; the reference keeps the true one */
  if (brk == 2) { /* BREAK=2: mel 40 shifted one bin inside the compact NEON table only (what the generated code reads) */
    for (int j = MLEN4[40] - 1; j > 0; j--) WC[MOFF[40] + j] = WC[MOFF[40] + j - 1];
    WC[MOFF[40]] = 0;
  }
  if (!strcmp(argv[1], "check")) {
    for (int k = 0; k < K_X_STFT; k++) if (want(k)) { lctx c; if (!setup(&c, k, F)) continue; run_core(&c); run_core(&c); check_out(kname[k], c.out, F); }
    return 0;
  }
  if (!strcmp(argv[1], "mt")) { mt_mode(argv[2], F); return 0; }
  if (!strcmp(argv[1], "e2e")) { e2e_mode(argv[2], F); return 0; }
  if (!strcmp(argv[1], "time")) {
    bench_prog pr[NK]; lctx cx[NK]; bench_stat st[NK]; int n = 0;
    for (int k = 0; k < NK; k++) if (want(k)) {
      if (!setup(&cx[n], k, F)) continue;
      run_core(&cx[n]); if (k < K_X_STFT && !check_out(kname[k], cx[n].out, F)) return 1;
      pr[n] = (bench_prog){kname[k], run_prog, NULL, &cx[n]}; n++;
    }
    bench_opts o = BENCH_OPTS_DEFAULT; o.units = F; o.unit_name = "frame";
    char title[96]; snprintf(title, sizeof title, "Whisper log-mel N=400 hop=160 %d mels F=%d", NM, F);
    bench_compare(title, pr, n, &o, st);
    return 0;
  }
  return 2;
}
