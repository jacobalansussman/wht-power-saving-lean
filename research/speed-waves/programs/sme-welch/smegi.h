// smegi.h -- ap-0 (variant of smeg.h with optional interleaved complex input (pass 1, LD2D) and output (last pass,
// MOVA + ST2D), i.e. numpy's complex128 layout with no separate conversion pass)
// smeg.h -- ap-0 (wave 13): batched complex fp64 FFT of N = 8^P points (P = 2..5) on the M4's SME matrix unit.
// Split complex, transforms contiguous (stride N), out of place, forward, unnormalised, natural order.
// Digits n = sum n_i 8^{P-i} (n_1 most significant), k = sum k_j 8^{j-1}. Pass j contracts n_j -> k_j as a complex
// 8x8 matrix (32 fp64 FMOPA/FMOPS per group of 64 points):
//   pass 1, group g = (n_2..n_{P-1}): tile[n_P][k1] = sum_{n_1} x[..] F1_g[k1][n_1], then * w_N^{n_P k1} in ZA;
//           rows -> S1[n_2..n_P][k1]
//   pass j>=2, group g = (n_{j+1}..n_P ; k_{j-1}..k_2): tile[k_j][k1] = sum_{n_j} Fj_{G_n}[k_j][n_j] S_{j-1}[..][k1];
//           all twiddles of pass j are folded into Fj (they depend only on the group); rows -> S_j (S_P = output)
// Two groups (4 tiles) accumulate while the previous job's 4 tiles drain; passes of P different transforms are
// interleaved so a job never reads what the job just before it wrote; L2 prefetch of the input PFD transforms ahead.
#include <arm_sme.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int P; long N, NG;           // NG = 8^{P-2} groups per pass
  double *Fr[6], *Fi[6];       // pass j (1..P): [G][m][k] with G in 0..(pass 1: NG, pass j: 8^{P-j}) - 1
  double *TWm, *TWi;           // [2][4][8]
  double *Sr[6][2], *Si[6][2]; // scratch S_1..S_{P-1}, double buffered
  int broken;
  long rs, st; const double *win, *mu;  // as-0 inter = 4: row stride, frame step, window, means [2t], [2t+1]
} smeg_plan;

static void smeg_w(long e, long n, double *re, double *im) { e %= n; if (e < 0) e += n; double th = 2.0 * M_PI * (double)e / (double)n; *re = cos(th); *im = -sin(th); }
static void *smeg_alloc(size_t b) { void *p = NULL; posix_memalign(&p, 128, b); memset(p, 0, b); return p; }
static smeg_plan *smeg_new(int P, int broken) {
  smeg_plan *pl = calloc(1, sizeof *pl); pl->P = P; long N = 1; for (int i = 0; i < P; i++) N *= 8; pl->N = N; pl->NG = N / 64; pl->broken = broken;
  long p8[8]; p8[0] = 1; for (int i = 1; i < 8; i++) p8[i] = p8[i - 1] * 8;
  // pass 1
  pl->Fr[1] = smeg_alloc(8 * pl->NG * 64); pl->Fi[1] = smeg_alloc(8 * pl->NG * 64);
  for (long g = 0; g < pl->NG; g++) for (int m = 0; m < 8; m++) for (int k = 0; k < 8; k++)
    smeg_w((long)k * (m * p8[P - 1] + 8 * g), N, &pl->Fr[1][(g * 8 + m) * 8 + k], &pl->Fi[1][(g * 8 + m) * 8 + k]);
  for (int j = 2; j <= P; j++) {
    long nG = p8[P - j];
    pl->Fr[j] = smeg_alloc(8 * nG * 64); pl->Fi[j] = smeg_alloc(8 * nG * 64);
    for (long G = 0; G < nG; G++) for (int m = 0; m < 8; m++) for (int k = 0; k < 8; k++)
      smeg_w((long)k * p8[j - 1] * (m * p8[P - j] + G), N, &pl->Fr[j][(G * 8 + m) * 8 + k], &pl->Fi[j][(G * 8 + m) * 8 + k]);
  }
  pl->TWm = smeg_alloc(8 * 64); pl->TWi = smeg_alloc(8 * 64);
  for (int h = 0; h < 2; h++) for (int i = 0; i < 4; i++) for (int k = 0; k < 8; k++) {
    int s = h + 2 * i; long e = ((long)s * k) % N; double th = 2.0 * M_PI * (double)e / (double)N, sn = sin(0.5 * th);
    pl->TWm[(h * 4 + i) * 8 + k] = -2.0 * sn * sn; pl->TWi[(h * 4 + i) * 8 + k] = -sin(th);
  }
  if (broken) pl->TWi[(1 * 4 + 0) * 8 + 3] = -pl->TWi[(1 * 4 + 0) * 8 + 3];
  for (int j = 1; j < P; j++) for (int b = 0; b < 2; b++) { pl->Sr[j][b] = smeg_alloc(8 * N); pl->Si[j][b] = smeg_alloc(8 * N); }
  return pl;
}

typedef struct { const double *ur, *ui, *vr, *vi; long us, vs; double *or_, *oi; long os; int tw; long t; int ii, io; const double *wp; double ma, mb; } ggrp;

#define PG svptrue_b64()
#define GACC2(A0, A1, A2, A3, g, m) do { \
    svfloat64_t ur0, ui0, ur1, ui1; \
    if ((g)[0].ii == 2) { svfloat64_t w0_ = svld1_f64(PG, (g)[0].wp + (m) * (g)[0].us), w1_ = svld1_f64(PG, (g)[1].wp + (m) * (g)[1].us); \
      ur0 = svmul_f64_x(PG, svsub_n_f64_x(PG, svld1_f64(PG, (g)[0].ur + (m) * (g)[0].us), (g)[0].ma), w0_); \
      ui0 = svmul_f64_x(PG, svsub_n_f64_x(PG, svld1_f64(PG, (g)[0].ui + (m) * (g)[0].us), (g)[0].mb), w0_); \
      ur1 = svmul_f64_x(PG, svsub_n_f64_x(PG, svld1_f64(PG, (g)[1].ur + (m) * (g)[1].us), (g)[1].ma), w1_); \
      ui1 = svmul_f64_x(PG, svsub_n_f64_x(PG, svld1_f64(PG, (g)[1].ui + (m) * (g)[1].us), (g)[1].mb), w1_); } \
    else if ((g)[0].ii) { svfloat64x2_t a_ = svld2_f64(PG, (g)[0].ur + (m) * (g)[0].us), b_ = svld2_f64(PG, (g)[1].ur + (m) * (g)[1].us); \
      ur0 = svget2_f64(a_, 0); ui0 = svget2_f64(a_, 1); ur1 = svget2_f64(b_, 0); ui1 = svget2_f64(b_, 1); } \
    else { ur0 = svld1_f64(PG, (g)[0].ur + (m) * (g)[0].us); ui0 = svld1_f64(PG, (g)[0].ui + (m) * (g)[0].us); \
      ur1 = svld1_f64(PG, (g)[1].ur + (m) * (g)[1].us); ui1 = svld1_f64(PG, (g)[1].ui + (m) * (g)[1].us); } \
    svfloat64_t vr0 = svld1_f64(PG, (g)[0].vr + (m) * (g)[0].vs), vi0 = svld1_f64(PG, (g)[0].vi + (m) * (g)[0].vs); \
    svfloat64_t vr1 = svld1_f64(PG, (g)[1].vr + (m) * (g)[1].vs), vi1 = svld1_f64(PG, (g)[1].vi + (m) * (g)[1].vs); \
    svmopa_za64_f64_m(A0, PG, PG, ur0, vr0); svmopa_za64_f64_m(A1, PG, PG, ur0, vi0); \
    svmopa_za64_f64_m(A2, PG, PG, ur1, vr1); svmopa_za64_f64_m(A3, PG, PG, ur1, vi1); \
    svmops_za64_f64_m(A0, PG, PG, ui0, vi0); svmopa_za64_f64_m(A1, PG, PG, ui0, vr0); \
    svmops_za64_f64_m(A2, PG, PG, ui1, vi1); svmopa_za64_f64_m(A3, PG, PG, ui1, vr1); } while (0)
#define GTW(DR, DI, h) do { \
    svcount_t pc_ = svptrue_c64(); \
    svfloat64x4_t R4 = svread_za64_f64_vg1x4((DR) + 8 * (h)), I4 = svread_za64_f64_vg1x4((DI) + 8 * (h)); \
    svfloat64x4_t Wm = svld1_f64_x4(pc_, twm + 32 * (h)), Wi = svld1_f64_x4(pc_, twi + 32 * (h)); \
    svmla_za64_f64_vg1x4((DR) + 8 * (h), R4, Wm); svmls_za64_f64_vg1x4((DR) + 8 * (h), I4, Wi); \
    svmla_za64_f64_vg1x4((DI) + 8 * (h), I4, Wm); svmla_za64_f64_vg1x4((DI) + 8 * (h), R4, Wi); } while (0)
#define GST(DR, DI, d, s) do { svst1_hor_za64(DR, s, PG, (d)->or_ + (s) * (d)->os); svst1_hor_za64(DI, s, PG, (d)->oi + (s) * (d)->os); } while (0)
#define GST4I(DR, DI, d, s0) do { svfloat64x4_t R4_ = svread_hor_za64_f64_vg4(DR, s0), I4_ = svread_hor_za64_f64_vg4(DI, s0); \
    svst2_f64(PG, (d)->or_ + ((s0) + 0) * (d)->os, svcreate2_f64(svget4_f64(R4_, 0), svget4_f64(I4_, 0))); \
    svst2_f64(PG, (d)->or_ + ((s0) + 1) * (d)->os, svcreate2_f64(svget4_f64(R4_, 1), svget4_f64(I4_, 1))); \
    svst2_f64(PG, (d)->or_ + ((s0) + 2) * (d)->os, svcreate2_f64(svget4_f64(R4_, 2), svget4_f64(I4_, 2))); \
    svst2_f64(PG, (d)->or_ + ((s0) + 3) * (d)->os, svcreate2_f64(svget4_f64(R4_, 3), svget4_f64(I4_, 3))); } while (0)
#define GAC1(d, s, q, R4_, I4_) do { double *a_ = (d)->or_ + (s) * (d)->os; svfloat64_t r_ = svget4_f64(R4_, q), i_ = svget4_f64(I4_, q); \
    svfloat64_t v_ = svld1_f64(PG, a_); v_ = svmla_f64_x(PG, v_, r_, r_); v_ = svmla_f64_x(PG, v_, i_, i_); svst1_f64(PG, a_, v_); } while (0)
#define GAC4I(DR, DI, d, s0) do { svfloat64x4_t R4_ = svread_hor_za64_f64_vg4(DR, s0), I4_ = svread_hor_za64_f64_vg4(DI, s0); \
    GAC1(d, (s0) + 0, 0, R4_, I4_); GAC1(d, (s0) + 1, 1, R4_, I4_); GAC1(d, (s0) + 2, 2, R4_, I4_); GAC1(d, (s0) + 3, 3, R4_, I4_); } while (0)
#define GOUT4(DR, DI, d, s0) do { if ((d)->io == 2) GAC4I(DR, DI, d, s0); else GST4I(DR, DI, d, s0); } while (0)
#define GJOBDEF(NAME, A0, A1, A2, A3, D0, D1, D2, D3, MASK) \
static inline __attribute__((always_inline)) void NAME(const ggrp *g, const ggrp *d, const double *twm, const double *twi) __arm_streaming __arm_inout("za") { \
  svzero_mask_za(MASK); \
  if (d && d[0].tw) { \
    GACC2(A0, A1, A2, A3, g, 0); GTW(D0, D1, 0); \
    GACC2(A0, A1, A2, A3, g, 1); GTW(D0, D1, 1); \
    GACC2(A0, A1, A2, A3, g, 2); GTW(D2, D3, 0); \
    GACC2(A0, A1, A2, A3, g, 3); GTW(D2, D3, 1); \
    GACC2(A0, A1, A2, A3, g, 4); GST(D0, D1, &d[0], 0); GST(D0, D1, &d[0], 1); GST(D2, D3, &d[1], 0); GST(D2, D3, &d[1], 1); \
    GACC2(A0, A1, A2, A3, g, 5); GST(D0, D1, &d[0], 2); GST(D0, D1, &d[0], 3); GST(D2, D3, &d[1], 2); GST(D2, D3, &d[1], 3); \
    GACC2(A0, A1, A2, A3, g, 6); GST(D0, D1, &d[0], 4); GST(D0, D1, &d[0], 5); GST(D2, D3, &d[1], 4); GST(D2, D3, &d[1], 5); \
    GACC2(A0, A1, A2, A3, g, 7); GST(D0, D1, &d[0], 6); GST(D0, D1, &d[0], 7); GST(D2, D3, &d[1], 6); GST(D2, D3, &d[1], 7); \
  } else if (d && d[0].io) { \
    GACC2(A0, A1, A2, A3, g, 0); GOUT4(D0, D1, &d[0], 0); GACC2(A0, A1, A2, A3, g, 1); \
    GACC2(A0, A1, A2, A3, g, 2); GOUT4(D0, D1, &d[0], 4); GACC2(A0, A1, A2, A3, g, 3); \
    GACC2(A0, A1, A2, A3, g, 4); GOUT4(D2, D3, &d[1], 0); GACC2(A0, A1, A2, A3, g, 5); \
    GACC2(A0, A1, A2, A3, g, 6); GOUT4(D2, D3, &d[1], 4); GACC2(A0, A1, A2, A3, g, 7); \
  } else if (d) { \
    _Pragma("clang loop unroll(full)") for (int m = 0; m < 8; m++) { GACC2(A0, A1, A2, A3, g, m); GST(D0, D1, &d[0], m); GST(D2, D3, &d[1], m); } \
  } else { \
    _Pragma("clang loop unroll(full)") for (int m = 0; m < 8; m++) GACC2(A0, A1, A2, A3, g, m); \
  } \
} \
static inline __attribute__((always_inline)) void NAME##_drain(const ggrp *d, const double *twm, const double *twi) __arm_streaming __arm_inout("za") { \
  if (d[0].tw) { GTW(D0, D1, 0); GTW(D0, D1, 1); GTW(D2, D3, 0); GTW(D2, D3, 1); } \
  if (d[0].io) { GOUT4(D0, D1, &d[0], 0); GOUT4(D0, D1, &d[0], 4); GOUT4(D2, D3, &d[1], 0); GOUT4(D2, D3, &d[1], 4); } \
  else for (int m = 0; m < 8; m++) { GST(D0, D1, &d[0], m); GST(D2, D3, &d[1], m); } \
}
GJOBDEF(gjobA, 0, 1, 2, 3, 4, 5, 6, 7, 0x0F)
GJOBDEF(gjobB, 4, 5, 6, 7, 0, 1, 2, 3, 0xF0)

static inline __attribute__((always_inline)) void smeg_groups(const smeg_plan *pl, ggrp *g, int j, long t, long p, const double *xr,
                                                              const double *xi, double *yr, double *yi, int inter) __arm_streaming {
  int P = pl->P; long N = pl->N, sP1 = N / 8;
  for (int q = 0; q < 2; q++) {
    long gi = 2 * p + q; ggrp *G = &g[q]; G->t = t; G->ii = 0; G->io = 0;
    if (j == 1 && inter == 4) {  // as-0: raw frames, detrend + window in pass 1
      G->ur = xr + t * pl->rs + 8 * gi; G->ui = G->ur + pl->st; G->us = sP1; G->ii = 2; G->wp = pl->win + 8 * gi; G->ma = pl->mu[2 * t]; G->mb = pl->mu[2 * t + 1];
      G->vr = pl->Fr[1] + gi * 64; G->vi = pl->Fi[1] + gi * 64; G->vs = 8;
      G->or_ = pl->Sr[1][t & 1] + 64 * gi; G->oi = pl->Si[1][t & 1] + 64 * gi; G->os = 8; G->tw = 1;
    } else if (j == 1 && (inter == 1 || inter == 2)) {
      G->ur = xr + 2 * (t * N + 8 * gi); G->ui = G->ur; G->us = 2 * sP1; G->ii = 1;
      G->vr = pl->Fr[1] + gi * 64; G->vi = pl->Fi[1] + gi * 64; G->vs = 8;
      G->or_ = pl->Sr[1][t & 1] + 64 * gi; G->oi = pl->Si[1][t & 1] + 64 * gi; G->os = 8; G->tw = 1;
    } else if (j == 1) {
      G->ur = xr + t * N + 8 * gi; G->ui = xi + t * N + 8 * gi; G->us = sP1;
      G->vr = pl->Fr[1] + gi * 64; G->vi = pl->Fi[1] + gi * 64; G->vs = 8;
      G->or_ = pl->Sr[1][t & 1] + 64 * gi; G->oi = pl->Si[1][t & 1] + 64 * gi; G->os = 8; G->tw = 1;
    } else {
      long kdig = 1; for (int i = 0; i < j - 2; i++) kdig *= 8;  // 8^{j-2}
      long Gn = gi / kdig, Gk = gi % kdig, pj = kdig * 64;         // 8^j
      G->ur = pl->Fr[j] + Gn * 64; G->ui = pl->Fi[j] + Gn * 64; G->us = 8;
      G->vr = pl->Sr[j - 1][t & 1] + 8 * gi; G->vi = pl->Si[j - 1][t & 1] + 8 * gi; G->vs = sP1;
      double *br = (j == P) ? yr + t * N : pl->Sr[j][t & 1], *bi = (j == P) ? yi + t * N : pl->Si[j][t & 1];
      G->or_ = br + Gn * pj + Gk * 8; G->oi = bi + Gn * pj + Gk * 8; G->os = kdig * 8; G->tw = 0;
      if (j == P && inter == 1) { G->or_ = yr + 2 * (t * N + Gk * 8); G->oi = G->or_; G->os = 2 * kdig * 8; G->io = 1; }
      if (j == P && (inter == 2 || inter == 3 || inter == 4)) { G->or_ = yr + Gk * 8; G->oi = G->or_; G->os = kdig * 8; G->io = 2; }  // as-0: acc[k] += |Z_k|^2
    }
  }
}
#ifndef GPFD
#define GPFD 2
#endif
// inter = 4 (as-0): pass 1 reads raw frames a = xr + t*rs, b = a + st, computes (x - mu) * win in streaming mode; last pass as 2
// inter = 3 (as-0): split input (xr, xi) as for 0, last pass accumulates as for 2
// inter = 2 (as-0): input interleaved as for 1; yr = acc[N] (double), acc[k] += |Y_t[k]|^2 summed over all B rows
// inter = 1: xr / yr are interleaved complex (re, im) arrays of B*N complex values; xi / yi ignored
__arm_locally_streaming __arm_new("za") void smeg_fft2(const smeg_plan *pl, const double *xr, const double *xi, double *yr, double *yi, long B, int inter) {
  ggrp pend[2], g[2]; int have = 0; long cnt = 0; int P = pl->P; long N = pl->N, npair = pl->NG / 2;
  const double *twm = pl->TWm, *twi = pl->TWi;
  long pf_per_job = (2 * N / 16 + npair * P - 1) / (npair * P);  // lines of 128 B per job so a round prefetches one transform
  for (long q = 0; q < GPFD && q < B; q++) for (long l = 0; l < 2 * N; l += 16) __builtin_prefetch(inter == 4 ? xr + q * pl->rs + l : (inter == 1 || inter == 2) ? xr + 2 * q * N + l : (l < N ? xr + q * N + l : xi + q * N + l - N), 0, 2);
  for (long r = 0; r < B + P - 1; r++) {
    long pfl = 0; long tp = r + GPFD;
    for (long p = 0; p < npair; p++)
      for (int j = 1; j <= P; j++) {
        if (tp < B) for (long u = 0; u < pf_per_job && pfl < 2 * N; u++, pfl += 16) {
          const double *a = inter == 4 ? xr + tp * pl->rs + pfl : (inter == 1 || inter == 2) ? xr + 2 * tp * N + pfl : (pfl < N ? xr + tp * N + pfl : xi + tp * N + (pfl - N)); __builtin_prefetch(a, 0, 2); }
        long t = r - (j - 1);
        if (t < 0 || t >= B) continue;
        smeg_groups(pl, g, j, t, p, xr, xi, yr, yi, inter);
        if (have && pend[0].t == t) { if (cnt & 1) gjobB_drain(pend, twm, twi); else gjobA_drain(pend, twm, twi); have = 0; }
        if ((cnt & 1) == 0) gjobA(g, have ? pend : NULL, twm, twi); else gjobB(g, have ? pend : NULL, twm, twi);
        pend[0] = g[0]; pend[1] = g[1]; have = 1; cnt++;
      }
  }
  if (have) { if (cnt & 1) gjobB_drain(pend, twm, twi); else gjobA_drain(pend, twm, twi); }
}
