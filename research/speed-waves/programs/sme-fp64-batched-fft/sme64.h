// sme64.h -- ap-0 (wave 13): batched complex fp64 FFT of 512 = 8 x 8 x 8 points on the M4's SME matrix unit.
// Split complex (re[], im[]), transforms contiguous (stride 512), out of place, forward, unnormalised, natural order.
// Plan: n = 64a + 8b + c, k = k1 + 8 k2 + 64 k3. Three radix-8 passes, each a complex 8x8 DFT matrix as 32 fp64
// FMOPA/FMOPS (8x8 tiles) per group of 64 points:
//   pass 1 (group b): tile[c][k1] = sum_a x[64a+8b+c] F1_b[k1][a], F1_b = diag(w64^{b k1}) F8; then x w512^{c k1}
//                    in ZA (multi-vector FMLA); rows -> S1[b][c][k1]
//   pass 2 (group c): tile[k2][k1] = sum_b F2_c[k2][b] S1[b][c][k1], F2_c = diag(w64^{c k2}) F8; rows -> S2[c][k2][k1]
//   pass 3 (group k2): tile[k3][k1] = sum_c F8[k3][c] S2[c][k2][k1]; rows -> y[k1 + 8 k2 + 64 k3] (natural order)
// Two groups (4 tiles) accumulate while the previous job's 4 tiles drain (twiddle + slice stores).
// Jobs of three transforms are interleaved (pass 1 of t, pass 2 of t-1, pass 3 of t-2) so no job reads what the
// job just before it wrote.
#include <arm_sme.h>
#include <math.h>
#include <string.h>

#define AL __attribute__((aligned(128)))
static double F1r[8][8][8] AL, F1i[8][8][8] AL;  // [b][a][k1]
static double F2r[8][8][8] AL, F2i[8][8][8] AL;  // [c][b][k2]
static double F3r[8][8] AL, F3i[8][8] AL;        // [c][k3]
static double TWm[2][4][8] AL, TWi[2][4][8] AL;  // vg1x4 group h: slices h, h+2, h+4, h+6: (Re w - 1), Im w, w = w512^{c k1}
static double S1r[2][512] AL, S1i[2][512] AL, S2r[2][512] AL, S2i[2][512] AL;
static int sme64_broken = 0;  // 1 = deliberately wrong sign of one twiddle (test must fail)

static void wpow(long e, long n, double *re, double *im, double *rem1) {  // w_n^e = exp(-2 pi i e / n)
  e %= n; if (e < 0) e += n;
  double th = 2.0 * M_PI * (double)e / (double)n;
  *re = cos(th); *im = -sin(th);
  if (rem1) { double s = sin(0.5 * th); *rem1 = -2.0 * s * s; }
}
static void sme64_init(void) {
  for (int b = 0; b < 8; b++) for (int a = 0; a < 8; a++) for (int k = 0; k < 8; k++) {
    // F1_b[k1][a] = w8^{a k1} w64^{b k1} = w64^{8 a k1 + b k1}
    wpow(8L * a * k + (long)b * k, 64, &F1r[b][a][k], &F1i[b][a][k], NULL);
    wpow(8L * a * k + (long)b * k, 64, &F2r[b][a][k], &F2i[b][a][k], NULL);  // same formula with (c, b, k2)
  }
  for (int c = 0; c < 8; c++) for (int k = 0; k < 8; k++) wpow((long)c * k, 8, &F3r[c][k], &F3i[c][k], NULL);
  for (int h = 0; h < 2; h++) for (int i = 0; i < 4; i++) for (int k = 0; k < 8; k++) {
    int s = h + 2 * i; double re;
    wpow((long)s * k, 512, &re, &TWi[h][i][k], &TWm[h][i][k]);
  }
  if (sme64_broken) TWi[1][0][3] = -TWi[1][0][3];
}

typedef struct { const double *ur, *ui, *vr, *vi; long us, vs; double *or_, *oi; long os; int tw; long t; } grp_t;

#define PG svptrue_b64()
#define ACC(TR, TI, g, m) do { \
    svfloat64_t ur_ = svld1_f64(PG, (g)->ur + (m) * (g)->us), ui_ = svld1_f64(PG, (g)->ui + (m) * (g)->us); \
    svfloat64_t vr_ = svld1_f64(PG, (g)->vr + (m) * (g)->vs), vi_ = svld1_f64(PG, (g)->vi + (m) * (g)->vs); \
    svmopa_za64_f64_m(TR, PG, PG, ur_, vr_); svmopa_za64_f64_m(TI, PG, PG, ur_, vi_); \
    svmops_za64_f64_m(TR, PG, PG, ui_, vi_); svmopa_za64_f64_m(TI, PG, PG, ui_, vr_); } while (0)
// interleave the two groups so 4 tiles rotate
#define ACC2(A0, A1, A2, A3, g, m) do { \
    svfloat64_t ur0 = svld1_f64(PG, (g)[0].ur + (m) * (g)[0].us), ui0 = svld1_f64(PG, (g)[0].ui + (m) * (g)[0].us); \
    svfloat64_t vr0 = svld1_f64(PG, (g)[0].vr + (m) * (g)[0].vs), vi0 = svld1_f64(PG, (g)[0].vi + (m) * (g)[0].vs); \
    svfloat64_t ur1 = svld1_f64(PG, (g)[1].ur + (m) * (g)[1].us), ui1 = svld1_f64(PG, (g)[1].ui + (m) * (g)[1].us); \
    svfloat64_t vr1 = svld1_f64(PG, (g)[1].vr + (m) * (g)[1].vs), vi1 = svld1_f64(PG, (g)[1].vi + (m) * (g)[1].vs); \
    svmopa_za64_f64_m(A0, PG, PG, ur0, vr0); svmopa_za64_f64_m(A1, PG, PG, ur0, vi0); \
    svmopa_za64_f64_m(A2, PG, PG, ur1, vr1); svmopa_za64_f64_m(A3, PG, PG, ur1, vi1); \
    svmops_za64_f64_m(A0, PG, PG, ui0, vi0); svmopa_za64_f64_m(A1, PG, PG, ui0, vr0); \
    svmops_za64_f64_m(A2, PG, PG, ui1, vi1); svmopa_za64_f64_m(A3, PG, PG, ui1, vr1); } while (0)
#define TW(DR, DI, h) do { \
    svcount_t pc_ = svptrue_c64(); \
    svfloat64x4_t R4 = svread_za64_f64_vg1x4((DR) + 8 * (h)), I4 = svread_za64_f64_vg1x4((DI) + 8 * (h)); \
    svfloat64x4_t Wm = svld1_f64_x4(pc_, &TWm[h][0][0]), Wi = svld1_f64_x4(pc_, &TWi[h][0][0]); \
    svmla_za64_f64_vg1x4((DR) + 8 * (h), R4, Wm); svmls_za64_f64_vg1x4((DR) + 8 * (h), I4, Wi); \
    svmla_za64_f64_vg1x4((DI) + 8 * (h), I4, Wm); svmla_za64_f64_vg1x4((DI) + 8 * (h), R4, Wi); } while (0)
#ifndef USE_STR
#define ST(DR, DI, d, s) do { svst1_hor_za64(DR, s, PG, (d)->or_ + (s) * (d)->os); svst1_hor_za64(DI, s, PG, (d)->oi + (s) * (d)->os); } while (0)
#else
#define ST(DR, DI, d, s) do { svstr_za(8 * (s) + (DR), (d)->or_ + (s) * (d)->os); svstr_za(8 * (s) + (DI), (d)->oi + (s) * (d)->os); } while (0)
#endif

// One job: accumulate groups g[0], g[1] into tiles (A0,A1),(A2,A3); drain d[0], d[1] from (D0,D1),(D2,D3).
#define JOBDEF(NAME, A0, A1, A2, A3, D0, D1, D2, D3, MASK) \
static inline __attribute__((always_inline)) void NAME(const grp_t *g, const grp_t *d) __arm_streaming __arm_inout("za") { \
  svzero_mask_za(MASK); \
  if (d && d[0].tw) { \
    ACC2(A0, A1, A2, A3, g, 0); TW(D0, D1, 0); \
    ACC2(A0, A1, A2, A3, g, 1); TW(D0, D1, 1); \
    ACC2(A0, A1, A2, A3, g, 2); TW(D2, D3, 0); \
    ACC2(A0, A1, A2, A3, g, 3); TW(D2, D3, 1); \
    ACC2(A0, A1, A2, A3, g, 4); ST(D0, D1, &d[0], 0); ST(D0, D1, &d[0], 1); ST(D2, D3, &d[1], 0); ST(D2, D3, &d[1], 1); \
    ACC2(A0, A1, A2, A3, g, 5); ST(D0, D1, &d[0], 2); ST(D0, D1, &d[0], 3); ST(D2, D3, &d[1], 2); ST(D2, D3, &d[1], 3); \
    ACC2(A0, A1, A2, A3, g, 6); ST(D0, D1, &d[0], 4); ST(D0, D1, &d[0], 5); ST(D2, D3, &d[1], 4); ST(D2, D3, &d[1], 5); \
    ACC2(A0, A1, A2, A3, g, 7); ST(D0, D1, &d[0], 6); ST(D0, D1, &d[0], 7); ST(D2, D3, &d[1], 6); ST(D2, D3, &d[1], 7); \
  } else if (d) { \
    _Pragma("clang loop unroll(full)") for (int m = 0; m < 8; m++) { ACC2(A0, A1, A2, A3, g, m); ST(D0, D1, &d[0], m); ST(D2, D3, &d[1], m); } \
  } else { \
    _Pragma("clang loop unroll(full)") for (int m = 0; m < 8; m++) ACC2(A0, A1, A2, A3, g, m); \
  } \
} \
static inline __attribute__((always_inline)) void NAME##_drain(const grp_t *d) __arm_streaming __arm_inout("za") { \
  if (d[0].tw) { TW(D0, D1, 0); TW(D0, D1, 1); TW(D2, D3, 0); TW(D2, D3, 1); } \
  for (int m = 0; m < 8; m++) { ST(D0, D1, &d[0], m); ST(D2, D3, &d[1], m); } \
}
JOBDEF(jobA, 0, 1, 2, 3, 4, 5, 6, 7, 0x0F)
JOBDEF(jobB, 4, 5, 6, 7, 0, 1, 2, 3, 0xF0)

static inline __attribute__((always_inline)) void mkgroups(grp_t *g, int pass, long t, int p, const double *xr, const double *xi,
                                                           double *yr, double *yi) __arm_streaming {
  for (int j = 0; j < 2; j++) {
    int q = 2 * p + j; grp_t *G = &g[j]; G->t = t;
    if (pass == 1) {
      G->ur = xr + t * 512 + 8 * q; G->ui = xi + t * 512 + 8 * q; G->us = 64;
      G->vr = &F1r[q][0][0]; G->vi = &F1i[q][0][0]; G->vs = 8;
      G->or_ = S1r[t & 1] + 64 * q; G->oi = S1i[t & 1] + 64 * q; G->os = 8; G->tw = 1;
    } else if (pass == 2) {
      G->ur = &F2r[q][0][0]; G->ui = &F2i[q][0][0]; G->us = 8;
      G->vr = S1r[t & 1] + 8 * q; G->vi = S1i[t & 1] + 8 * q; G->vs = 64;
      G->or_ = S2r[t & 1] + 64 * q; G->oi = S2i[t & 1] + 64 * q; G->os = 8; G->tw = 0;
    } else {
      G->ur = &F3r[0][0]; G->ui = &F3i[0][0]; G->us = 8;
      G->vr = S2r[t & 1] + 8 * q; G->vi = S2i[t & 1] + 8 * q; G->vs = 64;
      G->or_ = yr + t * 512 + 8 * q; G->oi = yi + t * 512 + 8 * q; G->os = 64; G->tw = 0;
    }
  }
}

__arm_locally_streaming __arm_new("za") void fft512_sme(const double *xr, const double *xi, double *yr, double *yi, long B) {
  grp_t pend[2], g[2]; int have = 0; long cnt = 0;
#ifndef PFD
#define PFD 0
#endif
  if (PFD) for (long q = 0; q < PFD && q < B; q++) for (int l = 0; l < 512; l += 16) { __builtin_prefetch(xr + q * 512 + l, 0, 2); __builtin_prefetch(xi + q * 512 + l, 0, 2); }
  for (long r = 0; r < B + 2; r++) {
    if (PFD && r + PFD < B) for (int l = 0; l < 512; l += 16) { __builtin_prefetch(xr + (r + PFD) * 512 + l, 0, 2); __builtin_prefetch(xi + (r + PFD) * 512 + l, 0, 2); }
    for (int p = 0; p < 4; p++)
      for (int pass = 1; pass <= 3; pass++) {
        long t = r - (pass - 1);
        if (t < 0 || t >= B) continue;
        mkgroups(g, pass, t, p, xr, xi, yr, yi);
        if (have && pend[0].t == t) {  // the next job depends on the pending one: drain first
          if (cnt & 1) jobB_drain(pend); else jobA_drain(pend);  // job cnt-1 even = jobA = tiles 0-3, drained by jobB_drain
          have = 0;
        }
        if ((cnt & 1) == 0) jobA(g, have ? pend : NULL); else jobB(g, have ? pend : NULL);
        pend[0] = g[0]; pend[1] = g[1]; have = 1; cnt++;
      }
  }
  if (have) { if (cnt & 1) jobB_drain(pend); else jobA_drain(pend); }
}
