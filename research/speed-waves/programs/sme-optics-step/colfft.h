// colfft.h -- as-1 (wave 14): FFTs along axis 0 (the COLUMNS) of a row-major complex128 array on the M4's SME unit.
// No transpose: a column FFT of a row-major array is a matrix product whose data operand is a contiguous piece of a
// ROW (8 consecutive columns), so every radix-R stage is "8x8 complex DFT matrix (twiddles folded in) times 8 rows".
// Stockham autosort, stage s (radix R, Ns = product of earlier radices, NR = N/R), group j = 0..NR-1, k = j mod Ns:
//   y[(j/Ns)*Ns*R + k + o*Ns][:] = sum_i M_k[o][i] x[j + i*NR][:],   M_k[o][i] = w_{Ns*R}^{i*(o*Ns + k)}, w_n = e^{-2 pi i/n}
// (formula checked against numpy in stockham_check.py). All twiddles live in the per-group matrices: no elementwise
// twiddle pass at all. ZA tiles: rows = output index o, columns = 8 consecutive array columns; 4 FMOPA/FMOPS per input
// row per 8 columns (32 per 8x8 complex block, as ap-0). Interleaved complex everywhere (LD2D in, MOVA + ST2D out).
// Two column blocks (16 columns) per tile set; two tile sets {0..3},{4..7}: one accumulates while the other drains.
// Radices: the small one (2 or 4) first (one matrix), then 8s.
#include <arm_sme.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int col_pfd = 0;  // software prefetch distance (groups) for the first stage, which reads the big array
typedef struct { long N; int S; int R[8]; long Ns[8]; double *Mr[8], *Mi[8]; } colplan;
static void *col_alloc(size_t b) { void *p = NULL; posix_memalign(&p, 128, b); memset(p, 0, b); return p; }
static colplan *col_new_o(long N, int broken, int small_last) {
  colplan *p = calloc(1, sizeof *p); p->N = N; long m = N; int S = 0;
  while (m % 8 == 0 && m > 1) m /= 8;
  if (m != 1 && !small_last) p->R[S++] = (int)m;  // 2 or 4 first
  for (long t = N / (m == 1 ? 1 : m); t > 1; t /= 8) p->R[S++] = 8;
  if (m != 1 && small_last) p->R[S++] = (int)m;
  p->S = S; long Ns = 1;
  for (int s = 0; s < S; s++) {
    int R = p->R[s]; p->Ns[s] = Ns; long NsR = Ns * R;
    p->Mr[s] = col_alloc(8 * 64 * Ns); p->Mi[s] = col_alloc(8 * 64 * Ns);
    for (long k = 0; k < Ns; k++) for (int i = 0; i < R; i++) for (int o = 0; o < R; o++) {
      long e = ((long)i * (o * Ns + k)) % NsR; double th = 2.0 * M_PI * (double)e / (double)NsR;
      p->Mr[s][(k * 8 + i) * 8 + o] = cos(th); p->Mi[s][(k * 8 + i) * 8 + o] = -sin(th);
    }
    Ns *= R;
  }
  if (broken) { int s = S - 1; long k = 1; p->Mi[s][(k * 8 + 1) * 8 + 2] = -p->Mi[s][(k * 8 + 1) * 8 + 2]; }  // one twiddle sign
  return p;
}
static colplan *col_new(long N, int broken) { return col_new_o(N, broken, 0); }

#ifndef PG
#define PG svptrue_b64()
#endif
#define CMOPA(T0, T1, T2, T3, i) do { \
    svfloat64x2_t a_ = svld2_f64(PG, ip[i] + c2), b_ = svld2_f64(PG, ip[i] + c2 + 16); \
    svfloat64_t xr_ = svget2_f64(a_, 0), xi_ = svget2_f64(a_, 1), yr_ = svget2_f64(b_, 0), yi_ = svget2_f64(b_, 1); \
    svmopa_za64_f64_m(T0, PG, PG, ur##i, xr_); svmopa_za64_f64_m(T1, PG, PG, ur##i, xi_); \
    svmopa_za64_f64_m(T2, PG, PG, ur##i, yr_); svmopa_za64_f64_m(T3, PG, PG, ur##i, yi_); \
    svmops_za64_f64_m(T0, PG, PG, ui##i, xi_); svmopa_za64_f64_m(T1, PG, PG, ui##i, xr_); \
    svmops_za64_f64_m(T2, PG, PG, ui##i, yi_); svmopa_za64_f64_m(T3, PG, PG, ui##i, yr_); } while (0)
#define LOADU(i) do { ur##i = svld1_f64(PG, Mr + (k * 8 + i) * 8); ui##i = svld1_f64(PG, Mi + (k * 8 + i) * 8); \
    ip[i] = in + 2 * (j + i * NR) * ins; op[i] = out + 2 * (base + i * Ns) * outs; } while (0)
#define LOAD8 do { LOADU(0); LOADU(1); LOADU(2); LOADU(3); LOADU(4); LOADU(5); LOADU(6); LOADU(7); } while (0)
#define LOAD4 do { LOADU(0); LOADU(1); LOADU(2); LOADU(3); } while (0)
#define LOAD2 do { LOADU(0); LOADU(1); } while (0)
#define CST4(TR, TI, s0, co) do { svfloat64x4_t R4_ = svread_hor_za64_f64_vg4(TR, s0), I4_ = svread_hor_za64_f64_vg4(TI, s0); \
    svst2_f64(PG, pop[(s0) + 0] + (co), svcreate2_f64(svget4_f64(R4_, 0), svget4_f64(I4_, 0))); \
    svst2_f64(PG, pop[(s0) + 1] + (co), svcreate2_f64(svget4_f64(R4_, 1), svget4_f64(I4_, 1))); \
    svst2_f64(PG, pop[(s0) + 2] + (co), svcreate2_f64(svget4_f64(R4_, 2), svget4_f64(I4_, 2))); \
    svst2_f64(PG, pop[(s0) + 3] + (co), svcreate2_f64(svget4_f64(R4_, 3), svget4_f64(I4_, 3))); } while (0)
#define CST2(TR, TI, co) do { svfloat64x2_t R2_ = svread_hor_za64_f64_vg2(TR, 0), I2_ = svread_hor_za64_f64_vg2(TI, 0); \
    svst2_f64(PG, pop[0] + (co), svcreate2_f64(svget2_f64(R2_, 0), svget2_f64(I2_, 0))); \
    svst2_f64(PG, pop[1] + (co), svcreate2_f64(svget2_f64(R2_, 1), svget2_f64(I2_, 1))); } while (0)

// body: accumulate into (A0..A3), drain previous set (D0..D3) interleaved
#define BODY8(A0, A1, A2, A3, D0, D1, D2, D3) do { \
    CMOPA(A0, A1, A2, A3, 0); CMOPA(A0, A1, A2, A3, 1); if (have) CST4(D0, D1, 0, pco); \
    CMOPA(A0, A1, A2, A3, 2); CMOPA(A0, A1, A2, A3, 3); if (have) CST4(D0, D1, 4, pco); \
    CMOPA(A0, A1, A2, A3, 4); CMOPA(A0, A1, A2, A3, 5); if (have) CST4(D2, D3, 0, pco + 16); \
    CMOPA(A0, A1, A2, A3, 6); CMOPA(A0, A1, A2, A3, 7); if (have) CST4(D2, D3, 4, pco + 16); } while (0)
#define DRAIN8(D0, D1, D2, D3) do { CST4(D0, D1, 0, pco); CST4(D0, D1, 4, pco); CST4(D2, D3, 0, pco + 16); CST4(D2, D3, 4, pco + 16); } while (0)
#define BODY4(A0, A1, A2, A3, D0, D1, D2, D3) do { \
    CMOPA(A0, A1, A2, A3, 0); CMOPA(A0, A1, A2, A3, 1); if (have) CST4(D0, D1, 0, pco); \
    CMOPA(A0, A1, A2, A3, 2); CMOPA(A0, A1, A2, A3, 3); if (have) CST4(D2, D3, 0, pco + 16); } while (0)
#define DRAIN4(D0, D1, D2, D3) do { CST4(D0, D1, 0, pco); CST4(D2, D3, 0, pco + 16); } while (0)
#define BODY2(A0, A1, A2, A3, D0, D1, D2, D3) do { \
    CMOPA(A0, A1, A2, A3, 0); if (have) CST2(D0, D1, pco); CMOPA(A0, A1, A2, A3, 1); if (have) CST2(D2, D3, pco + 16); } while (0)
#define DRAIN2(D0, D1, D2, D3) do { CST2(D0, D1, pco); CST2(D2, D3, pco + 16); } while (0)

// one stage: in rows (stride ins complex) -> out rows (stride outs complex); ncol complex columns (multiple of 16)
#define STAGEFN(NAME, R, BODY, DRAIN, LOAD) \
static void NAME(const colplan *p, int s, const double *in, long ins, double *out, long outs, long ncol) __arm_streaming __arm_inout("za") { \
  long N = p->N, NR = N / R, Ns = p->Ns[s]; const double *Mr = p->Mr[s], *Mi = p->Mi[s]; \
  double *pop[8]; long pco = 0; int have = 0, cur = 0; \
  for (long j = 0; j < NR; j++) { \
    long k = j % Ns, base = (j / Ns) * Ns * R + k; \
    svfloat64_t ur0, ur1, ur2, ur3, ur4, ur5, ur6, ur7, ui0, ui1, ui2, ui3, ui4, ui5, ui6, ui7; \
    const double *ip[8]; double *op[8]; LOAD; \
    if (s == 0 && col_pfd && j + col_pfd < NR) for (int i_ = 0; i_ < R; i_++) { \
      const double *q_ = in + 2 * (j + col_pfd + i_ * NR) * ins; for (long l_ = 0; l_ < 2 * ncol; l_ += 16) __builtin_prefetch(q_ + l_, 0, 3); } \
    for (long c2 = 0; c2 < 2 * ncol; c2 += 32) { \
      if (cur == 0) { svzero_mask_za(0x0F); BODY(0, 1, 2, 3, 4, 5, 6, 7); } \
      else { svzero_mask_za(0xF0); BODY(4, 5, 6, 7, 0, 1, 2, 3); } \
      _Pragma("clang loop unroll(full)") for (int i = 0; i < R; i++) pop[i] = op[i]; \
      pco = c2; have = 1; cur ^= 1; \
    } \
  } \
  if (have) { if (cur == 1) DRAIN(0, 1, 2, 3); else DRAIN(4, 5, 6, 7); } \
}
STAGEFN(col_stage8, 8, BODY8, DRAIN8, LOAD8)
STAGEFN(col_stage4, 4, BODY4, DRAIN4, LOAD4)
STAGEFN(col_stage2, 2, BODY2, DRAIN2, LOAD2)

// full column FFT of a strip of ncol columns: in (row stride ins) -> out (row stride outs); s0, s1: scratch N*ncol complex.
// out may equal in (the first stage reads everything before the last stage writes).
static void col_fft_r(const colplan *p, int sb, int se, const double *in, long ins, double *out, long outs, long ncol, double *s0, double *s1) __arm_streaming __arm_inout("za") {
  const double *src = in; long ss = ins;
  for (int s = sb; s < se; s++) {
    double *dst; long ds;
    if (s == se - 1) { dst = out; ds = outs; } else { dst = ((s - sb) & 1) ? s1 : s0; ds = ncol; }
    if (p->R[s] == 8) col_stage8(p, s, src, ss, dst, ds, ncol);
    else if (p->R[s] == 4) col_stage4(p, s, src, ss, dst, ds, ncol);
    else col_stage2(p, s, src, ss, dst, ds, ncol);
    src = dst; ss = ds;
  }
}
static void col_fft(const colplan *p, const double *in, long ins, double *out, long outs, long ncol, double *s0, double *s1) __arm_streaming __arm_inout("za") {
  col_fft_r(p, 0, p->S, in, ins, out, outs, ncol, s0, s1);
}
__arm_locally_streaming __arm_new("za") void col_fft_rls(const colplan *p, int sb, int se, const double *in, long ins, double *out, long outs, long ncol, double *s0, double *s1) {
  col_fft_r(p, sb, se, in, ins, out, outs, ncol, s0, s1);
}
__arm_locally_streaming __arm_new("za") void col_fft_ls(const colplan *p, const double *in, long ins, double *out, long outs, long ncol, double *s0, double *s1) {
  col_fft(p, in, ins, out, outs, ncol, s0, s1);
}
