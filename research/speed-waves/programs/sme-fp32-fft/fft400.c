/* fft400.c -- aj-2 (wave 11): batched 400-point complex FFT (fp32, forward, unscaled) on the SME2 matrix
 * unit. 400 is Whisper's n_fft; vDSP's DFT refuses it. Baselines: PFFFT z-order, FFTW 3.3.10 MEASURE
 * (both interleaved complex), and vDSP if its setup accepts 400.
 *
 * 400 = 25 x 16, n = 16 n1 + n2 (rows of 16 contiguous floats), k = k1 + 25 k2.
 *   pass A^T: tile(n2, k1) += X[n1,:] (x) F25[k1-block, n1]   k1 in two blocks: 0..15 (tiles 0 re, 1 im),
 *             16..24 (tiles 2 re, 3 im; rows 25..31 of F25 are zero) -> all four tiles accumulate.
 *   twiddle   Y'[k1][n2] = Y[k1][n2] w400^(k1 n2), in place in ZA (multi-vector FMLA on rows n2).
 *   stage     Y'^T rows (n2) stored to a 4 KB scratch with STR ZA (whole array vectors).
 *   pass B^T: tile(k2, k1-block) += F16[:, n2] (x) Y'^T[n2, k1-block]  (four tiles again)
 *   output    row k2 of tiles 0/1 = out[25 k2 + 0..15], of tiles 2/3 = out[25 k2 + 16..24] (9 lanes):
 *             natural order, split complex (re[400], im[400] per transform).
 * Not pipelined across transforms (each pass already uses the four tiles).
 * Usage: fft400 check | time same|separate B...
 */
#include "../timing-harness/bench.h"
#include <Accelerate/Accelerate.h>
#include <arm_sme.h>
#include <fftw3.h>
#include "pffft.h"

#define N 400
#ifndef PF400
#define PF400 8   /* L2 prefetch distance in transforms (0 = off) */
#endif
static float F25tab[25 * 64] __attribute__((aligned(128)));  /* per n1: Fr(k1 0..15), Fr(16..31), Fi(0..15), Fi(16..31) */
static float F16tab[16 * 32] __attribute__((aligned(128)));  /* per n2: Fr16(k2), Fi16(k2) */
static float W4tab[2 * 4 * 128] __attribute__((aligned(128))); /* block b, group j: Re-1 rows n2=j+4m (64), Im rows (64) */
static float S4[4 * 256] __attribute__((aligned(128)));

static void make_tables(void) {
  for (int n1 = 0; n1 < 25; n1++)
    for (int k = 0; k < 32; k++) {
      double th = 2 * M_PI * k * n1 / 25.0;
      F25tab[64 * n1 + k] = k < 25 ? (float)cos(th) : 0.f;
      F25tab[64 * n1 + 32 + k] = k < 25 ? (float)(-sin(th)) : 0.f;
    }
  for (int n2 = 0; n2 < 16; n2++)
    for (int k = 0; k < 16; k++) {
      double th = 2 * M_PI * k * n2 / 16.0;
      F16tab[32 * n2 + k] = (float)cos(th);
      F16tab[32 * n2 + 16 + k] = (float)(-sin(th));
    }
  for (int b = 0; b < 2; b++)
    for (int j = 0; j < 4; j++)
      for (int m = 0; m < 4; m++) {
        int n2 = j + 4 * m;
        for (int c = 0; c < 16; c++) {
          int k1 = 16 * b + c;
          double th = 2 * M_PI * k1 * n2 / 400.0;
          W4tab[512 * b + 128 * j + 16 * m + c] = k1 < 25 ? (float)(-2.0 * sin(th / 2) * sin(th / 2)) : 0.f;
          W4tab[512 * b + 128 * j + 64 + 16 * m + c] = k1 < 25 ? (float)(-sin(th)) : 0.f;
        }
      }
}

#define TWB(j, t0, wbase)                                              \
  {                                                                    \
    svfloat32x4_t a = svread_za32_f32_vg1x4(4 * (j) + (t0));           \
    svfloat32x4_t c = svread_za32_f32_vg1x4(4 * (j) + (t0) + 1);       \
    svfloat32x4_t wm = svld1_f32_x4(pc, (wbase) + 128 * (j));          \
    svfloat32x4_t ws = svld1_f32_x4(pc, (wbase) + 128 * (j) + 64);     \
    svmla_za32_f32_vg1x4(4 * (j) + (t0), a, wm);                       \
    svmls_za32_f32_vg1x4(4 * (j) + (t0), c, ws);                       \
    svmla_za32_f32_vg1x4(4 * (j) + (t0) + 1, a, ws);                   \
    svmla_za32_f32_vg1x4(4 * (j) + (t0) + 1, c, wm);                   \
  }

__arm_locally_streaming __arm_new("za") void fft400_sme(const float *re, const float *im, float *ore, float *oim, long B) {
  svbool_t pg = svptrue_b32(), p9 = svwhilelt_b32(0, 9);
  svcount_t pc = svptrue_c32();
  for (int q = 1; q < PF400 && q < B; q++) for (int l = 0; l < 13; l++) { __builtin_prefetch(re + N * q + 32 * l, 0, 2); __builtin_prefetch(im + N * q + 32 * l, 0, 2); }
  for (long b = 0; b < B; b++) {
    const float *xr = re + N * b, *xi = im + N * b;
    float *yr = ore + N * b, *yi = oim + N * b;
    if (PF400 && b + PF400 < B) for (int l = 0; l < 13; l++) { __builtin_prefetch(re + N * (b + PF400) + 32 * l, 0, 2); __builtin_prefetch(im + N * (b + PF400) + 32 * l, 0, 2); }
    svzero_mask_za(0xff);
    /* pass A^T */
    {
      svfloat32x4_t f = svld1_f32_x4(pc, F25tab);
      svfloat32_t ar = svld1_f32(pg, xr), ai = svld1_f32(pg, xi);
      svmopa_za32_f32_m(0, pg, pg, ar, svget4(f, 0)); svmopa_za32_f32_m(1, pg, pg, ai, svget4(f, 0));
      svmopa_za32_f32_m(2, pg, pg, ar, svget4(f, 1)); svmopa_za32_f32_m(3, pg, pg, ai, svget4(f, 1));
    }
    for (int n1 = 1; n1 < 25; n1++) {
      svfloat32x4_t f = svld1_f32_x4(pc, F25tab + 64 * n1);
      svfloat32_t ar = svld1_f32(pg, xr + 16 * n1), ai = svld1_f32(pg, xi + 16 * n1);
      svfloat32_t fra = svget4(f, 0), frb = svget4(f, 1), fia = svget4(f, 2), fib = svget4(f, 3);
      svmopa_za32_f32_m(0, pg, pg, ar, fra); svmopa_za32_f32_m(1, pg, pg, ar, fia);
      svmopa_za32_f32_m(2, pg, pg, ar, frb); svmopa_za32_f32_m(3, pg, pg, ar, fib);
      svmops_za32_f32_m(0, pg, pg, ai, fia); svmopa_za32_f32_m(1, pg, pg, ai, fra);
      svmops_za32_f32_m(2, pg, pg, ai, fib); svmopa_za32_f32_m(3, pg, pg, ai, frb);
    }
    /* twiddle in place: tiles (0,1) = k1 block 0, tiles (2,3) = block 1 */
    TWB(0, 0, W4tab) TWB(1, 0, W4tab) TWB(2, 0, W4tab) TWB(3, 0, W4tab)
    TWB(0, 2, W4tab + 512) TWB(1, 2, W4tab + 512) TWB(2, 2, W4tab + 512) TWB(3, 2, W4tab + 512)
    /* stage Y'^T rows */
    for (int n2 = 0; n2 < 16; n2++) {
      svstr_za(4 * n2 + 0, S4 + 16 * n2); svstr_za(4 * n2 + 1, S4 + 256 + 16 * n2);
      svstr_za(4 * n2 + 2, S4 + 512 + 16 * n2); svstr_za(4 * n2 + 3, S4 + 768 + 16 * n2);
    }
    svzero_mask_za(0xff);
    /* pass B^T */
    for (int n2 = 0; n2 < 16; n2++) {
      svfloat32x2_t f = svld1_f32_x2(pc, F16tab + 32 * n2);
      svfloat32_t fr = svget2(f, 0), fi = svget2(f, 1);
      svfloat32_t yra = svld1_f32(pg, S4 + 16 * n2), yia = svld1_f32(pg, S4 + 256 + 16 * n2);
      svfloat32_t yrb = svld1_f32(pg, S4 + 512 + 16 * n2), yib = svld1_f32(pg, S4 + 768 + 16 * n2);
      svmopa_za32_f32_m(0, pg, pg, fr, yra); svmopa_za32_f32_m(1, pg, pg, fr, yia);
      svmopa_za32_f32_m(2, pg, pg, fr, yrb); svmopa_za32_f32_m(3, pg, pg, fr, yib);
      if (n2 != 0 && n2 != 8) {
        svmops_za32_f32_m(0, pg, pg, fi, yia); svmopa_za32_f32_m(1, pg, pg, fi, yra);
        svmops_za32_f32_m(2, pg, pg, fi, yib); svmopa_za32_f32_m(3, pg, pg, fi, yrb);
      }
    }
    /* output, natural order */
    for (int k2 = 0; k2 < 16; k2++) {
      svst1_hor_za32(0, k2, pg, yr + 25 * k2); svst1_hor_za32(1, k2, pg, yi + 25 * k2);
      svst1_hor_za32(2, k2, p9, yr + 25 * k2 + 16); svst1_hor_za32(3, k2, p9, yi + 25 * k2 + 16);
    }
  }
}

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

enum { P_SME, P_PFFFTZ, P_PFFFTO, P_FFTW, P_VDFT, NP };
static const char *pname[NP] = {"SME2 25x16 (split)", "PFFFT z-order (interl.)", "PFFFT ordered (interl.)", "FFTW3f MEASURE (interl.)", "vDSP_DFT_Execute (split)"};
typedef struct {
  int kind; long B; bench_pool in, out; PFFFT_Setup *ps; float *work; fftwf_plan fp; vDSP_DFT_Setup ds;
} ctx;
static void prun(void *v, long i) {
  ctx *c = v;
  float *in = bench_pool_get(&c->in, i), *out = bench_pool_get(&c->out, i);
  long B = c->B, M = B * N;
  switch (c->kind) {
    case P_SME: fft400_sme(in, in + M, out, out + M, B); break;
    case P_PFFFTZ: for (long b = 0; b < B; b++) pffft_transform(c->ps, in + 2 * N * b, out + 2 * N * b, c->work, PFFFT_FORWARD); break;
    case P_PFFFTO: for (long b = 0; b < B; b++) pffft_transform_ordered(c->ps, in + 2 * N * b, out + 2 * N * b, c->work, PFFFT_FORWARD); break;
    case P_FFTW: for (long b = 0; b < B; b++) fftwf_execute_dft(c->fp, (fftwf_complex *)(in + 2 * N * b), (fftwf_complex *)(out + 2 * N * b)); break;
    case P_VDFT: for (long b = 0; b < B; b++) vDSP_DFT_Execute(c->ds, in + N * b, in + M + N * b, out + N * b, out + M + N * b); break;
  }
}
static int setup(ctx *c, int kind, long B, int mode) {
  memset(c, 0, sizeof *c); c->kind = kind; c->B = B;
  bench_pool_init(&c->in, sizeof(float) * 2 * N * B, mode);
  bench_pool_init(&c->out, sizeof(float) * 2 * N * B, mode);
  if (kind == P_PFFFTZ || kind == P_PFFFTO) { c->ps = pffft_new_setup(N, PFFFT_COMPLEX); c->work = pffft_aligned_malloc(sizeof(float) * 2 * N); if (!c->ps) return 0; }
  if (kind == P_FFTW) { c->fp = fftwf_plan_dft_1d(N, (fftwf_complex *)c->in.buf[0], (fftwf_complex *)c->out.buf[0], FFTW_FORWARD, FFTW_MEASURE); if (!c->fp) return 0; }
  if (kind == P_VDFT) { c->ds = vDSP_DFT_zop_CreateSetup(NULL, N, vDSP_DFT_FORWARD); if (!c->ds) return 0; }
  for (int p = 0; p < c->in.count; p++) bench_fill_f(c->in.buf[p], 2 * N * B);
  return 1;
}
static int check(int kind, long B) {
  ctx c;
  if (!setup(&c, kind, B, BENCH_SAME)) { printf("CHECK %-28s setup returned NULL (size %d not supported)\n", pname[kind], N); return 1; }
  float *in = c.in.buf[0], *out = c.out.buf[0]; long M = B * N;
  int inter = (kind == P_PFFFTZ || kind == P_PFFFTO || kind == P_FFTW);
  prun(&c, 0);
  double emax = 0, rmax = 0, Xr[N], Xi[N]; float ir[N], ii[N];
  float *tmp = bench_alloc(sizeof(float) * 2 * N);
  for (long b = 0; b < B; b++) {
    for (int n = 0; n < N; n++) {
      if (inter) { ir[n] = in[2 * (N * b + n)]; ii[n] = in[2 * (N * b + n) + 1]; } else { ir[n] = in[N * b + n]; ii[n] = in[M + N * b + n]; }
    }
    dft_ref(ir, ii, Xr, Xi);
    const float *gr, *gi; int st = 1;
    if (kind == P_PFFFTZ) { pffft_zreorder(c.ps, out + 2 * N * b, tmp, PFFFT_FORWARD); gr = tmp; gi = tmp + 1; st = 2; }
    else if (inter) { gr = out + 2 * N * b; gi = gr + 1; st = 2; }
    else { gr = out + N * b; gi = out + M + N * b; }
    for (int k = 0; k < N; k++) {
      double e = hypot(gr[st * k] - Xr[k], gi[st * k] - Xi[k]); if (!(e <= emax)) emax = e;
      double m = hypot(Xr[k], Xi[k]); if (m > rmax) rmax = m;
    }
  }
  int ok = emax / rmax < 1e-5;
  printf("CHECK %-28s B=%-4ld every transform vs fp64 DFT, natural order: max abs err %.3e  rel (to max |X|) %.3e  %s\n", pname[kind], B, emax, emax / rmax, ok ? "PASS" : "FAIL");
  return ok;
}
int main(int argc, char **argv) {
  bench_init(); make_tables();
  const char *wis = "$BASELINES/fftw/wisdom-m4-fftwf-measure.txt";
  fftwf_import_wisdom_from_filename(wis);
  if (argc < 2) return 2;
  if (!strcmp(argv[1], "check")) { int ok = 1; for (int k = 0; k < NP; k++) ok &= check(k, 16); ok &= check(P_SME, 1); ok &= check(P_SME, 3); return ok ? 0 : 1; }
  int mode = !strcmp(argv[2], "separate") ? BENCH_SEPARATE : BENCH_SAME;
  for (int a = 3; a < argc; a++) {
    long B = atol(argv[a]);
    ctx c[NP]; bench_prog p[NP]; bench_stat st[NP]; int np = 0;
    for (int k = 0; k < NP; k++) { if (!setup(&c[np], k, B, mode)) { printf("# %s: not supported at %d\n", pname[k], N); bench_pool_free(&c[np].in); bench_pool_free(&c[np].out); continue; } p[np] = (bench_prog){pname[k], prun, NULL, &c[np]}; np++; }
    char title[128]; snprintf(title, sizeof title, "cfft400 B=%ld %s", B, argv[2]);
    bench_opts o = BENCH_OPTS_DEFAULT; o.units = (double)N * B;
    bench_compare(title, p, np, &o, st);
  }
  return 0;
}
