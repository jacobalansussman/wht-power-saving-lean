/* fft_bench.c -- fp32 forward FFT baselines on the M4: Apple vDSP (fft_zop / DFT_Execute split,
 * DFT interleaved), FFTW 3.3.10 (NEON, FFTW_MEASURE) and PFFFT (NEON, ordered and z-domain order).
 * All out of place, input left unchanged (so repeated calls never overflow).
 * Usage: fft_bench c|r same|separate N...     (c = complex, r = real input)
 * Every program is checked against FFTW double precision (FFTW_ESTIMATE) before timing; FFTW double
 * itself is checked against a naive double DFT at N <= 4096.
 * Build: see baselines/README.md */
#include "../timing-harness/bench.h"
#include <Accelerate/Accelerate.h>
#include <fftw3.h>
#include "pffft.h"

enum { P_VFFT, P_VDFT, P_VDFTI, P_FFTW, P_PFFFT_ORD, P_PFFFT_Z, NPROG };
static const char *pname[NPROG] = {"vDSP_fft_zop(split)", "vDSP_DFT_Execute(split)", "vDSP_DFT_Interleaved",
                                   "FFTW3f-NEON MEASURE", "PFFFT ordered", "PFFFT z-order"};

typedef struct {
  int kind, real, n, log2n;
  bench_pool in, out;
  float *work;
  FFTSetup fs;
  vDSP_DFT_Setup ds;
  vDSP_DFT_Interleaved_Setup dsi;
  fftwf_plan fp;
  PFFFT_Setup *ps;
} fctx;

/* layout of one input buffer: complex: split = re[n] then im[n]; interleaved = (re,im)*n.
   real: split (vDSP even/odd) = re[n/2] im[n/2]; plain = x[n]. We size every buffer 2n floats. */
static void f_run(void *v, long i) {
  fctx *c = v;
  float *in = bench_pool_get(&c->in, i), *out = bench_pool_get(&c->out, i);
  int n = c->n;
  switch (c->kind) {
    case P_VFFT: {
      int h = c->real ? n / 2 : n;
      DSPSplitComplex a = {in, in + h}, b = {out, out + h};
      if (c->real) vDSP_fft_zrop(c->fs, &a, 1, &b, 1, c->log2n, kFFTDirection_Forward);
      else vDSP_fft_zop(c->fs, &a, 1, &b, 1, c->log2n, kFFTDirection_Forward);
    } break;
    case P_VDFT: {
      int h = c->real ? n / 2 : n;
      vDSP_DFT_Execute(c->ds, in, in + h, out, out + h);
    } break;
    case P_VDFTI: vDSP_DFT_Interleaved_Execute(c->dsi, (const DSPComplex *)in, (DSPComplex *)out); break;
    case P_FFTW:
      if (c->real) fftwf_execute_dft_r2c(c->fp, in, (fftwf_complex *)out);
      else fftwf_execute_dft(c->fp, (fftwf_complex *)in, (fftwf_complex *)out);
      break;
    case P_PFFFT_ORD: pffft_transform_ordered(c->ps, in, out, c->work, PFFFT_FORWARD); break;
    case P_PFFFT_Z: pffft_transform(c->ps, in, out, c->work, PFFFT_FORWARD); break;
  }
}

/* write the program's input buffer from a plain signal x (complex: 2n interleaved; real: n) */
static void put_input(fctx *c, float *buf, const double *x) {
  int n = c->n;
  int split = (c->kind == P_VFFT || c->kind == P_VDFT);
  if (c->real) {
    if (split) for (int k = 0; k < n / 2; k++) { buf[k] = x[2 * k]; buf[n / 2 + k] = x[2 * k + 1]; }
    else for (int k = 0; k < n; k++) buf[k] = x[k];
  } else {
    if (split) for (int k = 0; k < n; k++) { buf[k] = x[2 * k]; buf[n + k] = x[2 * k + 1]; }
    else for (int k = 0; k < 2 * n; k++) buf[k] = x[k];
  }
}
/* read the output into standard interleaved spectrum: complex n bins, real n/2+1 bins */
static void get_output(fctx *c, const float *o, double *X) {
  int n = c->n;
  int split = (c->kind == P_VFFT || c->kind == P_VDFT);
  if (c->kind == P_PFFFT_Z) { /* reorder a copy, for checking only */
    float *t = bench_alloc(sizeof(float) * 2 * n);
    pffft_zreorder(c->ps, o, t, PFFFT_FORWARD);
    fctx d = *c; d.kind = P_PFFFT_ORD;
    get_output(&d, t, X);
    free(t);
    return;
  }
  if (!c->real) {
    for (int k = 0; k < n; k++) {
      X[2 * k] = split ? o[k] : o[2 * k];
      X[2 * k + 1] = split ? o[n + k] : o[2 * k + 1];
    }
    return;
  }
  int h = n / 2;
  if (c->kind == P_FFTW) { for (int k = 0; k < 2 * (h + 1); k++) X[k] = o[k]; return; }
  double s = split ? 0.5 : 1.0; /* vDSP real FFTs return 2x the DFT */
  const float *re = split ? o : NULL, *im = split ? o + h : NULL;
  for (int k = 1; k < h; k++) {
    X[2 * k] = s * (split ? re[k] : o[2 * k]);
    X[2 * k + 1] = s * (split ? im[k] : o[2 * k + 1]);
  }
  X[0] = s * (split ? re[0] : o[0]); X[1] = 0;
  X[2 * h] = s * (split ? im[0] : o[1]); X[2 * h + 1] = 0;
}

static int setup(fctx *c, int kind, int real, int n, int mode) {
  memset(c, 0, sizeof *c);
  c->kind = kind; c->real = real; c->n = n;
  c->log2n = 0; while ((1 << c->log2n) < n) c->log2n++;
  bench_pool_init(&c->in, sizeof(float) * 2 * n, mode);
  bench_pool_init(&c->out, sizeof(float) * 2 * (n + 2), mode);
  float *i0 = c->in.buf[0], *o0 = c->out.buf[0];
  switch (kind) {
    case P_VFFT: c->fs = vDSP_create_fftsetup(c->log2n, kFFTRadix2); return c->fs != NULL;
    case P_VDFT:
      c->ds = real ? vDSP_DFT_zrop_CreateSetup(NULL, n, vDSP_DFT_FORWARD) : vDSP_DFT_zop_CreateSetup(NULL, n, vDSP_DFT_FORWARD);
      return c->ds != NULL;
    case P_VDFTI:
      if (real) return 0;
      c->dsi = vDSP_DFT_Interleaved_CreateSetup(NULL, n, vDSP_DFT_FORWARD, vDSP_DFT_Interleaved_ComplextoComplex);
      return c->dsi != NULL;
    case P_FFTW: {
      uint64_t t0 = bench_now_ns();
      c->fp = real ? fftwf_plan_dft_r2c_1d(n, i0, (fftwf_complex *)o0, FFTW_MEASURE)
                   : fftwf_plan_dft_1d(n, (fftwf_complex *)i0, (fftwf_complex *)o0, FFTW_FORWARD, FFTW_MEASURE);
      printf("# FFTW_MEASURE plan n=%d %s: %.3f s\n", n, real ? "r2c" : "c2c", (bench_now_ns() - t0) * 1e-9);
      return c->fp != NULL;
    }
    case P_PFFFT_ORD: case P_PFFFT_Z:
      if (!pffft_is_valid_size(n, real ? PFFFT_REAL : PFFFT_COMPLEX)) return 0;
      c->ps = pffft_new_setup(n, real ? PFFFT_REAL : PFFFT_COMPLEX);
      c->work = pffft_aligned_malloc(sizeof(float) * 2 * n);
      return c->ps != NULL;
  }
  return 0;
}

int main(int argc, char **argv) {
  bench_init();
  if (argc < 4) { fprintf(stderr, "usage: fft_bench c|r same|separate N...\n"); return 2; }
  int real = argv[1][0] == 'r';
  int mode = !strcmp(argv[2], "separate") ? BENCH_SEPARATE : BENCH_SAME;
  const char *wis = getenv("FFTW_WISDOM");
  if (wis) printf("# wisdom import %s: %d\n", wis, fftwf_import_wisdom_from_filename(wis));
  for (int a = 3; a < argc; a++) {
    int n = atoi(argv[a]);
    int nb = real ? n / 2 + 1 : n;
    double *x = malloc(sizeof(double) * 2 * n), *R = malloc(sizeof(double) * 2 * (n + 2)), *X = malloc(sizeof(double) * 2 * (n + 2));
    bench_fill_d(x, real ? n : 2 * n);
    for (int k = 0; k < (real ? n : 2 * n); k++) x[k] = (float)x[k]; /* reference sees the fp32 input */
    { /* reference: FFTW double, ESTIMATE */
      double *ti = fftw_malloc(sizeof(double) * 2 * n);
      memcpy(ti, x, sizeof(double) * (real ? n : 2 * n));
      fftw_plan pr = real ? fftw_plan_dft_r2c_1d(n, ti, (fftw_complex *)R, FFTW_ESTIMATE)
                          : fftw_plan_dft_1d(n, (fftw_complex *)ti, (fftw_complex *)R, FFTW_FORWARD, FFTW_ESTIMATE);
      memcpy(ti, x, sizeof(double) * (real ? n : 2 * n));
      fftw_execute(pr);
      fftw_destroy_plan(pr); fftw_free(ti);
      if (n <= 4096) { /* naive DFT check of the reference */
        double emax = 0, rmax = 0;
        for (int k = 0; k < nb; k++) {
          double sr = 0, si = 0;
          for (int j = 0; j < n; j++) {
            double ang = -2 * M_PI * (double)((long)j * k % n) / n, xr = real ? x[j] : x[2 * j], xi = real ? 0 : x[2 * j + 1];
            sr += xr * cos(ang) - xi * sin(ang); si += xr * sin(ang) + xi * cos(ang);
          }
          emax = fmax(emax, fmax(fabs(sr - R[2 * k]), fabs(si - R[2 * k + 1])));
          rmax = fmax(rmax, fmax(fabs(sr), fabs(si)));
        }
        printf("CHECK reference FFTW-double vs naive DFT n=%d: rel err %.3e %s\n", n, emax / rmax, emax / rmax < 1e-12 ? "PASS" : "FAIL");
      }
    }
    fctx c[NPROG];
    bench_prog p[NPROG];
    int np = 0;
    for (int kind = 0; kind < NPROG; kind++) {
      if (!setup(&c[np], kind, real, n, mode)) { printf("# skip %s n=%d (unsupported)\n", pname[kind], n); continue; }
      for (int b = 0; b < c[np].in.count; b++) put_input(&c[np], c[np].in.buf[b], x);
      f_run(&c[np], 0);
      get_output(&c[np], c[np].out.buf[0], X);
      char lab[96];
      snprintf(lab, sizeof lab, "%s %s n=%d", pname[kind], real ? "real" : "cplx", n);
      bench_check_d(lab, X, R, 2 * nb, 1e-5);
      p[np] = (bench_prog){pname[kind], f_run, NULL, &c[np]};
      np++;
    }
    bench_opts o = BENCH_OPTS_DEFAULT;
    o.units = n;
    bench_stat st[NPROG];
    char title[64];
    snprintf(title, sizeof title, "FFT fp32 %s n=%d %s", real ? "real" : "cplx", n, mode ? "separate" : "same");
    bench_compare(title, p, np, &o, st);
    for (int j = 0; j < np; j++) { bench_pool_free(&c[j].in); bench_pool_free(&c[j].out); }
    free(x); free(R); free(X);
  }
  const char *wout = getenv("FFTW_WISDOM_OUT");
  if (wout) printf("# wisdom export %s: %d\n", wout, fftwf_export_wisdom_to_filename(wout));
  return 0;
}
