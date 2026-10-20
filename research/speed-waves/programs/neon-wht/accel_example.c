/* accel_example.c -- small working examples of Apple Accelerate (vDSP) for the speed waves, each
 * checked against a plain loop: complex FFT (vDSP_fft_zop and vDSP_DFT_Execute), real FFT
 * (vDSP_fft_zrop: packed, 2x scaled), FIR / correlation (vDSP_conv), decimating FIR (vDSP_desamp).
 * Build: clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations accel_example.c -framework Accelerate -o accel_example */
#include <Accelerate/Accelerate.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double rnd(void) { return rand() / (double)RAND_MAX * 2 - 1; }
static void report(const char *what, double err, double scale) {
  printf("%-44s rel err %.2e %s\n", what, err / scale, err / scale < 1e-5 ? "PASS" : "FAIL");
}

int main(void) {
  enum { LOG2N = 10, N = 1 << LOG2N };
  static float re[N], im[N], ore[N], oim[N];
  static double xr[N], xi[N], Xr[N], Xi[N];
  for (int i = 0; i < N; i++) { re[i] = xr[i] = (float)rnd(); im[i] = xi[i] = (float)rnd(); }
  for (int k = 0; k < N; k++) { /* naive DFT reference */
    double sr = 0, si = 0;
    for (int j = 0; j < N; j++) {
      double a = -2 * M_PI * (double)((j * k) % N) / N;
      sr += xr[j] * cos(a) - xi[j] * sin(a);
      si += xr[j] * sin(a) + xi[j] * cos(a);
    }
    Xr[k] = sr; Xi[k] = si;
  }
  double mx = 0;
  for (int k = 0; k < N; k++) mx = fmax(mx, fmax(fabs(Xr[k]), fabs(Xi[k])));

  /* 1. complex FFT, classic API: setup once per max size; split complex layout */
  FFTSetup fs = vDSP_create_fftsetup(LOG2N, kFFTRadix2);
  DSPSplitComplex in = {re, im}, out = {ore, oim};
  vDSP_fft_zop(fs, &in, 1, &out, 1, LOG2N, kFFTDirection_Forward); /* unscaled forward */
  double e = 0;
  for (int k = 0; k < N; k++) e = fmax(e, fmax(fabs(ore[k] - Xr[k]), fabs(oim[k] - Xi[k])));
  report("vDSP_fft_zop complex N=1024", e, mx);

  /* 2. complex DFT API (any N = f*2^k, f in {1,3,5,15}); usually faster */
  vDSP_DFT_Setup ds = vDSP_DFT_zop_CreateSetup(NULL, N, vDSP_DFT_FORWARD);
  vDSP_DFT_Execute(ds, re, im, ore, oim);
  e = 0;
  for (int k = 0; k < N; k++) e = fmax(e, fmax(fabs(ore[k] - Xr[k]), fabs(oim[k] - Xi[k])));
  report("vDSP_DFT_Execute complex N=1024", e, mx);
  vDSP_DFT_DestroySetup(ds);

  /* 3. real FFT: input x[0..N) packed as even/odd into split (vDSP_ctoz); output packed:
        realp[0] = 2*X[0], imagp[0] = 2*X[N/2], (realp[k], imagp[k]) = 2*X[k] for 0<k<N/2 */
  static float x[N], er[N / 2], ei[N / 2], fr[N / 2], fi[N / 2];
  static double R[N / 2 + 1][2];
  for (int i = 0; i < N; i++) x[i] = (float)rnd();
  for (int k = 0; k <= N / 2; k++) {
    double sr = 0, si = 0;
    for (int j = 0; j < N; j++) { double a = -2 * M_PI * (double)((j * k) % N) / N; sr += x[j] * cos(a); si += x[j] * sin(a); }
    R[k][0] = sr; R[k][1] = si;
  }
  DSPSplitComplex ev = {er, ei}, fo = {fr, fi};
  vDSP_ctoz((const DSPComplex *)x, 2, &ev, 1, N / 2);
  vDSP_fft_zrop(fs, &ev, 1, &fo, 1, LOG2N, kFFTDirection_Forward);
  e = 0; mx = 0;
  for (int k = 1; k < N / 2; k++) {
    e = fmax(e, fmax(fabs(0.5 * fr[k] - R[k][0]), fabs(0.5 * fi[k] - R[k][1])));
    mx = fmax(mx, fabs(R[k][0]));
  }
  e = fmax(e, fmax(fabs(0.5 * fr[0] - R[0][0]), fabs(0.5 * fi[0] - R[N / 2][0])));
  report("vDSP_fft_zrop real N=1024 (x0.5, packed)", e, mx);
  vDSP_destroy_fftsetup(fs);

  /* 4. vDSP_conv: C[n] = sum_k A[n+k] * F[k*IF] (a CORRELATION). For a convolution pass a pointer
        to the LAST filter tap and stride -1. A must have N+P-1 samples. */
  enum { P = 31, NO = 1000 };
  static float A[NO + P - 1], F[P], C[NO];
  for (int i = 0; i < NO + P - 1; i++) A[i] = (float)rnd();
  for (int k = 0; k < P; k++) F[k] = (float)rnd();
  vDSP_conv(A, 1, F + P - 1, -1, C, 1, NO, P); /* convolution: C[n] = sum_k A[n+k] F[P-1-k] */
  e = 0; mx = 0;
  for (int n = 0; n < NO; n++) {
    double s = 0;
    for (int k = 0; k < P; k++) s += (double)A[n + k] * F[P - 1 - k];
    e = fmax(e, fabs(C[n] - s)); mx = fmax(mx, fabs(s));
  }
  report("vDSP_conv FIR 31 taps x 1000 outputs", e, mx);

  /* 5. vDSP_desamp: decimating FIR, C[n] = sum_p A[n*D + p] * F[p] (correlation form) */
  enum { D = 4, ND = 200 };
  static float Cd[ND];
  vDSP_desamp(A, D, F, Cd, ND, P);
  e = 0; mx = 0;
  for (int n = 0; n < ND; n++) {
    double s = 0;
    for (int p = 0; p < P; p++) s += (double)A[n * D + p] * F[p];
    e = fmax(e, fabs(Cd[n] - s)); mx = fmax(mx, fabs(s));
  }
  report("vDSP_desamp decimate by 4, 31 taps", e, mx);
  return 0;
}
