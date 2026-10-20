// ap-0 probe: complex 8x8 fp64 DFT matrix times an 8 x 64 complex panel (split re/im) on the SME unit.
// 8 ZA fp64 tiles = 4 complex 8x8 blocks in flight; 32 FMOPA per complex 8x8 block.
#include <arm_sme.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
static uint64_t now(void){ return clock_gettime_nsec_np(CLOCK_UPTIME_RAW); }
#define MO(t,a,b) svmopa_za64_f64_m(t, pg, pg, a, b)
#define MS(t,a,b) svmops_za64_f64_m(t, pg, pg, a, b)
// F: Fr[a][k] column a of F stored contiguously (8 doubles), Fi likewise. X: Xr[a][64], Xi[a][64]. out: 8 x 64 re, im
__arm_locally_streaming __arm_new("za") void panel(const double *Fr, const double *Fi, const double *Xr, const double *Xi,
     double *Or, double *Oi, long reps, int store) {
  svbool_t pg = svptrue_b64();
  for (long r = 0; r < reps; r++) {
    for (int h = 0; h < 64; h += 32) {
      svzero_za();
      for (int a = 0; a < 8; a++) {
        svfloat64_t fr = svld1_f64(pg, Fr + 8*a), fi = svld1_f64(pg, Fi + 8*a);
        const double *xr = Xr + 64*a + h, *xi = Xi + 64*a + h;
        svfloat64_t r0 = svld1_f64(pg, xr), r1 = svld1_f64(pg, xr+8), r2 = svld1_f64(pg, xr+16), r3 = svld1_f64(pg, xr+24);
        svfloat64_t i0 = svld1_f64(pg, xi), i1 = svld1_f64(pg, xi+8), i2 = svld1_f64(pg, xi+16), i3 = svld1_f64(pg, xi+24);
        MO(0,fr,r0); MO(1,fr,i0); MO(2,fr,r1); MO(3,fr,i1); MO(4,fr,r2); MO(5,fr,i2); MO(6,fr,r3); MO(7,fr,i3);
        MS(0,fi,i0); MO(1,fi,r0); MS(2,fi,i1); MO(3,fi,r1); MS(4,fi,i2); MO(5,fi,r2); MS(6,fi,i3); MO(7,fi,r3);
      }
      if (store) {
        for (int s = 0; s < 8; s++) {
          svst1_hor_za64(0, s, pg, Or + 64*s + h);    svst1_hor_za64(1, s, pg, Oi + 64*s + h);
          svst1_hor_za64(2, s, pg, Or + 64*s + h+8);  svst1_hor_za64(3, s, pg, Oi + 64*s + h+8);
          svst1_hor_za64(4, s, pg, Or + 64*s + h+16); svst1_hor_za64(5, s, pg, Oi + 64*s + h+16);
          svst1_hor_za64(6, s, pg, Or + 64*s + h+24); svst1_hor_za64(7, s, pg, Oi + 64*s + h+24);
        }
      }
    }
  }
}
// NEON reference rate: same complex product in plain C (auto-vectorised fp64, 2 lanes)
static void panel_c(const double *Fr, const double *Fi, const double *Xr, const double *Xi, double *Or, double *Oi) {
  for (int k = 0; k < 8; k++) {
    double sr[64] = {0}, si[64] = {0};
    for (int a = 0; a < 8; a++) { double fr = Fr[8*a+k], fi = Fi[8*a+k];
      for (int c = 0; c < 64; c++) { sr[c] += fr*Xr[64*a+c] - fi*Xi[64*a+c]; si[c] += fr*Xi[64*a+c] + fi*Xr[64*a+c]; } }
    for (int c = 0; c < 64; c++) { Or[64*k+c] = sr[c]; Oi[64*k+c] = si[c]; }
  }
}
int main(void) {
  pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
  static double Fr[64], Fi[64], Xr[512], Xi[512], Or[512], Oi[512], Rr[512], Ri[512];
  for (int a = 0; a < 8; a++) for (int k = 0; k < 8; k++) { Fr[8*a+k] = cos(-2*M_PI*a*k/8); Fi[8*a+k] = sin(-2*M_PI*a*k/8); }
  srand(1); for (int i = 0; i < 512; i++) { Xr[i] = rand()/(double)RAND_MAX - .5; Xi[i] = rand()/(double)RAND_MAX - .5; }
  panel(Fr, Fi, Xr, Xi, Or, Oi, 1, 1); panel_c(Fr, Fi, Xr, Xi, Rr, Ri);
  double e = 0; for (int i = 0; i < 512; i++) e = fmax(e, fmax(fabs(Or[i]-Rr[i]), fabs(Oi[i]-Ri[i])));
  printf("check max abs err %.3e %s\n", e, e < 1e-13 ? "PASS" : "FAIL");
  double macs = 32.0 * 64 * 8; // real MACs per panel (8 complex 8x8 blocks, 32 FMOPA x 64 MAC each... = 4 real products x 8x8x8 x 8 blocks)
  macs = 4.0 * 8 * 8 * 64;     // 4 real 8x8 products per complex, panel 8x8 times 8x64
  for (int store = 0; store < 2; store++) for (int rep = 0; rep < 5; rep++) {
    long reps = 200000; uint64_t t = now(); panel(Fr, Fi, Xr, Xi, Or, Oi, reps, store); double ns = (double)(now()-t)/reps;
    printf("SME panel store=%d: %.2f ns/panel  %.1f GMAC/s  %.4f ns per complex point per radix-8 pass\n", store, ns, macs/ns, ns/512);
  }
  for (int rep = 0; rep < 5; rep++) {
    long reps = 100000; uint64_t t = now(); for (long r = 0; r < reps; r++) { panel_c(Fr, Fi, Xr, Xi, Or, Oi); __asm__ volatile("" ::: "memory"); } double ns = (double)(now()-t)/reps;
    printf("NEON/C panel: %.2f ns/panel  %.1f GMAC/s  %.4f ns per complex point per pass\n", ns, macs/ns, ns/512);
  }
  return 0;
}
