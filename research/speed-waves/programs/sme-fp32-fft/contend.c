// Does vDSP's 256-point complex DFT use the P-cluster's SME unit? Time it with and without a second
// P-core thread that saturates the SME unit with FMOPAs. Our SME kernel is the positive control.
#define main main_unused
#include "fft256.c"
#undef main
#include <pthread.h>
static volatile int stop_hog = 0;
__arm_locally_streaming __arm_new("za") static void hog_sme(void) {
  svbool_t pg = svptrue_b32(); svfloat32_t a = svdup_f32(1e-3f), b = svdup_f32(1e-3f);
  while (!stop_hog) for (int i = 0; i < 100000; i++) {
    svmopa_za32_f32_m(0,pg,pg,a,b); svmopa_za32_f32_m(1,pg,pg,a,b); svmopa_za32_f32_m(2,pg,pg,a,b); svmopa_za32_f32_m(3,pg,pg,a,b); }
}
static void *hog(void *arg) { (void)arg; pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0); hog_sme(); return NULL; }
static double time_prog(ctx *c, long B) {
  double best = 1e30;
  for (int r = 0; r < 7; r++) { long it = 400; uint64_t t0 = bench_now_ns(); for (long i = 0; i < it; i++) prun(c, 0); double ns = (bench_now_ns() - t0) / (double)(it * B * N); if (ns < best) best = ns; }
  return best;
}
int main(void) {
  bench_init(); make_tables();
  long B = 64; int kinds[] = {P_SME5, P_VDFT, P_VFFT, P_PFFFTZ};
  ctx c[4]; for (int k = 0; k < 4; k++) setup(&c[k], kinds[k], B, BENCH_SAME);
  for (int round = 0; round < 2; round++) {
    double alone[4], with[4];
    for (int k = 0; k < 4; k++) alone[k] = time_prog(&c[k], B);
    stop_hog = 0; pthread_t th; pthread_create(&th, NULL, hog, NULL);
    for (volatile long w = 0; w < 20000000; w++);
    for (int k = 0; k < 4; k++) with[k] = time_prog(&c[k], B);
    stop_hog = 1; pthread_join(th, NULL);
    for (int k = 0; k < 4; k++) printf("round %d %-28s alone %.4f  with SME hog on another P-core %.4f ns/pt  ratio %.2f\n", round, pname[kinds[k]], alone[k], with[k], with[k] / alone[k]);
  }
}
