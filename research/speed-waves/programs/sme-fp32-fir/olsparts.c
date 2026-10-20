/* olsparts.c -- time the pieces of the vDSP overlap-save (forward DFT, zvmul, inverse DFT, ztoc, ctoz) */
#include <Accelerate/Accelerate.h>
#include "../timing-harness/bench.h"
typedef struct { int N; vDSP_DFT_Setup f, iv; float *a,*b,*c,*d,*Hr,*Hi,*buf; } C;
static void r_f(void *v, long i){ C*c=v; vDSP_DFT_Execute(c->f,c->a,c->b,c->c,c->d); }
static void r_i(void *v, long i){ C*c=v; vDSP_DFT_Execute(c->iv,c->c,c->d,c->a,c->b); }
static void r_m(void *v, long i){ C*c=v; DSPSplitComplex t={c->c,c->d},H={c->Hr,c->Hi}; vDSP_zvmul(&t,1,&H,1,&t,1,c->N/2,1); }
static void r_z(void *v, long i){ C*c=v; DSPSplitComplex s={c->a,c->b}; vDSP_ztoc(&s,1,(DSPComplex*)(c->buf+1),2,c->N/2-1); }
static void r_cz(void *v, long i){ C*c=v; DSPSplitComplex s={c->a,c->b}; vDSP_ctoz((DSPComplex*)c->buf,2,&s,1,c->N/2); }
static float *E,*F;
static void r_chain(void *v, long i){ C*c=v; vDSP_DFT_Execute(c->f,c->a,c->b,c->c,c->d);
  DSPSplitComplex t={c->c,c->d},H={c->Hr,c->Hi}; vDSP_zvmul(&t,1,&H,1,&t,1,c->N/2,1);
  vDSP_DFT_Execute(c->iv,c->c,c->d,E,F); DSPSplitComplex s={E,F}; vDSP_ztoc(&s,1,(DSPComplex*)(c->buf+1),2,c->N/2-1); }
static void r_chain2(void *v, long i){ C*c=v; vDSP_DFT_Execute(c->f,c->a,c->b,c->c,c->d);
  vDSP_DFT_Execute(c->iv,c->c,c->d,E,F); }
static float *G,*Hh;
#include <arm_neon.h>
static void cmul_neon(float *restrict ar, float *restrict ai, const float *restrict br, const float *restrict bi, float *restrict cr, float *restrict ci, int n){
  for(int k=0;k<n;k+=4){ float32x4_t xr=vld1q_f32(ar+k), xi=vld1q_f32(ai+k), hr=vld1q_f32(br+k), hi=vld1q_f32(bi+k);
    vst1q_f32(cr+k, vfmsq_f32(vmulq_f32(xr,hr), xi,hi)); vst1q_f32(ci+k, vfmaq_f32(vmulq_f32(xr,hi), xi,hr)); } }
static void r_chain3(void *v, long i){ C*c=v; vDSP_DFT_Execute(c->f,c->a,c->b,c->c,c->d);
  DSPSplitComplex t={c->c,c->d},H={c->Hr,c->Hi},o={G,Hh}; vDSP_zvmul(&t,1,&H,1,&o,1,c->N/2,1);
  vDSP_DFT_Execute(c->iv,G,Hh,E,F); DSPSplitComplex s={E,F}; vDSP_ztoc(&s,1,(DSPComplex*)(c->buf+1),2,c->N/2-1); }
static void r_chain4(void *v, long i){ C*c=v; vDSP_DFT_Execute(c->f,c->a,c->b,c->c,c->d);
  cmul_neon(c->c,c->d,c->Hr,c->Hi,G,Hh,c->N/2);
  vDSP_DFT_Execute(c->iv,G,Hh,E,F); DSPSplitComplex s={E,F}; vDSP_ztoc(&s,1,(DSPComplex*)(c->buf+1),2,c->N/2-1); }
static void r_chain5(void *v, long i){ C*c=v; vDSP_DFT_Execute(c->f,c->a,c->b,c->c,c->d);
  cmul_neon(c->c,c->d,c->Hr,c->Hi,G,Hh,c->N/2);
  vDSP_DFT_Execute(c->iv,G,Hh,E,F); }
int main(int argc,char**argv){
  bench_init(); int N=atoi(argv[1]);
  C c; c.N=N; c.f=vDSP_DFT_zrop_CreateSetup(NULL,N,vDSP_DFT_FORWARD); c.iv=vDSP_DFT_zrop_CreateSetup(c.f,N,vDSP_DFT_INVERSE);
  c.a=bench_alloc(2*N);c.b=bench_alloc(2*N);c.c=bench_alloc(2*N);c.d=bench_alloc(2*N);c.Hr=bench_alloc(2*N);c.Hi=bench_alloc(2*N);c.buf=bench_alloc(4*N+64);
  bench_fill_f(c.a,N/2);bench_fill_f(c.b,N/2);bench_fill_f(c.Hr,N/2);bench_fill_f(c.Hi,N/2);bench_fill_f(c.buf,N);
  E=bench_alloc(2*N);F=bench_alloc(2*N);G=bench_alloc(2*N);Hh=bench_alloc(2*N);
  bench_prog p[10]={{"chain_zvmul_oop",r_chain3,NULL,&c},{"chain_neoncmul",r_chain4,NULL,&c},{"chain_neoncmul_noztoc",r_chain5,NULL,&c},{"chain",r_chain,NULL,&c},{"fwd+inv",r_chain2,NULL,&c},{"fwdDFT",r_f,NULL,&c},{"invDFT",r_i,NULL,&c},{"zvmul",r_m,NULL,&c},{"ztoc",r_z,NULL,&c},{"ctoz",r_cz,NULL,&c}};
  bench_opts o=BENCH_OPTS_DEFAULT; o.units=N; o.print_reps=0; bench_stat st[10]; char t[64]; snprintf(t,64,"olsparts N=%d",N);
  bench_compare(t,p,10,&o,st); return 0; }
