/* ctr-window / lpmwu.c (own code). Weighted covering bound, one-sided accounting, for single-layer decoders.
   Pool: P line sets S, each needing >= 8 terms whose gather set meets S (exact 8-frames: every realisation;
   bases of R^8: rank-8 realisations, Gram up to gauge by #1212). A term with gather-label span W costs dim W and
   meets S iff S contains a line with label in W. For any distribution mu on the pool,
       cost >= 8 / z(mu),   z(mu) = max_W  mu{S : S meets L(W)} / dim W   (max over EVERY subspace W of F2^9).
   mu is tuned by multiplicative weights over a candidate set of W, then z(mu) is computed EXACTLY by enumerating
   all subspaces (RREF). usage: lpmwu pool.txt rounds iters eta  (pool format: n P F / labels / P lines of 8). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define MAXP 8192
#define PW (MAXP/64)
static int n, P, F, lab[112], lineof[512], pw;
static uint64_t sm[112][PW];
static double mu[MAXP], tab[MAXP/8][256];
static void buildtab(void){ int nb=(P+7)/8; for(int b=0;b<nb;b++){tab[b][0]=0; for(int v=1;v<256;v++){int low=__builtin_ctz(v); int s=b*8+low; tab[b][v]=tab[b][v&(v-1)]+(s<P?mu[s]:0);} } }
static double whit(const uint64_t *m){ double s=0; int nb=(P+7)/8; for(int b=0;b<nb;b++){ s+=tab[b][(m[b/8]>>((b%8)*8))&255]; } return s; }
/* candidates */
#define MAXC 20000
static int nc; static uint64_t cm[MAXC][PW]; static int cd[MAXC]; static uint64_t cl0[MAXC], cl1[MAXC];
static int havecand(uint64_t a,uint64_t b,int d){ for(int i=0;i<nc;i++) if(cl0[i]==a&&cl1[i]==b&&cd[i]<=d) return 1; return 0; }
/* full enumeration */
#define TOP 400
static double topr[TOP]; static uint64_t top0[TOP], top1[TOP]; static int topd[TOP]; static int ntop;
static double zmax; static int zd; static uint64_t z0,z1; static int bestbasis[9];
static int piv[9], rows[9];
static void evalW(int d){
  uint64_t l0=0,l1=0; int els=1<<d;
  for(int s=1;s<els;s++){int x=0;for(int r=0;r<d;r++)if(s>>r&1)x^=rows[r];int t=lineof[x];if(t>=0){if(t<64)l0|=1ULL<<t;else l1|=1ULL<<(t-64);}}
  if(!l0&&!l1)return;
  uint64_t m[PW]; memset(m,0,sizeof(uint64_t)*pw);
  for(int t=0;t<n;t++) if((t<64?(l0>>t):(l1>>(t-64)))&1) for(int i=0;i<pw;i++) m[i]|=sm[t][i];
  double r=whit(m)/d;
  if(r>zmax){zmax=r;zd=d;z0=l0;z1=l1;for(int q=0;q<d;q++)bestbasis[q]=rows[q];}
  /* keep top list */
  if(ntop<TOP){topr[ntop]=r;top0[ntop]=l0;top1[ntop]=l1;topd[ntop]=d;ntop++;}
  else { int mi=0; for(int i=1;i<TOP;i++) if(topr[i]<topr[mi]) mi=i; if(r>topr[mi]){topr[mi]=r;top0[mi]=l0;top1[mi]=l1;topd[mi]=d;} }
}
static void fill(int r,int d){ if(r==d){evalW(d);return;} int fb[9],nf=0; for(int b=0;b<piv[r];b++){int isp=0;for(int q=0;q<d;q++)if(piv[q]==b)isp=1;if(!isp)fb[nf++]=b;}
  for(int mask=0;mask<(1<<nf);mask++){int v=1<<piv[r];for(int i=0;i<nf;i++)if(mask>>i&1)v|=1<<fb[i];rows[r]=v;fill(r+1,d);} }
static void pivots(int r,int d,int below){ if(r==d){fill(0,d);return;} for(int p=below-1;p>=d-r-1;p--){piv[r]=p;pivots(r+1,d,p);} }
static void fullcheck(void){ buildtab(); zmax=0; ntop=0; for(int d=1;d<=9;d++) pivots(0,d,9); }
static void addcand(uint64_t l0,uint64_t l1,int d){ if(nc>=MAXC||havecand(l0,l1,d))return; cl0[nc]=l0;cl1[nc]=l1;cd[nc]=d; memset(cm[nc],0,sizeof(uint64_t)*pw);
  for(int t=0;t<n;t++) if((t<64?(l0>>t):(l1>>(t-64)))&1) for(int i=0;i<pw;i++) cm[nc][i]|=sm[t][i]; nc++; }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int rounds=atoi(argv[2]), iters=atoi(argv[3]); double eta=atof(argv[4]);
  if(fscanf(f,"%d %d %d",&n,&P,&F)!=3)return 1; pw=(P+63)/64; memset(lineof,-1,sizeof lineof);
  for(int i=0;i<n;i++){if(fscanf(f,"%d",&lab[i])!=1)return 1;lineof[lab[i]]=i;}
  for(int s=0;s<P;s++)for(int j=0;j<8;j++){int t;if(fscanf(f,"%d",&t)!=1)return 1;sm[t][s/64]|=1ULL<<(s%64);}
  double bestfloor=0; static double bestmu[MAXP];
  /* round 0: uniform mu, frames only and all sets */
  for(int s=0;s<P;s++)mu[s]=1.0/P;
  fullcheck(); printf("uniform over pool: z=%.6f floor=%.4f (worst W dim %d)\n",zmax,8/zmax,zd); fflush(stdout);
  if(8/zmax>bestfloor){bestfloor=8/zmax;memcpy(bestmu,mu,sizeof(double)*P);}
  for(int i=0;i<ntop;i++) addcand(top0[i],top1[i],topd[i]);
  static double w[MAXP], avg[MAXP];
  for(int rd=0;rd<rounds;rd++){
    for(int s=0;s<P;s++){w[s]=bestmu[s];avg[s]=0;}
    for(int it=0;it<iters;it++){
      double tot=0; for(int s=0;s<P;s++)tot+=w[s]; for(int s=0;s<P;s++){mu[s]=w[s]/tot; avg[s]+=mu[s];}
      buildtab(); int bc=0; double br=-1; for(int c=0;c<nc;c++){double r=whit(cm[c])/cd[c]; if(r>br){br=r;bc=c;}}
      double fac=exp(-eta); for(int s=0;s<P;s++) if((cm[bc][s/64]>>(s%64))&1) w[s]*=fac;
    }
    for(int s=0;s<P;s++)mu[s]=avg[s]/iters;
    fullcheck(); double fl=8/zmax;
    printf("round %d: candidates %d, exact z=%.6f -> floor %.4f (worst W dim %d, |L(W)| %d)\n",rd,nc,zmax,fl,zd,__builtin_popcountll(z0)+__builtin_popcountll(z1)); fflush(stdout);
    if(fl>bestfloor){bestfloor=fl;memcpy(bestmu,mu,sizeof(double)*P);}
    for(int i=0;i<ntop;i++) addcand(top0[i],top1[i],topd[i]);
  }
  /* final exact check of the best mu, report */
  memcpy(mu,bestmu,sizeof(double)*P); fullcheck();
  double fm=0; for(int s=0;s<F;s++)fm+=mu[s];
  printf("BEST: exact z=%.8f floor=%.4f per call; mass on frames %.4f, on bases %.4f; worst W dim %d basis:",zmax,8/zmax,fm,1-fm,zd);
  for(int q=0;q<zd;q++)printf(" %d",bestbasis[q]); printf("\n");
  FILE*o=fopen("lpmwu_mu.txt","w"); for(int s=0;s<P;s++)fprintf(o,"%.12g\n",mu[s]); fclose(o);
  return 0;
}
